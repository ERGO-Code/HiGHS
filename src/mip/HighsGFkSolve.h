/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsGFkLU.h
 * @brief linear system solve in GF(k) for mod-k cut separation
 * @author Leona Gottwald
 */

#ifndef HIGHS_GFk_SOLVE_H_
#define HIGHS_GFk_SOLVE_H_

#include <algorithm>
#include <cassert>
#include <queue>
#include <tuple>
#include <vector>

#include "lp_data/HConst.h"

// helper struct to compute the multipicative inverse by using fermats
// theorem and recursive repeated squaring.
// Under the assumption that k is a small prime and an 32bit int is enough
// to hold the number (k-1)^(k-2) good compilers should be able to optimize this
// code by inlining and unfolding the recursion to code that uses the fewest
// amount of integer multiplications and optimize the single integer division
// away as k is a compile time constant. Since the class was developed for the
// purpose of separating maximally violated mod-k cuts the assumption that k is
// a small constant prime number is not restrictive.

template <int k>
struct HighsGFk;

template <>
struct HighsGFk<2> {
  static constexpr unsigned int powk(unsigned int a) { return a * a; }
  static constexpr unsigned int inverse(unsigned int a) { return 1; }
};

template <>
struct HighsGFk<3> {
  static constexpr unsigned int powk(unsigned int a) { return a * a * a; }
  static constexpr unsigned int inverse(unsigned int a) { return a; }
};

template <int k>
struct HighsGFk {
  static constexpr unsigned int powk(unsigned int a) {
    return (k & 1) == 0 ? HighsGFk<2>::powk(HighsGFk<k / 2>::powk(a))
                        : HighsGFk<k - 1>::powk(a) * a;
  }

  static unsigned int inverse(unsigned int a) {
    return HighsGFk<k - 2>::powk(a) % k;
  }
};

class HighsGFkSolve {
  int numCol;
  int numRow;

  // triplet format
  std::vector<int> Arow;
  std::vector<int> Acol;
  std::vector<unsigned int> Avalue;

  // sizes of rows and columns
  std::vector<int> rowsize;
  std::vector<int> colsize;

  // linked list links for column based links for each nonzero
  std::vector<int> colhead;
  std::vector<int> Anext;
  std::vector<int> Aprev;

  // splay tree links for row based iteration and lookup
  std::vector<int> rowroot;
  std::vector<int> ARleft;
  std::vector<int> ARright;

  // right hand side vector
  std::vector<unsigned int> rhs;

  // column permutation for the factorization required for backwards solve
  std::vector<int> factorColPerm;
  std::vector<int> factorRowPerm;
  std::vector<int8_t> colBasisStatus;
  std::vector<int8_t> rowUsed;

  // working memory
  std::vector<int> iterstack;
  std::vector<int> rowpositions;
  std::vector<int> rowposColsizes;

  // priority queue to reuse free slots
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;

  void link(int pos);

  void unlink(int pos);

  void dropIfZero(int pos) {
    if (Avalue[pos] == 0) unlink(pos);
  }

  void storeRowPositions(int pos);

  void addNonzero(int row, int col, unsigned int val);

 public:
  // access to triplets and find nonzero function for unit test
  const std::vector<int>& getArow() const { return Arow; }
  const std::vector<int>& getAcol() const { return Acol; }
  const std::vector<unsigned>& getAvalue() const { return Avalue; }
  int numNonzeros() const { return int(Avalue.size() - freeslots.size()); }
  int findNonzero(int row, int col);

  template <unsigned int k, typename T>
  void fromCSC(const std::vector<T>& Aval, const std::vector<int>& Aindex,
               const std::vector<int>& Astart, int numRow) {
    Avalue.clear();
    Acol.clear();
    Arow.clear();

    freeslots = decltype(freeslots)();

    numCol = Astart.size() - 1;
    this->numRow = numRow;

    colhead.assign(numCol, -1);
    colsize.assign(numCol, 0);

    rhs.assign(numRow, 0);
    rowroot.assign(numRow, -1);
    rowsize.assign(numRow, 0);

    Avalue.reserve(Aval.size());
    Acol.reserve(Aval.size());
    Arow.reserve(Aval.size());

    for (int i = 0; i != numCol; ++i) {
      for (int j = Astart[i]; j != Astart[i + 1]; ++j) {
        assert(Aval[j] == (int64_t)Aval[j]);
        int64_t val = ((int64_t)Aval[j]) % k;
        if (val == 0) continue;

        if (val < 0) val += k;
        assert(val >= 0);

        Avalue.push_back(val);
        Acol.push_back(i);
        Arow.push_back(Aindex[j]);
      }
    }

    int nnz = Avalue.size();
    Anext.resize(nnz);
    Aprev.resize(nnz);
    ARleft.resize(nnz);
    ARright.resize(nnz);
    for (int pos = 0; pos != nnz; ++pos) link(pos);
  }

  template <unsigned int k, typename T>
  void setRhs(int row, T val) {
    rhs[row] = ((unsigned int)std::abs(val)) % k;
  }

  template <unsigned int k, typename ReportSolution>
  void solve(ReportSolution&& reportSolution) {
    auto cmpPrio = [](const std::pair<int, int>& a,
                      const std::pair<int, int>& b) {
      return a.first > b.first;
    };
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                        decltype(cmpPrio)>
        pqueue(cmpPrio);

    for (int i = 0; i != numCol; ++i) pqueue.emplace(colsize[i], i);

    int maxPivot = std::min(numRow, numCol);
    factorColPerm.clear();
    factorRowPerm.clear();
    factorColPerm.reserve(maxPivot);
    factorRowPerm.reserve(maxPivot);
    colBasisStatus.assign(numCol, 0);
    rowUsed.assign(numRow, 0);
    int numPivot = 0;

    while (!pqueue.empty()) {
      int pivotCol;
      int oldColSize;

      std::tie(oldColSize, pivotCol) = pqueue.top();
      pqueue.pop();

      if (colsize[pivotCol] == 0) continue;

      assert(colBasisStatus[pivotCol] == 0);

      if (colsize[pivotCol] != oldColSize) {
        pqueue.emplace(colsize[pivotCol], pivotCol);
        continue;
      }

      int pivot = -1;
      int pivotRow = -1;
      int pivotRowLen = HIGHS_CONST_I_INF;
      for (int coliter = colhead[pivotCol]; coliter != -1;
           coliter = Anext[coliter]) {
        int row = Arow[coliter];
        if (rowUsed[row]) continue;
        if (rowsize[row] < pivotRowLen) {
          pivotRowLen = rowsize[row];
          pivotRow = row;
          pivot = coliter;
        }
      }

      assert(pivot != -1);
      assert(Acol[pivot] == pivotCol);
      assert(Arow[pivot] == pivotRow);
      assert(Avalue[pivot] > 0);
      assert(Avalue[pivot] < k);

      unsigned int pivotInverse = HighsGFk<k>::inverse(Avalue[pivot]);
      assert((Avalue[pivot] * pivotInverse) % k == 1);

      rowpositions.clear();
      rowposColsizes.clear();
      storeRowPositions(rowroot[pivotRow]);
      assert(pivotRowLen == (int)rowpositions.size());
      int next;
      for (int coliter = colhead[pivotCol]; coliter != -1; coliter = next) {
        next = Anext[coliter];
        if (coliter == pivot) continue;

        assert(Acol[coliter] == pivotCol);
        assert(Arow[coliter] != pivotRow);
        assert(Avalue[coliter] != 0);

        int row = Arow[coliter];
        if (rowUsed[row]) continue;

        unsigned int pivotRowScale = pivotInverse * (k - Avalue[coliter]);

        rhs[row] = (rhs[row] + pivotRowScale * rhs[pivotRow]) % k;

        for (int pivotRowPos : rowpositions) {
          int nonzeroPos = findNonzero(Arow[coliter], Acol[pivotRowPos]);

          if (nonzeroPos == -1) {
            assert(Acol[pivotRowPos] != pivotCol);
            unsigned int val = (pivotRowScale * Avalue[pivotRowPos]) % k;
            if (val != 0) addNonzero(row, Acol[pivotRowPos], val);
          } else {
            Avalue[nonzeroPos] =
                (Avalue[nonzeroPos] + pivotRowScale * Avalue[pivotRowPos]) % k;
            assert(Acol[pivotRowPos] != pivotCol || Avalue[nonzeroPos] == 0);
            dropIfZero(nonzeroPos);
          }
        }
      }

      ++numPivot;
      factorColPerm.push_back(pivotCol);
      factorRowPerm.push_back(pivotRow);
      colBasisStatus[pivotCol] = 1;
      rowUsed[pivotRow] = 1;
      if (numPivot == maxPivot) break;

      for (int i = 0; i != pivotRowLen; ++i) {
        assert(Arow[rowpositions[i]] == pivotRow);
        int col = Acol[rowpositions[i]];
        int oldsize = rowposColsizes[i];

        // we only want to count rows that are not used so far, so we need to
        // decrease the counter for all columns in the pivot row by one
        --colsize[col];

        // the column size should never get negative
        assert(colsize[col] >= 0);

        if (colsize[col] == 0) continue;

        // the pivot column should occur in zero unused rows
        assert(col != pivotCol);

        // only reinsert the column if the size is smaller, as it otherwise
        // either is already at the correct position in the queue, or is at a
        // too good position and will be lazily reinserted when it is extracted
        // from the queue and the size does not match
        if (colsize[col] < oldsize) pqueue.emplace(colsize[col], col);
      }
    }

    // check if a solution exists by scanning the linearly dependent rows for
    // nonzero right hand sides
    for (int i = 0; i != numRow; ++i) {
      // if the row was used it is linearly independent
      if (rowUsed[i] == 1) continue;

      // if the row is linearly dependent, the right hand side must be zero,
      // otherwise no solution exists
      if (rhs[i] != 0) return;
    }

    // now iterate a subset of the basic solutions.
    // When a column leaves the basis we do not allow it to enter again so that
    // we iterate at most one solution for each nonbasic column

    std::vector<std::pair<int, unsigned int>> solution;
    solution.reserve(numCol);
    int numFactorRows = factorRowPerm.size();

    // create vector for swapping different columns into the basis
    // For each column we want to iterate one basic solution where the
    // column is basic
    std::vector<std::pair<int, int>> basisSwaps;
    assert(iterstack.empty());
    for (int i = numFactorRows - 1; i >= 0; --i) {
      int row = factorRowPerm[i];
      iterstack.push_back(rowroot[row]);

      while (!iterstack.empty()) {
        int rowpos = iterstack.back();
        iterstack.pop_back();
        assert(rowpos != -1);

        if (ARleft[rowpos] != -1) iterstack.push_back(ARleft[rowpos]);
        if (ARright[rowpos] != -1) iterstack.push_back(ARright[rowpos]);

        int col = Acol[rowpos];
        if (colBasisStatus[col] != 0) continue;

        colBasisStatus[col] = -1;
        basisSwaps.emplace_back(i, col);
      }
    }

    int basisSwapPos = 0;

    bool performedBasisSwap;
    do {
      performedBasisSwap = false;
      solution.clear();

      for (int i = numFactorRows - 1; i >= 0; --i) {
        int row = factorRowPerm[i];

        unsigned int solval = 0;

        for (const std::pair<int, unsigned int>& solentry : solution) {
          int pos = findNonzero(row, solentry.first);
          if (pos != -1) solval += Avalue[pos] * solentry.second;
        }

        solval = rhs[row] + k - (solval % k);

        int col = factorColPerm[i];
        int pos = findNonzero(row, col);
        assert(pos != -1);
        unsigned int colValInverse = HighsGFk<k>::inverse(Avalue[pos]);

        assert(solval >= 0);
        assert(colValInverse != 0);

        solval = (solval * colValInverse) % k;

        assert(solval >= 0 && solval < k);

        // only record nonzero solution values
        if (solval != 0) solution.emplace_back(col, solval);
      }

      reportSolution(solution);

      if (basisSwapPos < (int)basisSwaps.size()) {
        int basisIndex = basisSwaps[basisSwapPos].first;
        int enteringCol = basisSwaps[basisSwapPos].second;
        int leavingCol = factorColPerm[basisIndex];
        assert(colBasisStatus[leavingCol] == 1);
        factorColPerm[basisIndex] = enteringCol;
        colBasisStatus[enteringCol] = 1;
        colBasisStatus[leavingCol] = 0;
        performedBasisSwap = true;
        ++basisSwapPos;
      }
    } while (performedBasisSwap);
  }
};

#endif
