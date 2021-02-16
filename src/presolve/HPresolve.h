/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HAggregator.h
 * @brief
 * @author Leona Gottwald
 */
#ifndef PRESOLVE_HIGHS_PRESOLVE_H_
#define PRESOLVE_HAGGREGATOR_H_
#include <cassert>
#include <cmath>
#include <list>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"
#include "util/HighsMatrixSlice.h"

namespace presolve {

class HighsPostsolveStack;

class HPresolve {
  // pointer to model and options that where presolved
  HighsLp* model;
  const HighsOptions* options;

  // triplet storage
  std::vector<double> Avalue;
  std::vector<int> Arow;
  std::vector<int> Acol;

  // linked list links for column based links for each nonzero
  std::vector<int> colhead;
  std::vector<int> Anext;
  std::vector<int> Aprev;

  // splay tree links for row based iteration and nonzero lookup
  std::vector<int> rowroot;
  std::vector<int> ARleft;
  std::vector<int> ARright;

  // length of rows and columns
  std::vector<int> rowsize;
  std::vector<int> colsize;

  // vector to store the nonzero positions of a row
  std::vector<int> rowpositions;

  // vectors for remembering if a row implies the upper or lower bound of a
  // column
  std::vector<int> impliedLbRow;
  std::vector<int> impliedUbRow;

  // for each column the threshold of coefficient values for which a
  // substitution is considered numerically safe
  std::vector<double> col_numerics_threshold;

  // priority queue to reuse free slots from the front
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;

  // vectors holding implied row bounds
  std::vector<double> impliedRowLower;
  std::vector<double> impliedRowUpper;

  // set with the sizes and indices of equation rows sorted by the size and a
  // vector to access there iterator positions in the set by index for quick
  // removal
  std::set<std::pair<int, int>> equations;
  std::vector<std::set<std::pair<int, int>>::iterator> eqiters;

  // vectors storing singleton rows and columns
  std::vector<int> singletonRows;
  std::vector<int> singletonColumns;

  // flags to mark rows/columns as deleted
  std::vector<uint8_t> rowDeleted;
  std::vector<uint8_t> colDeleted;

  // counters for number of deleted rows and columns
  int numDeletedRows;
  int numDeletedCols;

  // private functions for different shared functionality and matrix
  // modification

  void link(int pos);

  void unlink(int pos);

  void dropIfZero(int pos);

  double getImpliedLb(int row, int col);

  double getImpliedUb(int row, int col);

  bool isImpliedFree(int col);

  bool isLowerImplied(int col);

  bool isUpperImplied(int col);

  bool impliedRowBoundsValid(int row) const;

  void invalidateImpliedRowBounds(int row);

  void computeImpliedRowBounds(int row);

  int countFillin(int row);

  bool checkFillin(HighsHashTable<int, int>& fillinCache, int row, int col);

  void substitute(HighsPostsolveStack& postsolveStack, int row, int col);

#ifndef NDEBUG
  void debugPrintRow(int row);

  void debugPrintSubMatrix(int row, int col);
#endif

  int findNonzero(int row, int col);

  void fromCSC(const std::vector<double>& Aval, const std::vector<int>& Aindex,
               const std::vector<int>& Astart);

  void fromCSR(const std::vector<double>& ARval,
               const std::vector<int>& ARindex,
               const std::vector<int>& ARstart);

  void toCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
             std::vector<int>& Astart);

  void toCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
             std::vector<int>& ARstart);

  void storeRow(int row);

  HighsTripletPositionSlice getStoredRow() const;

  HighsTripletListSlice getColumnVector(int col) const;

  HighsTripletTreeSlicePreOrder getRowVector(int row) const;

  HighsTripletTreeSliceInOrder getSortedRowVector(int row) const;

 public:
  HPresolve(std::vector<double>& rowLower, std::vector<double>& rowUpper,
            std::vector<double>& colCost, double& objOffset,
            std::vector<HighsVarType>& integrality,
            std::vector<double>& colLower, std::vector<double>& colUpper);

  void setInput(HighsLp& model_, HighsOptions& options_) {
    model = &model_;
    options = &options_;
  }

  int numNonzeros() const { return int(Avalue.size() - freeslots.size()); }

  void addNonzero(int row, int col, double val);

  void run(HighsPostsolveStack& postSolveStack);

  void substitute(int substcol, int staycol, double offset, double scale);

  void removeFixedCol(int col);

  void removeRow(int row);

  void aggregator(HighsPostsolveStack& postSolveStack);

  void removeRowSingletons(HighsPostsolveStack& postSolveStack);

  void removeForcingAndRedundantRows(HighsPostsolveStack& postSolveStack);

  int strengthenInequalities();

  void detectParallelRowsAndCols(HighsPostsolveStack& postsolveStack);
};

}  // namespace presolve
#endif