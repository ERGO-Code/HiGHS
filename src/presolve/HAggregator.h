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
#ifndef PRESOLVE_HAGGREGATOR_H_
#define PRESOLVE_HAGGREGATOR_H_
#include <cassert>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"

namespace presolve {

class HAggregator {
  // triplet storage
  std::vector<double> Avalue;
  std::vector<int> Arow;
  std::vector<int> Acol;

  // linked list links for column based links for each nonzero
  std::vector<int> colhead;
  std::vector<int> Anext;
  std::vector<int> Aprev;

  // splay tree links for row based iteration and lookup
  std::vector<int> rowroot;
  std::vector<int> ARleft;
  std::vector<int> ARright;

  std::vector<int> rowsize;
  std::vector<int> colsize;

  std::vector<int> iterstack;
  std::vector<int> rowpositions;
  std::unordered_map<int, int> fillinCache;
  std::vector<int> impliedLbRow;
  std::vector<int> impliedUbRow;
  std::vector<double> col_numerics_threshold;
  // priority queue to reuse free slots
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;

  // vectors holding row activities
  std::vector<HighsCDouble> minact;
  std::vector<HighsCDouble> maxact;
  std::vector<int> ninfmin;
  std::vector<int> ninfmax;

  struct ImpliedFreeVarReduction {
    int row;
    int col;
    int rowlen;
    int collen;
    int stackpos;
    double eqrhs;
    double colcost;
    double substcoef;
  };

 public:
  struct PostsolveStack {
    friend class HAggregator;

   private:
    std::vector<std::pair<int, double>> reductionValues;
    std::vector<ImpliedFreeVarReduction> reductionStack;

   public:
    void undo(HighsSolution& solution, HighsBasis& basis) const;

    void undo(std::vector<int>& colFlag, std::vector<int>& rowFlag,
              std::vector<double>& col_value, std::vector<double>& col_dual,
              std::vector<double>& row_dual,
              std::vector<HighsBasisStatus>& col_status,
              std::vector<HighsBasisStatus>& row_status) const;

    void undo(std::vector<double>& colvalue) const;

    void clear() {
      reductionStack.clear();
      reductionValues.clear();
    }

    void unsetFlags(std::vector<int>& rowFlag,
                    std::vector<int>& colFlag) const {
      for (const ImpliedFreeVarReduction& reduction : reductionStack) {
        rowFlag[reduction.row] = 0;
        colFlag[reduction.col] = 0;
      }
    }

    bool empty() const { return reductionStack.empty(); }
  };

 private:
  // set with equation rows and a vector to access there iterator positions in
  // the set by index
  std::set<std::pair<int, int>> equations;
  std::vector<std::set<std::pair<int, int>>::iterator> eqiters;

  // settings used for substitution behavior
  double drop_tolerance;
  double bound_tolerance;
  double markowitz_tol;
  int maxfillin;

  // references to row and column information. Row and objective information is
  // updated in the aggregator
  std::vector<double>& rowLower;
  std::vector<double>& rowUpper;
  std::vector<double>& colCost;
  double& objOffset;
  const std::vector<HighsVarType>& integrality;
  const std::vector<double>& colLower;
  const std::vector<double>& colUpper;

  void link(int pos);

  void unlink(int pos);

  void dropIfZero(int pos);

  void addNonzero(int row, int col, double val);

  double getImpliedLb(int row, int col);

  double getImpliedUb(int row, int col);

  bool isImpliedFree(int col);

  void computeActivities(int row);

  bool checkFillin(int row, int col);

  void substitute(PostsolveStack& postsolveStack, int row, int col);

#ifndef NDEBUG
  void debugPrintRow(int row);

  void debugPrintSubMatrix(int row, int col);
#endif

  template <typename Func>
  void loopRow(int row, Func&& func) {
    int current = rowroot[row];

    while (true) {
      while (current != -1) {
        iterstack.push_back(current);
        current = ARleft[current];
      }

      if (iterstack.empty()) return;

      if (func(iterstack.back())) {
        iterstack.clear();
        return;
      }

      current = ARright[iterstack.back()];
      iterstack.pop_back();
    }
  }

  int countFillin(int row);

  void storeRowPositions(int pos);

  int findNonzero(int row, int col);

 public:
  HAggregator(std::vector<double>& rowLower, std::vector<double>& rowUpper,
              std::vector<double>& colCost, double& objOffset,
              const std::vector<HighsVarType>& integrality,
              const std::vector<double>& colLower,
              const std::vector<double>& colUpper);

  void setMaxFillin(int maxfillin) { this->maxfillin = maxfillin; }

  void setDropTolerance(double drop_tolerance) {
    this->drop_tolerance = drop_tolerance;
  }

  void setBoundTolerance(double bound_tolerance) {
    this->bound_tolerance = bound_tolerance;
  }

  void setMarkowitzTolerance(double markowitz_tol) {
    this->markowitz_tol = markowitz_tol;
  }

  void fromCSC(const std::vector<double>& Aval, const std::vector<int>& Aindex,
               const std::vector<int>& Astart);

  void fromDynamicCSC(const std::vector<double>& Aval,
                      const std::vector<int>& Aindex,
                      const std::vector<int>& Astart,
                      const std::vector<int>& Aend,
                      const std::vector<int>& rowFlag,
                      const std::vector<int>& colFlag);

  void fromCSR(const std::vector<double>& ARval,
               const std::vector<int>& ARindex,
               const std::vector<int>& ARstart);

  void toCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
             std::vector<int>& Astart);

  void toCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
             std::vector<int>& ARstart);

  PostsolveStack run();
};

}  // namespace presolve
#endif