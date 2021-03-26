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
#define PRESOLVE_HIGHS_PRESOLVE_H_
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
#include "mip/HighsMipSolver.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsMatrixSlice.h"

namespace presolve {

class HighsPostsolveStack;

class HPresolve {
  // pointer to model and options that where presolved
  HighsLp* model;
  const HighsOptions* options;
  HighsMipSolver* mipsolver = nullptr;

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
  std::vector<int> rowsizeInteger;
  std::vector<int> rowsizeImplInt;
  std::vector<int> colsize;

  // vector to store the nonzero positions of a row
  std::vector<int> rowpositions;

  // stack to reuse free slots
  std::vector<int> freeslots;

  // vectors holding implied bounds on primal and dual variables as well as
  // their origins
  std::vector<double> implColLower;
  std::vector<double> implColUpper;
  std::vector<int> colLowerSource;
  std::vector<int> colUpperSource;
  std::vector<double> rowDualLower;
  std::vector<double> rowDualUpper;
  std::vector<double> implRowDualLower;
  std::vector<double> implRowDualUpper;
  std::vector<int> rowDualLowerSource;
  std::vector<int> rowDualUpperSource;

  // implied bounds on values of primal and dual rows computed from the bounds
  // of primal and dual variables
  HighsLinearSumBounds impliedRowBounds;
  HighsLinearSumBounds impliedDualRowBounds;

  std::vector<int> changedRowIndices;
  std::vector<uint8_t> changedRowFlag;
  std::vector<int> changedColIndices;
  std::vector<uint8_t> changedColFlag;

  std::vector<std::pair<int, int>> substitutionOpportunities;

  // set with the sizes and indices of equation rows sorted by the size and a
  // vector to access there iterator positions in the set by index for quick
  // removal
  std::set<std::pair<int, int>> equations;
  std::vector<std::set<std::pair<int, int>>::iterator> eqiters;

  bool shrinkProblemEnabled;
  size_t reductionLimit;

  // vectors storing singleton rows and columns
  std::vector<int> singletonRows;
  std::vector<int> singletonColumns;

  // flags to mark rows/columns as deleted
  std::vector<uint8_t> rowDeleted;
  std::vector<uint8_t> colDeleted;

  int probingContingent;
  int probingNumDelCol;
  int numProbed;

  // counters for number of deleted rows and columns
  int numDeletedRows;
  int numDeletedCols;

  // store old problem sizes to compute percentage redunctions in presolve loop
  int oldNumCol;
  int oldNumRow;

  enum class Result {
    Ok,
    PrimalInfeasible,
    DualInfeasible,
    Stopped,
  };

  // private functions for different shared functionality and matrix
  // modification

  void link(int pos);

  void unlink(int pos);

  void markChangedRow(int row);

  void markChangedCol(int col);

  double getMaxAbsColVal(int col) const;

  double getMaxAbsRowVal(int row) const;

  void updateColImpliedBounds(int row, int col, double val);

  void updateRowDualImpliedBounds(int row, int col, double val);

  bool rowCoefficientsIntegral(int row, double scale) const;

  bool isImpliedFree(int col) const;

  bool isDualImpliedFree(int row) const;

  bool isImpliedIntegral(int col);

  bool isImpliedInteger(int col);

  bool isLowerImplied(int col) const;

  bool isUpperImplied(int col) const;

  int countFillin(int row);

  bool checkFillin(HighsHashTable<int, int>& fillinCache, int row, int col);

#ifndef NDEBUG
  void debugPrintRow(HighsPostsolveStack& postSolveStack, int row);
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

  void markRowDeleted(int row);

  void markColDeleted(int col);

  void fixColToLower(HighsPostsolveStack& postsolveStack, int col);

  void fixColToUpper(HighsPostsolveStack& postsolveStack, int col);

  void fixColToZero(HighsPostsolveStack& postsolveStack, int col);

  void substitute(int row, int col, double rhs);

  void changeColUpper(int col, double newUpper);

  void changeColLower(int col, double newLower);

  void changeRowDualUpper(int row, double newUpper);

  void changeRowDualLower(int row, double newLower);

  void changeImplColUpper(int col, double newUpper, int originRow);

  void changeImplColLower(int col, double newLower, int originRow);

  void changeImplRowDualUpper(int row, double newUpper, int originCol);

  void changeImplRowDualLower(int row, double newLower, int originCol);

  Result applyConflictGraphSubstitutions(HighsPostsolveStack& postSolveStack);

  Result fastPresolveLoop(HighsPostsolveStack& postsolveStack);

  Result presolve(HighsPostsolveStack& postsolveStack);

  Result checkLimits(HighsPostsolveStack& postsolveStack);

  void storeCurrentProblemSize();

  double problemSizeReduction();

 public:
  // for LP presolve
  void setInput(HighsLp& model_, const HighsOptions& options_);

  // for MIP presolve
  void setInput(HighsMipSolver& mipsolver);

  void setReductionLimit(size_t reductionLimit) {
    this->reductionLimit = reductionLimit;
  }

  int numNonzeros() const { return int(Avalue.size() - freeslots.size()); }

  void shrinkProblem(HighsPostsolveStack& postSolveStack);

  void addToMatrix(int row, int col, double val);

  Result runProbing(HighsPostsolveStack& postSolveStack);

  Result doubletonEq(HighsPostsolveStack& postSolveStack, int row);

  Result singletonRow(HighsPostsolveStack& postSolveStack, int row);

  Result emptyCol(HighsPostsolveStack& postSolveStack, int col);

  Result singletonCol(HighsPostsolveStack& postSolveStack, int col);

  Result rowPresolve(HighsPostsolveStack& postSolveStack, int row);

  Result colPresolve(HighsPostsolveStack& postSolveStack, int col);

  Result solveOneRowComponent(HighsPostsolveStack& postsolveStack, int row);

  Result initialRowAndColPresolve(HighsPostsolveStack& postSolveStack);

  HighsModelStatus run(HighsPostsolveStack& postSolveStack);

  void computeIntermediateMatrix(std::vector<int>& flagRow,
                                 std::vector<int>& flagCol,
                                 size_t& numreductions);

  void substitute(int substcol, int staycol, double offset, double scale);

  void removeFixedCol(int col);

  void removeRow(int row);

  Result aggregator(HighsPostsolveStack& postSolveStack);

  Result removeRowSingletons(HighsPostsolveStack& postSolveStack);

  Result presolveColSingletons(HighsPostsolveStack& postSolveStack);

  Result presolveChangedRows(HighsPostsolveStack& postSolveStack);

  Result presolveChangedCols(HighsPostsolveStack& postSolveStack);

  Result removeDoubletonEquations(HighsPostsolveStack& postSolveStack);

  int strengthenInequalities();

  int detectImpliedIntegers();

  Result detectParallelRowsAndCols(HighsPostsolveStack& postsolveStack);

  Result sparsify(HighsPostsolveStack& postsolveStack);

  void setRelaxedImpliedBounds();

  static void debug(const HighsLp& lp, const HighsOptions& options);
};

}  // namespace presolve
#endif