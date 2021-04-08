/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <map>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "presolve/HAggregator.h"
#include "presolve/HPreData.h"
#include "presolve/PresolveAnalysis.h"
#include "test/DevKkt.h"

using std::list;
using std::string;

enum class HighsPostsolveStatus {
  NotPresolved = -1,
  ReducedSolutionEmpty,
  ReducedSolutionDimenionsError,
  SolutionRecovered,
  LpOrPresolveObjectMissing,
  BasisError,
  NoPostsolve
};

enum class HighsPresolveStatus {
  NotPresolved = -1,
  NotReduced,
  Infeasible,
  Unbounded,
  Empty,
  Reduced,
  ReducedToEmpty,
  Timeout,
  NullError,
  OptionsError,
};

namespace presolve {

enum class Presolver {
  kMainEmpty,
  kMainRowSingletons,
  kMainForcing,
  kMainColSingletons,
  kMainDoubletonEq,
  kMainDominatedCols,
  kMainSingletonsOnly,
  kMainMipDualFixing,
};

const std::map<Presolver, std::string> kPresolverNames{
    {Presolver::kMainEmpty, "Empty & fixed ()"},
    {Presolver::kMainRowSingletons, "Row singletons ()"},
    {Presolver::kMainForcing, "Forcing rows ()"},
    {Presolver::kMainColSingletons, "Col singletons ()"},
    {Presolver::kMainDoubletonEq, "Doubleton eq ()"},
    {Presolver::kMainDominatedCols, "Dominated Cols()"},
    {Presolver::kMainSingletonsOnly, "Singletons only()"},
    {Presolver::kMainMipDualFixing, "Dual fixing ()"}};

class Presolve : public HPreData {
 public:
  Presolve(HighsTimer& timer_ref) : timer(timer_ref) {}
  virtual ~Presolve() {}

  HighsPresolveStatus presolve();
  HighsPostsolveStatus postsolve(const HighsSolution& reduced_solution,
                                 const HighsBasis& reduced_basis,
                                 HighsSolution& recovered_solution,
                                 HighsBasis& recovered_basis);

  HighsPostsolveStatus primalPostsolve(
      const std::vector<double>& reduced_solution,
      HighsSolution& recovered_solution);

  void setNumericalTolerances();
  void load(const HighsLp& lp, bool mip = false);
  // todo: clear the public from below.
  string modelName;

  // Options
  std::vector<Presolver> order;

  struct AggregatorCall {
    HAggregator::PostsolveStack postsolveStack;
    std::vector<double> colCostAtCall;
    std::vector<double> ARvalueAtCall;
    std::vector<HighsInt> ARindexAtCall;
    std::vector<HighsInt> ARstartAtCall;
  };

  std::vector<AggregatorCall> aggregatorStack;

  HighsInt max_iterations = 0;

  void setTimeLimit(const double limit) {
    assert(limit < inf && limit > 0);
    timer.time_limit = limit;
  }

  HighsInt iPrint = 0;
  HighsLogOptions log_options;
  double objShift;

 private:
  HighsInt iKKTcheck = 0;
  HighsInt presolve(HighsInt print);

  const bool report_postsolve = false;

  void initializeVectors();
  void runAggregator();
  void runPropagator();
  void detectImpliedIntegers();
  void applyMipDualFixing();
  void setProblemStatus(const HighsInt s);
  void reportTimes();

  // new bounds on primal variables for implied free detection
  vector<double> implColLower;
  vector<double> implColUpper;
  vector<HighsInt> implColLowerRowIndex;
  vector<HighsInt> implColUpperRowIndex;

  vector<HighsInt> implRowDualLowerSingColRowIndex;
  vector<HighsInt> implRowDualUpperSingColRowIndex;

  // new bounds on row duals y_i
  vector<double> implRowDualLower;
  vector<double> implRowDualUpper;

  vector<double> implColDualLower;
  vector<double> implColDualUpper;
  vector<double> implRowValueLower;
  vector<double> implRowValueUpper;

  PresolveTimer timer;  // holds enum for main presolve rules

  enum stat {
    Unset = 0,
    Infeasible = 1,
    Unbounded = 2,
    Empty = 3,
    Optimal = 4,
    Reduced = 5,
    Timeout = 6,
  };

 private:
  bool mip;
  bool hasChange = true;
  HighsInt status = 0;  // 0 is unassigned, see enum stat

  list<HighsInt> singRow;  // singleton rows
  list<HighsInt> singCol;  // singleton columns

  // original data
 public:
  vector<double> colCostOriginal;

 private:
  vector<double> rowLowerOriginal;
  vector<double> rowUpperOriginal;
  vector<double> colLowerOriginal;
  vector<double> colUpperOriginal;

  // functions
  void setPrimalValue(const HighsInt j, const double value);
  void checkForChanges(HighsInt iteration);
  void resizeProblem();
  void resizeImpliedBounds();

  // easy transformations
  void removeFixedCol(HighsInt j);
  void removeEmpty();
  void removeFixed();
  void removeEmptyRow(HighsInt i);
  void removeEmptyColumn(HighsInt j);
  void removeRow(HighsInt i);
  void checkBoundsAreConsistent();

  // singleton rows
  void removeRowSingletons();
  HighsInt getSingRowElementIndexInAR(HighsInt i);
  HighsInt getSingColElementIndexInA(HighsInt j);

  // forcing constraints
  void removeForcingConstraints();
  pair<double, double> getImpliedRowBounds(HighsInt row);
  void setVariablesToBoundForForcingRow(const HighsInt row, const bool isLower);
  void dominatedConstraintProcedure(const HighsInt i, const double g,
                                    const double h);

  // doubleton equations
  void removeDoubletonEquations();
  pair<HighsInt, HighsInt> getXYDoubletonEquations(const HighsInt row);
  void processRowDoubletonEquation(const HighsInt row, const HighsInt x,
                                   const HighsInt y, const double akx,
                                   const double aky, const double b);
  pair<double, double> getNewBoundsDoubletonConstraint(HighsInt row,
                                                       HighsInt col, HighsInt j,
                                                       double aik, double aij);
  void UpdateMatrixCoeffDoubletonEquationXzero(
      const HighsInt i, const HighsInt x, const HighsInt y, const double aiy,
      const double akx, const double aky);
  void UpdateMatrixCoeffDoubletonEquationXnonZero(
      const HighsInt i, const HighsInt x, const HighsInt y, const double aiy,
      const double akx, const double aky);

  // column singletons
  void removeColumnSingletons();
  bool removeIfImpliedFree(HighsInt col, HighsInt i, HighsInt k);
  void removeFreeColumnSingleton(const HighsInt col, const HighsInt row,
                                 const HighsInt k);
  void removeZeroCostColumnSingleton(const HighsInt col, const HighsInt row,
                                     const HighsInt k);
  bool removeColumnSingletonInDoubletonInequality(const HighsInt col,
                                                  const HighsInt i,
                                                  const HighsInt k);
  void removeSecondColumnSingletonInDoubletonRow(const HighsInt j,
                                                 const HighsInt i);
  pair<double, double> getBoundsImpliedFree(double lowInit, double uppInit,
                                            const HighsInt col,
                                            const HighsInt i, const HighsInt k);
  void removeImpliedFreeColumn(const HighsInt col, const HighsInt i,
                               const HighsInt k);

  // dominated columns
  void removeDominatedColumns();
  void rowDualBoundsDominatedColumns();
  pair<double, double> getImpliedColumnBounds(HighsInt j);
  void removeIfWeaklyDominated(const HighsInt j, const double d,
                               const double e);

  //    void findDuplicateRows();
  //    void findDuplicateColumns();
  //    void removeDuplicateRows(HighsInt i, HighsInt k, double v);
  //    HighsInt makeCheckForDuplicateRows(HighsInt k, HighsInt i,
  //    vector<double>& coeff, vector<HighsInt>& colIndex, double v, HighsInt
  //    whichIsFirst); void removeDuplicateColumns(HighsInt j,HighsInt k, double
  //    v); bool checkDuplicateRows(HighsInt i, HighsInt k) ;
  //	  bool checkDuplicateColumns(HighsInt i, HighsInt k) ;

  // old or test
  // void updateRemovedColRow(HighsInt dim);
  // void updateRowsByNZ();
  void testAnAR(HighsInt post);

  void countRemovedRows(PresolveRule rule);
  void countRemovedCols(PresolveRule rule);

  double tol = 0.0000001;
  const double default_primal_feasiblility_tolerance = 1e-7;
  const double default_dual_feasiblility_tolerance = 1e-7;
  const double default_small_matrix_value = 1e-9;
  double inconsistent_bounds_tolerance;
  double fixed_column_tolerance;
  double doubleton_equation_bound_tolerance;
  double doubleton_inequality_bound_tolerance;
  double presolve_small_matrix_value;
  double empty_row_bound_tolerance;
  double dominated_column_tolerance;
  double weakly_dominated_column_tolerance;

  // postsolve
  bool noPostSolve = false;

  void addChange(const PresolveRule type, const HighsInt row,
                 const HighsInt col);
  void fillStackRowBounds(const HighsInt col);
  void setKKTcheckerData();

  void getBoundOnLByZj(const HighsInt row, const HighsInt j, double* lo,
                       double* up, const double colLow, const double colUpp);
  double getRowDualPost(const HighsInt row, const HighsInt col);
  double getColumnDualPost(const HighsInt col);
  void roundIntegerBounds(HighsInt col);
  string getDualsForcingRow(const HighsInt row, vector<HighsInt>& fRjs);
  void getDualsSingletonRow(const HighsInt row, const HighsInt col);
  void getDualsDoubletonEquation(const HighsInt row, const HighsInt col);
  void recordCounts(const string fileName);
  void trimA();

  void setBasisElement(const change c);

  // test basis matrix singularity
  //
  // public:
  //	vector<HighsInt> nbffull;
  //	vector<HighsInt> bindfull;
  //	void cmpNBF(HighsInt row, HighsInt col);
  //	void setNBFfullproblem(vector<HighsInt>& nbfFull, vector<HighsInt>&
  // bnFull); 	HighsInt testBasisMatrixSingularity();
  //

  // Dev presolve
  // April 2020
  void reportDevMainLoop();
  void reportDevMidMainLoop();
  PresolveStats stats;
  HighsInt runPresolvers(const std::vector<Presolver>& order);

  void checkKkt(const bool final = false);
  dev_kkt_check::State initState(const bool intermediate = false);

  void caseTwoSingletonsDoubletonInequality(const HighsInt row,
                                            const HighsInt x, const HighsInt y);

  // August 2020
  void removeSingletonsOnly();
  bool isKnapsack(const HighsInt col) const;
  void removeKnapsack(const HighsInt col);
};

}  // namespace presolve

#endif /* PRESOLVE_HPRESOLVE_H_ */
