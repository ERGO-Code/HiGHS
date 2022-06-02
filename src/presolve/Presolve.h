/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <map>
#include <stack>
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
  kNotPresolved = -1,
  kNoPrimalSolutionError,
  kSolutionRecovered,
  kBasisError
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

  HighsPostsolveStatus postsolve(const HighsSolution& reduced_solution,
                                 const HighsBasis& reduced_basis,
                                 HighsSolution& recovered_solution,
                                 HighsBasis& recovered_basis);

  void load(const HighsLp& lp, bool mip = false);
  // todo: clear the public from below.
  string modelName;

  // Options
  std::vector<Presolver> order;

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

  const bool report_postsolve = false;

  void detectImpliedIntegers();

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

  enum Stat {
    kUnset = 0,
    kInfeasible = 1,
    kUnboundedOrInfeasible = 2,
    kOptimal = 4,
    kReduced = 5,
    kTimeout = 6,
  };

 private:
  bool mip;
  bool hasChange = true;
  HighsInt status = Stat::kUnset;

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
  void resizeProblem();
  void resizeImpliedBounds();

  // easy transformations
  void removeFixedCol(HighsInt j);
  void removeEmptyRow(HighsInt i);

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

  void addChange(const PresolveRule type, const HighsInt row,
                 const HighsInt col);
  void fillStackRowBounds(const HighsInt col);

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
  PresolveStats stats;

  void checkKkt(const bool final = false);
  dev_kkt_check::State initState(const bool intermediate = false);

  // August 2020
  void removeSingletonsOnly();
};

}  // namespace presolve

#endif /* PRESOLVE_HPRESOLVE_H_ */
