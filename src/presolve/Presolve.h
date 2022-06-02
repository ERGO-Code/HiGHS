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

class Presolve : public HPreData {
 public:
  Presolve(HighsTimer& timer_ref) : timer(timer_ref) {}
  virtual ~Presolve() {}

  // todo: clear the public from below.

  // Options
  HighsLogOptions log_options;

 private:
  HighsInt iKKTcheck = 0;

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
