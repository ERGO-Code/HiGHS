/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
//#include <cassert>
//#include <iostream>

#include "io/HighsIO.h"
#include "simplex/HEkk.h"

//#include "lp_data/HConst.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
//#include "simplex/HFactorDebug.h"
//#include "simplex/SimplexTimer.h"
//#include "util/HighsRandom.h"
//#include "util/HighsUtils.h"

using std::cout;
using std::endl;

HighsStatus HEkk::init() { return HighsStatus::OK; }
HighsStatus HEkk::solve() {
  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "HEkk::solve called for LP with %d columns, %d rows and %d entries",
      simplex_lp.numCol_, simplex_lp.numRow_, simplex_lp.Astart_[simplex_lp.numCol_]);

  if (!simplex_lp_status.valid) {
    // Simplex LP is not valid so
    //
    // Set simplex options from HiGHS options.
    setSimplexOptions();
    // Initialise the simplex LP data
    //    initialiseSimplexLpDefinition();
    // Initialise the real and integer random vectors
    //    initialiseSimplexLpRandomVectors();
  }

  

  return HighsStatus::Error;
}

// Private methods

void HEkk::setSimplexOptions() {
  // Copy values of HighsOptions for the simplex solver
  //
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  simplex_info.simplex_strategy = options.simplex_strategy;
  simplex_info.dual_edge_weight_strategy =
      options.simplex_dual_edge_weight_strategy;
  simplex_info.price_strategy = options.simplex_price_strategy;
  simplex_info.dual_simplex_cost_perturbation_multiplier =
      options.dual_simplex_cost_perturbation_multiplier;
  simplex_info.factor_pivot_threshold = options.factor_pivot_threshold;
  simplex_info.update_limit = options.simplex_update_limit;

  // Set values of internal options
  simplex_info.store_squared_primal_infeasibility = true;
  // Option for analysing the LP solution
#ifdef HiGHSDEV
  bool useful_analysis = true;  // false;  //
  bool full_timing = false;
  // Options for reporting timing
  simplex_info.report_simplex_inner_clock = useful_analysis;  // true;
  simplex_info.report_simplex_outer_clock = full_timing;
  simplex_info.report_simplex_phases_clock = full_timing;
  simplex_info.report_HFactor_clock = useful_analysis;  // full_timing;//
  // Options for analysing the LP and simplex iterations
  simplex_info.analyse_lp = useful_analysis;  // false;  //
  simplex_info.analyse_iterations = useful_analysis;
  simplex_info.analyse_invert_form = useful_analysis;
  //  simplex_info.analyse_invert_condition = useful_analysis;
  simplex_info.analyse_invert_time = full_timing;
  simplex_info.analyse_rebuild_time = full_timing;
#endif
}

