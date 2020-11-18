/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkControl.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
//#include <cassert>
////#include <iostream>

#include "simplex/HEkk.h"

void HEkk::initialiseControl() {
  // Initialise the iteration count when control started. Need to
  // consider what to do if this isn't zero
  assert(iteration_count_==0);
  simplex_info_.control_iteration_count0 = iteration_count_;
  // Initialise the densities
  simplex_info_.col_aq_density = 0;
  simplex_info_.row_ep_density = 0;
  simplex_info_.row_ap_density = 0;
  simplex_info_.row_DSE_density = 0;
  simplex_info_.col_basic_feasibility_change_density = 0;
  simplex_info_.row_basic_feasibility_change_density = 0;
  simplex_info_.col_BFRT_density = 0;
  simplex_info_.primal_col_density = 0;
  // Set the row_dual_density to 1 since it's assumed all costs are at
  // least perturbed from zero, if not initially nonzero
  simplex_info_.dual_col_density = 1;
  simplex_info_.costly_DSE_frequency = 0;
  simplex_info_.num_costly_DSE_iteration = 0;
}

void HEkk::updateOperationResultDensity(
    const double local_density, double& density) {
  density = (1 - running_average_multiplier) * density +
            running_average_multiplier * local_density;
}

bool HEkk::switchToDevex() {
  bool switch_to_devex = false;
  // Firstly consider switching on the basis of NLA cost
  double costly_DSE_measure;
  double costly_DSE_measure_denominator;
  costly_DSE_measure_denominator =
      max(max(simplex_info_.row_ep_density, simplex_info_.col_aq_density), simplex_info_.row_ap_density);
  if (costly_DSE_measure_denominator > 0) {
    costly_DSE_measure = simplex_info_.row_DSE_density / costly_DSE_measure_denominator;
    costly_DSE_measure = costly_DSE_measure * costly_DSE_measure;
  } else {
    costly_DSE_measure = 0;
  }
  bool costly_DSE_iteration = costly_DSE_measure > costly_DSE_measure_limit &&
    simplex_info_.row_DSE_density > costly_DSE_minimum_density;
  simplex_info_.costly_DSE_frequency = (1 - running_average_multiplier) * simplex_info_.costly_DSE_frequency;
  if (costly_DSE_iteration) {
    simplex_info_.num_costly_DSE_iteration++;
    simplex_info_.costly_DSE_frequency += running_average_multiplier * 1.0;
    // What if non-dual iterations have been performed: need to think about this
    int lcNumIter = iteration_count_ - simplex_info_.control_iteration_count0;
    int numTot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
    // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations
    // have been costly
    switch_to_devex =
        simplex_info_.allow_dual_steepest_edge_to_devex_switch &&
        (simplex_info_.num_costly_DSE_iteration > lcNumIter * costly_DSE_fraction_num_costly_DSE_iteration_before_switch) &&
        (lcNumIter > costly_DSE_fraction_num_total_iteration_before_switch * numTot);

    if (switch_to_devex) {
      HighsLogMessage(
          options_.logfile, HighsMessageType::INFO,
          "Switch from DSE to Devex after %d costly DSE iterations of %d with "
          "densities C_Aq = %11.4g; R_Ep = %11.4g; R_Ap = "
          "%11.4g; DSE = %11.4g",
          simplex_info_.num_costly_DSE_iteration, lcNumIter, simplex_info_.col_aq_density, simplex_info_.row_ep_density,
          simplex_info_.row_ap_density, simplex_info_.row_DSE_density);
    }
  }
  /*
  if (!switch_to_devex) {
    // Secondly consider switching on the basis of weight accuracy
    double dse_weight_error_measure =
        average_log_low_dual_steepest_edge_weight_error +
        average_log_high_dual_steepest_edge_weight_error;
    double dse_weight_error_threshold =
        dual_steepest_edge_weight_log_error_threshold;
    switch_to_devex = allow_dual_steepest_edge_to_devex_switch &&
                      dse_weight_error_measure > dse_weight_error_threshold;
    if (switch_to_devex) {
      HighsLogMessage(logfile, HighsMessageType::INFO,
                      "Switch from DSE to Devex with log error measure of %g > "
                      "%g = threshold",
                      dse_weight_error_measure, dse_weight_error_threshold);
    }
  }
  */
  return switch_to_devex;
}

