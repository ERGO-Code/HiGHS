/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HEKK_H_
#define SIMPLEX_HEKK_H_

#include "HConfig.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"

/**
 * @brief Primal simplex solver for HiGHS
 *
 */

struct HighsEkkStatus {
  // Status of LP solved by the simplex method and its data
  bool valid = false;
  bool is_dualised = false;
  bool is_permuted = false;
  bool scaling_tried = false;
  bool has_basis = false;            // The simplex LP has a valid simplex basis
  bool has_matrix_col_wise = false;  // The HMatrix column-wise matrix is valid
  bool has_matrix_row_wise = false;  // The HMatrix row-wise matrix is valid
  bool has_factor_arrays =
      false;  // Has the arrays for the representation of B^{-1}
  bool has_dual_steepest_edge_weights = false;  // The DSE weights are known
  bool has_nonbasic_dual_values = false;  // The nonbasic dual values are known
  bool has_basic_primal_values = false;   // The basic primal values are known
  bool has_invert =
      false;  // The representation of B^{-1} corresponds to the current basis
  bool has_fresh_invert = false;  // The representation of B^{-1} corresponds to
                                  // the current basis and is fresh
  bool has_fresh_rebuild = false;  // The data are fresh from rebuild
  bool has_dual_objective_value =
      false;  // The dual objective function value is known
  bool has_primal_objective_value =
      false;                    // The dual objective function value is known
  bool has_dual_ray = false;    // A dual unbounded ray is known
  bool has_primal_ray = false;  // A primal unbounded ray is known
  SimplexSolutionStatus solution_status =
      SimplexSolutionStatus::UNSET;  // The solution status is UNSET
};

class HEkk {
 public:
  HEkk(HighsLp& lp, HighsOptions& options) : lp_(lp), options_(options) {}
  /**
   * @brief Solve a model instance
   */
  HighsStatus init();
  HighsStatus solve();

  const SimplexAlgorithm algorithm = SimplexAlgorithm::PRIMAL;

  HighsEkkStatus simplex_lp_status;
  HighsModelStatus model_status;

  HMatrix matrix;
  HFactor factor;

 private:
  // References:
  //
  // LP to be solved, HiGHS options to be used
  HighsLp& lp_;
  HighsOptions& options_;
};

#endif /* SIMPLEX_HEKK_H_ */
