/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsSolutionDebug.h"
#include "util/HighsUtils.h"

bool equalSolutionParams(const HighsSolutionParams& solution_params0,
                         const HighsSolutionParams& solution_params1) {
  bool equal = true;
  if (!equalSolutionObjectiveParams(solution_params0, solution_params1))
    equal = false;
  if (!equalSolutionStatusParams(solution_params0, solution_params1))
    equal = false;
  if (!equalSolutionInfeasibilityParams(solution_params0, solution_params1))
    equal = false;
  return equal;
}

bool equalSolutionObjectiveParams(const HighsSolutionParams& solution_params0,
                                  const HighsSolutionParams& solution_params1) {
  bool equal = true;
  double delta =
      highsRelativeDifference(solution_params0.objective_function_value,
                                solution_params1.objective_function_value);
  if (solution_params0.objective_function_value !=
      solution_params1.objective_function_value) {
    printf(
        "Solution params: objective_function_value %g != %g Difference = %g\n",
        solution_params0.objective_function_value,
        solution_params1.objective_function_value, delta);
    if (delta > 1e-12) equal = false;
  }
  return equal;
}

bool equalSolutionStatusParams(const HighsSolutionParams& solution_params0,
                               const HighsSolutionParams& solution_params1) {
  bool equal = true;
  if (solution_params0.primal_status != solution_params1.primal_status) {
    printf("Solution params: primal_status %d != %d\n",
           solution_params0.primal_status, solution_params1.primal_status);
    equal = false;
  }
  if (solution_params0.dual_status != solution_params1.dual_status) {
    printf("Solution params: dual_status %d != %d\n",
           solution_params0.dual_status, solution_params1.dual_status);
    equal = false;
  }
  return equal;
}

bool equalSolutionInfeasibilityParams(
    const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  double delta;
  bool equal = true;
  if (solution_params0.num_primal_infeasibilities !=
      solution_params1.num_primal_infeasibilities) {
    printf("Solution params: num_primal_infeasibilities %d != %d\n",
           solution_params0.num_primal_infeasibilities,
           solution_params1.num_primal_infeasibilities);
    equal = false;
  }

  delta =
      highsRelativeDifference(solution_params0.sum_primal_infeasibilities,
                                solution_params1.sum_primal_infeasibilities);
  if (solution_params0.sum_primal_infeasibilities !=
      solution_params1.sum_primal_infeasibilities) {
    printf(
        "Solution params: sum_primal_infeasibilities %g != %g Difference = "
        "%g\n",
        solution_params0.sum_primal_infeasibilities,
        solution_params1.sum_primal_infeasibilities, delta);
    if (delta > 1e-12) equal = false;
  }

  delta = highsRelativeDifference(solution_params0.max_primal_infeasibility,
                                    solution_params1.max_primal_infeasibility);
  if (solution_params0.max_primal_infeasibility !=
      solution_params1.max_primal_infeasibility) {
    printf(
        "Solution params: max_primal_infeasibility %g != %g Difference = %g\n",
        solution_params0.max_primal_infeasibility,
        solution_params1.max_primal_infeasibility, delta);
    if (delta > 1e-12) equal = false;
  }

  if (solution_params0.num_dual_infeasibilities !=
      solution_params1.num_dual_infeasibilities) {
    printf("Solution params: num_dual_infeasibilities %d != %d\n",
           solution_params0.num_dual_infeasibilities,
           solution_params1.num_dual_infeasibilities);
    equal = false;
  }

  delta = highsRelativeDifference(solution_params0.sum_dual_infeasibilities,
                                    solution_params1.sum_dual_infeasibilities);
  if (solution_params0.sum_dual_infeasibilities !=
      solution_params1.sum_dual_infeasibilities) {
    printf(
        "Solution params: sum_dual_infeasibilities %g != %g Difference = %g\n",
        solution_params0.sum_dual_infeasibilities,
        solution_params1.sum_dual_infeasibilities, delta);
    if (delta > 1e-12) equal = false;
  }

  delta = highsRelativeDifference(solution_params0.max_dual_infeasibility,
                                    solution_params1.max_dual_infeasibility);
  if (solution_params0.max_dual_infeasibility !=
      solution_params1.max_dual_infeasibility) {
    printf("Solution params: max_dual_infeasibility %g != %g Difference = %g\n",
           solution_params0.max_dual_infeasibility,
           solution_params1.max_dual_infeasibility, delta);
    if (delta > 1e-12) equal = false;
  }

  return equal;
}

