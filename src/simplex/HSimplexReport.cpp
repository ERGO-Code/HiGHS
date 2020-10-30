/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplex.h"
#include "simplex/HSimplexDebug.h"

void reportSimplexPhaseIterations(FILE* logfile, const int iteration_count,
                                  const HighsSimplexInfo& simplex_info,
                                  const SimplexAlgorithm algorithm,
                                  const bool initialise) {
  if (simplex_info.run_quiet) return;
  static int iteration_count0 = 0;
  static int dual_phase1_iteration_count0 = 0;
  static int dual_phase2_iteration_count0 = 0;
  static int primal_phase1_iteration_count0 = 0;
  static int primal_phase2_iteration_count0 = 0;
  static int primal_bound_swap0 = 0;
  if (initialise) {
    iteration_count0 = iteration_count;
    dual_phase1_iteration_count0 = simplex_info.dual_phase1_iteration_count;
    dual_phase2_iteration_count0 = simplex_info.dual_phase2_iteration_count;
    primal_phase1_iteration_count0 = simplex_info.primal_phase1_iteration_count;
    primal_phase2_iteration_count0 = simplex_info.primal_phase2_iteration_count;
    primal_bound_swap0 = simplex_info.primal_bound_swap;
    return;
  }
  const int delta_iteration_count = iteration_count - iteration_count0;
  const int delta_dual_phase1_iteration_count =
      simplex_info.dual_phase1_iteration_count - dual_phase1_iteration_count0;
  const int delta_dual_phase2_iteration_count =
      simplex_info.dual_phase2_iteration_count - dual_phase2_iteration_count0;
  const int delta_primal_phase1_iteration_count =
      simplex_info.primal_phase1_iteration_count -
      primal_phase1_iteration_count0;
  const int delta_primal_phase2_iteration_count =
      simplex_info.primal_phase2_iteration_count -
      primal_phase2_iteration_count0;
  const int delta_primal_bound_swap =
      simplex_info.primal_bound_swap - primal_bound_swap0;

  if (delta_dual_phase1_iteration_count + delta_dual_phase2_iteration_count +
          delta_primal_phase1_iteration_count +
          delta_primal_phase2_iteration_count !=
      delta_iteration_count) {
    printf("Iteration total error %d + %d + %d + %d != %d\n",
           delta_dual_phase1_iteration_count, delta_dual_phase2_iteration_count,
           delta_primal_phase1_iteration_count,
           delta_primal_phase2_iteration_count, delta_iteration_count);
  }
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    HighsLogMessage(
        logfile, HighsMessageType::INFO,
        "Primal simplex iterations [Ph1 %d; Ph2 %d; Swaps %d] Total %d",
        delta_primal_phase1_iteration_count,
        delta_primal_phase2_iteration_count, delta_primal_bound_swap,
        delta_iteration_count);
  } else {
    HighsLogMessage(logfile, HighsMessageType::INFO,
                    "Dual simplex iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d",
                    delta_dual_phase1_iteration_count,
                    delta_dual_phase2_iteration_count,
                    delta_primal_phase2_iteration_count, delta_iteration_count);
  }
}
