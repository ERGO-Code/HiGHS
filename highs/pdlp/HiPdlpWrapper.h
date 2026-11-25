/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.h
 * @brief
 */
#ifndef PDLP_HIPDLP_WRAPPER_H_
#define PDLP_HIPDLP_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "lp_data/HighsSolution.h"
HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object);

// A HiPDLP solve that continues an earlier one of the same LP: Highs retries
// the presolved LP this way when the postsolved solution fails the KKT check
// (Highs::retryHiPdlpAfterPostsolve). It stops on the lpKktCheck measures
// (hipdlp::StopCriterion::kLpKktCheck) whatever the first solve used.
struct HiPdlpRestart {
  const HighsSolution* start = nullptr;  // the earlier solution of this LP
  double omega = 0;  // its final primal weight; 0: not known
  // Multipliers of the tolerances of what lpKktCheck measures: primal
  // feasibility, dual feasibility, dual residual and objective (optimality)
  double primal_scale = 1, dual_scale = 1, residual_scale = 1;
  double objective_scale = 1;
};

// The KKT errors of the last lpKktCheck (HighsInfo) relative to their
// tolerances, one per tolerance of HiPDLP (see HiPdlpRestart)
struct HiPdlpKktViolations {
  double primal = 0;  // primal infeasibility and residual
  double dual = 0;
  double dual_residual = 0;
  double objective = 0;
  double max() const {
    return std::max({primal, dual, dual_residual, objective});
  }
};

HiPdlpKktViolations hiPdlpKktViolations(const HighsOptions& options,
                                        const HighsInfo& info);

// What Highs keeps of the HiPDLP solve of the presolved LP, to retry it when
// postsolve takes its KKT errors beyond what lpKktCheck accepts
struct HiPdlpWarmStart {
  bool valid = false;
  HighsSolution solution;  // of the presolved LP
  double omega = 0;        // the primal weight PDLP ended with; 0: unknown
  HiPdlpKktViolations violations;  // of `solution` on the presolved LP
};

// restart: null for an ordinary solve. final_omega (when not null) receives the
// primal weight PDLP ended with.
HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback,
                          const HiPdlpRestart* restart = nullptr,
                          double* final_omega = nullptr);

#endif
