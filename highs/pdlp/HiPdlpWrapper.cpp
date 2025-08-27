/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.cpp
 * @brief
 * @author Julian Hall
 */
#include "pdlp/HiPdlpWrapper.h"
HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object) {
  return solveLpHiPdlp(solver_object.options_, solver_object.timer_,
                       solver_object.lp_, solver_object.basis_,
                       solver_object.solution_, solver_object.model_status_,
                       solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  return HighsStatus::kError;
}

