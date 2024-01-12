/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/CupdlpWrapper.cpp
 * @brief
 * @author Julian Hall
 */
#include "pdlp/CupdlpWrapper.h"

HighsStatus solveLpCupdlp(HighsLpSolverObject& solver_object) {
  return solveLpCupdlp(solver_object.options_, solver_object.timer_, solver_object.lp_, 
		       solver_object.basis_, solver_object.solution_, 
		       solver_object.model_status_, solver_object.highs_info_,
		       solver_object.callback_);
}


HighsStatus solveLpCupdlp(const HighsOptions& options,
			  HighsTimer& timer,
			  const HighsLp& lp, 
			  HighsBasis& highs_basis,
			  HighsSolution& highs_solution,
			  HighsModelStatus& model_status,
			  HighsInfo& highs_info,
			  HighsCallback& callback) {
  // Indicate that there is no valid primal solution, dual solution or basis
  highs_basis.valid = false;
  highs_solution.value_valid = false;
  highs_solution.dual_valid = false;
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);
    assert(111==000);
  HighsStatus return_status = HighsStatus::kError;
  
  return return_status;
}

