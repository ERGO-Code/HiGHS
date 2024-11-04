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
#include <climits>

#include "mip/HighsMipSolverData.h"
#include "mip/feasibilityjump.hh"

void HighsMipSolverData::feasibilityJump() {

  const HighsLp* model = this->mipsolver.model_;
  const bool found_integer_feasible_solution = false;
  std::vector<double> col_value;
  double objective_function_value;
  if (found_integer_feasible_solution) {
    col_value.assign(model->num_col_, 0);
    objective_function_value = 0.0;
    addIncumbent(col_value, objective_function_value,
		 kSolutionSourceFeasibilityJump);
  }
}
