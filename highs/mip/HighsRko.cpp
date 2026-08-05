/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsRko.h"

// This method should return true if an integer feasible solution
// (returned as solution) has been found
bool rkoHeuristic(const HighsLp* lp, std::vector<double>& solution) {
  // The constraint matrix will be column-wise, with num_col+1 entries
  // in the start vector - simplifying loops over the columns. The
  // p_end_ vector will always be empty
  //
  // For coding methods, if you want to work row-wise, do this
  //
  // HighsSparseMatrix ar_matrix = lp.a_matrix_;
  //
  // ar_matrix.ensureRowwise();
  assert(lp->a_matrix_.isColwise());
  assert(lp->a_matrix_.start_.size() == static_cast<size_t>(lp->num_col_ + 1));
  assert(lp->a_matrix_.p_end_.empty());
  if (lp->mip_type_ == kMipTypeKnapsack) {
    printf("No RKO heuristic for a knapsack problem\n");
    return false;
  } else if (lp->mip_type_ == kMipTypeThlp) {
    printf("Calling the RKO heuristic for a THLP problem with n = %d\n",
           int(lp->thlp_data_.n));
    RKOConfig config;
    const bool dev_run = false;
    if (dev_run) {
      config.num_algorithms = 1;
      config.num_runs = 1;
    } else {
      config.num_algorithms = 8;
      config.num_runs = 4;
      config.enable_ipr = true;
    }     
    RKOOptimizer optimizer(config);
    const bool have_solution = optimizer.solveTHLP(lp->thlp_data_, solution);
    bool have_feasible_solution = false;
    if (have_solution) {
      const bool ok_solution_size = static_cast<HighsInt>(solution.size()) >= lp->num_col_;
      if (!ok_solution_size) {
	printf("RKO heuristic returns a solution with %d components, not %d\n",
	       int( solution.size()), int(lp->num_col_));
	return false;
      }
      const bool legal_solution =
	lp->thlp_data_.legalSolution(solution);
      HighsSparseMatrix ar_matrix = lp->a_matrix_;
      ar_matrix.ensureRowwise();
      const bool feasible_solution =
	lp->thlp_data_.feasibleSolution(lp->num_col_, lp->num_row_, 
					lp->col_lower_, lp->col_upper_, lp->col_names_,
					lp->row_lower_, lp->row_upper_, lp->row_names_,
					ar_matrix.start_, ar_matrix.index_,
					ar_matrix.value_, solution);
      have_feasible_solution = legal_solution && feasible_solution;
    }
      
    return have_feasible_solution;
  }
  return false;
}
