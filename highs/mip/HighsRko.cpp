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
  if (lp->mip_type_ != kMipTypeKnapsack) return false;
  printf("Calling the RKO heuristic for a knapsack problem with %d items\n",
         int(lp->num_col_));
  // solution is initialised to a vector of zeros, so this is
  // necessarily an integer feasible solutiuon of the MIP. Hence
  // return true
  return true;
}
