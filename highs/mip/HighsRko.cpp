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
  if (lp->mip_type_ != kMipTypeKnapsack) return false;
  printf("Calling the RKO heuristic for a knapsack problem with %d items\n",
         int(lp->num_col_));
  // solution is initialised to a vector of zeros, so this is
  // necessarily an integer feasible solutiuon of the MIP. Hence
  // return true
  return true;
}
