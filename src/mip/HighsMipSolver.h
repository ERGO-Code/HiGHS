/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsOptions.h"
#include "mip/SolveMip.h"

class HighsMipSolver : Highs {
 public:
  HighsMipSolver(const HighsOptions& options) : options_mip(options) {}

  HighsStatus runBnb();

 private:
  HighsStatus solveNode(Node& node);
  HighsStatus solveRootNode(Node& root);
  const HighsOptions options_mip;
};

#endif
