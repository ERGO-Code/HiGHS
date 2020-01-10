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

enum class HighsMipStatus {
  kOptimal,
  kTimeout,
  kError,
  kNodeOptimal,
  kNodeNotOptimal,
  kNodeError,
  kRootNodeNotOptimal,
  kRootNodeError,
  kUnderDevelopment,
  kTreeExhausted
};

class HighsMipSolver : Highs {
 public:
  HighsMipSolver(const HighsOptions& options, const HighsLp& lp)
      : options_mip_(options), mip_(lp) {
      }

  HighsMipStatus runMipSolver();

 private:
  HighsMipStatus solveRootNode();
  HighsMipStatus solveNode(Node& node);
  HighsMipStatus solveTree();
  
  Tree tree_;
  const HighsOptions options_mip_;
  const HighsLp mip_;
};

#endif
