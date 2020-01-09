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
  kRootNodeError,
  kRootNodeNotOptimal,
  kUnderDevelopment
};

class HighsMipSolver : Highs {
 public:
  HighsMipSolver(const HighsOptions& options, const HighsLp& lp)
      : options_mip_(options), mip_(lp) {
        root_.parent_id = -1;
      }

  HighsMipStatus runMipSolver();

 private:
  HighsMipStatus solveRootNode();
  HighsMipStatus solveNode(Node& node);
  HighsMipStatus solveTree(Node& root);
  
  Node root_;
  const HighsOptions options_mip_;
  const HighsLp mip_;
};

#endif
