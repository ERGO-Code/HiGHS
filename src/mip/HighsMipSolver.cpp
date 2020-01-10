/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolver.h"

#include "lp_data/HighsModelUtils.h"

// Branch-and-bound code below here:
// Solve a mixed integer problem using branch and bound.
HighsMipStatus HighsMipSolver::runMipSolver() {
  std::cout << "Warning: HiGHS MIP solver is under construction at the moment."
            << std::endl
            << "Running HiGHS MIP solver..." << std::endl;

  // Start timer.
  timer_.startRunHighsClock();
  double mip_solve_initial_time = timer_.readRunHighsClock();

  HighsMipStatus root_solve = solveRootNode();
  if (root_solve != HighsMipStatus::kNodeOptimal) return root_solve;

  // Start tree by making root node.
  // Highs ignores integrality constraints.
  Node root(-1, 0, 0);
  root.integer_variables = lp_.integrality_;
  root.col_lower_bound = lp_.colLower_;
  root.col_upper_bound = lp_.colUpper_;
  tree_.pushRootNode(root);

  HighsMipStatus tree_solve = solveTree();
  // todo: handle return value

  // Stop and read the HiGHS clock, then work out time for this call
  timer_.stopRunHighsClock();
  double mip_solve_final_time = timer_.readRunHighsClock();

  if (tree_.getBestSolution().size() > 0) {
    hmos_[0].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
    std::stringstream message;
    message << std::endl;
    message << "Optimal solution found.";
    message << std::endl;
    message << "Run status : "
            << highsModelStatusToString(hmos_[0].unscaled_model_status_)
            << std::endl;
    message << "Objective  : " << std::scientific << tree_.getBestObjective()
            << std::endl;
    message << "Time       : " << std::fixed << std::setprecision(3)
            << mip_solve_final_time - mip_solve_initial_time << std::endl;
    message << std::endl;

    HighsPrintMessage(options_mip_.output, options_mip_.message_level,
                      ML_MINIMAL, message.str().c_str());
    // todo: handle feasible vs optimal case once you have a timeout.
  } else {
    HighsPrintMessage(options_mip_.output, options_mip_.message_level,
                      ML_MINIMAL, "No feasible solution found.");
    // todo: handle infeasible vs timeout case once you have a timeout.
  }

  return HighsMipStatus::kUnderDevelopment;
}

HighsMipStatus HighsMipSolver::solveNode(Node& node) {
  // Apply changes to LP from node. For the moment only column bounds.
  changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0],
                   &node.col_upper_bound[0]);
  HighsStatus lp_solve_status = run();

  switch (lp_solve_status) {
    case HighsStatus::Warning:
      return HighsMipStatus::kNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kNodeError;
    default:
      break;
  }

  if (model_status_ != HighsModelStatus::OPTIMAL)
    return HighsMipStatus::kNodeNotOptimal;

  return HighsMipStatus::kNodeOptimal;
}

HighsMipStatus HighsMipSolver::solveRootNode() {
  // presolve off for the moment.
  options_.presolve = false;
  HighsStatus root_lp_solve_status = run();

  switch (root_lp_solve_status) {
    case HighsStatus::Warning:
      return HighsMipStatus::kRootNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kRootNodeError;
    default:
      break;
  }

  if (model_status_ != HighsModelStatus::OPTIMAL)
    return HighsMipStatus::kRootNodeNotOptimal;

  return HighsMipStatus::kNodeOptimal;
}

HighsMipStatus HighsMipSolver::solveTree() {
  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.
  tree_.branch(tree_.getRootNode());

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree_.empty()) {
    Node& node = tree_.next();
    HighsMipStatus node_solve_status = solveNode(node);
    if (node_solve_status != HighsMipStatus::kNodeOptimal) {
      // todo: handle case.
      std::cout << "Error or warning: Node " << node.id << " not solved to optimality" << std::endl; 
      continue;
    } 
    tree_.pop();
    tree_.branch(node);
  }
}