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
  
  // Load root node lp in highs and turn printing off.
  passModel(mip_);
  options_.message_level = 0;
  HighsMipStatus root_solve = solveRootNode();
  if (root_solve != HighsMipStatus::kNodeOptimal) return root_solve;

  // Start tree by making root node.
  // Highs ignores integrality constraints.
  Node root(-1, 0, 0);
  root.col_lower_bound = lp_.colLower_;
  root.col_upper_bound = lp_.colUpper_;
  root.integer_variables = lp_.integrality_;
  root.primal_solution = solution_.col_value;
  root.objective_value = info_.objective_function_value;

  // Add and solve children.
  HighsMipStatus tree_solve_status = solveTree(root);
  if (tree_solve_status != HighsMipStatus::kTreeExhausted) {
    std::cout << "Warning: tree not covered entirely." << std::endl;
    return tree_solve_status;
  }

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
  changeColsBounds(node.branch_col, node.branch_col, &node.col_lower_bound[node.branch_col],
                   &node.col_upper_bound[node.branch_col]);
  HighsStatus lp_solve_status = run();

  switch (lp_solve_status) {
    case HighsStatus::Warning:
      return HighsMipStatus::kNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kNodeError;
    default:
      break;
  }

  switch (model_status_) {
    case HighsModelStatus::OPTIMAL:
      node.primal_solution = solution_.col_value;
      node.objective_value = info_.objective_function_value;
      return HighsMipStatus::kNodeOptimal;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return HighsMipStatus::kNodeInfeasible;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return HighsMipStatus::kNodeUnbounded;
    default:
      break;
  }

  return HighsMipStatus::kNodeNotOptimal;
}

HighsMipStatus HighsMipSolver::solveRootNode() {
  // presolve off for the moment.
  options_.presolve = off_string;
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

HighsMipStatus HighsMipSolver::solveTree(Node& root) {
  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.
  tree_.branch(root);

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree_.empty()) {
    Node& node = tree_.next();
    HighsMipStatus node_solve_status = solveNode(node);
    switch (node_solve_status)
    {
    case HighsMipStatus::kNodeOptimal:
      std::cout << "Node " << node.id
                << " solved to optimality." << std::endl;
      tree_.pop();
      tree_.branch(node);
      break;
    case HighsMipStatus::kNodeUnbounded:
      return HighsMipStatus::kNodeUnbounded;
    case HighsMipStatus::kNodeInfeasible:
      std::cout << "Node " << node.id
                << " infeasible." << std::endl;
      tree_.pop();
      break;
    default:
      std::cout << "Error or warning: Node " << node.id
                << " not solved to optimality." << std::endl;
      break;
    }
  }

  return HighsMipStatus::kTreeExhausted;
}