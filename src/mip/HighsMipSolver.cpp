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
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_ALWAYS,
                    "Using branch and bound solver\n");

  // Need to start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  double mip_solve_initial_time = timer_.readRunHighsClock();

  // Start tree by making root node.
  std::unique_ptr<Node> root = std::unique_ptr<Node>(new Node(-1, 0, 0));

  root->integer_variables = lp_.integrality_;
  root->col_lower_bound = lp_.colLower_;
  root->col_upper_bound = lp_.colUpper_;

  call_status = solveRootNode(*(root.get()));
  return_status =
      interpretCallStatus(call_status, return_status, "solveRootNode");
  if (return_status == HighsStatus::Error)
    return HighsMipStatus::kRootNodeError;
  if (hmos_[0].scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    HighsPrintMessage(
        options_mip_.output, options_mip_.message_level, ML_ALWAYS,
        "Root note not solved to optimality. Status: %s\n",
        utilHighsModelStatusToString(hmos_[0].scaled_model_status_).c_str());
    call_status =
        highsStatusFromHighsModelStatus(hmos_[0].scaled_model_status_);
    return_status = interpretCallStatus(call_status, return_status);
    if (return_status == HighsStatus::Error)
      return HighsMipStatus::kRootNodeError;
    return HighsMipStatus::kRootNodeNotOptimal;
  }

  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.
  Tree tree(*(root.get()));
  tree.branch(options_mip_.output, options_mip_.message_level, *(root.get()));

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree.empty()) {
    Node& node = tree.next();
    call_status = solveNode(node);
    return_status =
        interpretCallStatus(call_status, return_status, "solveNode");
    if (return_status == HighsStatus::Error)
      return HighsMipStatus::kRootNodeError;
    tree.pop();

    if (hmos_[0].scaled_model_status_ == HighsModelStatus::PRIMAL_INFEASIBLE)
      continue;

    tree.branch(options_mip_.output, options_mip_.message_level, node);
  }

  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double mip_solve_final_time = timer_.readRunHighsClock();

  if (tree.getBestSolution().size() > 0) {
    hmos_[0].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
    std::stringstream message;
    message << std::endl;
    message << "Optimal solution found.";
    message << std::endl;
    message << "Run status : "
            << highsModelStatusToString(hmos_[0].unscaled_model_status_)
            << std::endl;
    message << "Objective  : " << std::scientific << tree.getBestObjective()
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

HighsStatus HighsMipSolver::solveNode(Node& node) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Apply column bounds from node to LP.
  const bool check_call = false;
  const bool call_changeColsBounds = true;
  if (call_changeColsBounds) {
    changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0],
                     &node.col_upper_bound[0]);
  } else {
    // Change the LP directly and invalidate the simplex information
    lp_.colLower_ = node.col_lower_bound;
    lp_.colUpper_ = node.col_upper_bound;
    hmos_[0].simplex_lp_status_.valid = false;
  }

  // Call warm start.
  //  HighsStatus status = run();
  // call works but simply calling run() should be enough and will call hot
  // start in the same way as a user would call it from the outside

  int iteration_count0;
  int iteration_count1;
  int solve0_iteration_count;
  int solve1_iteration_count;
  double solve0_objective_value;
  double solve1_objective_value;
  int solve0_status;
  int solve1_status;

  HighsSolutionParams& scaled_solution_params =
      hmos_[0].scaled_solution_params_;
  iteration_count0 = scaled_solution_params.simplex_iteration_count;

  call_status = run();
  return_status =
      interpretCallStatus(call_status, return_status, "solveLpSimplex");
  if (return_status == HighsStatus::Error) return return_status;

  iteration_count1 = scaled_solution_params.simplex_iteration_count;
  solve0_iteration_count = iteration_count1 - iteration_count0;
  solve0_objective_value = scaled_solution_params.objective_function_value;
  solve0_status = (int)hmos_[0].scaled_model_status_;
  printf("Solve0: Obj = %12g; Iter =%6d; Status =%2d\n", solve0_objective_value,
         solve0_iteration_count, solve0_status);

  if (check_call) {
    // Generate a fresh model object for the LP at this node
    hmos_[0].simplex_lp_status_.has_basis = false;
    hmos_[0].basis_.valid_ = false;
    iteration_count0 = scaled_solution_params.simplex_iteration_count;
    call_status = run();
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpSimplex");
    if (return_status == HighsStatus::Error) return return_status;
    iteration_count1 = scaled_solution_params.simplex_iteration_count;
    solve1_iteration_count = iteration_count1 - iteration_count0;
    solve1_objective_value = scaled_solution_params.objective_function_value;
    solve1_status = (int)hmos_[0].scaled_model_status_;
    printf("Solve1: Obj = %12g; Iter =%6d; Status =%2d\n",
           solve1_objective_value, solve1_iteration_count, solve1_status);
    double rlv_objective_value_difference =
        fabs(solve1_objective_value - solve0_objective_value) /
        max(1.0, fabs(solve1_objective_value));
    if (solve0_status != solve1_status) {
      // Look for unequal status
      printf(
          "!! NodeSolveInequality: Status difference: Status0=%2d; Status1=%2d "
          "!!\n",
          solve0_status, solve1_status);
    } else if (solve0_status != (int)HighsModelStatus::PRIMAL_INFEASIBLE) {
      // Unless infeasible, look for unequal objective
      if (rlv_objective_value_difference > 1e-12)
        printf(
            "!! NodeSolveInequality: Relative objective difference = %12g !!\n",
            rlv_objective_value_difference);
    }
  }

  // Set solution.
  if (hmos_[0].scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    node.primal_solution = hmos_[0].solution_.col_value;
    node.objective_value =
        hmos_[0].scaled_solution_params_.objective_function_value;
  }

  // Solve with a new hmo (replace with code above)
  // passModel(lp_);
  // lp_.colLower_ = node.col_lower_bound;
  // lp_.colUpper_ = node.col_upper_bound;

  // HighsStatus status = solveLpSimplex(hmos_[0]);

  // // Set solution.
  // if (status == HighsStatus::Optimal) {
  //   node.primal_solution = hmos_[0].solution_.col_value;
  //   node.objective_value =
  //   hmos_[0].scaled_solution_params_.objective_function_value;
  // }

  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(hmos_[0].scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}

HighsStatus HighsMipSolver::solveRootNode(Node& root) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // No presolve for the moment.
  // HighsStatus status = run();
  // call works but simply calling run() should be enough.
  call_status = run();
  return_status =
      interpretCallStatus(call_status, return_status, "solveLpSimplex");
  if (return_status == HighsStatus::Error) return return_status;

  if (hmos_[0].scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    root.primal_solution = hmos_[0].solution_.col_value;
    root.objective_value =
        hmos_[0].scaled_solution_params_.objective_function_value;
  }
  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(hmos_[0].scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}