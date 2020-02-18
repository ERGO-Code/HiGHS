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
  //printf("Writing out the MIP as MPS\n"); writeModel("mip.mps");
  options_.message_level = 0;
  HighsMipStatus root_solve = solveRootNode();
  if (root_solve != HighsMipStatus::kNodeOptimal) return root_solve;

  // Start tree by making root node.
  // Highs ignores integrality constraints.
  Node root(-1, 0.0, 0, 0);
  root.col_lower_bound = lp_.colLower_;
  root.col_upper_bound = lp_.colUpper_;
  root.integer_variables = lp_.integrality_;
  root.primal_solution = solution_.col_value;
  root.objective_value = info_.objective_function_value;

  //  writeSolutionForIntegerVariables(root);
  
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

void HighsMipSolver::writeSolutionForIntegerVariables(Node& node) {
  for (int iCol=0; iCol<lp_.numCol_; iCol++) {
    if (!lp_.integrality_[iCol]) continue;
#ifdef HiGHSDEV    
    printf("%2d [%10.4g, %10.4g, %10.4g]\n",
	   iCol, node.col_lower_bound[iCol], node.primal_solution[iCol], node.col_upper_bound[iCol]);
#endif
  }
}

HighsMipStatus HighsMipSolver::solveNode(Node& node, bool hotstart) {
  HighsStatus lp_solve_status = HighsStatus::Error;
  HighsModelStatus model_status = HighsModelStatus::NOTSET;

  // Force calls within run() to be silent by setting the HiGHS
  // logfile to NULL and the HiGHS message_level to zero.
  bool no_highs_log = true;
  // Use save_message_level and save_logfile to keep the original
  // message_level and logfile in case they are changed
  int save_message_level;
  FILE* save_logfile;

#ifdef HiGHSDEV
  // When full_highs_log is true, run() is verbose - for debugging
  bool full_highs_log = false;
  // Setting check_node_id forces full logging for a particular node
  const int check_node_id = HIGHS_CONST_I_INF;//517;//
#endif
  
#ifdef HiGHSDEV
  //  printf("SolveNode: Id = %d; ParentId = %d; BranchCol = %d\n", node.id, node.parent_id, node.branch_col);
  if (node.id == check_node_id) 
    // Switch on full logging for this node - {} used so VScode can
    // stop on this line
    full_highs_log = true;
  }
#endif
  if (hotstart) {
    // Apply changes to LP from node. For the moment only column bounds.
    // Get the original message_level and logfile in case they are set to something different for run() 
    getHighsOptionValue("message_level", save_message_level);
    save_logfile = options_.logfile;
#ifdef HiGHSDEV
    if (full_highs_log) {
      // Using full logging, so prevent "no logging"
      no_highs_log = false;
      setHighsOptionValue("message_level", 7);
      setHighsLogfile(stdout);
    }
#endif
    if (no_highs_log) {
      // Using no logging, so prevent it
      setHighsOptionValue("message_level", 0);
      setHighsLogfile(NULL);
    }
      
    changeColsBounds(0, mip_.numCol_ - 1, &node.col_lower_bound[0],
                    &node.col_upper_bound[0]);
    lp_solve_status = run();
    model_status = model_status_;
    if (no_highs_log
#ifdef HiGHSDEV
	|| full_highs_log
#endif
	) {
      // Reset the values of message_level and logfile 
      setHighsOptionValue("message_level", save_message_level);
      setHighsLogfile(save_logfile);
    }
#ifdef HiGHSDEV
    const bool check_hotstart = false;
    if (check_hotstart) {
      HighsModelStatus hotstart_model_status = model_status;
      double hotstart_objective;
      getHighsInfoValue("objective_function_value", hotstart_objective);
      
      Highs highs;
      highs.options_.message_level = 0;
      highs.options_.presolve = off_string;
      HighsLp lp_node = mip_;
      lp_node.colLower_ = node.col_lower_bound;
      lp_node.colUpper_ = node.col_upper_bound;
      highs.passModel(lp_node);
      if (no_highs_log) highs.setHighsLogfile(NULL);
      if (full_highs_log) {
	highs.setHighsOptionValue("message_level", 7);
      }
      lp_solve_status = highs.run();
      HighsModelStatus check_model_status = highs.model_status_;
      double check_objective;
      highs.getHighsInfoValue("objective_function_value", check_objective);
      if (check_model_status != hotstart_model_status) {
	// Check whether the model status is the same
	printf("SolveNode ERROR: %s = check_model_status != hotstart_model_status = %s\n",
	       highsModelStatusToString(check_model_status).c_str(),
	       highsModelStatusToString(hotstart_model_status).c_str());
      } else if (check_model_status == HighsModelStatus::OPTIMAL) {
	// Check that the same optimal objective value is found
	double dl_objective = fabs(check_objective-hotstart_objective)/max(1.0, fabs(check_objective));
	double tl_dl_objective = 1e-9;
	if (dl_objective > tl_dl_objective || full_highs_log) {
	  printf("SolveNode");
	  if (dl_objective > tl_dl_objective) {
	    printf(" ERROR");
	  }
	  printf(": Optimal objective difference = %g from (hotstart; check) = (%g; %g)\n", dl_objective, check_objective, hotstart_objective);
	}
      }
    }
#endif
  } else {
    // solve from scratch to test
    Highs highs;
    highs.options_.message_level = 0;
    HighsLp lp_node = mip_;
    lp_node.colLower_ = node.col_lower_bound;
    lp_node.colUpper_ = node.col_upper_bound;
    highs.passModel(lp_node);
    lp_solve_status = highs.run();
    model_status = highs.model_status_;
  }

  switch (lp_solve_status) {
    case HighsStatus::Warning:
      return HighsMipStatus::kNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kNodeError;
    default:
      break;
  }
  
  switch (model_status) {
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
  bool no_highs_log = true;
  int save_message_level = options_.message_level;
  FILE* save_logfile = options_.logfile;
  if (no_highs_log) {
    options_.message_level = 0;
    options_.logfile = NULL;
  }
  HighsStatus root_lp_solve_status = run();
  if (no_highs_log) {
    options_.message_level = save_message_level;
    options_.logfile = save_logfile;
  }

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
    if (node.parent_objective >= tree_.getBestObjective()) {
      // Don't solve if we can't better the best IFS
      tree_.pop();
      continue;
    }
    HighsMipStatus node_solve_status = solveNode(node);
    double best_objective;
    switch (node_solve_status)
    {
    case HighsMipStatus::kNodeOptimal:
#ifdef HiGHSDEV    
      printf("Node %9d (branch on %2d) optimal objective %10.4g: ", node.id, node.branch_col, node.objective_value);
#endif
      /*
      std::cout << "Node " << node.id
                << " solved to optimality." << std::endl;
      */
      tree_.pop();
      // Don't branch if we can't better the best IFS
      best_objective = tree_.getBestObjective();
      if (node.objective_value >= best_objective) {
#ifdef HiGHSDEV    
	printf("Don't branch since no better than best IFS of %10.4g\n", best_objective);
#endif
	break;
      }
      tree_.branch(node);
      //      if (node.branch_col > 47) writeSolutionForIntegerVariables(node);

      break;
    case HighsMipStatus::kNodeUnbounded:
      return HighsMipStatus::kNodeUnbounded;
    case HighsMipStatus::kNodeInfeasible:
#ifdef HiGHSDEV    
      printf("Node %9d (branch on %2d) infeasible\n", node.id, node.branch_col);
      /*
      std::cout << "Node " << node.id
                << " infeasible." << std::endl;
      */
#endif
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
