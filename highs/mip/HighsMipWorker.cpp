/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

#include "mip/HighsMipSolverData.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__,
                               HighsLpRelaxation& lprelax_)
    : mipsolver_(mipsolver__),
      mipdata_(*mipsolver_.mipdata_.get()),
      lprelaxation_(lprelax_),
      pseudocost_(mipsolver__),
      cutpool_(mipsolver_.numCol(), mipsolver_.options_mip_->mip_pool_age_limit,
               mipsolver_.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver_.options_mip_->mip_pool_age_limit,
                    mipsolver_.options_mip_->mip_pool_soft_limit) {
  // std::cout << mipdata_.domain.changedcolsflags_.size() << std::endl;
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));

  // Register cutpool and conflict pool in local search domain.
  // Add global cutpool.
  search_ptr_->getLocalDomain().addCutpool(mipsolver_.mipdata_->cutpool);
  search_ptr_->getLocalDomain().addConflictPool(
      mipsolver_.mipdata_->conflictPool);

  // cutpool_.matrix_.AheadNeg_.assign(mipsolver__.numCol(), -1);
  // cutpool_.matrix_.AheadPos_.assign(mipsolver__.numCol(), -1);

  // std::vector<HighsInt> AheadPos_;
  // std::vector<HighsInt> AheadNeg_;

  // add local cutpool
  search_ptr_->getLocalDomain().addCutpool(cutpool_);
  search_ptr_->getLocalDomain().addConflictPool(conflictpool_);

  // printf(
  //     "lprelax_ parameter address in constructor of mipworker %p, %d columns,
  //     " "and "
  //     "%d rows\n",
  //     (void*)&lprelax_, int(lprelax_.getLpSolver().getNumCol()),
  //     int(lprelax_.getLpSolver().getNumRow()));

  // printf(
  //     "lprelaxation_ address in constructor of mipworker %p, %d columns, and
  //     "
  //     "%d rows\n",
  //     (void*)&lprelaxation_, int(lprelaxation_.getLpSolver().getNumCol()),
  //     int(lprelaxation_.getLpSolver().getNumRow()));

  // HighsSearch has its own relaxation initialized no nullptr.

  search_ptr_->setLpRelaxation(&lprelaxation_);

  // printf(
  //     "Search has lp member in constructor of mipworker with address %p, %d "
  //     "columns, and %d rows\n",
  //     (void*)&search_ptr_->lp,
  //     int(search_ptr_->lp->getLpSolver().getNumCol()),
  //     int(search_ptr_->lp->getLpSolver().getNumRow()));
}

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

bool HighsMipWorker::addIncumbent(const std::vector<double>& sol, double solobj,
                                  const int solution_source,
                                  const bool print_display_line) {
  const bool execute_mip_solution_callback =
      !mipsolver_.submip &&
      (mipsolver_.callback_->user_callback
           ? mipsolver_.callback_->active[kCallbackMipSolution]
           : false);
  // Determine whether the potential new incumbent should be
  // transformed
  //
  // Happens if solobj improves on the upper bound or the MIP solution
  // callback is active

  // upper_bound from mipdata is solution_objective_ here

  const bool possibly_store_as_new_incumbent =
      solobj < solution.solution_objective_;

  const bool get_transformed_solution =
      possibly_store_as_new_incumbent || execute_mip_solution_callback;

  // Get the transformed objective and solution if required
  // todo:ig ???
  const double transformed_solobj =
      get_transformed_solution ? transformNewIntegerFeasibleSolution(
                                     sol, possibly_store_as_new_incumbent)
                               : 0;

  if (possibly_store_as_new_incumbent) {
    // #1463 use pre-computed transformed_solobj
    solobj = transformed_solobj;

    if (solobj >= solution.solution_objective_) return false;

    double prev_upper_bound = solution.solution_objective_;

    solution.solution_objective_ = solobj;

    bool bound_change = solution.solution_objective_ != prev_upper_bound;

    if (!mipsolver_.submip && bound_change)
      // todo:ig
      //   updatePrimalDualIntegral(lower_bound, lower_bound, prev_upper_bound,
      //                            upper_bound);

      solution.solution_ = sol;

    // double new_upper_limit = mipdata_.computeNewUpperLimit(solobj, 0.0, 0.0);

    // todo:ig
    // if (!mipsolver_.submip) saveReportMipSolution(new_upper_limit);

    // if (new_upper_limit < upper_limit) {
    //   ++numImprovingSols;
    //   upper_limit = new_upper_limit;
    //   optimality_limit =
    //       computeNewUpperLimit(solobj, mipsolver.options_mip_->mip_abs_gap,
    //                            mipsolver.options_mip_->mip_rel_gap);
    //   nodequeue.setOptimalityLimit(optimality_limit);
    //   debugSolution.newIncumbentFound();
    //   domain.propagate();
    //   if (!domain.infeasible())
    //   redcostfixing.propagateRootRedcost(mipsolver);

    //   // Two calls to printDisplayLine added for completeness,
    //   // ensuring that when the root node has an integer solution, a
    //   // logging line is issued

    //   if (domain.infeasible()) {
    //     pruned_treeweight = 1.0;
    //     nodequeue.clear();
    //     if (print_display_line)
    //       printDisplayLine(solution_source);  // Added for completeness
    //     return true;
    //   }
    //   cliquetable.extractObjCliques(mipsolver);
    //   if (domain.infeasible()) {
    //     pruned_treeweight = 1.0;
    //     nodequeue.clear();
    //     if (print_display_line)
    //       printDisplayLine(solution_source);  // Added for completeness
    //     return true;
    //   }
    //   pruned_treeweight += nodequeue.performBounding(upper_limit);
    //   printDisplayLine(solution_source);
    // }
  } else if (solution.solution_.empty())
    solution.solution_ = sol;

  return true;
}

double HighsMipWorker::transformNewIntegerFeasibleSolution(
    const std::vector<double>& sol,
    const bool possibly_store_as_new_incumbent) {
  //   HighsSolution solution;
  //   solution.col_value = sol;
  //   solution.value_valid = true;
  //   // Perform primal postsolve to get the original column values

  //   mipdata_.postSolveStack.undoPrimal(*mipsolver_.options_mip_, solution);
  //   // Determine the row values, as they aren't computed in primal
  //   // postsolve
  //   HighsInt first_check_row =
  //       -1;  // mipsolver.mipdata_->presolve.debugGetCheckRow();
  //   HighsStatus return_status =
  //       calculateRowValuesQuad(*mipsolver.orig_model_, solution,
  //       first_check_row);
  //   if (kAllowDeveloperAssert) assert(return_status == HighsStatus::kOk);
  //   bool allow_try_again = true;
  // try_again:

  //   // compute the objective value in the original space
  //   double bound_violation_ = 0;
  //   double row_violation_ = 0;
  //   double integrality_violation_ = 0;

  //   // Compute to quad precision the objective function value of the MIP
  //   // being solved - including the offset, and independent of objective
  //   // sense
  //   //
  //   HighsCDouble mipsolver_quad_precision_objective_value =
  //       mipsolver.orig_model_->offset_;
  //   if (kAllowDeveloperAssert)
  //     assert((HighsInt)solution.col_value.size() ==
  //            mipsolver.orig_model_->num_col_);
  //   HighsInt check_col = -1;
  //   HighsInt check_int = -1;
  //   HighsInt check_row = -1;
  //   const bool debug_report = false;
  //   for (HighsInt i = 0; i != mipsolver.orig_model_->num_col_; ++i) {
  //     const double value = solution.col_value[i];
  //     mipsolver_quad_precision_objective_value +=
  //         mipsolver.orig_model_->col_cost_[i] * value;

  //     if (mipsolver.orig_model_->integrality_[i] == HighsVarType::kInteger) {
  //       double integrality_infeasibility = fractionality(value);
  //       if (integrality_infeasibility >
  //           mipsolver.options_mip_->mip_feasibility_tolerance) {
  //         if (debug_report)
  //           printf("Col %d[%s] value %g has integrality infeasibility %g\n",
  //                  int(i), mipsolver.orig_model_->col_names_[i].c_str(),
  //                  value, integrality_infeasibility);
  //         check_int = i;
  //       }
  //       integrality_violation_ =
  //           std::max(integrality_infeasibility, integrality_violation_);
  //     }

  //     const double lower = mipsolver.orig_model_->col_lower_[i];
  //     const double upper = mipsolver.orig_model_->col_upper_[i];
  //     double primal_infeasibility = 0;
  //     if (value < lower - mipsolver.options_mip_->mip_feasibility_tolerance)
  //     {
  //       primal_infeasibility = lower - value;
  //     } else if (value >
  //                upper + mipsolver.options_mip_->mip_feasibility_tolerance) {
  //       primal_infeasibility = value - upper;
  //     } else
  //       continue;
  //     if (primal_infeasibility >
  //         mipsolver.options_mip_->primal_feasibility_tolerance) {
  //       if (debug_report)
  //         printf("Col %d[%s] [%g, %g, %g] has infeasibility %g\n", int(i),
  //                mipsolver.orig_model_->col_names_[i].c_str(), lower, value,
  //                upper, primal_infeasibility);
  //       check_col = i;
  //     }
  //     bound_violation_ = std::max(bound_violation_, primal_infeasibility);
  //   }

  //   for (HighsInt i = 0; i != mipsolver.orig_model_->num_row_; ++i) {
  //     const double value = solution.row_value[i];
  //     const double lower = mipsolver.orig_model_->row_lower_[i];
  //     const double upper = mipsolver.orig_model_->row_upper_[i];
  //     double primal_infeasibility;
  //     if (value < lower - mipsolver.options_mip_->mip_feasibility_tolerance)
  //     {
  //       primal_infeasibility = lower - value;
  //     } else if (value >
  //                upper + mipsolver.options_mip_->mip_feasibility_tolerance) {
  //       primal_infeasibility = value - upper;
  //     } else
  //       continue;
  //     if (primal_infeasibility >
  //         mipsolver.options_mip_->primal_feasibility_tolerance) {
  //       if (debug_report)
  //         printf("Row %d[%s] [%g, %g, %g] has infeasibility %g\n", int(i),
  //                mipsolver.orig_model_->row_names_[i].c_str(), lower, value,
  //                upper, primal_infeasibility);
  //       check_row = i;
  //     }
  //     row_violation_ = std::max(row_violation_, primal_infeasibility);
  //   }

  //   bool feasible =
  //       bound_violation_ <= mipsolver.options_mip_->mip_feasibility_tolerance
  //       && integrality_violation_ <=
  //           mipsolver.options_mip_->mip_feasibility_tolerance &&
  //       row_violation_ <= mipsolver.options_mip_->mip_feasibility_tolerance;

  //   if (!feasible && allow_try_again) {
  //     // printf(
  //     //     "trying to repair sol that is violated by %.12g bounds, %.12g "
  //     //     "integrality, %.12g rows\n",
  //     //     bound_violation_, integrality_violation_, row_violation_);
  //     HighsLp fixedModel = *mipsolver.orig_model_;
  //     fixedModel.integrality_.clear();
  //     for (HighsInt i = 0; i != mipsolver.orig_model_->num_col_; ++i) {
  //       if (mipsolver.orig_model_->integrality_[i] == HighsVarType::kInteger)
  //       {
  //         double solval = std::round(solution.col_value[i]);
  //         fixedModel.col_lower_[i] = std::max(fixedModel.col_lower_[i],
  //         solval); fixedModel.col_upper_[i] =
  //         std::min(fixedModel.col_upper_[i], solval);
  //       }
  //     }

  //     // this->total_repair_lp++;

  //     double time_available = std::max(
  //         mipsolver.options_mip_->time_limit - mipsolver.timer_.read(), 0.1);
  //     Highs tmpSolver;
  //     const bool debug_report = false;
  //     if (debug_report) {
  //       tmpSolver.setOptionValue("log_dev_level", 2);
  //       tmpSolver.setOptionValue("highs_analysis_level", 4);
  //     } else {
  //       tmpSolver.setOptionValue("output_flag", false);
  //     }
  //     // tmpSolver.setOptionValue("simplex_scale_strategy", 0);
  //     // tmpSolver.setOptionValue("presolve", "off");
  //     tmpSolver.setOptionValue("time_limit", time_available);
  //     tmpSolver.setOptionValue("primal_feasibility_tolerance",
  //                              mipsolver.options_mip_->mip_feasibility_tolerance);
  //     tmpSolver.passModel(std::move(fixedModel));

  //     // mipsolver.analysis_.mipTimerStart(kMipClockSimplexNoBasisSolveLp);

  //     tmpSolver.run();
  //     // mipsolver.analysis_.mipTimerStop(kMipClockSimplexNoBasisSolveLp);

  //     // this->total_repair_lp_iterations =
  //     //     tmpSolver.getInfo().simplex_iteration_count;
  //     if (tmpSolver.getInfo().primal_solution_status ==
  //     kSolutionStatusFeasible) {
  //     //   this->total_repair_lp_feasible++;
  //       solution = tmpSolver.getSolution();
  //       allow_try_again = false;
  //       goto try_again;
  //     }
  //   }

  //   // Get a double precision version of the objective function value of
  //   // the MIP being solved
  //   const double mipsolver_objective_value =
  //       double(mipsolver_quad_precision_objective_value);
  //   // Possible MIP solution callback
  //   if (!mipsolver.submip && feasible && mipsolver.callback_->user_callback
  //   &&
  //       mipsolver.callback_->active[kCallbackMipSolution]) {
  //     mipsolver.callback_->clearHighsCallbackDataOut();
  //     mipsolver.callback_->data_out.mip_solution = solution.col_value.data();
  //     // const bool interrupt = interruptFromCallbackWithData(
  //     //     kCallbackMipSolution, mipsolver_objective_value, "Feasible
  //     solution");
  //     // assert(!interrupt);
  //   }

  //   if (possibly_store_as_new_incumbent) {
  //     // Store the solution as incumbent in the original space if there
  //     // is no solution or if it is feasible
  //     if (feasible) {
  //       // if (!allow_try_again)
  //       //   printf("repaired solution with value %g\n",
  //       //   mipsolver_objective_value);
  //       // store
  //       mipsolver.row_violation_ = row_violation_;
  //       mipsolver.bound_violation_ = bound_violation_;
  //       mipsolver.integrality_violation_ = integrality_violation_;
  //       mipsolver.solution_ = std::move(solution.col_value);
  //       mipsolver.solution_objective_ = mipsolver_objective_value;
  //     } else {
  //       bool currentFeasible =
  //           mipsolver.solution_objective_ != kHighsInf &&
  //           mipsolver.bound_violation_ <=
  //               mipsolver.options_mip_->mip_feasibility_tolerance &&
  //           mipsolver.integrality_violation_ <=
  //               mipsolver.options_mip_->mip_feasibility_tolerance &&
  //           mipsolver.row_violation_ <=
  //               mipsolver.options_mip_->mip_feasibility_tolerance;
  //       //    check_col =
  //       37;//mipsolver.mipdata_->presolve.debugGetCheckCol();
  //       //    check_row =
  //       37;//mipsolver.mipdata_->presolve.debugGetCheckRow(); std::string
  //       check_col_data = ""; if (check_col >= 0) {
  //         check_col_data = " (col " + std::to_string(check_col);
  //         if (mipsolver.orig_model_->col_names_.size())
  //           check_col_data +=
  //               "[" + mipsolver.orig_model_->col_names_[check_col] + "]";
  //         check_col_data += ")";
  //       }
  //       std::string check_int_data = "";
  //       if (check_int >= 0) {
  //         check_int_data = " (col " + std::to_string(check_int);
  //         if (mipsolver.orig_model_->col_names_.size())
  //           check_int_data +=
  //               "[" + mipsolver.orig_model_->col_names_[check_int] + "]";
  //         check_int_data += ")";
  //       }
  //       std::string check_row_data = "";
  //       if (check_row >= 0) {
  //         check_row_data = " (row " + std::to_string(check_row);
  //         if (mipsolver.orig_model_->row_names_.size())
  //           check_row_data +=
  //               "[" + mipsolver.orig_model_->row_names_[check_row] + "]";
  //         check_row_data += ")";
  //       }
  //       highsLogUser(mipsolver.options_mip_->log_options,
  //       HighsLogType::kWarning,
  //                    "Solution with objective %g has untransformed
  //                    violations: " "bound = %.4g%s; integrality = %.4g%s; row
  //                    = %.4g%s\n", mipsolver_objective_value,
  //                    bound_violation_, check_col_data.c_str(),
  //                    integrality_violation_, check_int_data.c_str(),
  //                    row_violation_, check_row_data.c_str());

  //       const bool debug_repeat = false;  // true;//
  //       if (debug_repeat) {
  //         HighsSolution check_solution;
  //         check_solution.col_value = sol;
  //         check_solution.value_valid = true;
  //         postSolveStack.undoPrimal(*mipsolver.options_mip_, check_solution,
  //                                   check_col);
  //         fflush(stdout);
  //         if (kAllowDeveloperAssert) assert(111 == 999);
  //       }

  //       if (!currentFeasible) {
  //         // if the current incumbent is non existent or also not feasible we
  //         // still store the new one
  //         mipsolver.row_violation_ = row_violation_;
  //         mipsolver.bound_violation_ = bound_violation_;
  //         mipsolver.integrality_violation_ = integrality_violation_;
  //         mipsolver.solution_ = std::move(solution.col_value);
  //         mipsolver.solution_objective_ = mipsolver_objective_value;
  //       }

  //       // return infinity so that it is not used for bounding
  //       return kHighsInf;
  //     }
  //   }
  //   // return the objective value in the transformed space
  //   if (mipsolver.orig_model_->sense_ == ObjSense::kMaximize)
  //     return -double(mipsolver_quad_precision_objective_value +
  //                    mipsolver.model_->offset_);

  //   return double(mipsolver_quad_precision_objective_value -
  //                 mipsolver.model_->offset_);

  return 0;
}