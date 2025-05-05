/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__,
                               HighsLpRelaxation& lprelax_)
    : mipsolver_(mipsolver__),
      mipdata_(*mipsolver_.mipdata_.get()),
      lprelaxation_(lprelax_),
      pseudocost_(mipsolver__),
      cutpool_(mipsolver_.numCol(), mipsolver_.options_mip_->mip_pool_age_limit,
               mipsolver_.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver_.options_mip_->mip_pool_age_limit,
                    mipsolver_.options_mip_->mip_pool_soft_limit),
      upper_bound(kHighsInf) {
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
  const bool execute_mip_solution_callback = false;

  // Determine whether the potential new incumbent should be
  // transformed
  //
  // Happens if solobj improves on the upper bound or the MIP solution
  // callback is active
  const bool possibly_store_as_new_incumbent = solobj < upper_bound;
  const bool get_transformed_solution =
      possibly_store_as_new_incumbent || execute_mip_solution_callback;

  // Get the transformed objective and solution if required
  const double transformed_solobj =
      get_transformed_solution ? transformNewIntegerFeasibleSolution(
                                     sol, possibly_store_as_new_incumbent)
                               : 0;

  std::vector<double>& incumbent = solution_.solution_;

  if (possibly_store_as_new_incumbent) {
    solobj = transformed_solobj;
    if (solobj >= upper_bound) return false;

    double prev_upper_bound = upper_bound;

    upper_bound = solobj;

    bool bound_change = upper_bound != prev_upper_bound;
    // todo:
    // if (!mipsolver_.submip && bound_change)
    // updatePrimalDualIntegral(lower_bound, lower_bound, prev_upper_bound,
    //                          upper_bound);

    incumbent = sol;

    // todo:
    // double new_upper_limit = computeNewUpperLimit(solobj, 0.0, 0.0);

    // if (!is_user_solution && !mipsolver.submip)
    //   saveReportMipSolution(new_upper_limit);
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
  } else if (incumbent.empty())
    incumbent = sol;

  return true;
}

double HighsMipWorker::transformNewIntegerFeasibleSolution(
    const std::vector<double>& sol,
    const bool possibly_store_as_new_incumbent) {
  HighsSolution solution;
  solution.col_value = sol;
  solution.value_valid = true;

  // Perform primal postsolve to get the original column values
  mipsolver_.mipdata_->postSolveStack.undoPrimal(*mipsolver_.options_mip_,
                                                 solution);

  // Determine the row values, as they aren't computed in primal
  // postsolve
  HighsStatus return_status =
      calculateRowValuesQuad(*mipsolver_.orig_model_, solution);
  if (kAllowDeveloperAssert) assert(return_status == HighsStatus::kOk);
  bool allow_try_again = true;
try_again:

  // compute the objective value in the original space
  double bound_violation_ = 0;
  double row_violation_ = 0;
  double integrality_violation_ = 0;

  HighsCDouble mipsolver_quad_objective_value = 0;

  bool feasible = mipsolver_.solutionFeasible(
      mipsolver_.orig_model_, solution.col_value, &solution.row_value,
      bound_violation_, row_violation_, integrality_violation_,
      mipsolver_quad_objective_value);
  double mipsolver_objective_value = double(mipsolver_quad_objective_value);
  if (!feasible && allow_try_again) {
    // printf(
    //     "trying to repair sol that is violated by %.12g bounds, %.12g "
    //     "integrality, %.12g rows\n",
    //     bound_violation_, integrality_violation_, row_violation_);
    HighsLp fixedModel = *mipsolver_.orig_model_;
    fixedModel.integrality_.clear();
    for (HighsInt i = 0; i != mipsolver_.orig_model_->num_col_; ++i) {
      if (mipsolver_.orig_model_->integrality_[i] == HighsVarType::kInteger) {
        double solval = std::round(solution.col_value[i]);
        fixedModel.col_lower_[i] = std::max(fixedModel.col_lower_[i], solval);
        fixedModel.col_upper_[i] = std::min(fixedModel.col_upper_[i], solval);
      }
    }

    // todo:
    // this->total_repair_lp++;

    double time_available = std::max(
        mipsolver_.options_mip_->time_limit - mipsolver_.timer_.read(), 0.1);
    Highs tmpSolver;
    const bool debug_report = false;
    if (debug_report) {
      tmpSolver.setOptionValue("log_dev_level", 2);
      tmpSolver.setOptionValue("highs_analysis_level", 4);
    } else {
      tmpSolver.setOptionValue("output_flag", false);
    }
    // tmpSolver.setOptionValue("simplex_scale_strategy", 0);
    // tmpSolver.setOptionValue("presolve", kHighsOffString);
    tmpSolver.setOptionValue("time_limit", time_available);
    tmpSolver.setOptionValue(
        "primal_feasibility_tolerance",
        mipsolver_.options_mip_->mip_feasibility_tolerance);
    // check if only root presolve is allowed
    if (mipsolver_.options_mip_->mip_root_presolve_only)
      tmpSolver.setOptionValue("presolve", kHighsOffString);
    tmpSolver.passModel(std::move(fixedModel));

    mipsolver_.analysis_.mipTimerStart(kMipClockSimplexNoBasisSolveLp);
    tmpSolver.run();
    mipsolver_.analysis_.mipTimerStop(kMipClockSimplexNoBasisSolveLp);

    // todo:
    // this->total_repair_lp_iterations =
    //     tmpSolver.getInfo().simplex_iteration_count;
    if (tmpSolver.getInfo().primal_solution_status == kSolutionStatusFeasible) {
      // this->total_repair_lp_feasible++;
      solution = tmpSolver.getSolution();
      allow_try_again = false;
      goto try_again;
    }
  }

  // todo:
  // Possible MIP solution callback
  // if (!mipsolver.submip && feasible && mipsolver.callback_->user_callback &&
  //     mipsolver.callback_->active[kCallbackMipSolution]) {
  //   mipsolver.callback_->clearHighsCallbackDataOut();
  //   mipsolver.callback_->data_out.mip_solution = solution.col_value.data();
  //   const bool interrupt = interruptFromCallbackWithData(
  //       kCallbackMipSolution, mipsolver_objective_value, "Feasible
  //       solution");
  //   assert(!interrupt);
  // }

  if (possibly_store_as_new_incumbent) {
    // Store the solution as incumbent in the original space if there
    // is no solution or if it is feasible
    if (feasible) {
      // if (!allow_try_again)
      //   printf("repaired solution with value %g\n",
      //   mipsolver_objective_value);
      // store
      solution_.row_violation_ = row_violation_;
      solution_.bound_violation_ = bound_violation_;

      solution_.integrality_violation_ = integrality_violation_;
      solution_.solution_ = std::move(solution.col_value);
      solution_.solution_objective_ = mipsolver_objective_value;
    } else {
      bool currentFeasible =
          solution_.solution_objective_ != kHighsInf &&
          solution_.bound_violation_ <=
              mipsolver_.options_mip_->mip_feasibility_tolerance &&
          solution_.integrality_violation_ <=
              mipsolver_.options_mip_->mip_feasibility_tolerance &&
          solution_.row_violation_ <=
              mipsolver_.options_mip_->mip_feasibility_tolerance;

      highsLogUser(
          mipsolver_.options_mip_->log_options, HighsLogType::kWarning,
          "WORKER Solution with objective %g has untransformed violations: "
          "bound = %.4g; integrality = %.4g; row = %.4g\n",
          mipsolver_objective_value, bound_violation_, integrality_violation_,
          row_violation_);
      if (!currentFeasible) {
        // if the current incumbent is non existent or also not feasible we
        // still store the new one
        solution_.row_violation_ = row_violation_;
        solution_.bound_violation_ = bound_violation_;
        solution_.integrality_violation_ = integrality_violation_;
        solution_.solution_ = std::move(solution.col_value);
        solution_.solution_objective_ = mipsolver_objective_value;
      }

      // return infinity so that it is not used for bounding
      return kHighsInf;
    }
  }
  // return the objective value in the transformed space
  if (mipsolver_.orig_model_->sense_ == ObjSense::kMaximize)
    return -double(mipsolver_quad_objective_value + mipsolver_.model_->offset_);

  return double(mipsolver_quad_objective_value - mipsolver_.model_->offset_);
}