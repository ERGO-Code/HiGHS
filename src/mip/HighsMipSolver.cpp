/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolver.h"

#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "presolve/PresolveComponent.h"
#include "util/HighsCDouble.h"

#ifdef HIGHS_DEBUGSOL
std::vector<double> highsDebugSolution;
#endif

HighsMipSolver::HighsMipSolver(const HighsOptions& options, const HighsLp& lp,
                               bool submip)
    : options_mip_(&options), model_(&lp), submip(submip), rootbasis(nullptr) {}

HighsMipSolver::~HighsMipSolver() = default;

HighsPresolveStatus HighsMipSolver::runPresolve() {
  // todo: commented out parts or change Highs::runPresolve to operate on a
  // parameter LP rather than Highs::lp_. Not sure which approach is preferable.
  HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                    ML_MINIMAL, "\nrunning MIP presolve\n");
  const HighsLp& lp_ = *(model_);

  // Exit if the problem is empty or if presolve is set to off.
  if (options_mip_->presolve == off_string)
    return HighsPresolveStatus::NotPresolved;
  if (lp_.numCol_ == 0 && lp_.numRow_ == 0)
    return HighsPresolveStatus::NullError;

  // Clear info from previous runs if lp_ has been modified.
  if (presolve_.has_run_) presolve_.clear();
  double start_presolve = timer_.readRunHighsClock();

  // Set time limit.
  if (options_mip_->time_limit > 0 &&
      options_mip_->time_limit < HIGHS_CONST_INF) {
    double left = options_mip_->time_limit - start_presolve;
    if (left <= 0) {
      HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                        ML_VERBOSE,
                        "Time limit reached while reading in matrix\n");
      return HighsPresolveStatus::Timeout;
    }

    HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                      ML_VERBOSE,
                      "Time limit set: reading matrix took %.2g, presolve "
                      "time left: %.2g\n",
                      start_presolve, left);
    presolve_.options_.time_limit = left;
  }

  // Presolve.
  presolve_.init(lp_, timer_, true);

  if (options_mip_->time_limit > 0 &&
      options_mip_->time_limit < HIGHS_CONST_INF) {
    double current = timer_.readRunHighsClock();
    double time_init = current - start_presolve;
    double left = presolve_.options_.time_limit - time_init;
    if (left <= 0) {
      HighsPrintMessage(
          options_mip_->output, options_mip_->message_level, ML_VERBOSE,
          "Time limit reached while copying matrix into presolve.\n");
      return HighsPresolveStatus::Timeout;
    }

    HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                      ML_VERBOSE,
                      "Time limit set: copying matrix took %.2g, presolve "
                      "time left: %.2g\n",
                      time_init, left);
    presolve_.options_.time_limit = options_mip_->time_limit;
  }

  presolve_.data_.presolve_[0].message_level = options_mip_->message_level;
  presolve_.data_.presolve_[0].output = options_mip_->output;

  HighsPresolveStatus presolve_return_status = presolve_.run();

  // Handle max case.
  if (presolve_return_status == HighsPresolveStatus::Reduced &&
      lp_.sense_ == ObjSense::MAXIMIZE) {
    presolve_.negateReducedLpCost();
    presolve_.data_.reduced_lp_.sense_ = ObjSense::MAXIMIZE;
  }

  // Update reduction counts.
  switch (presolve_.presolve_status_) {
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_.getReducedProblem();
      presolve_.info_.n_cols_removed = lp_.numCol_ - reduced_lp.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_ - reduced_lp.numRow_;
      presolve_.info_.n_nnz_removed =
          (int)lp_.Avalue_.size() - (int)reduced_lp.Avalue_.size();
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      presolve_.info_.n_cols_removed = lp_.numCol_;
      presolve_.info_.n_rows_removed = lp_.numRow_;
      presolve_.info_.n_nnz_removed = (int)lp_.Avalue_.size();
      break;
    }
    default:
      break;
  }
  return presolve_return_status;
}

HighsPostsolveStatus HighsMipSolver::runPostsolve() {
  assert(presolve_.has_run_);
  if (!mipdata_) return HighsPostsolveStatus::ReducedSolutionEmpty;

  const std::vector<double>& incumbent = mipdata_->getSolution();
  bool solution_ok =
      presolve_.getReducedProblem().numCol_ == (int)incumbent.size();
  if (!solution_ok) return HighsPostsolveStatus::ReducedSolutionDimenionsError;

  if (presolve_.presolve_status_ != HighsPresolveStatus::Reduced &&
      presolve_.presolve_status_ != HighsPresolveStatus::ReducedToEmpty)
    return HighsPostsolveStatus::NoPostsolve;

  // Handle max case.
  // if (lp_.sense_ == ObjSense::MAXIMIZE)
  // presolve_.negateReducedLpColDuals(true);

  // Run postsolve
  HighsPostsolveStatus postsolve_status =
      presolve_.data_.presolve_[0].primalPostsolve(
          incumbent, presolve_.data_.recovered_solution_);

  return postsolve_status;

  // if (lp_.sense_ == ObjSense::MAXIMIZE)
  //   presolve_.negateReducedLpColDuals(false);

  // return HighsPostsolveStatus::SolutionRecovered;
}

void HighsMipSolver::run() {
  modelstatus_ = HighsModelStatus::NOTSET;
  // std::cout << options_mip_->presolve << std::endl;
  timer_.start(timer_.solve_clock);
  if (options_mip_->presolve != "off") {
    HighsPresolveStatus presolve_status = runPresolve();
    switch (presolve_status) {
      case HighsPresolveStatus::Reduced:
        reportPresolveReductions(*options_mip_, *model_,
                                 presolve_.getReducedProblem());
        model_ = &presolve_.getReducedProblem();
        break;
      case HighsPresolveStatus::Unbounded:
        modelstatus_ = HighsModelStatus::PRIMAL_UNBOUNDED;
        HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                          ML_MINIMAL,
                          "Presolve: Model detected to be unbounded\n");
        timer_.stop(timer_.solve_clock);
        return;
      case HighsPresolveStatus::Infeasible:
        modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
        HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                          ML_MINIMAL,
                          "Presolve: Model detected to be infeasible\n");
        timer_.stop(timer_.solve_clock);
        return;
      case HighsPresolveStatus::Timeout:
        modelstatus_ = HighsModelStatus::REACHED_TIME_LIMIT;
        HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                          ML_MINIMAL, "Time limit reached during presolve\n");
        timer_.stop(timer_.solve_clock);
        return;
      case HighsPresolveStatus::ReducedToEmpty:
        modelstatus_ = HighsModelStatus::OPTIMAL;
        reportPresolveReductions(*options_mip_, *model_, true);
        mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
        mipdata_->init();
        runPostsolve();
        timer_.stop(timer_.solve_clock);
        return;
      case HighsPresolveStatus::NotReduced:
        reportPresolveReductions(*options_mip_, *model_, false);
        break;
      default:
        assert(false);
    }
  }

  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  mipdata_->init();
  mipdata_->runSetup();
  if (modelstatus_ == HighsModelStatus::NOTSET) {
    mipdata_->evaluateRootNode();
  }
  if (mipdata_->nodequeue.empty()) {
    mipdata_->printDisplayLine();
    HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                      ML_MINIMAL, "\nSolving stopped with status: %s\n",
                      utilModelStatusToString(modelstatus_).c_str());
    bool haveSolution = mipdata_->upper_bound != HIGHS_CONST_INF;
    if (mipdata_->modelcleanup) {
      model_ = mipdata_->modelcleanup->origmodel;
      if (haveSolution)
        mipdata_->modelcleanup->recoverSolution(mipdata_->incumbent);
    }
    if (haveSolution) {
      if (options_mip_->presolve != "off")
        runPostsolve();
      else if (!mipdata_->getSolution().empty()) {
        presolve_.data_.recovered_solution_.col_value = mipdata_->getSolution();
        calculateRowValues(*model_, presolve_.data_.recovered_solution_);
      }
    }
    timer_.stop(timer_.solve_clock);
    return;
  }

  std::shared_ptr<const HighsBasis> basis;
  HighsSearch search{*this, mipdata_->pseudocost};
  HighsSeparation sepa;

  search.setLpRelaxation(&mipdata_->lp);
  sepa.setLpRelaxation(&mipdata_->lp);

  mipdata_->lower_bound = mipdata_->nodequeue.getBestLowerBound();
  search.installNode(mipdata_->nodequeue.popBestBoundNode());

  mipdata_->printDisplayLine();

  HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                    ML_MINIMAL, "\nstarting tree search\n");

  while (search.hasNode()) {
    // set iteration limit for each lp solve during the dive to 10 times the
    // average nodes

    int iterlimit = 100 * int(mipdata_->lp.getNumLpIterations() /
                              (double)std::max(size_t{1}, mipdata_->num_nodes));
    iterlimit = std::max(1000, iterlimit);

    mipdata_->lp.setIterationLimit(iterlimit);

    // perform the dive and put the open nodes to the queue
    size_t plungestart = mipdata_->num_nodes;
    bool limit_reached = false;
    while (true) {
      if (mipdata_->heuristic_lp_iterations <
          mipdata_->total_lp_iterations * mipdata_->heuristic_effort) {
        search.evaluateNode();
        if (!search.currentNodePruned()) search.heuristicSearchNew();
      }

      search.dive();
      ++mipdata_->num_leaves;

      search.flushStatistics();
      if (mipdata_->checkLimits()) {
        limit_reached = true;
        break;
      }

      if (!search.backtrack()) break;

      if (search.getCurrentEstimate() >= mipdata_->upper_limit) break;

      if (mipdata_->num_nodes - plungestart >=
          std::min(size_t{1000}, mipdata_->num_nodes / 10))
        break;

      if (mipdata_->dispfreq != 0) {
        if (mipdata_->num_leaves - mipdata_->last_displeave >=
            mipdata_->dispfreq)
          mipdata_->printDisplayLine();
      }

      // printf("continue plunging due to good esitmate\n");
    }
    search.openNodesToQueue(mipdata_->nodequeue);
    mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                     mipdata_->nodequeue.getBestLowerBound());

    if (limit_reached) break;

    if (mipdata_->dispfreq != 0) {
      if (mipdata_->num_leaves - mipdata_->last_displeave >= mipdata_->dispfreq)
        mipdata_->printDisplayLine();
    }

    // the search datastructure should have no installed node now
    assert(!search.hasNode());

    // propagate the global domain
    mipdata_->domain.propagate();

#ifdef HIGHS_DEBUGSOL
    assert(!mipdata_->domain.infeasible());
    for (int i = 0; i != numCol(); ++i) {
      assert(highsDebugSolution[i] + mipdata_->epsilon >=
             mipdata_->domain.colLower_[i]);
      assert(highsDebugSolution[i] - mipdata_->epsilon <=
             mipdata_->domain.colUpper_[i]);
    }
#endif

    // if global propagation detected infeasibility, stop here
    if (mipdata_->domain.infeasible()) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      break;
    }

    // if global propagation found bound changes, we update the local domain
    if (!mipdata_->domain.getChangedCols().empty()) {
      HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                        ML_MINIMAL, "added %d global bound changes\n",
                        (int)mipdata_->domain.getChangedCols().size());
      mipdata_->cliquetable.cleanupFixed(mipdata_->domain);
      mipdata_->domain.setDomainChangeStack(std::vector<HighsDomainChange>());
      search.resetLocalDomain();

      mipdata_->domain.clearChangedCols();
    }

    // remove the iteration limit when installing a new node
    mipdata_->lp.setIterationLimit();

    // loop to install the next node for the search
    while (!mipdata_->nodequeue.empty()) {
      // printf("popping node from nodequeue (length = %d)\n",
      // (int)nodequeue.size());
      assert(!search.hasNode());
      search.installNode(mipdata_->nodequeue.popBestNode());

      assert(search.hasNode());

      // set the current basis if available
      if (basis) {
        mipdata_->lp.setStoredBasis(basis);
        mipdata_->lp.recoverBasis();
      }

      // we evaluate the node directly here instead of performing a dive
      // because we first want to check if the node is not fathomed due to
      // new global information before we perform separation rounds for the node
      search.evaluateNode();

      // if the node was pruned we remove it from the search and install the
      // next node from the queue
      if (search.currentNodePruned()) {
        search.backtrack();
        ++mipdata_->num_leaves;
        ++mipdata_->num_nodes;
        search.flushStatistics();

        if (mipdata_->checkLimits()) {
          limit_reached = true;
          break;
        }

        mipdata_->lower_bound = std::min(
            mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound());

        if (mipdata_->dispfreq != 0) {
          if (mipdata_->num_leaves - mipdata_->last_displeave >=
              mipdata_->dispfreq)
            mipdata_->printDisplayLine();
        }
        continue;
      }

      // the node is still not fathomed, so perform separation
      sepa.separate(search.getLocalDomain());

      // after separation we store the new basis and proceed with the outer loop
      // to perform a dive from this node
      if (mipdata_->lp.getStatus() != HighsLpRelaxation::Status::Error &&
          mipdata_->lp.getStatus() != HighsLpRelaxation::Status::NotSet)
        mipdata_->lp.storeBasis();

      basis = mipdata_->lp.getStoredBasis();

      break;
    }

    if (limit_reached) break;
  }

  mipdata_->printDisplayLine();
  bool havesolution = mipdata_->upper_bound != HIGHS_CONST_INF;

  if (modelstatus_ == HighsModelStatus::NOTSET) {
    if (havesolution)
      modelstatus_ = HighsModelStatus::OPTIMAL;
    else
      modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
  }

  HighsPrintMessage(options_mip_->output, options_mip_->message_level,
                    ML_MINIMAL, "\nSolving stopped with status: %s\n",
                    utilModelStatusToString(modelstatus_).c_str());

  if (mipdata_->modelcleanup) {
    model_ = mipdata_->modelcleanup->origmodel;
    if (havesolution)
      mipdata_->modelcleanup->recoverSolution(mipdata_->incumbent);
  }

  if (havesolution) {
    if (options_mip_->presolve != "off")
      runPostsolve();
    else if (!mipdata_->getSolution().empty()) {
      presolve_.data_.recovered_solution_.col_value = mipdata_->getSolution();
      calculateRowValues(*model_, presolve_.data_.recovered_solution_);
    }
  }

  timer_.stop(timer_.solve_clock);

  assert(modelstatus_ != HighsModelStatus::NOTSET);
}
