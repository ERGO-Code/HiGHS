/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
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
#include "mip/HighsMipWorker.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolve.h"
#include "presolve/HighsPostsolveStack.h"
#include "presolve/PresolveComponent.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"

using std::fabs;

HighsMipSolver::HighsMipSolver(HighsCallback& callback,
                               const HighsOptions& options, const HighsLp& lp,
                               const HighsSolution& solution, bool submip,
                               HighsInt submip_level)
    : callback_(&callback),
      options_mip_(&options),
      model_(&lp),
      orig_model_(&lp),
      solution_objective_(kHighsInf),
      submip(submip),
      submip_level(submip_level),
      rootbasis(nullptr),
      pscostinit(nullptr),
      clqtableinit(nullptr),
      implicinit(nullptr) {
  timer_.setPrintfFlag(options_mip_->output_flag, options_mip_->log_to_console);
  assert(!submip || submip_level > 0);
  max_submip_level = 0;
  // Initialise empty terminator
  initialiseTerminator();
  assert(termination_status_ == HighsModelStatus::kNotset);
  if (solution.value_valid) {
#ifndef NDEBUG
    // MIP solver doesn't check row residuals, but they should be OK
    // so validate using assert
    bool valid, integral, feasible;
    assessLpPrimalSolution("For debugging: ", options, lp, solution, valid,
                           integral, feasible);
    assert(valid);
#endif
    // Initial solution can be infeasible, but need to set values for violation
    // and objective
    HighsCDouble quad_solution_objective_;
    solutionFeasible(orig_model_, solution.col_value, &solution.row_value,
                     bound_violation_, row_violation_, integrality_violation_,
                     quad_solution_objective_);
    solution_objective_ = double(quad_solution_objective_);
    solution_ = solution.col_value;
  }
}

HighsMipSolver::~HighsMipSolver() = default;

template <class F>
void HighsMipSolver::runTask(F&& f, highs::parallel::TaskGroup& tg,
                             bool parallel_lock, bool force_serial,
                             const std::vector<HighsInt>& indices) {
  if (indices.empty()) return;
  setParallelLock(parallel_lock);
  const bool spawn_tasks = !force_serial && indices.size() > 1 &&
                           !options_mip_->mip_search_simulate_concurrency;
  for (HighsInt i : indices) {
    if (spawn_tasks) {
      tg.spawn([&f, i] { f(i); });
    } else {
      f(i);
    }
  }
  if (spawn_tasks) {
    tg.taskWait();
  }
  setParallelLock(false);
}

void HighsMipSolver::run() {
  modelstatus_ = HighsModelStatus::kNotset;

  if (submip) {
    analysis_.analyse_mip_time = false;
  } else {
    analysis_.timer_ = &this->timer_;
    analysis_.sub_solver_call_time_ = &this->sub_solver_call_time_;
    analysis_.setup(*orig_model_, *options_mip_);
  }
  timer_.start();

  improving_solution_file_ = nullptr;
  if (!submip && options_mip_->mip_improving_solution_file != "")
    improving_solution_file_ =
        fopen(options_mip_->mip_improving_solution_file.c_str(), "w");

  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  analysis_.mipTimerStart(kMipClockPresolve);
  analysis_.mipTimerStart(kMipClockInit);
  mipdata_->init();
  analysis_.mipTimerStop(kMipClockInit);
#ifdef HIGHS_DEBUGSOL
  mipdata_->debugSolution.activate();
  bool debugSolActive = false;
  std::swap(mipdata_->debugSolution.debugSolActive, debugSolActive);
#endif
  analysis_.mipTimerStart(kMipClockRunPresolve);
  mipdata_->runMipPresolve(options_mip_->presolve_reduction_limit);
  analysis_.mipTimerStop(kMipClockRunPresolve);
  analysis_.mipTimerStop(kMipClockPresolve);

  if (analysis_.analyse_mip_time && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed presolve\n", timer_.read());
  // Identify whether time limit has been reached (in presolve)
  if (modelstatus_ == HighsModelStatus::kNotset &&
      timer_.read() >= options_mip_->time_limit)
    modelstatus_ = HighsModelStatus::kTimeLimit;

  if (modelstatus_ != HighsModelStatus::kNotset) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "Presolve: %s\n",
                 utilModelStatusToString(modelstatus_).c_str());
    if (modelstatus_ == HighsModelStatus::kOptimal) {
      mipdata_->lower_bound = 0;
      mipdata_->upper_bound = 0;
      mipdata_->transformNewIntegerFeasibleSolution(std::vector<double>());
      mipdata_->saveReportMipSolution();
    }
    cleanupSolve();
    return;
  }

  analysis_.mipTimerStart(kMipClockSolve);

  if (analysis_.analyse_mip_time && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - starting  setup\n", timer_.read());
  analysis_.mipTimerStart(kMipClockRunSetup);
#ifdef HIGHS_DEBUGSOL
  mipdata_->debugSolution.debugSolActive = debugSolActive;
#endif
  mipdata_->runSetup();
  analysis_.mipTimerStop(kMipClockRunSetup);
  if (analysis_.analyse_mip_time && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed setup\n", timer_.read());

  if (mipdata_->domain.infeasible()) {
    cleanupSolve();
    return;
  }
  // Initialise master worker.
  mipdata_->workers.emplace_back(*this, &mipdata_->lp, &mipdata_->domain,
                                 &mipdata_->cutpool, &mipdata_->conflictPool,
                                 &mipdata_->pseudocost);

  HighsMipWorker& master_worker = mipdata_->workers[0];

restart:
  if (modelstatus_ == HighsModelStatus::kNotset) {
    // Check limits have not been reached before evaluating root node
    if (mipdata_->checkLimits()) {
      cleanupSolve();
      return;
    }
    // Possibly query existence of an external solution
    if (!submip)
      mipdata_->queryExternalSolution(
          solution_objective_, kExternalMipSolutionQueryOriginAfterSetup);

    // Apply the trivial heuristics
    analysis_.mipTimerStart(kMipClockTrivialHeuristics);
    HighsModelStatus returned_model_status = mipdata_->trivialHeuristics();
    analysis_.mipTimerStop(kMipClockTrivialHeuristics);
    if (modelstatus_ == HighsModelStatus::kNotset &&
        returned_model_status == HighsModelStatus::kInfeasible) {
      // trivialHeuristics can spot trivial infeasibility, so act on it
      modelstatus_ = returned_model_status;
      cleanupSolve();
      return;
    }
    // Apply the feasibility jump heuristic (if enabled)
    if (options_mip_->mip_heuristic_run_feasibility_jump) {
      analysis_.mipTimerStart(kMipClockFeasibilityJump);
      HighsModelStatus returned_model_status = mipdata_->feasibilityJump();
      analysis_.mipTimerStop(kMipClockFeasibilityJump);
      if (modelstatus_ == HighsModelStatus::kNotset &&
          returned_model_status == HighsModelStatus::kInfeasible) {
        // feasibilityJump can spot trivial infeasibility, so act on it
        modelstatus_ = returned_model_status;
        cleanupSolve();
        return;
      }
    }
    // End of pre-root-node heuristics
    if (analysis_.analyse_mip_time && !submip)
      if (analysis_.analyse_mip_time & !submip)
        highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                     "MIP-Timing: %11.2g - starting evaluate root node\n",
                     timer_.read());
    analysis_.mipTimerStart(kMipClockEvaluateRootNode);

    mipdata_->evaluateRootNode(master_worker);

    analysis_.mipTimerStop(kMipClockEvaluateRootNode);
    if (this->terminate()) {
      modelstatus_ = this->terminationStatus();
      cleanupSolve();
      return;
    }
    // Sometimes the analytic centre calculation is not completed when
    // evaluateRootNode returns, so stop its clock if it's running
    if (analysis_.analyse_mip_time &&
        analysis_.mipTimerRunning(kMipClockIpxSolveLp))
      analysis_.mipTimerStop(kMipClockIpxSolveLp);
    if (analysis_.analyse_mip_time && !submip)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "MIP-Timing: %11.2g - completed evaluate root node\n",
                   timer_.read());
    // age 5 times to remove stored but never violated cuts after root
    // separation
    analysis_.mipTimerStart(kMipClockPerformAging0);
    mipdata_->cutpool.performAging();
    mipdata_->cutpool.performAging();
    mipdata_->cutpool.performAging();
    mipdata_->cutpool.performAging();
    mipdata_->cutpool.performAging();
    analysis_.mipTimerStop(kMipClockPerformAging0);
  }
  if (mipdata_->nodequeue.empty() || mipdata_->checkLimits()) {
    cleanupSolve();
    return;
  }

  mipdata_->updateLowerBound(mipdata_->nodequeue.getBestLowerBound());
  mipdata_->printDisplayLine();

  // Calculate maximum number of workers
  const HighsInt mip_search_concurrency = options_mip_->mip_search_concurrency;
  const HighsInt max_num_workers =
      highs::parallel::num_threads() == 1 || mip_search_concurrency <= 0 ||
              submip
          ? 1
          : mip_search_concurrency * highs::parallel::num_threads();
  HighsInt num_workers = 1;
  highs::parallel::TaskGroup tg;

  auto destroyOldWorkers = [&]() {
    if (mipdata_->workers.size() <= 1) return;
    while (mipdata_->domains.size() > 1) {
      mipdata_->domains.pop_back();
    }
    while (mipdata_->cutpools.size() > 1) {
      mipdata_->cutpools.pop_back();
    }
    while (mipdata_->conflictpools.size() > 1) {
      mipdata_->conflictpools.pop_back();
    }
    while (mipdata_->lps.size() > 1) {
      mipdata_->lps.pop_back();
    }
    // Global pseudo-cost not stored in pseudo-costs!
    while (!mipdata_->pseudocosts.empty()) {
      mipdata_->pseudocosts.pop_back();
    }
    while (mipdata_->workers.size() > 1) {
      mipdata_->workers.pop_back();
    }
  };

  auto createNewWorker = [&](HighsInt i) {
    mipdata_->domains.emplace_back(mipdata_->domain);
    mipdata_->lps.emplace_back(mipdata_->lp);
    mipdata_->cutpools.emplace_back(numCol(), options_mip_->mip_pool_age_limit,
                                    options_mip_->mip_pool_soft_limit, i + 1);
    mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                         options_mip_->mip_pool_soft_limit);
    mipdata_->domains.back().addCutpool(mipdata_->cutpools.back());
    assert(mipdata_->domains.back().getDomainChangeStack().empty());
    mipdata_->domains.back().addConflictPool(mipdata_->conflictpools.back());
    mipdata_->pseudocosts.emplace_back(*this);
    mipdata_->workers.emplace_back(
        *this, &mipdata_->lps.back(), &mipdata_->domains.back(),
        &mipdata_->cutpools.back(), &mipdata_->conflictpools.back(),
        &mipdata_->pseudocosts.back());
    mipdata_->lps.back().setMipWorker(mipdata_->workers.back());
    mipdata_->lp.notifyCutPoolsLpCopied(1);
    mipdata_->workers.back().randgen.initialise(options_mip_->random_seed +
                                                mipdata_->workers.size() - 1);
    mipdata_->workers.back().nodequeue.setNumCol(numCol());
    mipdata_->debugSolution.registerDomain(
        mipdata_->workers.back().search_ptr_->getLocalDomain());
  };

  // Use case: Change pointers in master worker to local copies of global info
  auto constructAdditionalWorkerData = [&](HighsMipWorker& worker) {
    assert(mipdata_->cutpools.size() == 1 &&
           mipdata_->conflictpools.size() == 1);
    assert(&worker == &mipdata_->workers[0]);
    mipdata_->cutpools.emplace_back(numCol(), options_mip_->mip_pool_age_limit,
                                    options_mip_->mip_pool_soft_limit, 1);
    worker.cutpool_ = &mipdata_->cutpools.back();
    mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                         options_mip_->mip_pool_soft_limit);
    worker.conflictpool_ = &mipdata_->conflictpools.back();
    mipdata_->domains.emplace_back(mipdata_->domain);
    worker.globaldom_ = &mipdata_->domains.back();
    worker.globaldom_->addCutpool(*worker.cutpool_);
    assert(worker.globaldom_->getDomainChangeStack().empty());
    worker.globaldom_->addConflictPool(*worker.conflictpool_);
    mipdata_->pseudocosts.emplace_back(*this);
    worker.pseudocost_ = &mipdata_->pseudocosts.back();
    worker.lp_->setMipWorker(worker);
    worker.resetSearch();
    worker.resetSepa();
    worker.nodequeue.clear();
    worker.nodequeue.setNumCol(numCol());
  };

  auto syncSolutions = [&]() -> void {
    // Note: Upper bound / limit of workers updated via addIncumbent
    for (HighsMipWorker& worker : mipdata_->workers) {
      for (auto& sol : worker.solutions_) {
        mipdata_->addIncumbent(std::get<0>(sol), std::get<1>(sol),
                               std::get<2>(sol));
      }
      worker.solutions_.clear();
    }
  };

  auto syncPools = [&](std::vector<HighsInt>& indices) -> void {
    if (!mipdata_->hasMultipleWorkers() || mipdata_->parallelLockActive())
      return;
    for (const HighsInt i : indices) {
      mipdata_->workers[i].conflictpool_->syncConflictPool(
          mipdata_->conflictPool);
      mipdata_->workers[i].cutpool_->syncCutPool(*this, mipdata_->cutpool);
    }
    mipdata_->cutpool.performAging();
    mipdata_->conflictPool.performAging();
  };

  auto syncGlobalDomain = [&](std::vector<HighsInt>& indices) -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    for (const HighsInt i : indices) {
      HighsMipWorker& worker = mipdata_->workers[i];
      const auto& domchgstack = worker.getGlobalDomain().getDomainChangeStack();
      for (const HighsDomainChange& domchg : domchgstack) {
        if ((domchg.boundtype == HighsBoundType::kLower &&
             domchg.boundval > mipdata_->domain.col_lower_[domchg.column]) ||
            (domchg.boundtype == HighsBoundType::kUpper &&
             domchg.boundval < mipdata_->domain.col_upper_[domchg.column])) {
          mipdata_->domain.changeBound(domchg,
                                       HighsDomain::Reason::unspecified());
        }
      }
    }
  };

  auto doResetWorkerDomain = [&](HighsInt i) {
    HighsMipWorker& worker = mipdata_->workers[i];
    for (const HighsDomainChange& domchg :
         mipdata_->domain.getDomainChangeStack()) {
      worker.getGlobalDomain().changeBound(domchg,
                                           HighsDomain::Reason::unspecified());
    }
#ifndef NDEBUG
    for (HighsInt col = 0; col < numCol(); ++col) {
      assert(mipdata_->domain.col_lower_[col] ==
             worker.globaldom_->col_lower_[col]);
      assert(mipdata_->domain.col_upper_[col] ==
             worker.globaldom_->col_upper_[col]);
    }
#endif
    worker.getGlobalDomain().setDomainChangeStack(
        std::vector<HighsDomainChange>());
    // Warning: Resetting local domain cannot be done in parallel (changes
    // propagationDomains of main pool)
    worker.search_ptr_->resetLocalDomain();
    worker.getGlobalDomain().clearChangedCols();
  };

  auto resetGlobalDomain = [&](bool force, bool resetWorkers) -> void {
    // if global propagation found bound changes, we update the domain
    if (!mipdata_->domain.getChangedCols().empty() || force) {
      analysis_.mipTimerStart(kMipClockUpdateLocalDomain);
      highsLogDev(options_mip_->log_options, HighsLogType::kInfo,
                  "added %" HIGHSINT_FORMAT " global bound changes\n",
                  (HighsInt)mipdata_->domain.getChangedCols().size());
      mipdata_->cliquetable.cleanupFixed(mipdata_->domain);
      if (mipdata_->hasMultipleWorkers() && resetWorkers) {
        // Sync worker domains here. cleanupFixed might have found extra changes
        std::vector<HighsInt> indices(num_workers);
        std::iota(indices.begin(), indices.end(), 0);
        runTask(doResetWorkerDomain, tg, false, true, indices);
      }
      for (const HighsInt col : mipdata_->domain.getChangedCols())
        mipdata_->implications.cleanupVarbounds(col);

      mipdata_->domain.setDomainChangeStack(std::vector<HighsDomainChange>());
      if (!mipdata_->hasMultipleWorkers())
        master_worker.search_ptr_->resetLocalDomain();
      mipdata_->domain.clearChangedCols();
      mipdata_->removeFixedIndices();
      analysis_.mipTimerStop(kMipClockUpdateLocalDomain);
    }
  };

  auto syncGlobalPseudoCost = [&]() -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    std::vector<HighsInt> nsamplesup = mipdata_->pseudocost.getNSamplesUp();
    std::vector<HighsInt> nsamplesdown = mipdata_->pseudocost.getNSamplesDown();
    std::vector<HighsInt> ninferencesup =
        mipdata_->pseudocost.getNInferencesUp();
    std::vector<HighsInt> ninferencesdown =
        mipdata_->pseudocost.getNInferencesDown();
    std::vector<HighsInt> ncutoffsup = mipdata_->pseudocost.getNCutoffsUp();
    std::vector<HighsInt> ncutoffsdown = mipdata_->pseudocost.getNCutoffsDown();
    for (HighsMipWorker& worker : mipdata_->workers) {
      mipdata_->pseudocost.flushPseudoCost(
          worker.getPseudocost(), nsamplesup, nsamplesdown, ninferencesup,
          ninferencesdown, ncutoffsup, ncutoffsdown);
    }
  };

  auto resetWorkerPseudoCosts = [&](std::vector<HighsInt>& indices) {
    if (!mipdata_->hasMultipleWorkers()) return;
    auto doResetWorkerPseudoCost = [&](HighsInt i) -> void {
      mipdata_->pseudocost.syncPseudoCost(mipdata_->workers[i].getPseudocost());
    };
    runTask(doResetWorkerPseudoCost, tg, false, false, indices);
  };

  destroyOldWorkers();
  master_worker.resetSearch();
  // master_worker.search_ptr_->resetLocalDomain();
  // TODO: This is only done to match seed from v1.12
  master_worker.resetSepa();
  master_worker.nodequeue.clear();
  master_worker.nodequeue.setNumCol(numCol());
  master_worker.upper_bound = mipdata_->upper_bound;
  master_worker.upper_limit = mipdata_->upper_limit;
  master_worker.optimality_limit = mipdata_->optimality_limit;

  HighsSearch& search = *master_worker.search_ptr_;
  mipdata_->debugSolution.registerDomain(search.getLocalDomain());

  analysis_.mipTimerStart(kMipClockSearch);
  search.installNode(mipdata_->nodequeue.popBestBoundNode());
  int64_t numStallNodes = 0;
  int64_t lastLbLeave = 0;
  int64_t numQueueLeaves = 0;
  HighsInt numHugeTreeEstim = 0;
  int64_t numNodesLastCheck = mipdata_->num_nodes;
  int64_t nextCheck = mipdata_->num_nodes;
  double treeweightLastCheck = 0.0;
  double upperLimLastCheck = mipdata_->upper_limit;
  double lowerBoundLastCheck = mipdata_->lower_bound;

  auto nodesInstalled = [&]() -> bool {
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.search_ptr_->hasNode()) return true;
    }
    return false;
  };

  auto infeasibleWorkerGlobalDomain = [&]() -> bool {
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.getGlobalDomain().infeasible()) return true;
    }
    return false;
  };

  auto getSearchIndicesWithNoNodes = [&](std::vector<HighsInt>& indices) {
    indices.clear();
    for (HighsInt i = 0;
         i != num_workers && i != mipdata_->nodequeue.numActiveNodes(); i++) {
      if (!mipdata_->workers[i].search_ptr_->hasNode()) {
        indices.emplace_back(i);
      }
    }
  };

  auto installNodes = [&](std::vector<HighsInt>& indices,
                          bool& limit_reached) -> void {
    for (const HighsInt i : indices) {
      if (numQueueLeaves - lastLbLeave >= 10) {
        mipdata_->workers[i].search_ptr_->installNode(
            mipdata_->nodequeue.popBestBoundNode());
        lastLbLeave = numQueueLeaves;
      } else {
        HighsInt bestBoundNodeStackSize =
            mipdata_->nodequeue.getBestBoundDomchgStackSize();
        double bestBoundNodeLb = mipdata_->nodequeue.getBestLowerBound();
        HighsNodeQueue::OpenNode nextNode(mipdata_->nodequeue.popBestNode());
        if (nextNode.lower_bound == bestBoundNodeLb &&
            (HighsInt)nextNode.domchgstack.size() == bestBoundNodeStackSize)
          lastLbLeave = numQueueLeaves;
        mipdata_->workers[i].search_ptr_->installNode(std::move(nextNode));
      }

      ++numQueueLeaves;

      if (mipdata_->workers[i].search_ptr_->getCurrentEstimate() >=
          mipdata_->upper_limit) {
        ++numStallNodes;
        if (options_mip_->mip_max_stall_nodes != kHighsIInf &&
            numStallNodes >= options_mip_->mip_max_stall_nodes) {
          limit_reached = true;
          modelstatus_ = HighsModelStatus::kSolutionLimit;
          break;
        }
      } else
        numStallNodes = 0;
    }
  };

  auto evaluateNodes = [&](std::vector<HighsInt>& indices) -> void {
    std::vector<HighsSearch::NodeResult> search_results(
        mipdata_->workers.size());
    auto doEvaluateNode = [&](HighsInt i) {
      search_results[i] = mipdata_->workers[i].search_ptr_->evaluateNode();
    };
    analysis_.mipTimerStart(kMipClockEvaluateNode1);
    runTask(doEvaluateNode, tg, true, false, indices);
    analysis_.mipTimerStop(kMipClockEvaluateNode1);
    analysis_.mipTimerStart(kMipClockCurrentNodeToQueue);
    for (HighsInt worker_id : indices) {
      if (search_results[worker_id] == HighsSearch::NodeResult::kSubOptimal) {
        mipdata_->workers[worker_id].search_ptr_->currentNodeToQueue(
            mipdata_->nodequeue);
      }
    }
    analysis_.mipTimerStop(kMipClockCurrentNodeToQueue);
  };

  auto handlePrunedNodes = [&](std::vector<HighsInt>& indices) -> bool {
    std::vector<uint8_t> infeasible(num_workers, 0);
    std::vector<uint8_t> flush(num_workers, 0);
    std::vector<uint8_t> prune(num_workers, 0);
    bool multiple_workers = num_workers > 1;
    auto doHandlePrunedNodes = [&](HighsInt i) {
      if (!mipdata_->workers[i].search_ptr_->currentNodePruned()) return;
      HighsDomain& globaldom = mipdata_->workers[i].getGlobalDomain();
      mipdata_->workers[i].search_ptr_->backtrack();
      flush[i] = 1;

      globaldom.propagate();
      if (!multiple_workers) {
        mipdata_->pruned_treeweight += mipdata_->nodequeue.pruneInfeasibleNodes(
            mipdata_->domain, mipdata_->feastol);
      }

      if (globaldom.infeasible()) {
        infeasible[i] = 1;
        if (!multiple_workers) {
          mipdata_->nodequeue.clear();
          mipdata_->pruned_treeweight = 1.0;
          mipdata_->updateLowerBound(
              std::min(kHighsInf, mipdata_->upper_bound));
        }
        return;
      }

      prune[i] = 1;

      if (multiple_workers || mipdata_->checkLimits()) {
        return;
      }

      mipdata_->updateLowerBound(std::min(
          mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
    };
    analysis_.mipTimerStart(kMipClockNodePrunedLoop);
    runTask(doHandlePrunedNodes, tg, true, false, indices);
    // Flush pruned nodes statistics that haven't yet been flushed
    for (HighsInt i : indices) {
      if (flush[i] == 1) {
        ++mipdata_->num_leaves;
        ++mipdata_->num_nodes;
        mipdata_->workers[i].search_ptr_->flushStatistics();
      }
    }
    // Remove search indices that need a new node
    HighsInt num_search_indices = static_cast<HighsInt>(indices.size());
    for (HighsInt i = num_search_indices - 1; i >= 0; i--) {
      if (prune[indices[i]] == 1) {
        num_search_indices--;
        std::swap(indices[i], indices[num_search_indices]);
      }
    }
    indices.resize(num_search_indices);

    for (uint8_t status : infeasible) {
      if (status == 1) {
        mipdata_->nodequeue.clear();
        mipdata_->pruned_treeweight = 1.0;
        mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
        analysis_.mipTimerStop(kMipClockNodePrunedLoop);
        return true;
      }
    }

    syncSolutions();
    // Handle case where all nodes have been pruned
    if (num_search_indices == 0) {
      mipdata_->updateLowerBound(std::min(
          mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
    }

    analysis_.mipTimerStop(kMipClockNodePrunedLoop);
    if (mipdata_->checkLimits()) {
      return true;
    }
    return false;
  };

  auto separateAndStoreBasis = [&](std::vector<HighsInt>& indices) -> bool {
    // the node is still not fathomed, so perform separation
    auto doSeparate = [&](HighsInt i) {
      mipdata_->workers[i].sepa_ptr_->separate(
          mipdata_->workers[i].search_ptr_->getLocalDomain());
    };
    if (options_mip_->mip_allow_cut_separation_at_nodes) {
      analysis_.mipTimerStart(kMipClockNodeSearchSeparation);
      runTask(doSeparate, tg, true, false, indices);
      // Age cutpools
      if (mipdata_->hasMultipleWorkers()) {
        for (HighsInt i : indices) {
          mipdata_->workers[i].cutpool_->performAging();
        }
      }
      analysis_.mipTimerStop(kMipClockNodeSearchSeparation);
    } else {
      for (HighsCutPool& cutpool : mipdata_->cutpools) {
        cutpool.performAging();
      }
    }

    auto syncSepaStats = [&](HighsMipWorker& worker) {
      mipdata_->cliquetable.getNumNeighbourhoodQueries() +=
          worker.sepa_stats.numNeighbourhoodQueries;
      worker.sepa_stats.numNeighbourhoodQueries = 0;
      mipdata_->sepa_lp_iterations += worker.sepa_stats.sepa_lp_iterations;
      mipdata_->total_lp_iterations += worker.sepa_stats.sepa_lp_iterations;
      worker.sepa_stats.sepa_lp_iterations = 0;
    };

    for (const HighsInt i : indices) {
      HighsMipWorker& worker = mipdata_->workers[i];
      syncSepaStats(worker);
      if (worker.getGlobalDomain().infeasible()) {
        worker.search_ptr_->cutoffNode();
        analysis_.mipTimerStart(kMipClockOpenNodesToQueue1);
        worker.search_ptr_->openNodesToQueue(mipdata_->nodequeue);
        analysis_.mipTimerStop(kMipClockOpenNodesToQueue1);
        mipdata_->nodequeue.clear();
        mipdata_->pruned_treeweight = 1.0;

        analysis_.mipTimerStart(kMipClockStoreBasis);
        mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
        return true;
      }
    }

    auto doStoreBasis = [&](HighsInt i) {
      // after separation we store the new basis and proceed with the outer loop
      // to perform a dive from this node
      if (mipdata_->workers[i].lp_->getStatus() !=
              HighsLpRelaxation::Status::kError &&
          mipdata_->workers[i].lp_->getStatus() !=
              HighsLpRelaxation::Status::kNotSet)
        mipdata_->workers[i].lp_->storeBasis();

      std::shared_ptr<const HighsBasis> basis =
          mipdata_->workers[i].lp_->getStoredBasis();
      if (!basis ||
          !isBasisConsistent(mipdata_->workers[i].lp_->getLp(), *basis)) {
        HighsBasis b = mipdata_->firstrootbasis;
        b.row_status.resize(mipdata_->workers[i].lp_->numRows(),
                            HighsBasisStatus::kBasic);
        basis = std::make_shared<const HighsBasis>(std::move(b));
        mipdata_->workers[i].lp_->setStoredBasis(basis);
      }
    };

    runTask(doStoreBasis, tg, false, false, indices);
    return false;
  };

  auto backtrackPlunge = [&](std::vector<HighsInt>& indices,
                             size_t plungestart) {
    HighsInt numPlungeNodes = mipdata_->num_nodes - plungestart;
    if (numPlungeNodes >=
        std::max(static_cast<double>(indices.size()) / 2, 1.0) * 100)
      return false;

    std::vector<uint8_t> backtracked(num_workers, 0);

    auto doBacktrackPlunge = [&](HighsInt i) {
      backtracked[i] = mipdata_->workers[i].search_ptr_->backtrackPlunge(
          mipdata_->hasMultipleWorkers() ? mipdata_->workers[i].nodequeue
                                         : mipdata_->nodequeue);
    };

    analysis_.mipTimerStart(kMipClockBacktrackPlunge);
    runTask(doBacktrackPlunge, tg, true, false, indices);
    analysis_.mipTimerStop(kMipClockBacktrackPlunge);

    // Remove search indices that were not backtracked
    HighsInt num_search_indices = static_cast<HighsInt>(indices.size());
    for (HighsInt i = num_search_indices - 1; i >= 0; i--) {
      if (backtracked[indices[i]] == 0) {
        num_search_indices--;
        std::swap(indices[i], indices[num_search_indices]);
      }
    }
    indices.resize(num_search_indices);
    if (num_search_indices == 0) return false;
#ifndef NDEBUG
    for (HighsInt i : indices) {
      assert(mipdata_->workers[i].search_ptr_->hasNode());
    }
#endif
    analysis_.mipTimerStart(kMipClockPerformAging2);
    for (HighsInt i : indices) {
      if (mipdata_->workers[i].conflictpool_->getNumConflicts() >
          options_mip_->mip_pool_soft_limit) {
        mipdata_->workers[i].conflictpool_->performAging();
      }
      mipdata_->workers[i].search_ptr_->flushStatistics();
    }
    analysis_.mipTimerStop(kMipClockPerformAging2);
    return true;
  };

  auto runHeuristics = [&](std::vector<HighsInt>& indices) -> void {
    std::vector<uint8_t> suboptimal(num_workers, 0);
    auto doRunHeuristics = [&](HighsInt i) -> void {
      HighsMipWorker& worker = mipdata_->workers[i];
      if (!mipdata_->parallelLockActive())
        analysis_.mipTimerStart(kMipClockDiveEvaluateNode, i);
      const HighsSearch::NodeResult evaluate_node_result =
          worker.search_ptr_->evaluateNode();
      if (!mipdata_->parallelLockActive())
        analysis_.mipTimerStop(kMipClockDiveEvaluateNode, i);

      if (evaluate_node_result == HighsSearch::NodeResult::kSubOptimal) {
        suboptimal[i] = 1;
        return;
      }

      if (!mipdata_->parallelLockActive())
        analysis_.mipTimerStart(kMipClockDivePrimalHeuristics, i);
      if (mipdata_->incumbent.empty()) {
        if (!mipdata_->parallelLockActive())
          analysis_.mipTimerStart(kMipClockDiveRandomizedRounding, i);
        mipdata_->heuristics.randomizedRounding(
            worker, worker.lp_->getLpSolver().getSolution().col_value);
        if (!mipdata_->parallelLockActive())
          analysis_.mipTimerStop(kMipClockDiveRandomizedRounding, i);
      }
      if (mipdata_->incumbent.empty()) {
        if (options_mip_->mip_heuristic_run_rens) {
          if (!mipdata_->parallelLockActive())
            analysis_.mipTimerStart(kMipClockDiveRens, i);
          mipdata_->heuristics.RENS(
              worker, worker.lp_->getLpSolver().getSolution().col_value);
          if (!mipdata_->parallelLockActive())
            analysis_.mipTimerStop(kMipClockDiveRens, i);
        }
      } else {
        if (options_mip_->mip_heuristic_run_rins) {
          if (!mipdata_->parallelLockActive())
            analysis_.mipTimerStart(kMipClockDiveRins, i);
          mipdata_->heuristics.RINS(
              worker, worker.lp_->getLpSolver().getSolution().col_value);
          if (!mipdata_->parallelLockActive())
            analysis_.mipTimerStop(kMipClockDiveRins, i);
        }
      }

      if (!mipdata_->parallelLockActive())
        analysis_.mipTimerStop(kMipClockDivePrimalHeuristics, i);
    };
    runTask(doRunHeuristics, tg, true, false, indices);
    for (const HighsInt i : indices) {
      if (suboptimal[i] == 0) {
        if (mipdata_->workers[i].search_ptr_->currentNodePruned()) {
          ++mipdata_->num_leaves;
          mipdata_->workers[i].search_ptr_->flushStatistics();
        }
        mipdata_->heuristics.flushStatistics(*this, mipdata_->workers[i]);
      }
    }
    // Remove search indices that have suboptimal status
    HighsInt num_search_indices = static_cast<HighsInt>(indices.size());
    for (HighsInt i = num_search_indices - 1; i >= 0; i--) {
      if (suboptimal[indices[i]] == 1) {
        num_search_indices--;
        std::swap(indices[i], indices[num_search_indices]);
      }
    }
    indices.resize(num_search_indices);
  };

  auto diveSearches = [&](std::vector<HighsInt>& indices) {
    analysis_.mipTimerStart(kMipClockTheDive);
    std::vector<HighsSearch::NodeResult> dive_results(
        mipdata_->workers.size(), HighsSearch::NodeResult::kBranched);

    // Create vector of non pruned indices
    std::vector<HighsInt> non_pruned_indices;
    for (HighsInt i : indices) {
      if (!mipdata_->workers[i].search_ptr_->currentNodePruned()) {
        non_pruned_indices.push_back(i);
      }
    }
    if (non_pruned_indices.empty()) return;

    auto doDiveSearch = [&](HighsInt i) {
      HighsMipWorker& worker = mipdata_->workers[i];
      if (!worker.search_ptr_->hasNode() ||
          worker.search_ptr_->currentNodePruned())
        return;
      dive_results[i] = worker.search_ptr_->dive();
    };
    runTask(doDiveSearch, tg, true, false, non_pruned_indices);
    analysis_.mipTimerStop(kMipClockTheDive);

    for (const HighsInt i : non_pruned_indices) {
      if (dive_results[i] != HighsSearch::NodeResult::kSubOptimal) {
        ++mipdata_->num_leaves;
        mipdata_->workers[i].search_ptr_->flushStatistics();
      }
    }

    // Remove search indices that have suboptimal status
    HighsInt num_search_indices = static_cast<HighsInt>(indices.size());
    for (HighsInt i = num_search_indices - 1; i >= 0; i--) {
      if (dive_results[indices[i]] == HighsSearch::NodeResult::kSubOptimal) {
        num_search_indices--;
        std::swap(indices[i], indices[num_search_indices]);
      }
    }
    indices.resize(num_search_indices);
  };

  // Search indices tracks which MIP workers were assigned nodes
  // Reduced search indices tracks which workers search haven't yet been pruned
  std::vector<HighsInt> search_indices(1, 0);
  std::vector<HighsInt> reduced_search_indices(1, 0);
  while (nodesInstalled()) {
    // Possibly query existence of an external solution
    if (!submip)
      mipdata_->queryExternalSolution(
          solution_objective_, kExternalMipSolutionQueryOriginBeforeDive);

    analysis_.mipTimerStart(kMipClockPerformAging1);
    // TODO: Is there a need to age local pools? They're essentially deleted.
    for (HighsConflictPool& conflict_pool : mipdata_->conflictpools) {
      conflict_pool.performAging();
    }
    analysis_.mipTimerStop(kMipClockPerformAging1);
    // set iteration limit for each lp solve during the dive to 10 times the
    // average nodes

    HighsInt iterlimit = 10 * std::max(mipdata_->lp.getAvgSolveIters(),
                                       mipdata_->avgrootlpiters);
    iterlimit = std::max({HighsInt{10000}, iterlimit,
                          HighsInt((3 * mipdata_->firstrootlpiters) / 2)});

    for (HighsLpRelaxation& lp : mipdata_->lps) {
      lp.setIterationLimit(iterlimit);
    }

    // perform the dive and put the open nodes to the queue
    size_t plungestart = mipdata_->num_nodes;
    bool limit_reached = false;

    bool considerHeuristics = true;
    analysis_.mipTimerStart(kMipClockDive);
    while (true) {
      // Possibly apply primal heuristics
      if (considerHeuristics && mipdata_->moreHeuristicsAllowed()) {
        runHeuristics(reduced_search_indices);
        if (reduced_search_indices.empty()) break;
      }

      considerHeuristics = false;

      if (infeasibleWorkerGlobalDomain()) break;
      syncSolutions();

      diveSearches(reduced_search_indices);
      syncSolutions();
      if (reduced_search_indices.empty()) break;

      if (mipdata_->checkLimits()) {
        limit_reached = true;
        break;
      }

      const bool backtrack_plunge =
          backtrackPlunge(reduced_search_indices, plungestart);
      if (!backtrack_plunge) break;

      mipdata_->printDisplayLine();
      // printf("continue plunging due to good estimate\n");
    }  // while (true)
    analysis_.mipTimerStop(kMipClockDive);

    analysis_.mipTimerStart(kMipClockOpenNodesToQueue0);
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.search_ptr_->hasNode()) {
        worker.search_ptr_->openNodesToQueue(mipdata_->nodequeue);
      }
      if (mipdata_->hasMultipleWorkers()) {
        // Remove nodes from worker node queues if backtrack plunged
        while (worker.nodequeue.numNodes() > 0) {
          HighsNodeQueue::OpenNode node =
              std::move(worker.nodequeue.popBestNode());
          mipdata_->nodequeue.emplaceNode(
              std::move(node.domchgstack), std::move(node.branchings),
              node.lower_bound, node.estimate, node.depth);
        }
      }
    }
    analysis_.mipTimerStop(kMipClockOpenNodesToQueue0);

    for (const HighsMipWorker& worker : mipdata_->workers) {
      worker.search_ptr_->flushStatistics();
    }

    if (limit_reached) {
      mipdata_->updateLowerBound(std::min(
          mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
      mipdata_->printDisplayLine();
      break;
    }

    // the search data structures should have no installed node now
    assert(!nodesInstalled());

    // propagate the global domain
    analysis_.mipTimerStart(kMipClockDomainPropgate);
    // sync global domain changes and cut + conflict pools from parallel dives
    syncPools(search_indices);
    syncGlobalDomain(search_indices);
    syncSolutions();
    mipdata_->domain.propagate();
    analysis_.mipTimerStop(kMipClockDomainPropgate);

    analysis_.mipTimerStart(kMipClockPruneInfeasibleNodes);
    mipdata_->pruned_treeweight += mipdata_->nodequeue.pruneInfeasibleNodes(
        mipdata_->domain, mipdata_->feastol);
    analysis_.mipTimerStop(kMipClockPruneInfeasibleNodes);

    // if global propagation detected infeasibility, stop here
    if (mipdata_->domain.infeasible()) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
      mipdata_->printDisplayLine();
      break;
    }

    mipdata_->updateLowerBound(std::min(
        mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
    mipdata_->printDisplayLine();
    if (mipdata_->nodequeue.empty()) break;

    // reset global domain and sync worker's global domains
    bool spawn_more_workers = num_workers < max_num_workers &&
                              mipdata_->nodequeue.numNodes() > num_workers;
    resetGlobalDomain(spawn_more_workers, mipdata_->hasMultipleWorkers());

    if (!submip && mipdata_->num_nodes >= nextCheck) {
      auto nTreeRestarts = mipdata_->numRestarts - mipdata_->numRestartsRoot;
      double currNodeEstim =
          numNodesLastCheck - mipdata_->num_nodes_before_run +
          (mipdata_->num_nodes - numNodesLastCheck) *
              double(1.0 - mipdata_->pruned_treeweight) /
              std::max(
                  double(mipdata_->pruned_treeweight - treeweightLastCheck),
                  mipdata_->epsilon);
      // printf(
      //     "nTreeRestarts: %d, numNodesThisRun: %ld, numNodesLastCheck: %ld,
      //     " "currNodeEstim: %g, " "prunedTreeWeightDelta: %g,
      //     numHugeTreeEstim: %d, numLeavesThisRun:
      //     "
      //     "%ld\n",
      //     nTreeRestarts, mipdata_->num_nodes -
      //     mipdata_->num_nodes_before_run, numNodesLastCheck -
      //     mipdata_->num_nodes_before_run, currNodeEstim, 100.0 *
      //     double(mipdata_->pruned_treeweight - treeweightLastCheck),
      //     numHugeTreeEstim,
      //     mipdata_->num_leaves - mipdata_->num_leaves_before_run);

      bool doRestart = false;

      double activeIntegerRatio =
          1.0 - mipdata_->percentageInactiveIntegers() / 100.0;
      activeIntegerRatio *= activeIntegerRatio;

      if (!doRestart) {
        double gapReduction = 1.0;
        if (mipdata_->upper_limit != kHighsInf) {
          double oldGap = upperLimLastCheck - lowerBoundLastCheck;
          double newGap = mipdata_->upper_limit - mipdata_->lower_bound;
          gapReduction = oldGap / newGap;
        }

        if (gapReduction < 1.0 + (0.05 / activeIntegerRatio) &&
            currNodeEstim >=
                activeIntegerRatio * 20 *
                    (mipdata_->num_nodes - mipdata_->num_nodes_before_run)) {
          nextCheck = mipdata_->num_nodes + 100;
          ++numHugeTreeEstim;
        } else {
          numHugeTreeEstim = 0;
          treeweightLastCheck = double(mipdata_->pruned_treeweight);
          numNodesLastCheck = mipdata_->num_nodes;
          upperLimLastCheck = mipdata_->upper_limit;
          lowerBoundLastCheck = mipdata_->lower_bound;
        }

        // Possibly prevent restart - necessary for debugging presolve
        // errors: see #1553
        if (options_mip_->mip_allow_restart) {
          int64_t minHugeTreeOffset =
              (mipdata_->num_leaves - mipdata_->num_leaves_before_run) / 1000;
          int64_t minHugeTreeEstim = HighsIntegers::nearestInteger(
              activeIntegerRatio * (10 + minHugeTreeOffset) *
              std::pow(1.5, nTreeRestarts));

          doRestart = numHugeTreeEstim >= minHugeTreeEstim;
        } else {
          doRestart = false;
        }
      } else {
        // count restart due to many fixings within the first 1000 nodes as
        // root restart
        ++mipdata_->numRestartsRoot;
      }

      if (doRestart) {
        highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                     "\nRestarting search from the root node\n");
        mipdata_->performRestart();
        analysis_.mipTimerStop(kMipClockSearch);
        goto restart;
      }
    }  // if (!submip && mipdata_->num_nodes >= nextCheck))

    // remove the iteration limit when installing a new node
    // mipdata_->lp.setIterationLimit();

    // Create new workers if there's sufficient nodes
    if (spawn_more_workers) {
      HighsInt new_max_num_workers =
          std::min(static_cast<HighsInt>(mipdata_->nodequeue.numNodes()),
                   max_num_workers);
      mipdata_->pseudocost.removeChanged();
      if (num_workers == 1) {
        constructAdditionalWorkerData(master_worker);
      }
      for (HighsInt i = num_workers; i != new_max_num_workers; i++) {
        createNewWorker(i);
        num_workers++;
      }
    }

    // loop to install the next node for the search
    double this_node_search_time = -analysis_.mipTimerRead(kMipClockNodeSearch);
    analysis_.mipTimerStart(kMipClockNodeSearch);

    while (!mipdata_->nodequeue.empty()) {
      // Update global pseudo-cost with worker information
      syncGlobalPseudoCost();

      // Get new candidate worker search indices
      getSearchIndicesWithNoNodes(search_indices);
      reduced_search_indices = search_indices;

      // Only update worker's pseudo-costs that have been assigned a node
      resetWorkerPseudoCosts(search_indices);

      installNodes(search_indices, limit_reached);
      if (limit_reached) break;

      // we evaluate the node directly here instead of performing a dive
      // because we first want to check if the node is not fathomed due to
      // new global information before we perform separation rounds for the node
      evaluateNodes(search_indices);

      // if the node was pruned we remove it from the search
      // Warning: Overloading limit_reached with an infeasible status here.
      limit_reached = handlePrunedNodes(reduced_search_indices);
      if (limit_reached) break;
      if (reduced_search_indices.empty()) {
        if (mipdata_->hasMultipleWorkers()) {
          syncGlobalDomain(search_indices);
        }
        resetGlobalDomain(false, mipdata_->hasMultipleWorkers());
        continue;
      }

      bool infeasible = separateAndStoreBasis(reduced_search_indices);
      if (infeasible) break;
      syncSolutions();
      break;
    }  // while(!mipdata_->nodequeue.empty())
    analysis_.mipTimerStop(kMipClockNodeSearch);
    if (analysis_.analyse_mip_time) {
      this_node_search_time += analysis_.mipTimerRead(kMipClockNodeSearch);
      analysis_.node_search_time.push_back(this_node_search_time);
    }
    if (limit_reached) {
      break;
    }
  }  // while(search.hasNode())
  syncSolutions();
  analysis_.mipTimerStop(kMipClockSearch);

  cleanupSolve();
}

void HighsMipSolver::cleanupSolve() {
  for (HighsMipWorker& worker : mipdata_->workers) {
    assert(worker.solutions_.empty());
  }
  if (mipdata_->terminatorActive()) {
    if (mipdata_->terminatorTerminated()) {
      // Indicate that this instance has been interrupted
      modelstatus_ = HighsModelStatus::kHighsInterrupt;
    } else if (!submip) {
      // When sub-MIPs call cleanupSolve(), they generally don't have
      // a termination criterion for the whole MIP solver
      //
      // Possibly allow sub-MIPs to terminate if the time limit is
      // reached
      //
      // No other instance has terminated, so terminate
      mipdata_->terminatorTerminate();
    }
  }

  // Force a final logging line
  mipdata_->printDisplayLine(kSolutionSourceCleanup);
  // Stop the solve clock - which won't be running if presolve
  // determines the model status
  if (analysis_.mipTimerRunning(kMipClockSolve))
    analysis_.mipTimerStop(kMipClockSolve);

  // Need to complete the calculation of P-D integral, checking for NO
  // gap change
  mipdata_->updatePrimalDualIntegral(
      mipdata_->lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
      mipdata_->upper_bound, false);
  analysis_.mipTimerStart(kMipClockPostsolve);

  bool havesolution = solution_objective_ != kHighsInf;
  bool feasible;
  if (havesolution)
    feasible =
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance;
  else
    feasible = false;

  dual_bound_ = mipdata_->lower_bound;
  if (mipdata_->objectiveFunction.isIntegral()) {
    double rounded_lower_bound =
        std::ceil(mipdata_->lower_bound *
                      mipdata_->objectiveFunction.integralScale() -
                  mipdata_->feastol) /
        mipdata_->objectiveFunction.integralScale();
    dual_bound_ = std::max(dual_bound_, rounded_lower_bound);
  }
  dual_bound_ += model_->offset_;
  primal_bound_ = mipdata_->upper_bound + model_->offset_;
  node_count_ = mipdata_->num_nodes;
  total_lp_iterations_ = mipdata_->total_lp_iterations;
  dual_bound_ = std::min(dual_bound_, primal_bound_);
  primal_dual_integral_ = mipdata_->primal_dual_integral.value;

  // adjust objective sense in case of maximization problem
  if (orig_model_->sense_ == ObjSense::kMaximize) {
    dual_bound_ = -dual_bound_;
    primal_bound_ = -primal_bound_;
  }

  if (modelstatus_ == HighsModelStatus::kNotset ||
      modelstatus_ == HighsModelStatus::kInfeasible) {
    if (feasible && havesolution)
      modelstatus_ = HighsModelStatus::kOptimal;
    else
      modelstatus_ = HighsModelStatus::kInfeasible;
  }

  analysis_.mipTimerStop(kMipClockPostsolve);
  timer_.stop();

  std::string solutionstatus = "-";

  if (havesolution) {
    // Surely this definition of feasible is unnecessary
    bool lc_feasible =
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance;
    assert(feasible == lc_feasible);
    solutionstatus = feasible ? "feasible" : "infeasible";
  }

  gap_ = fabs(primal_bound_ - dual_bound_);
  if (primal_bound_ == 0.0)
    gap_ = dual_bound_ == 0.0 ? 0.0 : kHighsInf;
  else if (primal_bound_ != kHighsInf)
    gap_ = fabs(primal_bound_ - dual_bound_) / fabs(primal_bound_);
  else
    gap_ = kHighsInf;

  std::array<char, 128> gapString =
      getGapString(gap_, primal_bound_, options_mip_);

  bool timeless_log = options_mip_->timeless_log;
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "\nSolving report\n");
  if (this->orig_model_->model_name_.length())
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Model             %s\n",
                 this->orig_model_->model_name_.c_str());
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Status            %s\n"
               "  Primal bound      %.12g\n"
               "  Dual bound        %.12g\n"
               "  Gap               %s\n",
               utilModelStatusToString(modelstatus_).c_str(), primal_bound_,
               dual_bound_, gapString.data());
  if (!timeless_log)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  P-D integral      %.12g\n",
                 mipdata_->primal_dual_integral.value);
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Solution status   %s\n", solutionstatus.c_str());
  if (solutionstatus != "-")
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "                    %.12g (objective)\n"
                 "                    %.12g (bound viol.)\n"
                 "                    %.12g (int. viol.)\n"
                 "                    %.12g (row viol.)\n",
                 solution_objective_, bound_violation_, integrality_violation_,
                 row_violation_);
  if (!timeless_log) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Timing            %.2f\n", timer_.read());
    if (analysis_.analyse_mip_time)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "                    %.2f (presolve)\n"
                   "                    %.2f (solve)\n"
                   "                    %.2f (postsolve)\n",
                   analysis_.mipTimerRead(kMipClockPresolve),
                   analysis_.mipTimerRead(kMipClockSolve),
                   analysis_.mipTimerRead(kMipClockPostsolve));
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Max sub-MIP depth %d\n"
               "  Nodes             %llu\n",
               int(max_submip_level), (long long unsigned)mipdata_->num_nodes);
  if (mipdata_->total_repair_lp) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Repair LPs        %llu (%llu feasible; %llu iterations)\n",
                 (long long unsigned)mipdata_->total_repair_lp,
                 (long long unsigned)mipdata_->total_repair_lp_feasible,
                 (long long unsigned)mipdata_->total_repair_lp_iterations);
  } else {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Repair LPs        0\n");
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  LP iterations     %llu\n",
               (long long unsigned)mipdata_->total_lp_iterations);
  if (mipdata_->total_lp_iterations)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "                    %llu (strong br.)\n"
                 "                    %llu (separation)\n"
                 "                    %llu (heuristics)\n",
                 (long long unsigned)mipdata_->sb_lp_iterations,
                 (long long unsigned)mipdata_->sepa_lp_iterations,
                 (long long unsigned)mipdata_->heuristic_lp_iterations);

  if (!timeless_log) analysis_.reportMipTimer();

  analysis_.checkSubSolverCallTime(sub_solver_call_time_);

  assert(modelstatus_ != HighsModelStatus::kNotset);

  if (improving_solution_file_ != nullptr) fclose(improving_solution_file_);
}

// Only called in Highs::runPresolve
void HighsMipSolver::runMipPresolve(const HighsInt presolve_reduction_limit) {
  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  mipdata_->init();
  mipdata_->runMipPresolve(presolve_reduction_limit);
}

const HighsLp& HighsMipSolver::getPresolvedModel() const {
  return mipdata_->presolvedModel;
}

HighsPresolveStatus HighsMipSolver::getPresolveStatus() const {
  return mipdata_->presolve_status;
}

presolve::HighsPostsolveStack HighsMipSolver::getPostsolveStack() const {
  return mipdata_->postSolveStack;
}

void HighsMipSolver::callbackGetCutPool() const {
  assert(callback_->user_callback);
  assert(callback_->callbackActive(kCallbackMipGetCutPool));
  callback_->clearHighsCallbackOutput();
  HighsCallbackOutput& data_out = callback_->data_out;

  HighsSparseMatrix cut_matrix;
  mipdata_->lp.getCutPool(data_out.cutpool_num_col, data_out.cutpool_num_cut,
                          data_out.cutpool_lower, data_out.cutpool_upper,
                          cut_matrix);

  // take ownership
  data_out.cutpool_start = std::move(cut_matrix.start_);
  data_out.cutpool_index = std::move(cut_matrix.index_);
  data_out.cutpool_value = std::move(cut_matrix.value_);

  const bool interrupt = mipdata_->interruptFromCallbackWithData(
      kCallbackMipGetCutPool, solution_objective_, "MIP cut pool");
  assert(!interrupt);
}

std::array<char, 128> getGapString(const double gap_,
                                   const double primal_bound_,
                                   const HighsOptions* options_mip_) {
  std::array<char, 128> gapString = {};
  if (gap_ == kHighsInf)
    std::strcpy(gapString.data(), "inf");
  else {
    double printTol = std::max(std::min(1e-2, 1e-1 * gap_), 1e-6);
    auto gapValString = highsDoubleToString(100.0 * gap_, printTol);
    double gapTol = options_mip_->mip_rel_gap;

    if (options_mip_->mip_abs_gap > options_mip_->mip_feasibility_tolerance) {
      gapTol = primal_bound_ == 0.0
                   ? kHighsInf
                   : std::max(gapTol,
                              options_mip_->mip_abs_gap / fabs(primal_bound_));
    }

    if (gapTol == 0.0)
      std::snprintf(gapString.data(), gapString.size(), "%s%%",
                    gapValString.data());
    else if (gapTol != kHighsInf) {
      printTol = std::max(std::min(1e-2, 1e-1 * gapTol), 1e-6);
      auto gapTolString = highsDoubleToString(100.0 * gapTol, printTol);
      std::snprintf(gapString.data(), gapString.size(),
                    "%s%% (tolerance: %s%%)", gapValString.data(),
                    gapTolString.data());
    } else
      std::snprintf(gapString.data(), gapString.size(), "%s%% (tolerance: inf)",
                    gapValString.data());
  }

  return gapString;
}

bool HighsMipSolver::solutionFeasible(const HighsLp* lp,
                                      const std::vector<double>& col_value,
                                      const std::vector<double>* pass_row_value,
                                      double& bound_violation,
                                      double& row_violation,
                                      double& integrality_violation,
                                      HighsCDouble& obj) const {
  bound_violation = 0;
  row_violation = 0;
  integrality_violation = 0;
  const double mip_feasibility_tolerance =
      options_mip_->mip_feasibility_tolerance;

  obj = lp->offset_;

  if (kAllowDeveloperAssert)
    assert(col_value.size() == static_cast<size_t>(lp->num_col_));
  for (HighsInt i = 0; i != lp->num_col_; ++i) {
    const double value = col_value[i];
    obj += lp->col_cost_[i] * value;

    if (lp->integrality_[i] == HighsVarType::kInteger) {
      integrality_violation =
          std::max(fractionality(value), integrality_violation);
    }

    const double lower = lp->col_lower_[i];
    const double upper = lp->col_upper_[i];
    double primal_infeasibility;
    if (value < lower - mip_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + mip_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    } else
      continue;

    bound_violation = std::max(bound_violation, primal_infeasibility);
  }

  // Check row feasibility if there are a positive number of rows.
  //
  // If there are no rows and pass_row_value is nullptr, then
  // row_value_p is also nullptr since row_value is not resized
  if (lp->num_row_ > 0) {
    std::vector<double> row_value;
    if (pass_row_value) {
      if (kAllowDeveloperAssert)
        assert((*pass_row_value).size() == static_cast<size_t>(lp->num_row_));
    } else {
      calculateRowValuesQuad(*lp, col_value, row_value);
    }
    const double* row_value_p =
        pass_row_value ? (*pass_row_value).data() : row_value.data();
    assert(row_value_p);

    for (HighsInt i = 0; i != lp->num_row_; ++i) {
      const double value = row_value_p[i];
      const double lower = lp->row_lower_[i];
      const double upper = lp->row_upper_[i];

      double primal_infeasibility;
      if (value < lower - mip_feasibility_tolerance) {
        primal_infeasibility = lower - value;
      } else if (value > upper + mip_feasibility_tolerance) {
        primal_infeasibility = value - upper;
      } else
        continue;

      row_violation = std::max(row_violation, primal_infeasibility);
    }
  }

  const bool feasible = bound_violation <= mip_feasibility_tolerance &&
                        integrality_violation <= mip_feasibility_tolerance &&
                        row_violation <= mip_feasibility_tolerance;
  return feasible;
}

std::vector<HighsModelStatus> HighsMipSolver::initialiseTerminatorRecord(
    HighsInt num_instance) const {
  std::vector<HighsModelStatus> record(num_instance, HighsModelStatus::kNotset);
  return record;
}

void HighsMipSolver::initialiseTerminator(HighsInt num_instance_,
                                          HighsInt my_instance_,
                                          HighsModelStatus* record_) {
  this->termination_status_ = HighsModelStatus::kNotset;
  this->terminator_.clear();
  this->terminator_.initialise(num_instance_, my_instance_, record_);
}

void HighsMipSolver::initialiseTerminator(const HighsMipSolver& mip_solver) {
  this->terminator_.clear();
  if (!mip_solver.mipdata_->terminatorActive()) return;
  assert(mip_solver.mipdata_->terminatorConcurrency() > 0);
  this->initialiseTerminator(mip_solver.mipdata_->terminatorConcurrency(),
                             mip_solver.mipdata_->terminatorMyInstance(),
                             mip_solver.terminator_.record);
}

void HighsMipSolver::setParallelLock(bool lock) const {
  if (!mipdata_->hasMultipleWorkers()) return;
  mipdata_->parallel_lock = lock;
  for (HighsConflictPool& conflictpool : mipdata_->conflictpools) {
    conflictpool.setAgeLock(lock);
  }
}
