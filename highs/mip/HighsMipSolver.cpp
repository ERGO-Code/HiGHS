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
  assert(!submip || submip_level > 0);
  max_submip_level = 0;
  if (solution.value_valid) {
    // MIP solver doesn't check row residuals, but they should be OK
    // so validate using assert
#ifndef NDEBUG
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

void HighsMipSolver::run() {
  const bool debug_logging = false;  // true;
  modelstatus_ = HighsModelStatus::kNotset;

  if (submip) {
    analysis_.analyse_mip_time = false;
  } else {
    analysis_.timer_ = &this->timer_;
    analysis_.setup(*orig_model_, *options_mip_);
  }
  timer_.start();

  improving_solution_file_ = nullptr;
  if (!submip && options_mip_->mip_improving_solution_file != "")
    improving_solution_file_ =
        fopen(options_mip_->mip_improving_solution_file.c_str(), "w");

  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));

  // todo:ig mipdata_. initialize worker

  analysis_.mipTimerStart(kMipClockPresolve);
  analysis_.mipTimerStart(kMipClockInit);

  mipdata_->init();

  analysis_.mipTimerStop(kMipClockInit);
  analysis_.mipTimerStart(kMipClockRunPresolve);
  mipdata_->runPresolve(options_mip_->presolve_reduction_limit);
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

  mipdata_->runSetup();

  analysis_.mipTimerStop(kMipClockRunSetup);
  if (analysis_.analyse_mip_time && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed setup\n", timer_.read());

  // Initialize master worker.
  // Now the worker lives in mipdata.
  // The master worker is used in evaluateRootNode.
  if (mipdata_->domain.infeasible()) {
    cleanupSolve();
    return;
  }
  mipdata_->workers.emplace_back(*this, &mipdata_->lp, &mipdata_->domain,
                                 &mipdata_->cutpool, &mipdata_->conflictPool);

  HighsMipWorker& master_worker = mipdata_->workers.at(0);

restart:
  if (modelstatus_ == HighsModelStatus::kNotset) {
    // Check limits have not been reached before evaluating root node
    if (mipdata_->checkLimits()) {
      cleanupSolve();
      return;
    }
    // Possibly look for primal solution from the user
    if (!submip && callback_->user_callback &&
        callback_->active[kCallbackMipUserSolution])
      mipdata_->callbackUserSolution(solution_objective_,
                                     kUserMipSolutionCallbackOriginAfterSetup);

    if (options_mip_->mip_allow_feasibility_jump) {
      // Apply the feasibility jump before evaluating the root node
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
      const bool bailout_after_feasibility_jump = false;
      if (bailout_after_feasibility_jump) {
        highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                     "HighsMipSolver: Bailing out after Feasibility Jump with "
                     "model status = %s\n",
                     utilModelStatusToString(returned_model_status).c_str());
        modelstatus_ = HighsModelStatus::kInterrupt;
        cleanupSolve();
        return;
      }
    }
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
    if (analysis_.analyse_mip_time && !submip)
      if (analysis_.analyse_mip_time & !submip)
        highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                     "MIP-Timing: %11.2g - starting evaluate root node\n",
                     timer_.read());
    analysis_.mipTimerStart(kMipClockEvaluateRootNode);

    mipdata_->evaluateRootNode(master_worker);

    analysis_.mipTimerStop(kMipClockEvaluateRootNode);
    // Sometimes the analytic centre calculation is not completed when
    // evaluateRootNode returns, so stop its clock if it's running
    if (analysis_.analyse_mip_time &&
        analysis_.mipTimerRunning(kMipClockIpmSolveLp))
      analysis_.mipTimerStop(kMipClockIpmSolveLp);
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

  printf(
      "MIPSOLVER mipdata lp deque member  with address %p, %d "
      "columns, and %d rows\n",
      (void*)&mipdata_->lps.at(0),
      int(mipdata_->lps.at(0).getLpSolver().getNumCol()),
      int(mipdata_->lps.at(0).getLpSolver().getNumRow()));

  printf("Passed to search atm:  \n");

  printf(
      "MIPSOLVER mipdata lp ref with address %p, %d "
      "columns, and %d rows\n",
      (void*)&mipdata_->lp, int(mipdata_->lp.getLpSolver().getNumCol()),
      int(mipdata_->lp.getLpSolver().getNumRow()));

  printf(
      "master_worker lprelaxation_ member  with address %p, %d "
      "columns, and %d rows\n",
      (void*)&master_worker.lprelaxation_,
      int(master_worker.lprelaxation_->getLpSolver().getNumCol()),
      int(mipdata_->lps.at(0).getLpSolver().getNumRow()));

  std::shared_ptr<const HighsBasis> basis;

  double prev_lower_bound = mipdata_->lower_bound;

  mipdata_->lower_bound = mipdata_->nodequeue.getBestLowerBound();

  bool bound_change = mipdata_->lower_bound != prev_lower_bound;
  if (!submip && bound_change)
    mipdata_->updatePrimalDualIntegral(prev_lower_bound, mipdata_->lower_bound,
                                       mipdata_->upper_bound,
                                       mipdata_->upper_bound);

  mipdata_->printDisplayLine();
  int64_t num_nodes = mipdata_->nodequeue.numNodes();
  if (num_nodes > 1) {
    // Should be exactly one node on the queue?
    if (debug_logging)
      printf(
          "HighsMipSolver::run() popping node from nodequeue with %d > 1 "
          "nodes\n",
          HighsInt(num_nodes));
    assert(num_nodes == 1);
  }

  // Initialize worker relaxations and mipworkers
  const HighsInt mip_search_concurrency = options_mip_->mip_search_concurrency;
  const HighsInt num_workers =
      highs::parallel::num_threads() == 1 || mip_search_concurrency <= 1
          ? 1
          : mip_search_concurrency * highs::parallel::num_threads();
  highs::parallel::TaskGroup tg;

  auto destroyOldWorkers = [&]() {
    if (mipdata_->workers.size() <= 1) return;
    while (mipdata_->cutpools.size() > 1) {
      mipdata_->cutpools.pop_back();
    }
    while (mipdata_->conflictpools.size() > 1) {
      mipdata_->conflictpools.pop_back();
    }
    while (mipdata_->domains.size() > 1) {
      mipdata_->domains.pop_back();
    }
    while (mipdata_->lps.size() > 1) {
      mipdata_->lps.pop_back();
    }
    while (mipdata_->workers.size() > 1) {
      mipdata_->workers.pop_back();
    }
  };

  auto constructAdditionalWorkerData = [&](HighsMipWorker& worker) {
    // A use case: Change pointer in master worker to local copies of global
    // info
    assert(mipdata_->cutpools.size() == 1 &&
           mipdata_->conflictpools.size() == 1);
    assert(&worker == &mipdata_->workers.at(0));
    mipdata_->cutpools.emplace_back(numCol(), options_mip_->mip_pool_age_limit,
                                    options_mip_->mip_pool_soft_limit, 1);
    worker.cutpool_ = &mipdata_->cutpools.back();
    mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                         options_mip_->mip_pool_soft_limit);
    worker.conflictpool_ = &mipdata_->conflictpools.back();
    mipdata_->domains.emplace_back(mipdata_->domain);
    worker.globaldom_ = &mipdata_->domains.back();
    worker.globaldom_->addCutpool(*worker.cutpool_);
    worker.globaldom_->addConflictPool(*worker.conflictpool_);
    worker.resetSearch();
    worker.lprelaxation_->setMipWorker(worker);
  };

  auto createNewWorker = [&](HighsInt i) {
    mipdata_->domains.emplace_back(mipdata_->domain);
    mipdata_->lps.emplace_back(mipdata_->lp);
    mipdata_->cutpools.emplace_back(numCol(), options_mip_->mip_pool_age_limit,
                                    options_mip_->mip_pool_soft_limit, i + 1);
    mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                         options_mip_->mip_pool_soft_limit);
    mipdata_->workers.emplace_back(
        *this, &mipdata_->lps.back(), &mipdata_->domains.back(),
        &mipdata_->cutpools.back(), &mipdata_->conflictpools.back());
    mipdata_->lps.back().setMipWorker(mipdata_->workers.back());
    mipdata_->lp.notifyCutPoolsLpCopied(1);
  };

  auto resetGlobalDomain = [&](bool force = false) -> void {
    // if global propagation found bound changes, we update the domain
    if (!mipdata_->domain.getChangedCols().empty() || force) {
      analysis_.mipTimerStart(kMipClockUpdateLocalDomain);
      highsLogDev(options_mip_->log_options, HighsLogType::kInfo,
                  "added %" HIGHSINT_FORMAT " global bound changes\n",
                  (HighsInt)mipdata_->domain.getChangedCols().size());
      mipdata_->cliquetable.cleanupFixed(mipdata_->domain);
      for (HighsInt col : mipdata_->domain.getChangedCols())
        mipdata_->implications.cleanupVarbounds(col);

      mipdata_->domain.setDomainChangeStack(std::vector<HighsDomainChange>());
      mipdata_->domain.clearChangedCols();
      mipdata_->removeFixedIndices();
      analysis_.mipTimerStop(kMipClockUpdateLocalDomain);
    }
  };

  // TODO: Should we be propagating this first?
  if (num_workers > 1) resetGlobalDomain(true);
  destroyOldWorkers();
  if (num_workers > 1) constructAdditionalWorkerData(master_worker);
  master_worker.upper_bound = mipdata_->upper_bound;
  assert(master_worker.solutions_.empty());
  master_worker.solutions_.clear();
  for (HighsInt i = 1; i != num_workers; ++i) {
    createNewWorker(i);
  }

  master_worker.resetSepa();
  HighsSearch& search = *master_worker.search_ptr_;
  mipdata_->debugSolution.registerDomain(search.getLocalDomain());
  HighsSeparation& sepa = *master_worker.sepa_ptr_;

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

  // Lambda for combining limit_reached across searches
  auto limitReached = [&]() -> bool {
    bool limit_reached = false;
    for (HighsInt iSearch = 0; iSearch < mip_search_concurrency; iSearch++)
      limit_reached = limit_reached ||
                      mipdata_->workers[iSearch].search_ptr_->limit_reached_;
    return limit_reached;
  };

  // Lambda checking whether to break out of search
  auto breakSearch = [&]() -> bool {
    bool break_search = false;
    for (HighsInt iSearch = 0; iSearch < mip_search_concurrency; iSearch++)
      break_search =
          break_search || mipdata_->workers[iSearch].search_ptr_->break_search_;
    return break_search;
  };

  auto setParallelLock = [&](bool lock) -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    mipdata_->parallel_lock = lock;
    for (HighsConflictPool& conflictpool : mipdata_->conflictpools) {
      conflictpool.setAgeLock(lock);
    }
  };

  auto syncSolutions = [&]() -> void {
    for (HighsMipWorker& worker : mipdata_->workers) {
      for (auto& sol : worker.solutions_) {
        mipdata_->addIncumbent(std::get<0>(sol), std::get<1>(sol),
                               std::get<2>(sol));
      }
      worker.solutions_.clear();
    }
    // Pass the new upper bound information back to the worker
    for (HighsMipWorker& worker : mipdata_->workers) {
      assert(mipdata_->upper_bound <= worker.upper_bound);
      worker.upper_bound = mipdata_->upper_bound;
    }
  };

  auto syncPools = [&]() -> void {
    if (!mipdata_->hasMultipleWorkers() || mipdata_->parallelLockActive())
      return;
    for (HighsInt i = 1; i < mipdata_->conflictpools.size(); ++i) {
      mipdata_->conflictpools[i].syncConflictPool(mipdata_->conflictPool);
    }
    for (HighsInt i = 1; i < mipdata_->cutpools.size(); ++i) {
      mipdata_->cutpools[i].syncCutPool(*this, mipdata_->cutpool);
    }
    mipdata_->cutpool.performAging();
  };

  auto syncGlobalDomain = [&]() -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    for (HighsMipWorker& worker : mipdata_->workers) {
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

  auto resetWorkerDomains = [&]() -> void {
    // 1. Backtrack to global domain for all local global domains
    // 2. Push all changes from the true global domain
    // 3. Clear changedCols and domChgStack, and reset local search domain for
    // all workers
    // TODO MT: Is it simpler to just copy the domain each time
    if (mipdata_->hasMultipleWorkers()) {
      for (HighsMipWorker& worker : mipdata_->workers) {
        for (const HighsDomainChange& domchg :
             mipdata_->domain.getDomainChangeStack()) {
          worker.getGlobalDomain().changeBound(
              domchg, HighsDomain::Reason::unspecified());
             }
        worker.getGlobalDomain().setDomainChangeStack(
            std::vector<HighsDomainChange>());
        worker.getGlobalDomain().clearChangedCols();
        worker.search_ptr_->resetLocalDomain();
      }
    }
  };

  auto nodesRemaining = [&]() -> bool {
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.search_ptr_->hasNode()) return true;
    }
    return false;
  };

  auto infeasibleGlobalDomain = [&]() -> bool {
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.getGlobalDomain().infeasible()) return true;
    }
    return false;
  };

  auto getSearchIndicesWithNoNodes = [&]() -> std::vector<HighsInt> {
    std::vector<HighsInt> search_indices;
    for (HighsInt i = 0; i < mipdata_->workers.size(); i++) {
      if (!mipdata_->workers[i].search_ptr_->hasNode()) {
        search_indices.emplace_back(i);
      }
    }
    if (search_indices.size() > mipdata_->nodequeue.numActiveNodes()) {
      search_indices.resize(mipdata_->nodequeue.numActiveNodes());
    }
    return search_indices;
  };

  auto getSearchIndicesWithNodes = [&]() -> std::vector<HighsInt> {
    std::vector<HighsInt> search_indices;
    for (HighsInt i = 0; i < mipdata_->workers.size(); i++) {
      if (mipdata_->workers[i].search_ptr_->hasNode()) {
        search_indices.emplace_back(i);
      }
    }
    return search_indices;
  };

  auto installNodes = [&](std::vector<HighsInt>& search_indices,
                          bool& limit_reached) -> void {
    for (HighsInt index : search_indices) {
      if (numQueueLeaves - lastLbLeave >= 10) {
        mipdata_->workers[index].search_ptr_->installNode(
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
        mipdata_->workers[index].search_ptr_->installNode(std::move(nextNode));
      }

      ++numQueueLeaves;

      if (mipdata_->workers[index].search_ptr_->getCurrentEstimate() >=
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

  auto evaluateNodes = [&](std::vector<HighsInt>& search_indices) -> void {
    std::vector<HighsSearch::NodeResult> search_results(search_indices.size());
    analysis_.mipTimerStart(kMipClockEvaluateNode1);
    setParallelLock(true);
    for (HighsInt i = 0; i != search_indices.size(); i++) {
      if (mipdata_->parallelLockActive() && search_indices.size() > 1) {
        tg.spawn([&, i]() {
          search_results[i] =
              mipdata_->workers[search_indices[i]].search_ptr_->evaluateNode();
        });
      } else {
        search_results[i] =
            mipdata_->workers[search_indices[i]].search_ptr_->evaluateNode();
      }
    }
    if (mipdata_->parallelLockActive()) tg.taskWait();
    setParallelLock(false);
    analysis_.mipTimerStop(kMipClockEvaluateNode1);
    for (HighsInt i = 0; i != search_indices.size(); i++) {
      if (search_results[i] == HighsSearch::NodeResult::kSubOptimal) {
        analysis_.mipTimerStart(kMipClockCurrentNodeToQueue);
        mipdata_->workers[search_indices[i]].search_ptr_->currentNodeToQueue(
            mipdata_->nodequeue);
        analysis_.mipTimerStop(kMipClockCurrentNodeToQueue);
      }
    }
  };

  auto doHandlePrunedNodes = [&](HighsInt index, bool thread_safe, bool& flush,
                                 bool& infeasible) {
    HighsDomain& globaldom = mipdata_->workers[index].getGlobalDomain();
    mipdata_->workers[index].search_ptr_->backtrack();
    if (!thread_safe) {
      ++mipdata_->num_leaves;
      ++mipdata_->num_nodes;
      mipdata_->workers[index].search_ptr_->flushStatistics();
    } else {
      flush = true;
    }

    globaldom.propagate();
    if (!thread_safe) {
      mipdata_->pruned_treeweight += mipdata_->nodequeue.pruneInfeasibleNodes(
          mipdata_->domain, mipdata_->feastol);
    }

    if (globaldom.infeasible()) {
      infeasible = true;
      if (!thread_safe) {
        mipdata_->nodequeue.clear();
        mipdata_->pruned_treeweight = 1.0;

        double prev_lower_bound = mipdata_->lower_bound;

        mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);

        bool bound_change = mipdata_->lower_bound != prev_lower_bound;
        if (!submip && bound_change)
          mipdata_->updatePrimalDualIntegral(
              prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
              mipdata_->upper_bound);
      }
      return;
    }

    if (!thread_safe && mipdata_->checkLimits()) {
      return;
    }

    double prev_lower_bound = mipdata_->lower_bound;

    if (!thread_safe) {
      mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                       mipdata_->nodequeue.getBestLowerBound());
    }

    bool bound_change = mipdata_->lower_bound != prev_lower_bound;
    if (!submip && bound_change)
      mipdata_->updatePrimalDualIntegral(
          prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
          mipdata_->upper_bound);

    if (!thread_safe) {
      assert(index == 0);
      resetGlobalDomain();
    }

    analysis_.mipTimerStop(kMipClockNodePrunedLoop);
  };

  auto handlePrunedNodes =
      [&](std::vector<HighsInt>& search_indices) -> std::pair<bool, bool> {
    analysis_.mipTimerStart(kMipClockNodePrunedLoop);
    // If flush then change statistics for all searches where this was the case
    // If infeasible then global domain is infeasible and stop the solve
    // If limit_reached then return something appropriate
    // In multi-thread case now check limits again after everything has been
    // flushed
    std::deque<bool> infeasible(search_indices.size(), false);
    std::deque<bool> flush(search_indices.size(), false);
    std::vector<bool> prune(search_indices.size(), false);
    setParallelLock(true);
    for (HighsInt i = 0; i < search_indices.size(); i++) {
      if (!mipdata_->workers[search_indices[i]]
               .search_ptr_->currentNodePruned())
        continue;
      if (mipdata_->parallelLockActive() && search_indices.size() > 1) {
        tg.spawn([&, i]() {
          doHandlePrunedNodes(search_indices[i], mipdata_->parallelLockActive(),
                              flush[i], infeasible[i]);
        });
      } else {
        doHandlePrunedNodes(search_indices[i], mipdata_->parallelLockActive(),
                            flush[i], infeasible[i]);
      }
      // This search object is "finished" and needs a new node
      prune[i] = true;
    }
    if (mipdata_->parallelLockActive()) {
      tg.taskWait();
      for (HighsInt i = 0; i < search_indices.size() && i < flush.size(); i++) {
        if (flush[i]) {
          ++mipdata_->num_leaves;
          ++mipdata_->num_nodes;
          mipdata_->workers[search_indices[i]].search_ptr_->flushStatistics();
        }
      }
    }
    setParallelLock(false);

    // Remove search indices that need a new node
    HighsInt num_search_indices = static_cast<HighsInt>(search_indices.size());
    for (HighsInt i = num_search_indices - 1; i >= 0; i--) {
      if (prune[i]) {
        num_search_indices--;
        std::swap(search_indices[i], search_indices[num_search_indices]);
      }
    }
    search_indices.resize(num_search_indices);

    for (bool status : infeasible) {
      if (status) {
        mipdata_->nodequeue.clear();
        mipdata_->pruned_treeweight = 1.0;

        double prev_lower_bound = mipdata_->lower_bound;

        mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);

        bool bound_change = mipdata_->lower_bound != prev_lower_bound;
        if (!submip && bound_change)
          mipdata_->updatePrimalDualIntegral(
              prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
              mipdata_->upper_bound);
        analysis_.mipTimerStop(kMipClockNodePrunedLoop);
        return std::make_pair(false, true);
      }
    }

    // Handle case where all nodes have been pruned (and lb hasn't been updated
    // due to parallelism)
    syncSolutions();
    if (mipdata_->hasMultipleWorkers() && num_search_indices == 0) {
      double prev_lower_bound = mipdata_->lower_bound;
      mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                       mipdata_->nodequeue.getBestLowerBound());
      if (!submip && (mipdata_->lower_bound != prev_lower_bound))
        mipdata_->updatePrimalDualIntegral(
            prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
            mipdata_->upper_bound);
    }

    if (mipdata_->checkLimits()) {
      analysis_.mipTimerStop(kMipClockNodePrunedLoop);
      return std::make_pair(true, false);
    }
    analysis_.mipTimerStop(kMipClockNodePrunedLoop);
    return std::make_pair(false, false);
  };

  auto separateAndStoreBasis =
      [&](std::vector<HighsInt>& search_indices) -> bool {
    // the node is still not fathomed, so perform separation
    analysis_.mipTimerStart(kMipClockNodeSearchSeparation);
    setParallelLock(true);
    for (HighsInt i : search_indices) {
      if (mipdata_->parallelLockActive()) {
        tg.spawn([&, i]() {
          mipdata_->workers[i].sepa_ptr_->separate(
              mipdata_->workers[i].search_ptr_->getLocalDomain());
        });
      } else {
        mipdata_->workers[i].sepa_ptr_->separate(
            mipdata_->workers[i].search_ptr_->getLocalDomain());
      }
    }
    analysis_.mipTimerStop(kMipClockNodeSearchSeparation);
    if (mipdata_->parallelLockActive()) tg.taskWait();
    setParallelLock(false);

    for (HighsInt i : search_indices) {
      if (mipdata_->workers[i].getGlobalDomain().infeasible()) {
        mipdata_->workers[i].search_ptr_->cutoffNode();
        analysis_.mipTimerStart(kMipClockOpenNodesToQueue1);
        mipdata_->workers[i].search_ptr_->openNodesToQueue(mipdata_->nodequeue);
        analysis_.mipTimerStop(kMipClockOpenNodesToQueue1);
        mipdata_->nodequeue.clear();
        mipdata_->pruned_treeweight = 1.0;

        analysis_.mipTimerStart(kMipClockStoreBasis);
        double prev_lower_bound = mipdata_->lower_bound;

        mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);

        bool bound_change = mipdata_->lower_bound != prev_lower_bound;
        if (!submip && bound_change)
          mipdata_->updatePrimalDualIntegral(
              prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
              mipdata_->upper_bound);
        return true;
      }
      // after separation we store the new basis and proceed with the outer loop
      // to perform a dive from this node
      if (mipdata_->workers[i].lprelaxation_->getStatus() !=
              HighsLpRelaxation::Status::kError &&
          mipdata_->workers[i].lprelaxation_->getStatus() !=
              HighsLpRelaxation::Status::kNotSet)
        mipdata_->workers[i].lprelaxation_->storeBasis();

      basis = mipdata_->workers[i].lprelaxation_->getStoredBasis();
      if (!basis || !isBasisConsistent(
                        mipdata_->workers[i].lprelaxation_->getLp(), *basis)) {
        HighsBasis b = mipdata_->firstrootbasis;
        b.row_status.resize(mipdata_->workers[i].lprelaxation_->numRows(),
                            HighsBasisStatus::kBasic);
        basis = std::make_shared<const HighsBasis>(std::move(b));
        mipdata_->workers[i].lprelaxation_->setStoredBasis(basis);
      }
    }
    return false;
  };

  auto doRunHeuristics = [&](HighsMipWorker& worker) -> void {
    bool clocks = !mipdata_->parallelLockActive();
    if (clocks) analysis_.mipTimerStart(kMipClockDiveEvaluateNode);
    const HighsSearch::NodeResult evaluate_node_result =
        worker.search_ptr_->evaluateNode();
    if (clocks) analysis_.mipTimerStop(kMipClockDiveEvaluateNode);

    if (evaluate_node_result == HighsSearch::NodeResult::kSubOptimal) return;

    if (worker.search_ptr_->currentNodePruned()) {
      if (clocks) {
        ++mipdata_->num_leaves;
        search.flushStatistics();
      }
    } else {
      if (clocks) analysis_.mipTimerStart(kMipClockDivePrimalHeuristics);
      // TODO MT: Make trivial heuristics work locally
      if (mipdata_->incumbent.empty() && clocks) {
        analysis_.mipTimerStart(kMipClockDiveRandomizedRounding);
        mipdata_->heuristics.randomizedRounding(
            worker,
            worker.lprelaxation_->getLpSolver().getSolution().col_value);
        analysis_.mipTimerStop(kMipClockDiveRandomizedRounding);
      }

      if (mipdata_->incumbent.empty()) {
        if (options_mip_->mip_heuristic_run_rens) {
          if (clocks) analysis_.mipTimerStart(kMipClockDiveRens);
          mipdata_->heuristics.RENS(
              worker,
              worker.lprelaxation_->getLpSolver().getSolution().col_value);
          if (clocks) analysis_.mipTimerStop(kMipClockDiveRens);
        }
      } else {
        if (options_mip_->mip_heuristic_run_rins) {
          if (clocks) analysis_.mipTimerStart(kMipClockDiveRins);
          mipdata_->heuristics.RINS(
              worker,
              worker.lprelaxation_->getLpSolver().getSolution().col_value);
          if (clocks) analysis_.mipTimerStop(kMipClockDiveRins);
        }
      }

      if (clocks) mipdata_->heuristics.flushStatistics(master_worker);
      if (clocks) analysis_.mipTimerStop(kMipClockDivePrimalHeuristics);
    }
  };

  auto runHeuristics = [&]() -> void {
    setParallelLock(true);
    std::vector<HighsInt> search_indices = getSearchIndicesWithNodes();
    for (HighsInt i : search_indices) {
      if (mipdata_->parallelLockActive()) {
        tg.spawn([&, i]() { doRunHeuristics(mipdata_->workers[i]); });
      } else {
        doRunHeuristics(mipdata_->workers[i]);
      }
    }
    if (mipdata_->parallelLockActive()) {
      tg.taskWait();
      for (const HighsInt i : search_indices) {
        if (mipdata_->workers[i].search_ptr_->currentNodePruned()) {
          ++mipdata_->num_leaves;
          search.flushStatistics();
        } else {
          mipdata_->heuristics.flushStatistics(mipdata_->workers[i]);
        }
      }
    }
    setParallelLock(false);
  };

  auto diveAllSearches = [&]() -> bool {
    std::vector<double> dive_times(mipdata_->workers.size(),
                                   -analysis_.mipTimerRead(kMipClockTheDive));
    analysis_.mipTimerStart(kMipClockTheDive);
    std::vector<HighsSearch::NodeResult> dive_results(
        mipdata_->workers.size(), HighsSearch::NodeResult::kBranched);
    setParallelLock(true);
    if (mipdata_->workers.size() > 1) {
      for (int i = 0; i < mipdata_->workers.size(); i++) {
        tg.spawn([&, i]() {
          if (!mipdata_->workers[i].search_ptr_->hasNode() ||
              mipdata_->workers[i].search_ptr_->currentNodePruned()) {
            dive_times[i] = -1;
          } else {
            dive_results[i] = mipdata_->workers[i].search_ptr_->dive();
            dive_times[i] += analysis_.mipTimerRead(kMipClockNodeSearch);
          }
        });
      }
      tg.taskWait();
    } else {
      if (!search.currentNodePruned()) {
        dive_results[0] = search.dive();
        dive_times[0] += analysis_.mipTimerRead(kMipClockNodeSearch);
      }
    }
    analysis_.mipTimerStop(kMipClockTheDive);
    setParallelLock(false);
    bool suboptimal = false;
    for (int i = 0; i < mipdata_->workers.size(); i++) {
      if (dive_times[i] != -1) {
        analysis_.dive_time.push_back(dive_times[i]);
        if (dive_results[i] == HighsSearch::NodeResult::kSubOptimal) {
          suboptimal = true;
        } else {
          ++mipdata_->num_leaves;
          mipdata_->workers[i].search_ptr_->flushStatistics();
        }
      }
    }
    return suboptimal;
  };

  while (nodesRemaining()) {
    // Possibly look for primal solution from the user
    if (!submip && callback_->user_callback &&
        callback_->active[kCallbackMipUserSolution])
      mipdata_->callbackUserSolution(solution_objective_,
                                     kUserMipSolutionCallbackOriginBeforeDive);

    analysis_.mipTimerStart(kMipClockPerformAging1);
    for (HighsConflictPool& conflictpool : mipdata_->conflictpools) {
      conflictpool.performAging();
    }
    analysis_.mipTimerStop(kMipClockPerformAging1);
    // set iteration limit for each lp solve during the dive to 10 times the
    // average nodes

    HighsInt iterlimit = 10 * std::max(mipdata_->lp.getAvgSolveIters(),
                                       mipdata_->avgrootlpiters);
    iterlimit = std::max({HighsInt{10000}, iterlimit,
                          HighsInt((3 * mipdata_->firstrootlpiters) / 2)});

    mipdata_->lp.setIterationLimit(iterlimit);
    for (HighsLpRelaxation& lp : mipdata_->lps) {
      lp.setIterationLimit(iterlimit);
    }

    // perform the dive and put the open nodes to the queue
    size_t plungestart = mipdata_->num_nodes;
    bool limit_reached = false;

    // atm heuristics in the dive break lseu debug64
    // bool considerHeuristics = true;
    bool considerHeuristics = false;

    analysis_.mipTimerStart(kMipClockDive);
    while (true) {
      // Possibly apply primal heuristics
      if (considerHeuristics && mipdata_->moreHeuristicsAllowed()) {
        runHeuristics();
      }

      considerHeuristics = false;

      syncSolutions();
      if (infeasibleGlobalDomain()) break;

      bool suboptimal = diveAllSearches();
      if (suboptimal) break;

      syncSolutions();
      if (mipdata_->checkLimits()) {
        limit_reached = true;
        break;
      }

      if (!mipdata_->hasMultipleWorkers()) {
        HighsInt numPlungeNodes = mipdata_->num_nodes - plungestart;
        if (numPlungeNodes >= 100) break;

        analysis_.mipTimerStart(kMipClockBacktrackPlunge);
        const bool backtrack_plunge =
            search.backtrackPlunge(mipdata_->nodequeue);
        analysis_.mipTimerStop(kMipClockBacktrackPlunge);
        if (!backtrack_plunge) break;
      }

      if (!mipdata_->hasMultipleWorkers()) assert(search.hasNode());

      analysis_.mipTimerStart(kMipClockPerformAging2);
      for (HighsConflictPool& conflictpool : mipdata_->conflictpools) {
        if (conflictpool.getNumConflicts() >
            options_mip_->mip_pool_soft_limit) {
          conflictpool.performAging();
        }
      }
      analysis_.mipTimerStop(kMipClockPerformAging2);

      for (HighsMipWorker& worker : mipdata_->workers) {
        worker.search_ptr_->flushStatistics();
      }
      mipdata_->printDisplayLine();
      if (mipdata_->hasMultipleWorkers()) break;
      // printf("continue plunging due to good estimate\n");
    }  // while (true)
    analysis_.mipTimerStop(kMipClockDive);

    analysis_.mipTimerStart(kMipClockOpenNodesToQueue0);
    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.search_ptr_->hasNode()) {
        worker.search_ptr_->openNodesToQueue(mipdata_->nodequeue);
      }
    }
    analysis_.mipTimerStop(kMipClockOpenNodesToQueue0);

    for (HighsMipWorker& worker : mipdata_->workers) {
      worker.search_ptr_->flushStatistics();
    }

    syncSolutions();

    if (limit_reached) {
      double prev_lower_bound = mipdata_->lower_bound;

      mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                       mipdata_->nodequeue.getBestLowerBound());

      bool bound_change = mipdata_->lower_bound != prev_lower_bound;
      if (!submip && bound_change)
        mipdata_->updatePrimalDualIntegral(
            prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
            mipdata_->upper_bound);
      mipdata_->printDisplayLine();
      break;
    }

    // the search datastructure should have no installed node now
    assert(!nodesRemaining());

    // propagate the global domain
    analysis_.mipTimerStart(kMipClockDomainPropgate);
    // sync global domain changes from parallel dives
    syncGlobalDomain();
    syncPools();
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

      double prev_lower_bound = mipdata_->lower_bound;

      mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);

      bool bound_change = mipdata_->lower_bound != prev_lower_bound;
      if (!submip && bound_change)
        mipdata_->updatePrimalDualIntegral(
            prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
            mipdata_->upper_bound);
      mipdata_->printDisplayLine();
      break;
    }

    // set local global domains of all workers to copy changes of global
    if (mipdata_->hasMultipleWorkers()) resetWorkerDomains();

    double prev_lower_bound = mipdata_->lower_bound;

    mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                     mipdata_->nodequeue.getBestLowerBound());
    bool bound_change = mipdata_->lower_bound != prev_lower_bound;
    if (!submip && bound_change)
      mipdata_->updatePrimalDualIntegral(
          prev_lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
          mipdata_->upper_bound);
    mipdata_->printDisplayLine();
    if (mipdata_->nodequeue.empty()) break;

    // flush all changes made to the global domain
    resetGlobalDomain();

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

    // loop to install the next node for the search
    double this_node_search_time = -analysis_.mipTimerRead(kMipClockNodeSearch);
    analysis_.mipTimerStart(kMipClockNodeSearch);

    while (!mipdata_->nodequeue.empty()) {
      // printf("popping node from nodequeue (length = %" HIGHSINT_FORMAT ")\n",
      // (HighsInt)nodequeue.size());
      std::vector<HighsInt> search_indices = getSearchIndicesWithNoNodes();

      installNodes(search_indices, limit_reached);
      if (limit_reached) break;

      // we evaluate the node directly here instead of performing a dive
      // because we first want to check if the node is not fathomed due to
      // new global information before we perform separation rounds for the node
      evaluateNodes(search_indices);

      // if the node was pruned we remove it from the search and install the
      // next node from the queue
      std::pair<bool, bool> limit_or_infeas = handlePrunedNodes(search_indices);
      if (limit_or_infeas.first) limit_reached = true;
      if (limit_or_infeas.first || limit_or_infeas.second) break;
      // TODO MT: If everything was pruned then do a global sync!
      if (search_indices.empty()) continue;

      bool infeasible = separateAndStoreBasis(search_indices);
      if (infeasible) break;
      syncSolutions();
      break;
    }  // while(!mipdata_->nodequeue.empty())
    analysis_.mipTimerStop(kMipClockNodeSearch);
    if (analysis_.analyse_mip_time) {
      this_node_search_time += analysis_.mipTimerRead(kMipClockNodeSearch);
      analysis_.node_search_time.push_back(this_node_search_time);
    }
    if (limit_reached) break;
  }  // while(search.hasNode())
  analysis_.mipTimerStop(kMipClockSearch);

  cleanupSolve();
}

void HighsMipSolver::cleanupSolve() {
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
    bool feasible =
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance;
    solutionstatus = feasible ? "feasible" : "infeasible";
  }

  gap_ = fabs(primal_bound_ - dual_bound_);
  if (primal_bound_ == 0.0)
    gap_ = dual_bound_ == 0.0 ? 0.0 : kHighsInf;
  else if (primal_bound_ != kHighsInf)
    gap_ = fabs(primal_bound_ - dual_bound_) / fabs(primal_bound_);
  else
    gap_ = kHighsInf;

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
  if (!timeless_log)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Timing            %.2f (total)\n"
                 "                    %.2f (presolve)\n"
                 "                    %.2f (solve)\n"
                 "                    %.2f (postsolve)\n",
                 timer_.read(), analysis_.mipTimerRead(kMipClockPresolve),
                 analysis_.mipTimerRead(kMipClockSolve),
                 analysis_.mipTimerRead(kMipClockPostsolve));
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Max sub-MIP depth %d\n"
               "  Nodes             %llu\n"
               "  Repair LPs        %llu (%llu feasible; %llu iterations)\n"
               "  LP iterations     %llu (total)\n"
               "                    %llu (strong br.)\n"
               "                    %llu (separation)\n"
               "                    %llu (heuristics)\n",
               int(max_submip_level), (long long unsigned)mipdata_->num_nodes,
               (long long unsigned)mipdata_->total_repair_lp,
               (long long unsigned)mipdata_->total_repair_lp_feasible,
               (long long unsigned)mipdata_->total_repair_lp_iterations,
               (long long unsigned)mipdata_->total_lp_iterations,
               (long long unsigned)mipdata_->sb_lp_iterations,
               (long long unsigned)mipdata_->sepa_lp_iterations,
               (long long unsigned)mipdata_->heuristic_lp_iterations);

  if (!timeless_log) analysis_.reportMipTimer();

  assert(modelstatus_ != HighsModelStatus::kNotset);
}

// Only called in Highs::runPresolve
void HighsMipSolver::runPresolve(const HighsInt presolve_reduction_limit) {
  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  mipdata_->init();
  mipdata_->runPresolve(presolve_reduction_limit);
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
  HighsCallbackDataOut& data_out = callback_->data_out;

  std::vector<double> cut_lower;
  std::vector<double> cut_upper;
  HighsSparseMatrix cut_matrix;

  mipdata_->lp.getCutPool(data_out.cutpool_num_col, data_out.cutpool_num_cut,
                          cut_lower, cut_upper, cut_matrix);

  data_out.cutpool_num_nz = cut_matrix.numNz();
  data_out.cutpool_start = cut_matrix.start_.data();
  data_out.cutpool_index = cut_matrix.index_.data();
  data_out.cutpool_value = cut_matrix.value_.data();
  data_out.cutpool_lower = cut_lower.data();
  data_out.cutpool_upper = cut_upper.data();
  callback_->user_callback(kCallbackMipGetCutPool, "MIP cut pool",
                           &callback_->data_out, &callback_->data_in,
                           callback_->user_callback_data);
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
