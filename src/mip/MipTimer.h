/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/MipTimer.h
 * @brief Indices of mip iClocks
 */
#ifndef MIP_MIPTIMER_H_
#define MIP_MIPTIMER_H_

// Clocks for profiling the MIP dual mip solver
enum iClockMip {
  kMipClockTotal = 0,
  kMipClockPresolve,
  kMipClockSolve,
  kMipClockPostsolve,
  // Level 1
  kMipClockInit,
  kMipClockRunPresolve,
  kMipClockRunSetup,
  kMipClockTrivialHeuristics,
  kMipClockEvaluateRootNode,
  kMipClockPerformAging0,
  kMipClockSearch,
  // Search
  kMipClockProbingPresolve,
  kMipClockPerformAging1,
  kMipClockDive,
  kMipClockOpenNodesToQueue0,
  kMipClockDomainPropgate,
  kMipClockPruneInfeasibleNodes,
  kMipClockUpdateLocalDomain,
  kMipClockNodeSearch,

  // Evaluate root node
  kMipClockStartSymmetryDetection,
  kMipClockStartAnalyticCentreComputation,
  kMipClockEvaluateRootLp,
  kMipClockSeparateLpCuts,
  kMipClockRandomizedRounding,
  kMipClockPerformRestart,
  kMipClockRootSeparation,
  kMipClockFinishAnalyticCentreComputation,
  kMipClockRootCentralRounding,
  kMipClockRootSeparationRound0,
  kMipClockRootHeuristicsReducedCost,
  kMipClockRootSeparationRound1,
  kMipClockRootHeuristicsRens,
  kMipClockRootSeparationRound2,
  kMipClockRootFeasibilityPump,
  kMipClockRootSeparationRound3,
  //  kMipClock@,
  //  kMipClock@,
  //  kMipClock@,
  //  kMipClock@,

  kMipClockEvaluateRootNode0,
  kMipClockEvaluateRootNode1,
  kMipClockEvaluateRootNode2,

  // Dive
  kMipClockDiveEvaluateNode,
  kMipClockDivePrimalHeuristics,
  kMipClockTheDive,
  kMipClockBacktrackPlunge,
  kMipClockPerformAging2,

  // Dive primal heuristics
  kMipClockDiveRandomizedRounding,
  kMipClockDiveRens,
  kMipClockDiveRins,

  // NodeSearch
  kMipClockCurrentNodeToQueue,
  kMipClockSearchBacktrack,
  kMipClockNodePrunedLoop,
  kMipClockOpenNodesToQueue1,
  kMipClockEvaluateNode1,
  kMipClockNodeSearchSeparation,
  kMipClockStoreBasis,
  //  kMipClock@,

  // Separation
  kMipClockRootSeparationRound,
  kMipClockRootSeparationFinishAnalyticCentreComputation,
  kMipClockRootSeparationCentralRounding,
  kMipClockRootSeparationEvaluateRootLp,

  // LP solves
  kMipClockSimplexBasisSolveLp,
  kMipClockSimplexNoBasisSolveLp,
  kMipClockIpmSolveLp,

  // Sub-MIP solves
  kMipClockSubMipSolve,

  kMipClockProbingImplications,

  kNumMipClock  //!< Number of MIP clocks
};

const double tolerance_percent_report = 0.1;

class MipTimer {
 public:
  void initialiseMipClocks(HighsTimerClock& mip_timer_clock) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    clock.resize(kNumMipClock);
    clock[kMipClockTotal] = 0;
    clock[kMipClockPresolve] = timer_pointer->clock_def("MIP presolve");
    clock[kMipClockSolve] = timer_pointer->clock_def("MIP solve");
    clock[kMipClockPostsolve] = timer_pointer->clock_def("MIP postsolve");

    // Sometimes the analytic centre clock isn't stopped - because it
    // runs on a separate thread. Although it would be good to
    // understand this better, for now don't assert that this clock
    // has stopped in HighsTimer.h. This is done with a hard-coded
    // clock ID that needs to equal clock[kMipClockIpmSolveLp]
    //
    // Define the clocks for evaluating the LPs first, so that
    // clock[kMipClockIpmSolveLp] isn't changed by inserting new
    // clocks
    clock[kMipClockSimplexBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - simplex basis");
    clock[kMipClockSimplexNoBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - simplex no basis");
    clock[kMipClockIpmSolveLp] = timer_pointer->clock_def("Solve LP: IPM");
    assert(clock[kMipClockIpmSolveLp] == 9);

    // Level 1 - Should correspond to kMipClockTotal
    clock[kMipClockInit] = timer_pointer->clock_def("Initialise");
    clock[kMipClockRunPresolve] = timer_pointer->clock_def("Run presolve");
    clock[kMipClockRunSetup] = timer_pointer->clock_def("Run setup");
    clock[kMipClockTrivialHeuristics] =
        timer_pointer->clock_def("Trivial heuristics");
    clock[kMipClockEvaluateRootNode] =
        timer_pointer->clock_def("Evaluate root node");
    clock[kMipClockPerformAging0] = timer_pointer->clock_def("Perform aging 0");
    clock[kMipClockSearch] = timer_pointer->clock_def("Search");
    // kMipClockPostsolve

    // Evaluate root node
    clock[kMipClockStartSymmetryDetection] =
        timer_pointer->clock_def("Start symmetry detection");
    clock[kMipClockStartAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - start");
    clock[kMipClockEvaluateRootLp] =
        timer_pointer->clock_def("Evaluate root LP");
    clock[kMipClockSeparateLpCuts] =
        timer_pointer->clock_def("Separate LP cuts");
    clock[kMipClockRandomizedRounding] =
        timer_pointer->clock_def("Randomized rounding");
    clock[kMipClockPerformRestart] =
        timer_pointer->clock_def("Perform restart");
    clock[kMipClockRootSeparation] =
        timer_pointer->clock_def("Root separation");
    clock[kMipClockFinishAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - finish");
    clock[kMipClockRootCentralRounding] =
        timer_pointer->clock_def("Root central rounding");
    clock[kMipClockRootSeparationRound0] =
        timer_pointer->clock_def("Root separation round 0");
    clock[kMipClockRootHeuristicsReducedCost] =
        timer_pointer->clock_def("Root heuristics reduced cost");
    clock[kMipClockRootSeparationRound1] =
        timer_pointer->clock_def("Root separation round 1");
    clock[kMipClockRootHeuristicsRens] =
        timer_pointer->clock_def("Root heuristics RENS");
    clock[kMipClockRootSeparationRound2] =
        timer_pointer->clock_def("Root separation round 2");
    clock[kMipClockRootFeasibilityPump] =
        timer_pointer->clock_def("Root feasibility pump");
    clock[kMipClockRootSeparationRound3] =
        timer_pointer->clock_def("Root separation round 3");
    //    clock[kMipClock@] = timer_pointer->clock_def("@");

    clock[kMipClockEvaluateRootNode0] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode0");
    clock[kMipClockEvaluateRootNode1] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode1");
    clock[kMipClockEvaluateRootNode2] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode2");

    // Separation
    clock[kMipClockRootSeparationRound] =
        timer_pointer->clock_def("Separation");
    clock[kMipClockRootSeparationFinishAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - finish");
    clock[kMipClockRootSeparationCentralRounding] =
        timer_pointer->clock_def("Central rounding");
    clock[kMipClockRootSeparationEvaluateRootLp] =
        timer_pointer->clock_def("Evaluate root LP");

    // Presolve - Should correspond to kMipClockRunPresolve
    clock[kMipClockProbingPresolve] =
        timer_pointer->clock_def("Probing - presolve");

    // Search - Should correspond to kMipClockSearch
    clock[kMipClockPerformAging1] = timer_pointer->clock_def("Perform aging 1");
    clock[kMipClockDive] = timer_pointer->clock_def("Dive");
    clock[kMipClockOpenNodesToQueue0] =
        timer_pointer->clock_def("Open nodes to queue 0");
    clock[kMipClockDomainPropgate] =
        timer_pointer->clock_def("Domain propagate");
    clock[kMipClockPruneInfeasibleNodes] =
        timer_pointer->clock_def("Prune infeasible nodes");
    clock[kMipClockUpdateLocalDomain] =
        timer_pointer->clock_def("Update local domain");
    clock[kMipClockNodeSearch] = timer_pointer->clock_def("Node search");
    //    clock[kMipClock@] = timer_pointer->clock_def("@");

    // Dive - Should correspond to kMipClockDive
    clock[kMipClockDiveEvaluateNode] =
        timer_pointer->clock_def("Evaluate node");
    clock[kMipClockDivePrimalHeuristics] =
        timer_pointer->clock_def("Dive primal heuristics");
    clock[kMipClockTheDive] = timer_pointer->clock_def("The dive");
    clock[kMipClockBacktrackPlunge] =
        timer_pointer->clock_def("Backtrack plunge");
    clock[kMipClockPerformAging2] = timer_pointer->clock_def("Perform aging 2");

    // Primal heuristics - Should correspond to kMipDiveClockPrimalHeuristics
    clock[kMipClockDiveRandomizedRounding] =
        timer_pointer->clock_def("Dive Randomized rounding");
    clock[kMipClockDiveRens] = timer_pointer->clock_def("Dive RENS");
    clock[kMipClockDiveRins] = timer_pointer->clock_def("Dive RINS");

    // Node search
    clock[kMipClockCurrentNodeToQueue] =
        timer_pointer->clock_def("Current node to queue");
    clock[kMipClockSearchBacktrack] =
        timer_pointer->clock_def("Search backtrack");
    clock[kMipClockNodePrunedLoop] =
        timer_pointer->clock_def("Pruned loop search");
    clock[kMipClockOpenNodesToQueue1] =
        timer_pointer->clock_def("Open nodes to queue 1");
    clock[kMipClockEvaluateNode1] = timer_pointer->clock_def("Evaluate node 1");
    clock[kMipClockNodeSearchSeparation] =
        timer_pointer->clock_def("Node search separation");
    clock[kMipClockStoreBasis] = timer_pointer->clock_def("Store basis");
    //    clock[] = timer_pointer->clock_def("");

    // Sub-MIP clock
    clock[kMipClockSubMipSolve] = timer_pointer->clock_def("Sub-MIP solves");

    clock[kMipClockProbingImplications] =
        timer_pointer->clock_def("Probing - implications");
    //    clock[] = timer_pointer->clock_def("");
  };

  bool reportMipClockList(const char* grepStamp,
                          const std::vector<HighsInt> mip_clock_list,
                          const HighsTimerClock& mip_timer_clock,
                          const HighsInt kMipClockIdeal = kMipClockTotal,
                          const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    HighsInt mip_clock_list_size = mip_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(mip_clock_list_size);
    for (HighsInt en = 0; en < mip_clock_list_size; en++) {
      clockList[en] = clock[mip_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kMipClockIdeal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void csvMipClockList(const std::string grep_query,
                       const std::string model_name,
                       const std::vector<HighsInt> mip_clock_list,
                       const HighsTimerClock& mip_timer_clock,
                       const HighsInt kMipClockIdeal, const bool header,
                       const bool end_line) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kMipClockIdeal]];
    if (ideal_sum_time < 1e-2) return;
    const HighsInt num_clock = mip_clock_list.size();
    if (header) {
      printf("grep_%s,model,ideal", grep_query.c_str());
      for (HighsInt iX = 0; iX < num_clock; iX++) {
        HighsInt iclock = clock[mip_clock_list[iX]];
        printf(",%s", timer_pointer->clock_names[iclock].c_str());
      }
      printf(",Unaccounted");
      if (end_line) printf("\n");
      return;
    }
    double sum_time = 0;
    printf("grep_%s,%s,%11.4g", grep_query.c_str(), model_name.c_str(),
           ideal_sum_time);
    for (HighsInt iX = 0; iX < num_clock; iX++) {
      HighsInt iclock = clock[mip_clock_list[iX]];
      double time = timer_pointer->read(iclock);
      sum_time += time;
      printf(",%11.4g", time);
    }
    printf(",%11.4g", ideal_sum_time - sum_time);
    if (end_line) printf("\n");
  }

  void reportMipCoreClock(const HighsTimerClock& mip_timer_clock) {
    //    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    const std::vector<HighsInt> mip_clock_list{
        kMipClockPresolve, kMipClockSolve, kMipClockPostsolve};
    reportMipClockList("MipCore_", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);
  };

  void reportMipLevel1Clock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockInit,
                                               kMipClockRunPresolve,
                                               kMipClockRunSetup,
                                               kMipClockTrivialHeuristics,
                                               kMipClockEvaluateRootNode,
                                               kMipClockPerformAging0,
                                               kMipClockSearch,
                                               kMipClockPostsolve};
    reportMipClockList("MipLevl1", mip_clock_list, mip_timer_clock,
                       kMipClockTotal, tolerance_percent_report);
  };

  void reportMipSolveLpClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockSimplexBasisSolveLp,
                                               kMipClockSimplexNoBasisSolveLp,
                                               kMipClockIpmSolveLp};
    reportMipClockList("MipSlvLp", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);  //, tolerance_percent_report);
  };

  void reportMipSubMipSolveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockSubMipSolve};
    reportMipClockList("MipSlvLp", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);  //, tolerance_percent_report);
  };

  void reportMipPresolveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockProbingPresolve};
    reportMipClockList("MipPrslv", mip_clock_list, mip_timer_clock,
                       kMipClockRunPresolve, tolerance_percent_report);
  };

  void reportAltEvaluateRootNodeClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockEvaluateRootNode0,
                                               kMipClockEvaluateRootNode1,
                                               kMipClockEvaluateRootNode2};
    reportMipClockList(
        "AltEvaluateRootNode", mip_clock_list, mip_timer_clock,
        kMipClockEvaluateRootNode);  //, tolerance_percent_report);
  };

  void reportMipEvaluateRootNodeClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockStartSymmetryDetection,
        kMipClockStartAnalyticCentreComputation,
        kMipClockEvaluateRootLp,
        kMipClockSeparateLpCuts,
        kMipClockRandomizedRounding,
        kMipClockPerformRestart,
        kMipClockRootSeparation,
        kMipClockFinishAnalyticCentreComputation,
        kMipClockRootCentralRounding,
        kMipClockRootSeparationRound0,
        kMipClockRootHeuristicsReducedCost,
        kMipClockRootSeparationRound1,
        kMipClockRootHeuristicsRens,
        kMipClockRootSeparationRound2,
        kMipClockRootFeasibilityPump,
        kMipClockRootSeparationRound3
        //	kMipClock@,
        //	kMipClock@
    };
    reportMipClockList(
        "MipEvaluateRootNode", mip_clock_list, mip_timer_clock,
        kMipClockEvaluateRootNode);  //, tolerance_percent_report);
  };

  void reportMipSeparationClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockRootSeparationRound,
        kMipClockRootSeparationFinishAnalyticCentreComputation,
        kMipClockRootSeparationCentralRounding,
        kMipClockRootSeparationEvaluateRootLp};
    reportMipClockList("MipRootSeparation", mip_clock_list, mip_timer_clock,
                       kMipClockRootSeparation);  //, tolerance_percent_report);
  };

  void reportMipSearchClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockPerformAging1,        kMipClockDive,
        kMipClockOpenNodesToQueue0,    kMipClockDomainPropgate,
        kMipClockPruneInfeasibleNodes, kMipClockUpdateLocalDomain,
        kMipClockNodeSearch,
        //	kMipClock@
    };
    reportMipClockList("MipSerch", mip_clock_list, mip_timer_clock,
                       kMipClockSearch, tolerance_percent_report);
  };

  void reportMipDiveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockDiveEvaluateNode, kMipClockDivePrimalHeuristics,
        kMipClockTheDive, kMipClockBacktrackPlunge, kMipClockPerformAging2};
    reportMipClockList("MipDive_", mip_clock_list, mip_timer_clock,
                       kMipClockDive, tolerance_percent_report);
  };

  void reportMipDivePrimalHeuristicsClock(
      const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockDiveRandomizedRounding, kMipClockDiveRens, kMipClockDiveRins};
    reportMipClockList("MipDivePrimalHeuristics", mip_clock_list,
                       mip_timer_clock, kMipClockDivePrimalHeuristics,
                       tolerance_percent_report);
  };

  void reportMipNodeSearchClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockCurrentNodeToQueue, kMipClockNodePrunedLoop,
        //      kMipClockSearchBacktrack,
        kMipClockOpenNodesToQueue1, kMipClockEvaluateNode1,
        kMipClockNodeSearchSeparation};  //, kMipClockStoreBasis};
    reportMipClockList("MipNodeSearch", mip_clock_list, mip_timer_clock,
                       kMipClockNodeSearch);  //, tolerance_percent_report);
  };

  void csvMipClock(const std::string model_name,
                   const HighsTimerClock& mip_timer_clock, const bool header,
                   const bool end_line) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockRunPresolve, kMipClockEvaluateRootNode,
        kMipClockDivePrimalHeuristics, kMipClockTheDive, kMipClockNodeSearch};
    csvMipClockList("csvMIP", model_name, mip_clock_list, mip_timer_clock,
                    kMipClockTotal, header, end_line);
  };

  void csvEvaluateRootNodeClock(const std::string model_name,
                                const HighsTimerClock& mip_timer_clock,
                                const bool header, const bool end_line) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockStartSymmetryDetection,
        kMipClockStartAnalyticCentreComputation,
        kMipClockEvaluateRootLp,
        kMipClockSeparateLpCuts,
        kMipClockRandomizedRounding,
        kMipClockPerformRestart,
        kMipClockRootSeparation,
        kMipClockFinishAnalyticCentreComputation,
        kMipClockRootCentralRounding,
        kMipClockRootSeparationRound0,
        kMipClockRootHeuristicsReducedCost,
        kMipClockRootSeparationRound1,
        kMipClockRootHeuristicsRens,
        kMipClockRootSeparationRound2,
        kMipClockRootFeasibilityPump,
        kMipClockRootSeparationRound3};
    csvMipClockList("csvRootNode", model_name, mip_clock_list, mip_timer_clock,
                    kMipClockEvaluateRootNode, header, end_line);
  };
};

#endif /* MIP_MIPTIMER_H_ */
