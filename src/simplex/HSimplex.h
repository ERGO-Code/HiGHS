/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "HConfig.h"
#include "HighsModelObject.h"
#include <cassert>
#include <vector>
#include <cstring> // For strcmp

enum SIMPLEX_STRATEGY {
  SIMPLEX_STRATEGY_DUAL_PLAIN = 0,
  SIMPLEX_STRATEGY_DUAL_TASKS,
  SIMPLEX_STRATEGY_DUAL_MULTI,
  SIMPLEX_STRATEGY_PRIMAL
};
  
enum SIMPLEX_CRASH_STRATEGY {
  SIMPLEX_CRASH_STRATEGY_OFF = 0,
  SIMPLEX_CRASH_STRATEGY_DF,
  SIMPLEX_CRASH_STRATEGY_LTSSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF,
  SIMPLEX_CRASH_STRATEGY_BIXBY,
  SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS,
  SIMPLEX_CRASH_STRATEGY_BASIC,
  SIMPLEX_CRASH_STRATEGY_TEST_SING
};

enum SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY {
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG = 0,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH
};

enum SIMPLEX_PRICE_STRATEGY {
  SIMPLEX_PRICE_STRATEGY_COL = 0,
  SIMPLEX_PRICE_STRATEGY_ROW,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_ULTRA
};

/**
 * @brief Class for simplex utilities
 */
class HSimplex {
 public:

  int crash_strategy(const char *crashMode) {
#ifdef HiGHSDEV
    printf("crashMode = %s\n", crashMode);
#endif
    int crash_strategy;
    if (strcmp(crashMode, "off") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_OFF;
    else if (strcmp(crashMode, "ltssf") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_DF;
    else if (strcmp(crashMode, "ltssf1") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_LTSSF_K;
    else if (strcmp(crashMode, "ltssf2") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_LTSSF_PRI;
    else if (strcmp(crashMode, "ltssf3") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_LTSF_K;
    else if (strcmp(crashMode, "ltssf4") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_LTSF_PRI;
    else if (strcmp(crashMode, "ltssf5") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_LTSF;
    else if (strcmp(crashMode, "ltssf6") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_BIXBY;
    else if (strcmp(crashMode, "ltssf7") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS;
    else if (strcmp(crashMode, "bs") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_BASIC;
#ifdef HiGHSDEV
    else if (strcmp(crashMode, "tssing") == 0)
      crash_strategy = SIMPLEX_CRASH_STRATEGY_TEST_SING;
#endif
    else {
      printf("crash_strategy unrecognised crashMode = %s - using No crash\n", crashMode);
      crash_strategy = SIMPLEX_CRASH_STRATEGY_OFF;
    }
    return crash_strategy;
  }
  
  int dual_edge_weight_strategy(const char *edWtMode) {
#ifdef HiGHSDEV
    printf("edWtMode = %s\n", edWtMode);
#endif
    int dual_edge_weight_strategy;
    if (strcmp(edWtMode, "dan") == 0)
      dual_edge_weight_strategy = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG;
    else if (strcmp(edWtMode, "dvx") == 0)
      dual_edge_weight_strategy = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX;
    else if (strcmp(edWtMode, "dse") == 0)
      dual_edge_weight_strategy = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE;
    else if (strcmp(edWtMode, "dse0") == 0)
      dual_edge_weight_strategy = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL;
    else if (strcmp(edWtMode, "dse2dvx") == 0)
      dual_edge_weight_strategy = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH;
    else {
      printf("edWt_strategy unrecognised edWtMode = %s - using using DSE with possible switch to Devex\n", edWtMode);
      dual_edge_weight_strategy =SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH;
    }
    return dual_edge_weight_strategy;
  }

  int price_strategy(const char *priceMode) {
#ifdef HiGHSDEV
    printf("priceMode = %s\n", priceMode);
#endif
    int price_strategy;
    if (strcmp(priceMode, "col") == 0)
      price_strategy = SIMPLEX_PRICE_STRATEGY_COL;
    else if (strcmp(priceMode, "row") == 0)
      price_strategy = SIMPLEX_PRICE_STRATEGY_ROW;
    else if (strcmp(priceMode, "rowsw") == 0)
      price_strategy = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH;
    else if (strcmp(priceMode, "rowswcolsw") == 0)
      price_strategy = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH;
    else if (strcmp(priceMode, "ultra") == 0)
      price_strategy = SIMPLEX_PRICE_STRATEGY_ROW_ULTRA;
    else {
      printf("price_strategy unrecognised priceMode = %s - using row Price with switch or colump price switch\n", priceMode);
      price_strategy = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH;
    }
    return price_strategy;
  }

  void options(
	       HighsModelObject *ptr_highs_model,
	       const HighsOptions& opt  //!< HiGHS options
	       ) {
#ifdef HiGHSDEV
    printf("Calling HSimplex::options\n");
#endif
    HighsSimplexInfo &simplex_info_ = ptr_highs_model->simplex_info_;

    // Deduce values of options from HighsOptions strings
#ifdef HiGHSDEV
    printf("presolveMode = %s\n", opt.presolveMode.c_str());
#endif

    if (opt.pami) {
      simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_DUAL_MULTI;
    } else if (opt.sip) {
      simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_DUAL_TASKS;
    } else {
      simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_DUAL_PLAIN;
    }

    simplex_info_.crash_strategy = crash_strategy(opt.crashMode.c_str());
#ifdef HiGHSDEV
    printf("crashMode =    %s: strategy = %d\n", opt.crashMode.c_str(), simplex_info_.crash_strategy);
#endif  

    simplex_info_.dual_edge_weight_strategy = dual_edge_weight_strategy(opt.edWtMode.c_str());
#ifdef HiGHSDEV
    printf("edWtMode =     %s: strategy = %d\n", opt.edWtMode.c_str(), simplex_info_.dual_edge_weight_strategy);
#endif

    simplex_info_.price_strategy = price_strategy(opt.priceMode.c_str());
#ifdef HiGHSDEV
    printf("priceMode =    %s: strategy = %d\n", opt.priceMode.c_str(), simplex_info_.price_strategy);
#endif

    // Copy values of HighsOptions for the simplex solver
    simplex_info_.primal_feasibility_tolerance = opt.primal_feasibility_tolerance;
    simplex_info_.dual_feasibility_tolerance = opt.dual_feasibility_tolerance;
    simplex_info_.perturb_costs = opt.perturb_costs_simplex;
    simplex_info_.iteration_limit = opt.iteration_limit_simplex;
    simplex_info_.dual_objective_value_upper_bound = opt.dual_objective_value_upper_bound;

    // Set values of internal options

    // Options for reporting timing
    simplex_info_.reportSimplexInnerClock = false;
    simplex_info_.reportSimplexOuterClock = false;
#ifdef HiGHSDEV
    simplex_info_.reportSimplexPhasesClock = false;
    // Option for analysing simplex iterations
    simplex_info_.analyseLp = true;//false;
    simplex_info_.analyseSimplexIterations = true;//false;
    simplex_info_.analyseLpSolution = true;//false;
    simplex_info_.analyseInvertTime = true;//false;
    simplex_info_.analyseRebuildTime = true;//false;
#endif
    
  }

  void computeDualObjectiveValue(
				 HighsModelObject highs_model_object,
				 int phase = 2) {
    HighsLp &lp_ = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
    
    simplex_info_.dualObjectiveValue = 0;
    const int numTot = lp_.numCol_ + lp_.numRow_;
    for (int i = 0; i < numTot; i++) {
      if (highs_model_object.basis_.nonbasicFlag_[i]) {
	simplex_info_.dualObjectiveValue += simplex_info_.workValue_[i] * simplex_info_.workDual_[i];
      }
    }
    if (phase != 1) {
      simplex_info_.dualObjectiveValue *= highs_model_object.scale_.cost_;
      simplex_info_.dualObjectiveValue -= lp_.offset_;
    }
    // Now have dual objective value
    highs_model_object.haveDualObjectiveValue = 1;
  }

  void initialiseHighsModelObjectRandomVectors(
				  HighsModelObject &highs_model
				  ) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
    const int numCol = highs_model.solver_lp_.numCol_;
    // Instantiate and (re-)initialise the random number generator
    HighsRandom random;
    random.initialise();
    // Generate a random permutation of the column indices
    simplex_info_.numColPermutation_.resize(numCol);
    int *numColPermutation = &simplex_info_.numColPermutation_[0];
    for (int i = 0; i < numCol; i++) numColPermutation[i] = i;
    for (int i = numCol - 1; i >= 1; i--) {
      int j = random.integer() % (i + 1);
      std::swap(numColPermutation[i], numColPermutation[j]);
    }
    // Generate a vector of random reals numbers 
    const int numTot = highs_model.solver_lp_.numCol_ + highs_model.solver_lp_.numRow_;
    simplex_info_.numTotRandomValue_.resize(numTot);
    double *numTotRandomValue = &simplex_info_.numTotRandomValue_[0];
    for (int i = 0; i < numTot; i++) {
      numTotRandomValue[i] = random.fraction();
    }
  }
};

#endif // SIMPLEX_HSIMPLEX_H_
