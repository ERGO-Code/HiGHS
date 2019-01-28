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
#include "SimplexConst.h" // For simplex strategy constants
#include "HighsModelObject.h"
#include <cassert>
#include <vector>
#include <cstring> // For strcmp

/**
 * @brief Class for simplex utilities
 */
class HSimplex {
 public:

  void options(
	       HighsModelObject & highs_model_object, //!< Model object in which simplex options are to be set
	       const HighsOptions& opt                //!< HiGHS options
	       ) {
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;

    // Copy values of HighsOptions for the simplex solver
    // TODO: Get this right with proper simplex strategy
    simplex_info_.simplex_strategy = opt.simplex_strategy;
    simplex_info_.crash_strategy = opt.simplex_crash_strategy;
    simplex_info_.dual_edge_weight_strategy = opt.simplex_dual_edge_weight_strategy;
    simplex_info_.price_strategy = opt.simplex_price_strategy;
    simplex_info_.primal_feasibility_tolerance = opt.primal_feasibility_tolerance;
    simplex_info_.dual_feasibility_tolerance = opt.dual_feasibility_tolerance;
    simplex_info_.dual_objective_value_upper_bound = opt.dual_objective_value_upper_bound;
    simplex_info_.perturb_costs = opt.simplex_perturb_costs;
    simplex_info_.iteration_limit = opt.simplex_iteration_limit;
    simplex_info_.highs_run_time_limit = opt.highs_run_time_limit;

    // Set values of internal options

    // Options for reporting timing
    simplex_info_.reportSimplexInnerClock = false;
    simplex_info_.reportSimplexOuterClock = false;
#ifdef HiGHSDEV
    simplex_info_.reportSimplexPhasesClock = false;
    // Option for analysing simplex iterations
    simplex_info_.analyseLp = false;
    simplex_info_.analyseSimplexIterations = true;//false
    simplex_info_.analyseLpSolution = false;
    simplex_info_.analyseInvertTime = false;
    simplex_info_.analyseRebuildTime = false;
#endif
    
  }

  void computeDualObjectiveValue(
				 HighsModelObject &highs_model_object,
				 int phase = 2
				 ) {
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
