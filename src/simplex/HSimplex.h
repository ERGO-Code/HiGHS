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

/**
 * @brief Class for simplex utilities
 */
class HSimplex {
 public:

  enum SIMPLEX_CRASH_STRATEGY {
    simplex_crash_strategy_off,
    simplex_crash_strategy_df,
    simplex_crash_strategy_ltssf_k,
    simplex_crash_strategy_ltssf_pri,
    simplex_crash_strategy_ltsf_k,
    simplex_crash_strategy_ltsf_pri,
    simplex_crash_strategy_ltsf,
    simplex_crash_strategy_bixby,
    simplex_crash_strategy_bixby_no_nonzero_col_costs,
    simplex_crash_strategy_basic,
    simplex_crash_strategy_test_sing
  };

  enum SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY {
    simplex_dual_edge_weight_strategy_dantzig,
    simplex_dual_edge_weight_strategy_devex,
    simplex_dual_edge_weight_strategy_steepest_edge,
    simplex_dual_edge_weight_strategy_steepest_edge_unit_initial,
    simplex_dual_edge_weight_strategy_steepest_edge_to_devex_switch
  };

  enum SIMPLEX_PRICE_STRATEGY {
    simplex_price_strategy_col,
    simplex_price_strategy_row,
    simplex_price_strategy_row_switch,
    simplex_price_strategy_row_switch_col_switch,
    simplex_price_strategy_row_ultra
  };

  int crash_strategy(const char *crashMode) {
    printf("crashMode = %s\n", crashMode);
    int crashStrategy;
    if (strcmp(crashMode, "off") == 0)
      crashStrategy = simplex_crash_strategy_off;
    else if (strcmp(crashMode, "ltssf") == 0)
      crashStrategy = simplex_crash_strategy_df;
    else if (strcmp(crashMode, "ltssf1") == 0)
      crashStrategy = simplex_crash_strategy_ltssf_k;
    else if (strcmp(crashMode, "ltssf2") == 0)
      crashStrategy = simplex_crash_strategy_ltssf_pri;
    else if (strcmp(crashMode, "ltssf3") == 0)
      crashStrategy = simplex_crash_strategy_ltsf_k;
    else if (strcmp(crashMode, "ltssf4") == 0)
      crashStrategy = simplex_crash_strategy_ltsf_pri;
    else if (strcmp(crashMode, "ltssf5") == 0)
      crashStrategy = simplex_crash_strategy_ltsf;
    else if (strcmp(crashMode, "ltssf6") == 0)
      crashStrategy = simplex_crash_strategy_bixby;
    else if (strcmp(crashMode, "ltssf7") == 0)
      crashStrategy = simplex_crash_strategy_bixby_no_nonzero_col_costs;
    else if (strcmp(crashMode, "bs") == 0)
      crashStrategy = simplex_crash_strategy_basic;
#ifdef HiGHSDEV
    else if (strcmp(crashMode, "tssing") == 0)
      crashStrategy = simplex_crash_strategy_test_sing;
#endif
    else {
      printf("crash_strategy unrecognised crashMode = %s - using No crash\n", crashMode);
      crashStrategy = simplex_crash_strategy_off;
    }
    return crashStrategy;
  }
  
  int dual_edge_weight_strategy(const char *edWtMode) {
    printf("edWtMode = %s\n", edWtMode);
    int dualEdgeWeightStrategy;
    if (strcmp(edWtMode, "dan") == 0)
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_dantzig;
    else if (strcmp(edWtMode, "dvx") == 0)
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_devex;
    else if (strcmp(edWtMode, "sed") == 0)
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_steepest_edge;
    else if (strcmp(edWtMode, "sed0") == 0)
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_steepest_edge_unit_initial;
    else if (strcmp(edWtMode, "dse2dvx") == 0)
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_steepest_edge_to_devex_switch;
    else {
      printf("edWt_strategy unrecognised edWtMode = %s - using using DSE with possible switch to Devex\n", edWtMode);
      dualEdgeWeightStrategy = simplex_dual_edge_weight_strategy_steepest_edge_to_devex_switch;
    }
    return dualEdgeWeightStrategy;
  }

  int price_strategy(const char *priceMode) {
    printf("priceMode = %s\n", priceMode);
    int priceStrategy;
    if (strcmp(priceMode, "col") == 0)
      priceStrategy = simplex_price_strategy_col;
    else if (strcmp(priceMode, "row") == 0)
      priceStrategy = simplex_price_strategy_row;
    else if (strcmp(priceMode, "rowsw") == 0)
      priceStrategy = simplex_price_strategy_row_switch;
    else if (strcmp(priceMode, "rowswcolsw") == 0)
      priceStrategy = simplex_price_strategy_row_switch_col_switch;
    else if (strcmp(priceMode, "ultra") == 0)
      priceStrategy = simplex_price_strategy_row_ultra;
    else {
      printf("price_strategy unrecognised priceMode = %s - using row Price with switch or colump price switch\n", priceMode);
      priceStrategy = simplex_price_strategy_row_switch_col_switch;
    }
    return priceStrategy;
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
    printf("presolveMode = %s\n", opt.presolveMode.c_str());

    simplex_info_.crashStrategy = crash_strategy(opt.crashMode.c_str());
    printf("crashMode =    %s: Strategy = %d\n", opt.crashMode.c_str(), simplex_info_.crashStrategy);
  
    simplex_info_.dualEdgeWeightStrategy = dual_edge_weight_strategy(opt.edWtMode.c_str());
    printf("edWtMode =     %s: Strategy = %d\n", opt.edWtMode.c_str(), simplex_info_.dualEdgeWeightStrategy);

    simplex_info_.priceStrategy = price_strategy(opt.priceMode.c_str());
    printf("priceMode =    %s: Strategy = %d\n", opt.priceMode.c_str(), simplex_info_.priceStrategy);

    //    if (strcmp(EdWt_ArgV, "Dan") == 0)

    // Copy values of HighsOptions for the simplex solver
    simplex_info_.primalFeasibilityTolerance = opt.primalFeasibilityTolerance;
    simplex_info_.dualFeasibilityTolerance = opt.dualFeasibilityTolerance;
    simplex_info_.perturbCosts = opt.perturbCostsSimplex;
    simplex_info_.iterationLimit = opt.iterationLimitSimplex;
    simplex_info_.dualObjectiveValueUpperBound = opt.dualObjectiveValueUpperBound;

    // Set values of internal options
    simplex_info_.pamiCutoff = 0.95;

    // Options for reporting timing
    simplex_info_.reportSimplexInnerClock = false;
    simplex_info_.reportSimplexOuterClock = false;
#ifdef HiGHSDEV
    simplex_info_.reportSimplexPhasesClock = false;
    // Option for analysing simplex iterations
    simplex_info_.analyseLp = false;
    simplex_info_.analyseSimplexIterations = false;
    simplex_info_.analyseLpSolution = false;
    simplex_info_.analyseInvertTime = false;
    simplex_info_.analyseRebuildTime = false;
#endif
    
  }

  void computeDualObjectiveValue(
				 HighsModelObject *ptr_highs_model,
				 int phase = 2) {
    HighsLp &lp_ = ptr_highs_model->solver_lp_;
    HighsSimplexInfo &simplex_info_ = ptr_highs_model->simplex_info_;
    
    simplex_info_.dualObjectiveValue = 0;
    const int numTot = lp_.numCol_ + lp_.numRow_;
    for (int i = 0; i < numTot; i++) {
      if (ptr_highs_model->basis_.nonbasicFlag_[i]) {
	simplex_info_.dualObjectiveValue += simplex_info_.workValue_[i] * simplex_info_.workDual_[i];
      }
    }
    if (phase != 1) {
      simplex_info_.dualObjectiveValue *= ptr_highs_model->scale_.cost_;
      simplex_info_.dualObjectiveValue -= lp_.offset_;
    }
    // Now have dual objective value
    ptr_highs_model->haveDualObjectiveValue = 1;
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

//void simplexAllocateWorkAndBaseArrays(std::string& highs_model) {HighsModelObject obj;}

/*
void simplexAllocateWorkAndBaseArrays(HighsModelObject &highs_model) {
  // Allocate bounds and solution spaces
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  const int numTot = lp_->numCol_ + lp_->numRow_;
  simplex_->workCost_.resize(numTot);
  simplex_->workDual_.resize(numTot);
  // Was workShift.assign(numTot, 0); but shift is populated by call to
  // initCost()
  simplex_->workShift_.resize(numTot);

  simplex_->workLower_.resize(numTot);
  simplex_->workUpper_.resize(numTot);
  simplex_->workRange_.resize(numTot);
  simplex_->workValue_.resize(numTot);

  simplex_->baseLower_.resize(lp_->numRow_);
  simplex_->baseUpper_.resize(lp_->numRow_);
  simplex_->baseValue_.resize(lp_->numRow_); 
}

void simplexInitPh2ColCost(HighsModelObject &highs_model, int firstcol, int lastcol) {
  // Copy the Phase 2 cost and zero the shift
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    simplex_->workCost_[var] = lp_->sense_ * lp_->colCost_[col];
    simplex_->workShift_[var] = 0.;
  }
}

void simplexInitPh2RowCost(HighsModelObject &highs_model, int firstrow, int lastrow) {
  // Zero the cost and shift
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp_->numCol_ + row;
    simplex_->workCost_[var] = 0;
    simplex_->workShift_[var] = 0.;
  }
}

void simplexInitCost(HighsModelObject &highs_model, int perturb) {
  // Copy the cost
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  simplexInitPh2ColCost(highs_model, 0, lp_->numCol_ - 1);
  simplexInitPh2RowCost(highs_model, 0, lp_->numRow_ - 1);
  // See if we want to skip perturbation
  problemPerturbed = 0;
  if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0) return;
  problemPerturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < lp_->numCol_; i++) bigc = max(bigc, fabs(simplex_->workCost_[i]));
  if (bigc > 100) bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  const int numTot = lp_->numCol_ + lp_->numRow_;
  for (int i = 0; i < numTot; i++) boxedRate += (simplex_->workRange_[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) bigc = min(bigc, 1.0);
  if (bigc < 1) {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < lp_->numCol_; i++) {
    double lower = lp_->colLower_[i];
    double upper = lp_->colUpper_[i];
    double xpert = (fabs(simplex_->workCost_[i]) + 1) * base * (1 + numTotRandomValue[i]);
    if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper == HIGHS_CONST_INF) {  // Lower
      simplex_->workCost_[i] += xpert;
    } else if (lower == -HIGHS_CONST_INF) {  // Upper
      simplex_->workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex_->workCost_[i] += (simplex_->workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }

  for (int i = lp_->numCol_; i < numTot; i++) {
    simplex_->workCost_[i] += (0.5 - numTotRandomValue[i]) * 1e-12;
  }
}

void simplexInitPh2ColBound(HighsModelObject &highs_model, int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  assert(firstcol >= 0);
  assert(lastcol < lp_->numCol_);
  for (int col = firstcol; col <= lastcol; col++) {
    simplex_->workLower_[col] = lp_->colLower_[col];
    simplex_->workUpper_[col] = lp_->colUpper_[col];
    simplex_->workRange_[col] = simplex_->workUpper_[col] - simplex_->workLower_[col];
  }
}

void simplexInitPh2RowBound(HighsModelObject &highs_model, int firstrow, int lastrow) {
  // Copy bounds and compute ranges
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  assert(firstrow >= 0);
  assert(lastrow < lp_->numRow_);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp_->numCol_ + row;
    simplex_->workLower_[var] = -lp_->rowUpper_[row];
    simplex_->workUpper_[var] = -lp_->rowLower_[row];
    simplex_->workRange_[var] = simplex_->workUpper_[var] - simplex_->workLower_[var];
  }
}

void SimplexInitBound(HighsModelObject &highs_model, int phase) {
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  SimplexInitPh2ColBound(highs_model, 0, lp_->numCol_ - 1);
  SimplexInitPh2RowBound(highs_model, 0, lp_->numRow_ - 1);
  if (phase == 2) return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HIGHS_CONST_INF;
  const int numTot = lp_->numCol_ + lp_->numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_->workLower_[i] == -inf && simplex_->workUpper_[i] == inf) {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= lp_->numCol_) continue;
      simplex_->workLower_[i] = -1000, simplex_->workUpper_[i] = 1000;  // FREE
    } else if (simplex_->workLower_[i] == -inf) {
      simplex_->workLower_[i] = -1, simplex_->workUpper_[i] = 0;  // UPPER
    } else if (simplex_->workUpper_[i] == inf) {
      simplex_->workLower_[i] = 0, simplex_->workUpper_[i] = 1;  // LOWER
    } else {
      simplex_->workLower_[i] = 0, simplex_->workUpper_[i] = 0;  // BOXED or FIXED
    }
    simplex_->workRange_[i] = simplex_->workUpper_[i] - simplex_->workLower_[i];
  }
}

void simplexInitValueFromNonbasic(HighsModelObject &highs_model, int firstvar, int lastvar) {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  assert(firstvar >= 0);
  const int numTot = lp_->numCol_ + lp_->numRow_;
  assert(lastvar < numTot);
  // double dl_pr_act, norm_dl_pr_act;
  // norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++) {
    if (basis_->nonbasicFlag_[var]) {
      // Nonbasic variable
      // double prev_pr_act = simplex_->workValue_[var];
      if (simplex_->workLower_[var] == simplex_->workUpper_[var]) {
        // Fixed
        simplex_->workValue_[var] = simplex_->workLower_[var];
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      } else if (!highs_isInfinity(-simplex_->workLower_[var])) {
        // Finite lower bound so boxed or lower
        if (!highs_isInfinity(simplex_->workUpper_[var])) {
          // Finite upper bound so boxed
          if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP) {
            // Set at lower
            simplex_->workValue_[var] = simplex_->workLower_[var];
          } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN) {
            // Set at upper
            simplex_->workValue_[var] = simplex_->workUpper_[var];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
            simplex_->workValue_[var] = simplex_->workLower_[var];
          }
        } else {
          // Lower
          simplex_->workValue_[var] = simplex_->workLower_[var];
          basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
        }
      } else if (!highs_isInfinity(simplex_->workUpper_[var])) {
        // Upper
        simplex_->workValue_[var] = simplex_->workUpper_[var];
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        simplex_->workValue_[var] = 0;
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      }
      // dl_pr_act = simplex_->workValue_[var] - prev_pr_act;
      // norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
      //      %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var, simplex_->workLower_[var],
      // simplex_->workValue_[var], simplex_->workUpper_[var], simplex_->workDual_[var], dl_pr_act);
    } else {
      // Basic variable
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
    }
  }
  //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
  //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is
  //  %g\n", norm_dl_pr_act);
}

void SimplexInitValue(HighsModelObject &highs_model) {
  HighsLp *lp_ = &highs_model.solver_lp_;
  const int numTot = lp_->numCol_ + lp_->numRow_;
  SimplexInitValueFromNonbasic(highs_model, 0, numTot - 1);
}

void simplexPopulateWorkArrays(HighsModelObject &highs_model) {
  // Initialize the values
  simplexInitCost(highs_model);
  simplexInitBound(highs_model);
  simplexInitValue(highs_model);
}

void simplexInitWithLogicalBasis(HighsModelObject &highs_model) {
  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  HighsLp *lp_ = &highs_model.solver_lp_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  for (int row = 0; row < lp_->numRow_; row++) basis_->basicIndex_[row] = lp_->numCol_ + row;
  for (int col = 0; col < lp_->numCol_; col++) basis_->nonbasicFlag_[col] = 1;
  numBasicLogicals = lp_->numRow_;

  simplexAllocateWorkAndBaseArrays(HighsModelObject &highs_model);
  simplexPopulateWorkArrays(HighsModelObject &highs_model);

  // Deduce the consequences of a new basis
  //  mlFg_Update(mlFg_action_NewBasis);
}
*/
#endif // SIMPLEX_HSIMPLEX_H_
