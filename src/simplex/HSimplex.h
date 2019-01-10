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
//#include <vector>
//#include <iostream>

/**
 * @brief Class for simplex utilities
 */
class HSimplex {
 public:

  void options(
	       HighsModelObject *ptr_highs_model,
	       const HighsOptions& opt  //!< HiGHS options
	       ) {
    printf("Calling HSimplex::options\n");
    
    HighsSimplexInfo &simplex_info_ = ptr_highs_model->simplex_info_;
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
    HighsLp &lp_ = ptr_highs_model->lp_scaled_;
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

};

//void simplexAllocateWorkAndBaseArrays(std::string& highs_model) {HighsModelObject obj;}

/*
void simplexAllocateWorkAndBaseArrays(HighsModelObject &highs_model) {
  // Allocate bounds and solution spaces
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    simplex_->workCost_[var] = lp_->sense_ * lp_->colCost_[col];
    simplex_->workShift_[var] = 0.;
  }
}

void simplexInitPh2RowCost(HighsModelObject &highs_model, int firstrow, int lastrow) {
  // Zero the cost and shift
  HighsLp *lp_ = &highs_model.lp_scaled_;
  HighsSimplexInfo *simplex_info_ = &highs_model.simplex_info_;
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp_->numCol_ + row;
    simplex_->workCost_[var] = 0;
    simplex_->workShift_[var] = 0.;
  }
}

void simplexInitCost(HighsModelObject &highs_model, int perturb) {
  // Copy the cost
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
    double xpert = (fabs(simplex_->workCost_[i]) + 1) * base * (1 + colRandomValue[i]);
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
    simplex_->workCost_[i] += (0.5 - colRandomValue[i]) * 1e-12;
  }
}

void simplexInitPh2ColBound(HighsModelObject &highs_model, int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
  HighsLp *lp_ = &highs_model.lp_scaled_;
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
