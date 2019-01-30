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
#include "HighsIO.h"
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
    
    simplex_info_.transpose_solver_lp = opt.transpose_solver_lp;
    simplex_info_.scale_solver_lp = opt.scale_solver_lp;
    simplex_info_.permute_solver_lp = opt.permute_solver_lp;
    simplex_info_.tighten_solver_lp = opt.tighten_solver_lp;
    
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
  
  void initialiseSolverLpRandomVectors(
				       HighsModelObject &highs_model
				       ) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
    const int numCol = highs_model.solver_lp_.numCol_;
    const int numTot = highs_model.solver_lp_.numCol_ + highs_model.solver_lp_.numRow_;
    // Instantiate and (re-)initialise the random number generator
    HighsRandom random;
    random.initialise();
    // Generate the random vectors in the same order as hsol to
    // maintain repeatable performance
    //
    // Generate a random permutation of all the indices
    simplex_info_.numTotPermutation_.resize(numTot);
    int *numTotPermutation = &simplex_info_.numTotPermutation_[0];
    for (int i = 0; i < numTot; i++) numTotPermutation[i] = i;
    for (int i = numTot - 1; i >= 1; i--) {
      int j = random.integer() % (i + 1);
      std::swap(numTotPermutation[i], numTotPermutation[j]);
    }
    
    // Generate a vector of random reals 
    simplex_info_.numTotRandomValue_.resize(numTot);
    double *numTotRandomValue = &simplex_info_.numTotRandomValue_[0];
    for (int i = 0; i < numTot; i++) {
      numTotRandomValue[i] = random.fraction();
    }
    
    //
    // Generate a random permutation of the column indices
    simplex_info_.numColPermutation_.resize(numCol);
    int *numColPermutation = &simplex_info_.numColPermutation_[0];
    for (int i = 0; i < numCol; i++) numColPermutation[i] = i;
    for (int i = numCol - 1; i >= 1; i--) {
      int j = random.integer() % (i + 1);
      std::swap(numColPermutation[i], numColPermutation[j]);
    }
    
  }
  
  // TRANSPOSE:
  
  void transpose_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called transpose_solver_lp: simplex_info_.transposed_solver_lp = %d\n", simplex_info_.transposed_solver_lp);
#endif
    if (simplex_info_.transposed_solver_lp) return;
    HighsLp& primal_lp = highs_model.lp_;
    
    int primalNumCol = primal_lp.numCol_;
    int primalNumRow = primal_lp.numRow_;
    
    int transposeCancelled = 0;
    if (1.0 * primalNumCol / primalNumRow > 0.2) {
      //        cout << "transpose-cancelled-by-ratio" << endl;
      transposeCancelled = 1;
      return;
    }
    
    vector<int>& primalAstart = primal_lp.Astart_;
    vector<int>& primalAindex = primal_lp.Aindex_;
    vector<double>& primalAvalue = primal_lp.Avalue_;
    vector<double>& primalColCost = primal_lp.colCost_;
    vector<double>& primalColLower = primal_lp.colLower_;
    vector<double>& primalColUpper = primal_lp.colUpper_;
    vector<double>& primalRowLower = primal_lp.rowLower_;
    vector<double>& primalRowUpper = primal_lp.rowUpper_;
    
    HighsLp& dual_lp = highs_model.solver_lp_;
    /*
      vector<int>& dualAstart = dual_lp.Astart_;
      vector<int>& dualAindex = dual_lp.Aindex_;
      vector<double>& dualAvalue = dual_lp.Avalue_;
      vector<double>& dualColCost = dual_lp.colCost_;
      vector<double>& dualColLower = dual_lp.colLower_;
      vector<double>& dualColUpper = dual_lp.colUpper_;
      vector<double>& dualRowLower = dual_lp.rowLower_;
      vector<double>& dualRowUpper = dual_lp.rowUpper_;
    */
    
    // Convert primal cost to dual bound
    const double inf = HIGHS_CONST_INF;
    vector<double> dualRowLower(primalNumCol);
    vector<double> dualRowUpper(primalNumCol);
    for (int j = 0; j < primalNumCol; j++) {
      double lower = primalColLower[j];
      double upper = primalColUpper[j];
      
      /*
       * Primal      Dual
       * Free        row = c
       * x > 0       row < c
       * x < 0       row > c
       * x = 0       row free
       * other       cancel
       */
      
      if (lower == -inf && upper == inf) {
	dualRowLower[j] = primalColCost[j];
	dualRowUpper[j] = primalColCost[j];
      } else if (lower == 0 && upper == inf) {
	dualRowLower[j] = -inf;
	dualRowUpper[j] = primalColCost[j];
      } else if (lower == -inf && upper == 0) {
	dualRowLower[j] = primalColCost[j];
	dualRowUpper[j] = +inf;
      } else if (lower == 0 && upper == 0) {
	dualRowLower[j] = -inf;
	dualRowUpper[j] = +inf;
      } else {
	transposeCancelled = 1;
	break;
      }
    }
    
    // Check flag
    if (transposeCancelled == 1) {
      //        cout << "transpose-cancelled-by-column" << endl;
      return;
    }
    
    // Convert primal row bound to dual variable cost
    vector<double> dualColLower(primalNumRow);
    vector<double> dualColUpper(primalNumRow);
    vector<double> dualCost(primalNumRow);
    for (int i = 0; i < primalNumRow; i++) {
      double lower = primalRowLower[i];
      double upper = primalRowUpper[i];
      
      /*
       * Primal      Dual
       * row = b     Free
       * row < b     y < 0
       * row > b     y > 0
       * row free    y = 0
       * other       cancel
       */
      
      if (lower == upper) {
	dualColLower[i] = -inf;
	dualColUpper[i] = +inf;
	dualCost[i] = -lower;
      } else if (lower == -inf && upper != inf) {
	dualColLower[i] = -inf;
	dualColUpper[i] = 0;
	dualCost[i] = -upper;
      } else if (lower != -inf && upper == inf) {
	dualColLower[i] = 0;
	dualColUpper[i] = +inf;
	dualCost[i] = -lower;
      } else if (lower == -inf && upper == inf) {
	dualColLower[i] = 0;
	dualColUpper[i] = 0;
	dualCost[i] = 0;
      } else {
	transposeCancelled = 1;
	break;
      }
    }
    
    // Check flag
    if (transposeCancelled == 1) {
      //        cout << "transpose-cancelled-by-row" << endl;
      return;
    }
    
    // We can now really transpose things
    vector<int> iwork(primalNumRow, 0);
    vector<int> ARstart(primalNumRow + 1, 0);
    int AcountX = primalAindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++) iwork[primalAindex[k]]++;
    for (int i = 1; i <= primalNumRow; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < primalNumRow; i++) iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < primalNumCol; iCol++) {
      for (int k = primalAstart[iCol]; k < primalAstart[iCol + 1]; k++) {
	int iRow = primalAindex[k];
	int iPut = iwork[iRow]++;
	ARindex[iPut] = iCol;
	ARvalue[iPut] = primalAvalue[k];
      }
    }
    
    // Transpose the problem!
    std::swap(primalNumRow, primalNumCol);
    dual_lp.Astart_.swap(ARstart);
    dual_lp.Aindex_.swap(ARindex);
    dual_lp.Avalue_.swap(ARvalue);
    dual_lp.colLower_.swap(dualColLower);
    dual_lp.colUpper_.swap(dualColUpper);
    dual_lp.rowLower_.swap(dualRowLower);
    dual_lp.rowUpper_.swap(dualRowUpper);
    dual_lp.colCost_.swap(dualCost);
    //    cout << "problem-transposed" << endl;
    // Deduce the consequences of transposing the LP
    //  mlFg_Update(mlFg_action_TransposeLP);
    simplex_info_.transposed_solver_lp = true;
  }
  
  // SCALING:
  // Limits on scaling factors
  const double minAlwScale = 1 / 1024.0;
  const double maxAlwScale = 1024.0;
  const double maxAlwCostScale = maxAlwScale;
  const double minAlwColScale = minAlwScale;
  const double maxAlwColScale = maxAlwScale;
  const double minAlwRowScale = minAlwScale;
  const double maxAlwRowScale = maxAlwScale;
  
#ifdef HiGHSDEV
  // Information on large costs
  const double tlLargeCo = 1e5;
  int numLargeCo;
  vector<int> largeCostFlag;
  double largeCostScale;
#endif
  
  void scaleHighsModelInit(HighsModelObject &highs_model) {
    highs_model.scale_.col_.assign(highs_model.solver_lp_.numCol_, 1);
    highs_model.scale_.row_.assign(highs_model.solver_lp_.numRow_, 1);
    highs_model.scale_.cost_ = 1;
#ifdef HiGHSDEV
    //  largeCostScale = 1;
#endif
  }
  
  void scaleCosts(HighsModelObject &highs_model) {
    // Scale the costs by no less than minAlwCostScale
    double costScale = highs_model.scale_.cost_;
    double maxNzCost = 0;
    for (int iCol = 0; iCol < highs_model.solver_lp_.numCol_; iCol++) {
      if (highs_model.solver_lp_.colCost_[iCol]) {
	maxNzCost = max(fabs(highs_model.solver_lp_.colCost_[iCol]), maxNzCost);
      }
    }
    // Scaling the costs up effectively increases the dual tolerance to
    // which the problem is solved - so, if the max cost is small the
    // scaling factor pushes it up by a power of 2 so it's close to 1
    // Scaling the costs down effectively decreases the dual tolerance
    // to which the problem is solved - so this can't be done too much
    costScale = 1;
    const double ln2 = log(2.0);
    // Scale the costs if the max cost is positive and outside the range [1/16,
    // 16]
    if ((maxNzCost > 0) && ((maxNzCost < (1.0 / 16)) || (maxNzCost > 16))) {
      costScale = maxNzCost;
      costScale = pow(2.0, floor(log(costScale) / ln2 + 0.5));
      costScale = min(costScale, maxAlwCostScale);
    }
#ifdef HiGHSDEV
    HighsPrintMessage(ML_MINIMAL, "MaxNzCost = %11.4g: scaling all costs by %11.4g\ngrep_CostScale,%g,%g\n",
		      maxNzCost, costScale, maxNzCost, costScale);
#endif
    if (costScale == 1) return;
    // Scale the costs (and record of maxNzCost) by costScale, being at most
    // maxAlwCostScale
    for (int iCol = 0; iCol < highs_model.solver_lp_.numCol_; iCol++) {
      highs_model.solver_lp_.colCost_[iCol] /= costScale;
    }
    maxNzCost /= costScale;
    
#ifdef HiGHSDEV
    bool alwLargeCostScaling = false;
    /*
      if (alwLargeCostScaling && (numLargeCo > 0)) {
      // Scale any large costs by largeCostScale, being at most (a further)
      // maxAlwCostScale
      largeCostScale = maxNzCost;
      largeCostScale = pow(2.0, floor(log(largeCostScale) / ln2 + 0.5));
      largeCostScale = min(largeCostScale, maxAlwCostScale);
      printf(
      "   Scaling all |cost| > %11.4g by %11.4g\ngrep_LargeCostScale,%g,%g\n",
      tlLargeCo, largeCostScale, tlLargeCo, largeCostScale);
      for (int iCol = 0; iCol < highs_model.solver_lp_.numCol_; iCol++) {
      if (largeCostFlag[iCol]) {
      highs_model.solver_lp_.colCost_[iCol] /= largeCostScale;
      }
      }
      }
    */
    HighsPrintMessage(ML_MINIMAL, "After cost scaling\n");
    //  utils.util_analyseVectorValues("Column costs", highs_model.solver_lp_.numCol_, highs_model.solver_lp_.colCost_, false);
#endif
  }
  
  void scale_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called scale_solver_lp: simplex_info_.scaled_solver_lp = %d\n", simplex_info_.scaled_solver_lp);
#endif
    if (simplex_info_.scaled_solver_lp) return;
    // Scale the LP highs_model.solver_lp_, assuming all data are in place
    // Reset all scaling to 1
    HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
    HighsTimer &timer = highs_model.timer_;
    timer.start(timer.scaleClock);
    scaleHighsModelInit(highs_model);
    int numCol = highs_model.solver_lp_.numCol_;
    int numRow = highs_model.solver_lp_.numRow_;
    double *colScale = &highs_model.scale_.col_[0];
    double *rowScale = &highs_model.scale_.row_[0];
    int *Astart = &highs_model.solver_lp_.Astart_[0];
    int *Aindex = &highs_model.solver_lp_.Aindex_[0];
    double *Avalue = &highs_model.solver_lp_.Avalue_[0];
    double *colCost = &highs_model.solver_lp_.colCost_[0];
    double *colLower = &highs_model.solver_lp_.colLower_[0];
    double *colUpper = &highs_model.solver_lp_.colUpper_[0];
    double *rowLower = &highs_model.solver_lp_.rowLower_[0];
    double *rowUpper = &highs_model.solver_lp_.rowUpper_[0];
    
    // Allow a switch to/from the original scaling rules
    bool originalScaling = true;
    bool alwCostScaling = true;
    if (originalScaling) alwCostScaling = false;
    
    
    // Find out range of matrix values and skip matrix scaling if all
    // |values| are in [0.2, 5]
    const double inf = HIGHS_CONST_INF;
    double min0 = inf, max0 = 0;
    for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
      double value = fabs(Avalue[k]);
      min0 = min(min0, value);
      max0 = max(max0, value);
    }
    bool noScaling = min0 >= 0.2 && max0 <= 5;
    //   noScaling = false; printf("!!!! FORCE SCALING !!!!\n");
    if (noScaling) {
      // No matrix scaling, but possible cost scaling
#ifdef HiGHSDEV
      HighsPrintMessage(ML_MINIMAL, "grep_Scaling,%s,Obj,0,Row,1,1,Col,1,1,0\n", highs_model.lp_.model_name_.c_str());
#endif
      // Possibly scale the costs
      if (!originalScaling && alwCostScaling) scaleCosts(highs_model);
      timer.stop(timer.scaleClock);
      simplex_info_.scaled_solver_lp = true;
      return;
    }
    // See if we want to include cost include if minimum nonzero cost is less than
    // 0.1
    double minNzCost = inf;
    for (int i = 0; i < numCol; i++) {
      if (colCost[i]) minNzCost = min(fabs(colCost[i]), minNzCost);
    }
    bool includeCost = false;
    //  if (originalScaling)
    includeCost = minNzCost < 0.1;
    
    // Search up to 6 times
    vector<double> rowMin(numRow, inf);
    vector<double> rowMax(numRow, 1 / inf);
    for (int search_count = 0; search_count < 6; search_count++) {
      // Find column scale, prepare row data
      for (int iCol = 0; iCol < numCol; iCol++) {
	// For column scale (find)
	double colMin = inf;
	double colMax = 1 / inf;
	double myCost = fabs(colCost[iCol]);
	if (includeCost && myCost != 0)
	  colMin = min(colMin, myCost), colMax = max(colMax, myCost);
	for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	  double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
	  colMin = min(colMin, value), colMax = max(colMax, value);
	}
	colScale[iCol] = 1 / sqrt(colMin * colMax);
	if (!originalScaling) {
	  // Ensure that column scale factor is not excessively large or small
	  colScale[iCol] =
            min(max(minAlwColScale, colScale[iCol]), maxAlwColScale);
	}
	// For row scale (only collect)
	for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	  int iRow = Aindex[k];
	  double value = fabs(Avalue[k]) * colScale[iCol];
	  rowMin[iRow] = min(rowMin[iRow], value);
	  rowMax[iRow] = max(rowMax[iRow], value);
	}
      }
      
      // For row scale (find)
      for (int iRow = 0; iRow < numRow; iRow++) {
	rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
	if (!originalScaling) {
	  // Ensure that row scale factor is not excessively large or small
	  rowScale[iRow] = min(max(minAlwRowScale, rowScale[iRow]), maxAlwRowScale);
	}
      }
      rowMin.assign(numRow, inf);
      rowMax.assign(numRow, 1 / inf);
    }
    
    // Make it numerical better
    // Also determine the max and min row and column scaling factors
    double minColScale = inf;
    double maxColScale = 1 / inf;
    double minRowScale = inf;
    double maxRowScale = 1 / inf;
    const double ln2 = log(2.0);
    for (int iCol = 0; iCol < numCol; iCol++) {
      colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
      minColScale = min(colScale[iCol], minColScale);
      maxColScale = max(colScale[iCol], maxColScale);
    }
    for (int iRow = 0; iRow < numRow; iRow++) {
      rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));
      minRowScale = min(rowScale[iRow], minRowScale);
      maxRowScale = max(rowScale[iRow], maxRowScale);
    }
#ifdef HiGHSDEV
    bool excessScaling =
      (minColScale < minAlwColScale) || (maxColScale > maxAlwColScale) ||
      (minRowScale < minAlwRowScale) || (maxRowScale > maxAlwRowScale);
    
    HighsPrintMessage(ML_MINIMAL, "grep_Scaling,%s,%d,%d,Obj,%g,%d,Row,%g,%g,Col,%g,%g,%d\n",
		      highs_model.lp_.model_name_.c_str(), originalScaling, alwCostScaling, minNzCost,
		      includeCost, minColScale, maxColScale, minRowScale, maxRowScale,
		      excessScaling);
#endif
    
    // Apply scaling to matrix and bounds
    for (int iCol = 0; iCol < numCol; iCol++)
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
	Avalue[k] *= (colScale[iCol] * rowScale[Aindex[k]]);
    
    for (int iCol = 0; iCol < numCol; iCol++) {
      colLower[iCol] /= colLower[iCol] == -inf ? 1 : colScale[iCol];
      colUpper[iCol] /= colUpper[iCol] == +inf ? 1 : colScale[iCol];
      colCost[iCol] *= colScale[iCol];
    }
    for (int iRow = 0; iRow < numRow; iRow++) {
      rowLower[iRow] *= rowLower[iRow] == -inf ? 1 : rowScale[iRow];
      rowUpper[iRow] *= rowUpper[iRow] == +inf ? 1 : rowScale[iRow];
    }
    // Deduce the consequences of scaling the LP
    //  mlFg_Update(mlFg_action_ScaleLP);
#ifdef HiGHSDEV
    // Analyse the scaled LP
    //  if (simplex_info.analyse_lp) {
    //    util_analyseLp(highs_model.solver_lp_, "Scaled");
    //  }
    //  if (mlFg_scaledLP) {
    //  utils.util_analyseVectorValues("Column scaling factors", numCol, colScale, false);
    //  utils.util_analyseVectorValues("Row scaling factors", numRow, rowScale, false);
    //  }
#endif
    // Possibly scale the costs
    if (!originalScaling && alwCostScaling) scaleCosts(highs_model);
    simplex_info_.scaled_solver_lp = true;
    timer.stop(timer.scaleClock);
  }
  
  // PERMUTE:
  
  void permute_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called permute_solver_lp: simplex_info_.permuted_solver_lp = %d\n", simplex_info_.permuted_solver_lp);
#endif
    if (simplex_info_.permuted_solver_lp) return;
    //  HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
    HSimplex simplex_method_;
    simplex_method_.initialiseSolverLpRandomVectors(highs_model);
    
    int numCol = highs_model.solver_lp_.numCol_;
    vector<int>& numColPermutation = highs_model.simplex_info_.numColPermutation_;
    vector<int>& Astart = highs_model.solver_lp_.Astart_;
    vector<int>& Aindex = highs_model.solver_lp_.Aindex_;
    vector<double>& Avalue = highs_model.solver_lp_.Avalue_;
    vector<double>& colCost = highs_model.solver_lp_.colCost_;
    vector<double>& colLower = highs_model.solver_lp_.colLower_;
    vector<double>& colUpper = highs_model.solver_lp_.colUpper_;
    vector<double>& colScale = highs_model.scale_.col_;
    
    // 2. Duplicate the original data to copy from
    vector<int> saveAstart = highs_model.solver_lp_.Astart_;
    vector<int> saveAindex = highs_model.solver_lp_.Aindex_;
    vector<double> saveAvalue = highs_model.solver_lp_.Avalue_;
    vector<double> saveColCost = highs_model.solver_lp_.colCost_;
    vector<double> saveColLower = highs_model.solver_lp_.colLower_;
    vector<double> saveColUpper = highs_model.solver_lp_.colUpper_;
    vector<double> saveColScale = highs_model.scale_.col_;
    
    // 3. Generate the permuted matrix and corresponding vectors of column data
    int countX = 0;
    for (int i = 0; i < numCol; i++) {
      int fromCol = numColPermutation[i];
      Astart[i] = countX;
      for (int k = saveAstart[fromCol]; k < saveAstart[fromCol + 1]; k++) {
	Aindex[countX] = saveAindex[k];
	Avalue[countX] = saveAvalue[k];
	countX++;
      }
      colCost[i] = saveColCost[fromCol];
      colLower[i] = saveColLower[fromCol];
      colUpper[i] = saveColUpper[fromCol];
      colScale[i] = saveColScale[fromCol];
    }
    assert(Astart[numCol] == countX);
    // Deduce the consequences of shuffling the LP
    //  mlFg_Update(mlFg_action_ShuffleLP);
    simplex_info_.permuted_solver_lp = true;
  }
  
  // TIGHTEN:
  
  void tighten_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called tighten_solver_lp: simplex_info_.tightened_solver_lp = %d\n", simplex_info_.tightened_solver_lp);
#endif
    if (simplex_info_.tightened_solver_lp) return;
    HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
    
    int numCol = highs_model.solver_lp_.numCol_;
    int numRow = highs_model.solver_lp_.numRow_;
    vector<int>& Astart = highs_model.solver_lp_.Astart_;
    vector<int>& Aindex = highs_model.solver_lp_.Aindex_;
    vector<double>& Avalue = highs_model.solver_lp_.Avalue_;
    vector<double>& colCost = highs_model.solver_lp_.colCost_;
    vector<double>& colLower = highs_model.solver_lp_.colLower_;
    vector<double>& colUpper = highs_model.solver_lp_.colUpper_;
    vector<double>& rowLower = highs_model.solver_lp_.rowLower_;
    vector<double>& rowUpper = highs_model.solver_lp_.rowUpper_;
    
    
    vector<int> iwork(numRow, 0);
    vector<int> ARstart(numRow + 1, 0);
    int AcountX = Aindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++) iwork[Aindex[k]]++;
    for (int i = 1; i <= numRow; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < numRow; i++) iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < numCol; iCol++) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	int iRow = Aindex[k];
	int iPut = iwork[iRow]++;
	ARindex[iPut] = iCol;
	ARvalue[iPut] = Avalue[k];
      }
    }
    
    // Save column bounds
    vector<double> colLower_0 = highs_model.solver_lp_.colLower_;
    vector<double> colUpper_0 = highs_model.solver_lp_.colUpper_;
    
    double big_B = 1e10;
    int iPass = 0;
    for (;;) {
      int numberChanged = 0;
      for (int iRow = 0; iRow < numRow; iRow++) {
	// SKIP free rows
	if (rowLower[iRow] < -big_B && rowUpper[iRow] > big_B) continue;
	
	// possible row
	int ninfU = 0;
	int ninfL = 0;
	double xmaxU = 0.0;
	double xminL = 0.0;
	int myStart = ARstart[iRow];
	int myEnd = ARstart[iRow + 1];
	// Compute possible lower and upper ranges
	
	for (int k = myStart; k < myEnd; ++k) {
	  int iCol = ARindex[k];
	  double value = ARvalue[k];
	  double upper = value > 0 ? colUpper[iCol] : -colLower[iCol];
	  double lower = value > 0 ? colLower[iCol] : -colUpper[iCol];
	  value = fabs(value);
	  if (upper < big_B)
	    xmaxU += upper * value;
	  else
	    ++ninfU;
	  if (lower > -big_B)
	    xminL += lower * value;
	  else
	    ++ninfL;
	}
	
	// Build in a margin of error
	xmaxU += 1.0e-8 * fabs(xmaxU);
	xminL -= 1.0e-8 * fabs(xminL);
	
	double xminLmargin = (fabs(xminL) > 1.0e8) ? 1e-12 * fabs(xminL) : 0;
	double xmaxUmargin = (fabs(xmaxU) > 1.0e8) ? 1e-12 * fabs(xmaxU) : 0;
	
	// Skip redundant row : also need to consider U < L  case
	double comp_U = xmaxU + ninfU * 1.0e31;
	double comp_L = xminL - ninfL * 1.0e31;
	if (comp_U <= rowUpper[iRow] + 1e-7 && comp_L >= rowLower[iRow] - 1e-7)
	  continue;
	
	double row_L = rowLower[iRow];
	double row_U = rowUpper[iRow];
	
	// Now see if we can tighten column bounds
	for (int k = myStart; k < myEnd; ++k) {
	  double value = ARvalue[k];
	  int iCol = ARindex[k];
	  double col_L = colLower[iCol];
	  double col_U = colUpper[iCol];
	  double new_L = -HIGHS_CONST_INF;
	  double new_U = +HIGHS_CONST_INF;
	  
	  if (value > 0.0) {
	    if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_U > +big_B))
	      new_L = (row_L - xmaxU) / value + (1 - ninfU) * col_U - xmaxUmargin;
	    if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_L < -big_B))
	      new_U = (row_U - xminL) / value + (1 - ninfL) * col_L + xminLmargin;
	  } else {
	    if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_L < -big_B))
	      new_U = (row_L - xmaxU) / value + (1 - ninfU) * col_L + xmaxUmargin;
	    if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_U > +big_B))
	      new_L = (row_U - xminL) / value + (1 - ninfL) * col_U - xminLmargin;
	  }
	  
	  if (new_U < col_U - 1.0e-12 && new_U < big_B) {
	    colUpper[iCol] = max(new_U, col_L);
	    numberChanged++;
	  }
	  if (new_L > col_L + 1.0e-12 && new_L > -big_B) {
	    colLower[iCol] = min(new_L, col_U);
	    numberChanged++;
	  }
	}
      }
      
      if (numberChanged == 0) break;
      iPass++;
      if (iPass > 10) break;
    }
    
    double useTolerance = 1.0e-3;
    for (int iCol = 0; iCol < numCol; iCol++) {
      if (colUpper_0[iCol] > colLower_0[iCol] + useTolerance) {
	const double relax = 100.0 * useTolerance;
	if (colUpper[iCol] - colLower[iCol] < useTolerance + 1.0e-8) {
	  colLower[iCol] = max(colLower_0[iCol], colLower[iCol] - relax);
	  colUpper[iCol] = min(colUpper_0[iCol], colUpper[iCol] + relax);
	} else {
	  if (colUpper[iCol] < colUpper_0[iCol]) {
	    colUpper[iCol] = min(colUpper[iCol] + relax, colUpper_0[iCol]);
	  }
	  if (colLower[iCol] > colLower_0[iCol]) {
	    colLower[iCol] = min(colLower[iCol] - relax, colLower_0[iCol]);
	  }
	}
      }
    }
    simplex_info_.tightened_solver_lp = true;
  }
  
};
#endif // SIMPLEX_HSIMPLEX_H_
