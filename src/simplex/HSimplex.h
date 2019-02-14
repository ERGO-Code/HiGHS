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
#include "HighsUtils.h"
#include "HighsLpUtils.h"
#include "HighsModelObject.h"
#include "SimplexTimer.h"

#include <cassert>
#include <vector>
#include <cstring> // For strcmp

/**
 * @brief Class for simplex utilities
 */
class HSimplex {
 public:
  
  /*
  // Increment iteration count (here!) and (possibly) store the pivots for
  // debugging NLA
  void record_pivots(int columnIn, int columnOut, double alpha) {
    // NB This is where the iteration count is updated!
    if (columnIn >= 0) simplex_info_.iteration_count++;
#ifdef HiGHSDEV
    historyColumnIn.push_back(columnIn);
    historyColumnOut.push_back(columnOut);
    historyAlpha.push_back(alpha);
#endif
  }
#ifdef HiGHSDEV
  // Store and write out the pivots for debugging NLA
  void writePivots(const char* suffix) {
    string filename = "z-" + solver_lp_->model_name_ + "-" + suffix;
    ofstream output(filename.c_str());
    int count = historyColumnIn.size();
    double current_run_highs_time = timer_->readRunHighsClock();
    output << solver_lp_->model_name_ << " " << count << "\t" << current_run_highs_time << endl;
    output << setprecision(12);
    for (int i = 0; i < count; i++) {
      output << historyColumnIn[i] << "\t";
      output << historyColumnOut[i] << "\t";
      output << historyAlpha[i] << endl;
    }
    output.close();
  }
#endif
  */
  void clear_solver_lp_data(
	     HighsModelObject & highs_model_object //!< Model object in which data for LP to be solved is to be cleared
	     ) {
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
    simplex_info_.solver_lp_has_matrix_col_wise = false;
    simplex_info_.solver_lp_has_matrix_row_wise = false;
    simplex_info_.solver_lp_has_dual_steepest_edge_weights = false;
    simplex_info_.solver_lp_has_nonbasic_dual_values = false;
    simplex_info_.solver_lp_has_basic_primal_values = false;
    simplex_info_.solver_lp_has_invert = false;
    simplex_info_.solver_lp_has_fresh_invert = false;
    simplex_info_.solver_lp_has_fresh_rebuild = false;
    simplex_info_.solver_lp_has_dual_objective_value = false;
  }
  
  void clear_solver_lp(
	     HighsModelObject & highs_model_object //!< Model object in which LP to be solved is to be cleared
	     ) {
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
    // Once the solver LP has its own basis
    //    highs_model_object.solver_basis_.valid_ = false;
    simplex_info_.solver_lp_is_transposed = false;
    simplex_info_.solver_lp_is_scaled = false;
    simplex_info_.solver_lp_is_permuted = false;
    simplex_info_.solver_lp_is_tightened = false;
    clear_solver_lp_data(highs_model_object);
  }
  
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
    simplex_info_.update_limit = opt.simplex_update_limit;
    simplex_info_.highs_run_time_limit = opt.highs_run_time_limit;
    
    simplex_info_.transpose_solver_lp = opt.transpose_solver_lp;
    simplex_info_.scale_solver_lp = opt.scale_solver_lp;
    simplex_info_.permute_solver_lp = opt.permute_solver_lp;
    simplex_info_.tighten_solver_lp = opt.tighten_solver_lp;
    
    // Set values of internal options
    
    // Options for reporting timing
    simplex_info_.report_simplex_inner_clock = true;//false;
    simplex_info_.report_simplex_outer_clock = false;
#ifdef HiGHSDEV
    simplex_info_.report_simplex_phases_clock = true;//false;
    // Option for analysing simplex iterations
    simplex_info_.analyseLp = true;//false;
    simplex_info_.analyseSimplexIterations = true;//false
    simplex_info_.analyseLpSolution = true;//false;
    simplex_info_.analyse_invert_time = false;
    simplex_info_.analyseRebuildTime = false;
#endif
    
  }
  
  void update_solver_lp_status_flags(
				 HighsModelObject &highs_model_object,
				 LpAction action
				 ) {
    
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
    switch (action) {
    case LpAction::TRANSPOSE:
#ifdef HIGHSDEV
      printf(" LpAction::TRANSPOSE\n");
#endif
      simplex_info_.solver_lp_is_transposed = true;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::SCALE:
#ifdef HIGHSDEV
      printf(" LpAction::SCALE\n");
#endif
      simplex_info_.solver_lp_is_scaled = true;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::PERMUTE:
#ifdef HIGHSDEV
      printf(" LpAction::PERMUTE\n");
#endif
      simplex_info_.solver_lp_is_permuted = true;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::TIGHTEN:
#ifdef HIGHSDEV
      printf(" LpAction::TIGHTEN\n");
#endif
      simplex_info_.solver_lp_is_tightened = true;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::NEW_COSTS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COSTS\n");
#endif
      //      initCost();
      simplex_info_.solver_lp_has_nonbasic_dual_values = false;
      simplex_info_.solver_lp_has_fresh_rebuild = false;
      simplex_info_.solver_lp_has_dual_objective_value = false;
      break;
    case LpAction::NEW_BOUNDS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BOUNDS\n");
#endif
      //      simplex_info_.solver_lp_ = true;
      //     initBound();
      //     initValue();
      simplex_info_.solver_lp_has_basic_primal_values = false;
      simplex_info_.solver_lp_has_fresh_rebuild = false;
      simplex_info_.solver_lp_has_dual_objective_value = false;
      break;
    case LpAction::NEW_BASIS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BASIS\n");
#endif
      highs_model_object.basis_.valid_ = true;
      //      highs_model_object.solver_basis_.valid_ = false;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::NEW_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COLS\n");
#endif
      highs_model_object.basis_.valid_ = true;
      //      highs_model_object.solver_basis_.valid_ = false;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::NEW_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_ROWS\n");
#endif
      highs_model_object.basis_.valid_ = true;
      //      highs_model_object.solver_basis_.valid_ = false;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::DEL_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_COLS\n");
#endif
      highs_model_object.basis_.valid_ = false;
      //      highs_model_object.solver_basis_.valid_ = false;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::DEL_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS\n");
#endif
      highs_model_object.basis_.valid_ = false;
      //      highs_model_object.solver_basis_.valid_ = false;
      clear_solver_lp_data(highs_model_object);
      break;
    case LpAction::DEL_ROWS_BASIS_OK:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS_BASIS_OK\n");
#endif
      //      simplex_info_.solver_lp_ = true;
      break;
    default:
#ifdef HIGHSDEV
      printf(" Unrecognised LpAction::%d\n", (int) action);
#endif
      break;
    }
  }

void report_basis(HighsModelObject &highs_model_object) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsBasis &basis = highs_model_object.basis_;
  
  printf("\nReporting current basis: solver_lp.numCol_ = %d; solver_lp.numRow_ = %d\n", solver_lp.numCol_,
         solver_lp.numRow_);
  if (solver_lp.numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < solver_lp.numCol_; col++) {
    int var = col;
    if (basis.nonbasicFlag_[var])
      printf("%6d %6d        %6d %6d\n", var, col, basis.nonbasicFlag_[var],
             basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d\n", var, col, basis.nonbasicFlag_[var]);
  }
  if (solver_lp.numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < solver_lp.numRow_; row++) {
    int var = solver_lp.numCol_ + row;
    if (basis.nonbasicFlag_[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basis.basicIndex_[row],
             basis.nonbasicFlag_[var], basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basis.basicIndex_[row], basis.nonbasicFlag_[var]);
  }
}


  void report_solver_lp_status_flags(
				 HighsModelObject &highs_model_object
				 ) {
    HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
    printf("\nReporting solver status and flags:\n\n");
    printf("  is_transposed =                  %d\n", simplex_info_.solver_lp_is_transposed);
    printf("  is_scaled =                      %d\n", simplex_info_.solver_lp_is_scaled);
    printf("  is_permuted =                    %d\n", simplex_info_.solver_lp_is_permuted);
    printf("  is_tightened =                   %d\n", simplex_info_.solver_lp_is_tightened);
    printf("  has_matrix_col_wise =            %d\n", simplex_info_.solver_lp_has_matrix_col_wise);
    printf("  has_matrix_row_wise =            %d\n", simplex_info_.solver_lp_has_matrix_row_wise);
    printf("  has_dual_steepest_edge_weights = %d\n", simplex_info_.solver_lp_has_dual_steepest_edge_weights);
    printf("  has_nonbasic_dual_values =       %d\n", simplex_info_.solver_lp_has_nonbasic_dual_values);
    printf("  has_basic_primal_values =        %d\n", simplex_info_.solver_lp_has_basic_primal_values);
    printf("  has_invert =                     %d\n", simplex_info_.solver_lp_has_invert);
    printf("  has_fresh_invert =               %d\n", simplex_info_.solver_lp_has_fresh_invert);
    printf("  has_fresh_rebuild =              %d\n", simplex_info_.solver_lp_has_fresh_rebuild);
    printf("  has_dual_objective_value =       %d\n", simplex_info_.solver_lp_has_dual_objective_value);
  }
  void compute_dual_objective_value(
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
    simplex_info_.solver_lp_has_dual_objective_value = true;
  }
  
  void initialise_solver_lp_random_vectors(
				       HighsModelObject &highs_model
				       ) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
    const int numCol = highs_model.solver_lp_.numCol_;
    const int numTot = highs_model.solver_lp_.numCol_ + highs_model.solver_lp_.numRow_;
    // Instantiate and (re-)initialise the random number generator
    HighsRandom &random = highs_model.random_;
    random.initialise();
    //
    // Generate a random permutation of the column indices
    simplex_info_.numColPermutation_.resize(numCol);
    int *numColPermutation = &simplex_info_.numColPermutation_[0];
    for (int i = 0; i < numCol; i++) numColPermutation[i] = i;
    for (int i = numCol - 1; i >= 1; i--) {
      int j = random.integer() % (i + 1);
      std::swap(numColPermutation[i], numColPermutation[j]);
    }
    
    // Re-initialise the random number generator and generate the
    // random vectors in the same order as hsol to maintain repeatable
    // performance
    random.initialise();
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
    
  }
  
  // TRANSPOSE:
  
  void transpose_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called transpose_solver_lp: simplex_info_.solver_lp_is_transposed = %d\n", simplex_info_.solver_lp_is_transposed);
#endif
    if (simplex_info_.solver_lp_is_transposed) return;
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
    update_solver_lp_status_flags(highs_model, LpAction::TRANSPOSE);
      //    simplex_info_.solver_lp_is_transposed = true;
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
    printf("Called scale_solver_lp: simplex_info_.solver_lp_is_scaled = %d\n", simplex_info_.solver_lp_is_scaled);
#endif
    if (simplex_info_.solver_lp_is_scaled) return;
    // Scale the LP highs_model.solver_lp_, assuming all data are in place
    // Reset all scaling to 1
    HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
    HighsTimer &timer = highs_model.timer_;
    timer.start(timer.scale_clock);
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
      timer.stop(timer.scale_clock);
      update_solver_lp_status_flags(highs_model, LpAction::SCALE);
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
    update_solver_lp_status_flags(highs_model, LpAction::SCALE);
    // Possibly scale the costs
    if (!originalScaling && alwCostScaling) scaleCosts(highs_model);
    timer.stop(timer.scale_clock);
  }
  
  // PERMUTE:
  
  void permute_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called permute_solver_lp: simplex_info_.solver_lp_is_permuted = %d\n", simplex_info_.solver_lp_is_permuted);
#endif
    if (simplex_info_.solver_lp_is_permuted) return;
    //  HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
    HSimplex simplex_method_;
    simplex_method_.initialise_solver_lp_random_vectors(highs_model);
    
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
    // Deduce the consequences of permuting the LP
    update_solver_lp_status_flags(highs_model, LpAction::PERMUTE); 
  }
  
  // TIGHTEN:
  
  void tighten_solver_lp(HighsModelObject &highs_model) {
    HighsSimplexInfo &simplex_info_ = highs_model.simplex_info_;
#ifdef HiGHSDEV
    printf("Called tighten_solver_lp: simplex_info_.solver_lp_is_tightened = %d\n", simplex_info_.solver_lp_is_tightened);
#endif
    if (simplex_info_.solver_lp_is_tightened) return;
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
    simplex_info_.solver_lp_is_tightened = true;
  }

  void initialise_basic_index(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;

    int num_basic_variables = 0;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int var = 0; var < numTot; var++) {
      if (!basis.nonbasicFlag_[var]) {
	assert(num_basic_variables < solver_lp.numRow_);
	basis.basicIndex_[num_basic_variables] = var;
	num_basic_variables++;
      }
    }
    assert(num_basic_variables = solver_lp.numRow_ - 1);
  }

  void allocate_work_and_base_arrays(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Allocate bounds and solution spaces
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    simplex_info.workCost_.resize(numTot);
    simplex_info.workDual_.resize(numTot);
    simplex_info.workShift_.resize(numTot);
    
    simplex_info.workLower_.resize(numTot);
    simplex_info.workUpper_.resize(numTot);
    simplex_info.workRange_.resize(numTot);
    simplex_info.workValue_.resize(numTot);
    
    simplex_info.baseLower_.resize(solver_lp.numRow_);
    simplex_info.baseUpper_.resize(solver_lp.numRow_);
    simplex_info.baseValue_.resize(solver_lp.numRow_);
  }

  void initialise_from_nonbasic(HighsModelObject &highs_model_object) {
    // Initialise basicIndex from nonbasic* then allocate and populate
    // (where possible) work* arrays and allocate basis* arrays
    initialise_basic_index(highs_model_object);
    allocate_work_and_base_arrays(highs_model_object);
    populate_work_arrays(highs_model_object);
    
    // Deduce the consequences of a new basis
    update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  }
  
  void replace_from_nonbasic(HighsModelObject &highs_model_object) {
    // Initialise basicIndex using nonbasic* then populate (where possible)
    // work* arrays
    initialise_basic_index(highs_model_object);
    populate_work_arrays(highs_model_object);
    
    // Deduce the consequences of a new basis
    update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  }
  
  void initialise_with_logical_basis(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Initialise with a logical basis then allocate and populate (where
    // possible) work* arrays and allocate basis* arrays
    
    for (int row = 0; row < solver_lp.numRow_; row++) basis.basicIndex_[row] = solver_lp.numCol_ + row;
    for (int col = 0; col < solver_lp.numCol_; col++) basis.nonbasicFlag_[col] = 1;
    simplex_info.num_basic_logicals = solver_lp.numRow_;
    
    allocate_work_and_base_arrays(highs_model_object);
    populate_work_arrays(highs_model_object);
    
    // Deduce the consequences of a new basis
    update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  }

  void initialise_value_from_nonbasic(HighsModelObject &highs_model_object, int firstvar, int lastvar) {
    // Initialise workValue and nonbasicMove from nonbasicFlag and
    // bounds, except for boxed variables when nonbasicMove is used to
    // set workValue=workLower/workUpper
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    assert(firstvar >= 0);
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    assert(lastvar < numTot);
    // double dl_pr_act, norm_dl_pr_act;
    // norm_dl_pr_act = 0.0;
    for (int var = firstvar; var <= lastvar; var++) {
      if (basis.nonbasicFlag_[var]) {
	// Nonbasic variable
	// double prev_pr_act = simplex_info.workValue_[var];
	if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
	  // Fixed
	  simplex_info.workValue_[var] = simplex_info.workLower_[var];
	  basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
	} else if (!highs_isInfinity(-simplex_info.workLower_[var])) {
	  // Finite lower bound so boxed or lower
	  if (!highs_isInfinity(simplex_info.workUpper_[var])) {
	    // Finite upper bound so boxed
	    if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
	      // Set at lower
	      simplex_info.workValue_[var] = simplex_info.workLower_[var];
	    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
	      // Set at upper
	      simplex_info.workValue_[var] = simplex_info.workUpper_[var];
	    } else {
	      // Invalid nonbasicMove: correct and set value at lower
	      basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
	      simplex_info.workValue_[var] = simplex_info.workLower_[var];
	    }
	  } else {
	    // Lower
	    simplex_info.workValue_[var] = simplex_info.workLower_[var];
	    basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
	  }
	} else if (!highs_isInfinity(simplex_info.workUpper_[var])) {
	  // Upper
	  simplex_info.workValue_[var] = simplex_info.workUpper_[var];
	  basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
	} else {
	  // FREE
	  simplex_info.workValue_[var] = 0;
	  basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
	}
	// dl_pr_act = simplex_info.workValue_[var] - prev_pr_act;
	// norm_dl_pr_act += dl_pr_act*dl_pr_act;
	//      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
	//      %8g; %8g] Du = %8g; DlPr = %8g\n",
	//					var, simplex_info.workLower_[var],
	// simplex_info.workValue_[var], simplex_info.workUpper_[var], simplex_info.workDual_[var], dl_pr_act);
      } else {
	// Basic variable
	basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      }
    }
    //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
    //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is
    //  %g\n", norm_dl_pr_act);
  }

  void initialise_value(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    initialise_value_from_nonbasic(highs_model_object, 0, numTot - 1);
  }

  void initialise_phase2_col_bound(HighsModelObject &highs_model_object, int firstcol, int lastcol) {
    // Copy bounds and compute ranges
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    assert(firstcol >= 0);
    assert(lastcol < solver_lp.numCol_);
    for (int col = firstcol; col <= lastcol; col++) {
      simplex_info.workLower_[col] = solver_lp.colLower_[col];
      simplex_info.workUpper_[col] = solver_lp.colUpper_[col];
      simplex_info.workRange_[col] = simplex_info.workUpper_[col] - simplex_info.workLower_[col];
    }
  }

  void initialise_phase2_row_bound(HighsModelObject &highs_model_object, int firstrow, int lastrow) {
    // Copy bounds and compute ranges
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    assert(firstrow >= 0);
    assert(lastrow < solver_lp.numRow_);
    for (int row = firstrow; row <= lastrow; row++) {
      int var = solver_lp.numCol_ + row;
      simplex_info.workLower_[var] = -solver_lp.rowUpper_[row];
      simplex_info.workUpper_[var] = -solver_lp.rowLower_[row];
      simplex_info.workRange_[var] = simplex_info.workUpper_[var] - simplex_info.workLower_[var];
    }
  }

  void initialise_bound(HighsModelObject &highs_model_object, int phase = 2) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
    // necessary to compute Phase 1 bounds
    initialise_phase2_col_bound(highs_model_object, 0, solver_lp.numCol_ - 1);
    initialise_phase2_row_bound(highs_model_object, 0, solver_lp.numRow_ - 1);
    if (phase == 2) return;

    // In Phase 1: change to dual phase 1 bound
    const double inf = HIGHS_CONST_INF;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) {
      if (simplex_info.workLower_[i] == -inf && simplex_info.workUpper_[i] == inf) {
	// Won't change for row variables: they should never become
	// non basic
	if (i >= solver_lp.numCol_) continue;
	simplex_info.workLower_[i] = -1000, simplex_info.workUpper_[i] = 1000;  // FREE
      } else if (simplex_info.workLower_[i] == -inf) {
	simplex_info.workLower_[i] = -1, simplex_info.workUpper_[i] = 0;  // UPPER
      } else if (simplex_info.workUpper_[i] == inf) {
	simplex_info.workLower_[i] = 0, simplex_info.workUpper_[i] = 1;  // LOWER
      } else {
	simplex_info.workLower_[i] = 0, simplex_info.workUpper_[i] = 0;  // BOXED or FIXED
      }
      simplex_info.workRange_[i] = simplex_info.workUpper_[i] - simplex_info.workLower_[i];
    }
  }

  void initialise_phase2_col_cost(HighsModelObject &highs_model_object, int firstcol, int lastcol) {
    // Copy the Phase 2 cost and zero the shift
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    for (int col = firstcol; col <= lastcol; col++) {
      int var = col;
      simplex_info.workCost_[var] = solver_lp.sense_ * solver_lp.colCost_[col];
      simplex_info.workShift_[var] = 0.;
    }
  }
  
  void initialise_phase2_row_cost(HighsModelObject &highs_model_object, int firstrow, int lastrow) {
    // Zero the cost and shift
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    for (int row = firstrow; row <= lastrow; row++) {
      int var = solver_lp.numCol_ + row;
      simplex_info.workCost_[var] = 0;
      simplex_info.workShift_[var] = 0.;
    }
  }

  void initialise_cost(HighsModelObject &highs_model_object, int perturb = 0) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Copy the cost
    initialise_phase2_col_cost(highs_model_object, 0, solver_lp.numCol_ - 1);
    initialise_phase2_row_cost(highs_model_object, 0, solver_lp.numRow_ - 1);
    // See if we want to skip perturbation
    simplex_info.costs_perturbed = 0;
    if (perturb == 0 || simplex_info.perturb_costs == 0) return;
    simplex_info.costs_perturbed = 1;

    // Perturb the original costs, scale down if is too big
    double bigc = 0;
    for (int i = 0; i < solver_lp.numCol_; i++) bigc = max(bigc, fabs(simplex_info.workCost_[i]));
    if (bigc > 100) bigc = sqrt(sqrt(bigc));

    // If there's few boxed variables, we will just use Simple perturbation
    double boxedRate = 0;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) boxedRate += (simplex_info.workRange_[i] < 1e30);
    boxedRate /= numTot;
    if (boxedRate < 0.01) bigc = min(bigc, 1.0);
    if (bigc < 1) {
      //        bigc = sqrt(bigc);
    }

    // Determine the perturbation base
    double base = 5e-7 * bigc;

    // Now do the perturbation
    for (int i = 0; i < solver_lp.numCol_; i++) {
      double lower = solver_lp.colLower_[i];
      double upper = solver_lp.colUpper_[i];
      double xpert = (fabs(simplex_info.workCost_[i]) + 1) * base * (1 + simplex_info.numTotRandomValue_[i]);
      if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF) {
	// Free - no perturb
      } else if (upper == HIGHS_CONST_INF) {  // Lower
	simplex_info.workCost_[i] += xpert;
      } else if (lower == -HIGHS_CONST_INF) {  // Upper
	simplex_info.workCost_[i] += -xpert;
      } else if (lower != upper) {  // Boxed
	simplex_info.workCost_[i] += (simplex_info.workCost_[i] >= 0) ? xpert : -xpert;
      } else {
	// Fixed - no perturb
      }
    }
    
    for (int i = solver_lp.numCol_; i < numTot; i++) {
      simplex_info.workCost_[i] += (0.5 - simplex_info.numTotRandomValue_[i]) * 1e-12;
    }
  }

  // Get the nonbasicMove value for a particular variable - may not be used
  int get_one_nonbasicMove(HighsModelObject &highs_model_object, int var) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    //  printf("Calling get_one_nonbasicMove with var = %2d; numTot = %2d\n", var, numTot); 
    assert(var >= 0);
    assert(var < numTot);
    if (!highs_isInfinity(-simplex_info.workLower_[var])) {
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
	// Finite lower and upper bounds so nonbasic move depends on whether they
	// are equal
	if (simplex_info.workLower_[var] == simplex_info.workUpper_[var])
	  // Fixed variable so nonbasic move is zero
	  return NONBASIC_MOVE_ZE;
	// Boxed variable so nonbasic move is up (from lower bound)
	return NONBASIC_MOVE_UP;
      } else
	// Finite lower bound and infinite upper bound so nonbasic move is up
	// (from lower bound)
	return NONBASIC_MOVE_UP;
    } else
      // Infinite lower bound so nonbasic move depends on whether the upper
      // bound is finite
      if (!highs_isInfinity(simplex_info.workUpper_[var]))
	// Finite upper bound so nonbasic move is down (from upper bound)
	return NONBASIC_MOVE_DN;
    // Infinite upper bound so free variable: nonbasic move is zero
    return NONBASIC_MOVE_ZE;
  }

  void populate_work_arrays(HighsModelObject &highs_model_object) {
    // Initialize the values
    initialise_cost(highs_model_object);
    initialise_bound(highs_model_object);
    initialise_value(highs_model_object);
  }

  void replace_with_logical_basis(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Replace basis with a logical basis then populate (where possible)
    // work* arrays
    for (int row = 0; row < solver_lp.numRow_; row++) {
      int var = solver_lp.numCol_ + row;
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.basicIndex_[row] = var;
    }
    for (int col = 0; col < solver_lp.numCol_; col++) {
      basis.nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
    }
    simplex_info.num_basic_logicals = solver_lp.numRow_;
    
    populate_work_arrays(highs_model_object);

    // Deduce the consequences of a new basis
    update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
   
  }

  void replace_with_new_basis(HighsModelObject &highs_model_object, const int *XbasicIndex) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    // Replace basis with a new basis then populate (where possible)
    // work* arrays
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int var = 0; var < numTot; var++) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    }
    simplex_info.num_basic_logicals = 0;
    for (int row = 0; row < solver_lp.numRow_; row++) {
      int var = XbasicIndex[row];
      if (var >= solver_lp.numCol_) simplex_info.num_basic_logicals++;
      basis.basicIndex_[row] = var;
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    }
    
    populate_work_arrays(highs_model_object);

    // Deduce the consequences of a new basis
    update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  }

  void extend_with_logical_basis(HighsModelObject &highs_model_object, 
				 int firstcol, int lastcol, int firstrow, int lastrow) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  // Add nonbasic structurals and basic slacks according to model bounds.
  //
  // NB Assumes that the basis data structures and work vectors on
  // entry are assigned for columns 0..firstcol-1 and rows
  // 0..firstrow-1 and that they constitute a valid basis. Thus they
  // correspond to "firstcol" number of columns and "firstrow" number
  // of rows. Also assumes that solver_lp_->numCol_ and solver_lp_->numRow_ have already been
  // updated to correspond to any additional columns and rows. This is
  // necessary so that generic methods can be used to assign model
  // data to arrays dimensioned 0..numTot
  //
  // Null intervals firstcol...lastcol and firstrow...lastrow are
  // permitted, but this is achieved by setting the "last" to be less
  // than "first" since the latter is used to indicate what's
  // currently in the data structure.

  assert(firstcol >= 0);
  assert(firstrow >= 0);

  // printf("Called extendWithLogicalBasis:\n   solver_lp.numCol_ =   %d\n   firstcol =
  // %d\n   lastcol =  %d\n   solver_lp.numRow_ =   %d\n   firstrow = %d\n   lastrow =
  // %d\n", solver_lp.numCol_, firstcol, lastcol, solver_lp.numRow_, firstrow, lastrow);
  // Determine the numbers of columns and rows to be added

  int numAddCol = max(lastcol - firstcol + 1, 0);
  int numAddRow = max(lastrow - firstrow + 1, 0);
  int numAddTot = numAddCol + numAddRow;
  if (numAddTot == 0) return;

  // Determine the numbers of columns and rows before and after this method

  int local_oldNumCol = firstcol;
  int local_oldNumRow = firstrow;
  int local_oldNumTot = local_oldNumCol + local_oldNumRow;

  int local_newNumCol = max(local_oldNumCol, lastcol + 1);
  int local_newNumRow = max(local_oldNumRow, lastrow + 1);
  int local_newNumTot = local_newNumCol + local_newNumRow;

  const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
#ifdef SCIPDEV
  printf("extendWithLogicalBasis\n");
  printf("solver_lp.numCol_/Row/Tot = %d/%d/%d\n", solver_lp.numCol_, solver_lp.numRow_, numTot);
  printf("local_newNumCol/Row/Tot = %d/%d/%d\n", local_newNumCol,
         local_newNumRow, local_newNumTot);
  cout << flush;
#endif
  // ToDo: Replace references to local_newNum* by references to num* from here
  // on
  assert(local_newNumCol == solver_lp.numCol_);
  assert(local_newNumRow == solver_lp.numRow_);
  assert(local_newNumTot == numTot);

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  bool basisOK = nonbasic_flag_basic_index_ok(highs_model_object, local_oldNumCol, local_oldNumRow);
  if (!basisOK)
    printf("extend_with_logical_basis: basisOK = %d\n", basisOK);
  assert(basisOK);
#endif

  //  Resize if necessary

  if (solver_lp.numRow_ > local_oldNumRow) {
    basis.basicIndex_.resize(solver_lp.numRow_);

    simplex_info.baseLower_.resize(solver_lp.numRow_);
    simplex_info.baseUpper_.resize(solver_lp.numRow_);
    simplex_info.baseValue_.resize(solver_lp.numRow_);
  }
  if (numTot > local_oldNumTot) {
    basis.nonbasicFlag_.resize(numTot);
    basis.nonbasicMove_.resize(numTot);

    simplex_info.workCost_.resize(numTot);
    simplex_info.workDual_.resize(numTot);
    simplex_info.workShift_.resize(numTot);

    simplex_info.workLower_.resize(numTot);
    simplex_info.workUpper_.resize(numTot);
    simplex_info.workRange_.resize(numTot);
    simplex_info.workValue_.resize(numTot);
  }

  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary

  int rowShift = solver_lp.numCol_ - local_oldNumCol;
  if (rowShift > 0) {
    // printf("Shifting row data by %d using row=%d..0\n", rowShift,
    // local_oldNumRow-1);cout << flush;
    for (int row = local_oldNumRow - 1; row >= 0; row--) {
      basis.basicIndex_[row] += rowShift;
      basis.nonbasicFlag_[solver_lp.numCol_ + row] = basis.nonbasicFlag_[local_oldNumCol + row];
      basis.nonbasicMove_[solver_lp.numCol_ + row] = basis.nonbasicMove_[local_oldNumCol + row];

      simplex_info.workCost_[solver_lp.numCol_ + row] = simplex_info.workCost_[local_oldNumCol + row];
      simplex_info.workDual_[solver_lp.numCol_ + row] = simplex_info.workDual_[local_oldNumCol + row];
      simplex_info.workShift_[solver_lp.numCol_ + row] = simplex_info.workShift_[local_oldNumCol + row];

      simplex_info.workLower_[solver_lp.numCol_ + row] = simplex_info.workLower_[local_oldNumCol + row];
      simplex_info.workUpper_[solver_lp.numCol_ + row] = simplex_info.workUpper_[local_oldNumCol + row];
      simplex_info.workRange_[solver_lp.numCol_ + row] = simplex_info.workRange_[local_oldNumCol + row];
      simplex_info.workValue_[solver_lp.numCol_ + row] = simplex_info.workValue_[local_oldNumCol + row];

      // printf("Setting basicIndex[%2d] = %2d; basis.nonbasicFlag_[%2d] = %2d;
      // basis.nonbasicMove_[%2d] = %2d\n",
      //      row, basicIndex[row],
      //      solver_lp.numCol_+row, basis.nonbasicFlag_[local_oldNumCol+row],
      //      solver_lp.numCol_+row, basis.nonbasicMove_[local_oldNumCol+row]);cout << flush;
    }
  }
  // report_basis(highs_model_object);
  // printf("After possibly shifting row data\n");
  // Make any new columns nonbasic
  //  printf("Make any new cols nonbasic: %d %d %d\n", solver_lp.numCol_, firstcol,
  //  lastcol);
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    //    printf("Setting basis.nonbasicFlag_[%2d] = NONBASIC_FLAG_TRUE\n", var);
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    //    printf("Calling get_one_nonbasicMove(%2d)\n", var);
    //    basis.nonbasicMove_[var] = get_one_nonbasicMove(highs_model_object, var);
  }
  // Make any new rows basic
  //  printf("Make any new rows basic: %d %d %d\n", solver_lp.numRow_, firstrow, lastrow);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = solver_lp.numCol_ + row;
    //    printf("Setting basis.nonbasicFlag_[%2d] = NONBASIC_FLAG_FALSE; Setting
    //    basicIndex[%2d] = %2d\n", var, row, var);
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    basis.basicIndex_[row] = var;
  }

  // Initialise costs for the new columns and rows
  printf("init_Phase2_col_cost(firstcol, lastcol);\n");
  printf("init_Phase2_row_cost(firstrow, lastrow);\n");

  // Initialise bounds for the new columns and rows
  printf("init_Phase2_col_bound(firstcol, lastcol);\n");
  printf("init_Phase2_row_bound(firstrow, lastrow);\n");

  // Initialise values (and nonbasicMove) for the new columns
  printf("Call init_value_from_nonbasic(firstcol, lastcol);\n");

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  basisOK = nonbasic_flag_basic_index_ok(highs_model_object, solver_lp.numCol_, solver_lp.numRow_);
  assert(basisOK);
#endif

  simplex_info.num_basic_logicals += numAddRow;

  //  report_basis(highs_model_object);

  // Deduce the consequences of adding new columns and/or rows
  if (numAddCol) update_solver_lp_status_flags(highs_model_object, LpAction::NEW_COLS);
  if (numAddRow) update_solver_lp_status_flags(highs_model_object, LpAction::NEW_ROWS);
  }

  void setup_num_basic_logicals(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    simplex_info.num_basic_logicals = 0;
    for (int i = 0; i < solver_lp.numRow_; i++)
      if (basis.basicIndex_[i] >= solver_lp.numCol_)
	simplex_info.num_basic_logicals += 1;
#ifdef HiGHSDEV
    printf("Determined num_basic_logicals = %d of %d\n", simplex_info.num_basic_logicals, solver_lp.numRow_);
#endif
  }

  void setup_for_solve(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    int solver_num_row = solver_lp.numRow_;
    int solver_num_col = solver_lp.numCol_;
    if (solver_num_row == 0) return;

    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HMatrix &matrix = highs_model_object.matrix_;
    HFactor &factor = highs_model_object.factor_;
#ifdef HiGHSDEV
    report_solver_lp_status_flags(highs_model_object);
#endif
    bool basis_valid = highs_model_object.basis_.valid_;
#ifdef HiGHSDEV
    printf("In setup_for_solve: basis_valid = %d \n", basis_valid);
#endif
    if (basis_valid) {
    // Model has a basis so just count the number of basic logicals
      setup_num_basic_logicals(highs_model_object);
    } else {
      // Model has no basis: set up a logical basis then populate (where
      // possible) work* arrays
      replace_with_logical_basis(highs_model_object);
      printf("Called replaceWithLogicalBasis\n");
  }

  if (!(simplex_info.solver_lp_has_matrix_col_wise && simplex_info.solver_lp_has_matrix_row_wise)) {
    // Make a copy of col-wise matrix for HMatrix and create its row-wise matrix
    if (simplex_info.num_basic_logicals == solver_num_row) {
      matrix.setup_lgBs(solver_num_col, solver_num_row,
			  &solver_lp.Astart_[0],
			  &solver_lp.Aindex_[0],
			  &solver_lp.Avalue_[0]);
      //      printf("Called matrix_->setup_lgBs\n");cout<<flush;
    } else {
      matrix.setup(solver_num_col, solver_num_row,
		     &solver_lp.Astart_[0],
		     &solver_lp.Aindex_[0],
		     &solver_lp.Avalue_[0],
		     &basis.nonbasicFlag_[0]);
      //      printf("Called matrix_->setup\n");cout<<flush;
    }
    // Indicate that there is a colum-wise and row-wise copy of the
    // matrix: can't be done in matrix_->setup_lgBs
    //    simplex_info.solver_lp_has_matrix_col_wise = true;
    //    simplex_info.solver_lp_has_matrix_row_wise = true;
  }

    // TODO Put something in to skip factor_->setup
    // Initialise factor arrays, passing &basis.basicIndex_[0] so that its
    // address can be copied to the internal Factor pointer
    factor.setup(solver_num_col, solver_num_row,
		   &solver_lp.Astart_[0],
		   &solver_lp.Aindex_[0],
		   &solver_lp.Avalue_[0],
		   &basis.basicIndex_[0]);
    // Indicate that the model has factor arrays: can't be done in factor.setup
    //simplex_info.solver_lp_has_factor_arrays = true;
  }

  bool nonbasic_flag_basic_index_ok(HighsModelObject &highs_model_object, int XnumCol, int XnumRow) {
    HighsBasis &basis = highs_model_object.basis_;
    assert(XnumCol >= 0);
    assert(XnumRow >= 0);
    //  printf("Called nonbasic_flag_basic_index_ok(%d, %d)\n", XnumCol, XnumRow);
    int XnumTot = XnumCol + XnumRow;
    int num_basic_variables = 0;
    if (XnumTot > 0) {
      for (int var = 0; var < XnumTot; var++) if (!basis.nonbasicFlag_[var]) num_basic_variables++;
    }
    assert(num_basic_variables == XnumRow);
    if (num_basic_variables != XnumRow) return false;
    if (XnumRow > 0) {
      for (int row = 0; row < XnumRow; row++) {
	int flag = basis.nonbasicFlag_[basis.basicIndex_[row]];
	assert(!flag);
	if (flag) return false;
      }
    }
    return true;
  }

  bool work_arrays_ok(HighsModelObject &highs_model_object, int phase) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    //  printf("Called work_arrays_ok(%d)\n", phase);cout << flush;
    bool ok = true;
    // Only check phase 2 bounds: others will have been set by solve() so can be
    // trusted
    if (phase == 2) {
      for (int col = 0; col < solver_lp.numCol_; ++col) {
	int var = col;
	if (!highs_isInfinity(-simplex_info.workLower_[var])) {
	  ok = simplex_info.workLower_[var] == solver_lp.colLower_[col];
	  if (!ok) {
	    printf("For col %d, simplex_info.workLower_ should be %g but is %g\n", col,
		   solver_lp.colLower_[col], simplex_info.workLower_[var]);
	    return ok;
	  }
	}
	if (!highs_isInfinity(simplex_info.workUpper_[var])) {
	  ok = simplex_info.workUpper_[var] == solver_lp.colUpper_[col];
	  if (!ok) {
	    printf("For col %d, simplex_info.workUpper_ should be %g but is %g\n", col,
		   solver_lp.colUpper_[col], simplex_info.workUpper_[var]);
	    return ok;
	  }
	}
      }
      for (int row = 0; row < solver_lp.numRow_; ++row) {
	int var = solver_lp.numCol_ + row;
	if (!highs_isInfinity(-simplex_info.workLower_[var])) {
	  ok = simplex_info.workLower_[var] == -solver_lp.rowUpper_[row];
	  if (!ok) {
	    printf("For row %d, simplex_info.workLower_ should be %g but is %g\n", row,
		   -solver_lp.rowUpper_[row], simplex_info.workLower_[var]);
	    return ok;
	  }
	}
	if (!highs_isInfinity(simplex_info.workUpper_[var])) {
	  ok = simplex_info.workUpper_[var] == -solver_lp.rowLower_[row];
	  if (!ok) {
	    printf("For row %d, simplex_info.workUpper_ should be %g but is %g\n", row,
		   -solver_lp.rowLower_[row], simplex_info.workUpper_[var]);
	    return ok;
	  }
	}
      }
    }
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int var = 0; var < numTot; ++var) {
      ok = simplex_info.workRange_[var] == (simplex_info.workUpper_[var] - simplex_info.workLower_[var]);
      if (!ok) {
	printf("For variable %d, simplex_info.workRange_ should be %g = %g - %g but is %g\n",
	       var, simplex_info.workUpper_[var] - simplex_info.workLower_[var], simplex_info.workUpper_[var],
	       simplex_info.workLower_[var], simplex_info.workRange_[var]);
	return ok;
      }
    }
    // Don't check perturbed costs: these will have been set by solve() so can be
    // trusted
    if (!simplex_info.costs_perturbed) {
      for (int col = 0; col < solver_lp.numCol_; ++col) {
	int var = col;
	ok = simplex_info.workCost_[var] == solver_lp.sense_ * solver_lp.colCost_[col];
	if (!ok) {
	  printf("For col %d, simplex_info.workLower_ should be %g but is %g\n", col,
		 solver_lp.colLower_[col], simplex_info.workCost_[var]);
	  return ok;
	}
      }
      for (int row = 0; row < solver_lp.numRow_; ++row) {
	int var = solver_lp.numCol_ + row;
	ok = simplex_info.workCost_[var] == 0.;
	if (!ok) {
	  printf("For row %d, simplex_info.workCost_ should be zero but is %g\n", row,
		 simplex_info.workCost_[var]);
	  return ok;
	}
      }
    }
    // ok must be true if we reach here
    assert(ok);
    return ok;
  }

  bool one_nonbasic_move_vs_work_arrays_ok(HighsModelObject &highs_model_object, int var) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
  const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
  //  printf("Calling oneNonbasicMoveVsWorkArrays_ok with var = %2d; numTot =
  //  %2d\n Bounds [%11g, %11g] nonbasicMove = %d\n",
  //	 var, numTot, simplex_info.workLower_[var], simplex_info.workUpper_[var], basis.nonbasicMove_[var]);
  // cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  // Make sure we're not checking a basic variable
  if (!basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-simplex_info.workLower_[var])) {
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
        // Fixed variable
        ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          printf(
              "Fixed variable %d (solver_lp.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
              "move should be zero but is %d\n",
              var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var], simplex_info.workUpper_[var],
              basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
        if (!ok) {
          printf(
              "Fixed variable %d (solver_lp.numCol_ = %d) so simplex_info.work value should be %g but "
              "is %g\n",
              var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          printf(
              "Boxed variable %d (solver_lp.numCol_ = %d) [%11g, %11g, %11g] range %g so "
              "nonbasic move should be up/down but is  %d\n",
              var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_info.workUpper_[var] - simplex_info.workLower_[var], basis.nonbasicMove_[var]);
          return ok;
        }
        if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (solver_lp.numCol_ = %d) with NONBASIC_MOVE_UP so work "
                "value should be %g but is %g\n",
                var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (solver_lp.numCol_ = %d) with NONBASIC_MOVE_DN so work "
                "value should be %g but is %g\n",
                var, solver_lp.numCol_, simplex_info.workUpper_[var], simplex_info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (solver_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d\n",
            var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var], simplex_info.workUpper_[var],
            NONBASIC_MOVE_UP, basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (solver_lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (solver_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d\n",
            var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var], simplex_info.workUpper_[var],
            basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (solver_lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, solver_lp.numCol_, simplex_info.workUpper_[var], simplex_info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        printf(
            "Free variable %d (solver_lp.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
            "move should be zero but is  %d\n",
            var, solver_lp.numCol_, simplex_info.workLower_[var], simplex_info.workValue_[var], simplex_info.workUpper_[var],
            basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == 0.0;
      if (!ok) {
        printf(
            "Free variable %d (solver_lp.numCol_ = %d) so work value should be zero but "
            "is %g\n",
            var, solver_lp.numCol_, simplex_info.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
    
  }

  bool all_nonbasic_move_vs_work_arrays_ok(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    //    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    bool ok;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int var = 0; var < numTot; ++var) {
      printf("NonbasicMoveVsWorkArrays: var = %2d; basis.nonbasicFlag_[var] = %2d\n",
	     var, basis.nonbasicFlag_[var]);
      if (!basis.nonbasicFlag_[var]) continue;
      ok = one_nonbasic_move_vs_work_arrays_ok(highs_model_object, var);
      if (!ok) {printf("Error in NonbasicMoveVsWorkArrays for nonbasic variable %d\n", var);
	assert(ok);
	return ok;
      }
    }
    // ok must be true if we reach here
    assert(ok);
    return ok;
  }

  bool ok_to_solve(HighsModelObject &highs_model_object, int level, int phase) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    //  printf("Called ok_to_solve(%1d, %1d)\n", level, phase);
    bool ok;
    // Level 0: Minimal check - just look at flags. This means we trust them!
    ok =
      basis.valid_ &&
      simplex_info.solver_lp_has_matrix_col_wise &&
      simplex_info.solver_lp_has_matrix_row_wise &&
      //    simplex_info.solver_lp_has_factor_arrays &&
      simplex_info.solver_lp_has_dual_steepest_edge_weights &&
      simplex_info.solver_lp_has_invert;
    // TODO: Eliminate the following line ASAP!!!
    ok = true;
    if (!ok) {
      if (!basis.valid_)
	printf("Not OK to solve since basis.valid_ = %d\n", basis.valid_);
      if (!simplex_info.solver_lp_has_matrix_col_wise)
	printf("Not OK to solve since simplex_info.solver_lp_has_matrix_col_wise = %d\n",
	       simplex_info.solver_lp_has_matrix_col_wise);
      if (!simplex_info.solver_lp_has_matrix_row_wise)
	printf("Not OK to solve since simplex_info.solver_lp_has_matrix_row_wise = %d\n",
	       simplex_info.solver_lp_has_matrix_row_wise);
      //    if (!simplex_info.solver_lp_has_factor_arrays)
      //      printf("Not OK to solve since simplex_info.solver_lp_has_factor_arrays = %d\n",
      //             simplex_info.solver_lp_has_factor_arrays);
      if (!simplex_info.solver_lp_has_dual_steepest_edge_weights)
	printf("Not OK to solve since simplex_info.solver_lp_has_dual_steepest_edge_weights = %d\n",
	       simplex_info.solver_lp_has_dual_steepest_edge_weights);
      if (!simplex_info.solver_lp_has_invert)
	printf("Not OK to solve since simplex_info.solver_lp_has_invert = %d\n",
	       simplex_info.solver_lp_has_invert); 
    }
    assert(ok);
    if (level <= 0) return ok;
    // Level 1: Basis and data check
    ok = nonbasic_flag_basic_index_ok(highs_model_object, solver_lp.numCol_, solver_lp.numRow_);
    if (!ok) {
      printf("Error in nonbasicFlag and basicIndex\n"); 
      assert(ok);
      return ok;
    }
    ok = work_arrays_ok(highs_model_object, phase);
    if (!ok) {
      printf("Error in workArrays\n"); 
      assert(ok);
      return ok;
    }
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int var = 0; var < numTot; ++var) {
      if (basis.nonbasicFlag_[var]) {
	// Nonbasic variable
	ok = one_nonbasic_move_vs_work_arrays_ok(highs_model_object, var);
	if (!ok) {
	  printf("Error in nonbasicMoveVsWorkArrays for variable %d of %d\n", var, numTot); 
	  assert(ok);
	  return ok;
	}
      }
    }
    if (level <= 1) return ok;
    printf("OKtoSolve(%1d) not implemented\n", level); 
    return ok;
  }


  void flip_bound(HighsModelObject &highs_model_object, int iCol) {
    int *nonbasicMove = &highs_model_object.basis_.nonbasicMove_[0];
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
    simplex_info.workValue_[iCol] = move == 1 ? simplex_info.workLower_[iCol] : simplex_info.workUpper_[iCol];
  }
  /*
  int handle_rank_deficiency(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HFactor &factor = highs_model_object.factor_;
    HighsBasis &basis = highs_model_object.basis_;
    int rankDeficiency = factor.rankDeficiency;
    const int *noPvC = factor.getNoPvC();
    printf("Returned %d = factor.build();\n", rankDeficiency);
    fflush(stdout);
    vector<int> basicRows;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    basicRows.resize(numTot);
    //    printf("Before - basis.basicIndex_:"); for (int iRow=0; iRow<solver_lp.numRow_; iRow++)
    //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
    for (int iRow = 0; iRow < solver_lp.numRow_; iRow++) basicRows[basis.basicIndex_[iRow]] = iRow;
    for (int k = 0; k < rankDeficiency; k++) {
      //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
      //      k, noPvC[k]);fflush(stdout);
      int columnIn = solver_lp.numCol_ + factor.noPvR[k];
      int columnOut = noPvC[k];
      int rowOut = basicRows[columnOut];
      //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
      //      %11.4g]\n", columnIn, columnOut, rowOut, simplex_info.workLower_[columnOut],
      //      simplex_info.workUpper_[columnOut]);
      if (basis.basicIndex_[rowOut] != columnOut) {
	printf("%d = basis.basicIndex_[rowOut] != noPvC[k] = %d\n", basis.basicIndex_[rowOut],
	       columnOut);
	fflush(stdout);
      }
      int sourceOut = setSourceOutFmBd(columnOut);
      updatePivots(columnIn, rowOut, sourceOut);
      updateMatrix(columnIn, columnOut);
    }
    //    printf("After  - basis.basicIndex_:"); for (int iRow=0; iRow<solver_lp.numRow_; iRow++)
    //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
#ifdef HiGHSDEV
    factor.checkInvert();
#endif
    return 0;
  }
  */
  int compute_factor(HighsModelObject &highs_model_object) {
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HMatrix &matrix = highs_model_object.matrix_;
    HFactor &factor = highs_model_object.factor_;
#ifdef HiGHSDEV
    HighsTimer &timer = highs_model_object.timer_;
    double tt0 = 0;
    int iClock = simplex_info.clock_[InvertClock];
    if (simplex_info.analyse_invert_time) tt0 = timer.clock_time[iClock];
#endif
    // TODO Understand why handling noPvC and noPvR in what seem to be
    // different ways ends up equivalent.
    int rankDeficiency = factor.build();
    if (rankDeficiency) {
      //    handle_rank_deficiency();
      //    simplex_info.solution_status = SimplexSolutionStatus::SINGULAR;
#ifdef HiGHSDEV
      //    writePivots("failed");
#endif
      //      return rankDeficiency;
    }
    //    printf("INVERT: After %d iterations and %d updates\n", simplex_info.iteration_count,
    //    simplex_info.update_count);
    simplex_info.update_count = 0;
    
#ifdef HiGHSDEV
    if (simplex_info.analyse_invert_time) {
      int iClock = simplex_info.clock_[InvertClock];
      simplex_info.total_inverts = timer.clock_num_call[iClock];
      simplex_info.total_invert_time = timer.clock_time[iClock];
      double invertTime = simplex_info.total_invert_time - tt0;
      printf(
	     "           INVERT  %4d     on iteration %9d: INVERT  time = %11.4g; "
	     "Total INVERT  time = %11.4g\n",
	     simplex_info.total_inverts,
	     simplex_info.iteration_count, invertTime, simplex_info.total_invert_time);
    }
#endif
    
    // Now have a representation of B^{-1}, and it is fresh!
    simplex_info.solver_lp_has_invert = true;
    simplex_info.solver_lp_has_fresh_invert = true;
    return 0;
  }
  
  void compute_primal(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HMatrix &matrix = highs_model_object.matrix_;
    HFactor &factor = highs_model_object.factor_;
    // Setup a local buffer for the values of basic variables
    HVector buffer;
    buffer.setup(solver_lp.numRow_);
    buffer.clear();
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) {
      if (basis.nonbasicFlag_[i] && simplex_info.workValue_[i] != 0) {
	matrix.collect_aj(buffer, i, simplex_info.workValue_[i]);
      }
    }
    factor.ftran(buffer, 1);
    
    for (int i = 0; i < solver_lp.numRow_; i++) {
      int iCol = basis.basicIndex_[i];
      simplex_info.baseValue_[i] = -buffer.array[i];
      simplex_info.baseLower_[i] = simplex_info.workLower_[iCol];
      simplex_info.baseUpper_[i] = simplex_info.workUpper_[iCol];
    }
    // Now have basic primals
    simplex_info.solver_lp_has_basic_primal_values = true;
  }
  
  void compute_dual(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HMatrix &matrix = highs_model_object.matrix_;
    HFactor &factor = highs_model_object.factor_;
    bool an_compute_dual_norm2 = false;
    double btran_rhs_norm2;
    double btran_sol_norm2;
    double work_dual_norm2;
    
    // Create a local buffer for the pi vector
    HVector buffer;
    buffer.setup(solver_lp.numRow_);
    buffer.clear(); 
    for (int iRow = 0; iRow < solver_lp.numRow_; iRow++) {
      buffer.index[iRow] = iRow;
      buffer.array[iRow] =
        simplex_info.workCost_[basis.basicIndex_[iRow]] + simplex_info.workShift_[basis.basicIndex_[iRow]];
    }
    buffer.count = solver_lp.numRow_;
    if (an_compute_dual_norm2) {
      btran_rhs_norm2 = buffer.norm2();
      btran_rhs_norm2 = sqrt(btran_rhs_norm2);
    }
    //  printf("compute_dual: Before BTRAN\n");cout<<flush;
    factor.btran(buffer, 1);
    //  printf("compute_dual: After  BTRAN\n");cout<<flush;
    if (an_compute_dual_norm2) {
      btran_sol_norm2 = buffer.norm2();
      btran_sol_norm2 = sqrt(btran_sol_norm2);
    }
    
    // Create a local buffer for the values of reduced costs
    HVector bufferLong;
    bufferLong.setup(solver_lp.numCol_);
    bufferLong.clear();
    matrix.price_by_col(bufferLong, buffer);
    for (int i = 0; i < solver_lp.numCol_; i++) {
      simplex_info.workDual_[i] = simplex_info.workCost_[i] - bufferLong.array[i];
    }
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = solver_lp.numCol_; i < numTot; i++) {
      simplex_info.workDual_[i] = simplex_info.workCost_[i] - buffer.array[i - solver_lp.numCol_];
    }
    
    if (an_compute_dual_norm2) {
      work_dual_norm2 = 0;
      for (int i = 0; i < numTot; i++)
	work_dual_norm2 += simplex_info.workDual_[i] * simplex_info.workDual_[i];
      work_dual_norm2 = sqrt(work_dual_norm2);
      //  printf("compute_dual: B.pi=c_B has ||c_B||=%11.4g; ||pi||=%11.4g;
      //  ||pi^TA-c||=%11.4g\n", btran_rhs_norm2, btran_sol_norm2, work_dual_norm2);
      double current_dual_feasibility_tolerance = simplex_info.dual_feasibility_tolerance;
      double new_dual_feasibility_tolerance = work_dual_norm2 / 1e16;
      if (new_dual_feasibility_tolerance > 1e-1) {
	printf(
	       "Seriously: do you expect to solve an LP with ||pi^TA-c||=%11.4g?\n",
	       work_dual_norm2);
      } else if (new_dual_feasibility_tolerance > 10 * current_dual_feasibility_tolerance) {
	printf(
	       "||pi^TA-c|| = %12g so solving with dual_feasibility_tolerance = %12g\n",
	       work_dual_norm2, new_dual_feasibility_tolerance);
	simplex_info.dual_feasibility_tolerance = new_dual_feasibility_tolerance;
      }
    }
    
    // Now have nonbasic duals
    simplex_info.solver_lp_has_nonbasic_dual_values = true;
  }
  
  void correct_dual(HighsModelObject &highs_model_object, int* free_infeasibility_count) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsRandom &random = highs_model_object.random_;
    const double tau_d = simplex_info.dual_feasibility_tolerance;
    const double inf = HIGHS_CONST_INF;
    int workCount = 0;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) {
      if (basis.nonbasicFlag_[i]) {
	if (simplex_info.workLower_[i] == -inf && simplex_info.workUpper_[i] == inf) {
	  // FREE variable
	  workCount += (fabs(simplex_info.workDual_[i]) >= tau_d);
	} else if (basis.nonbasicMove_[i] * simplex_info.workDual_[i] <= -tau_d) {
	  if (simplex_info.workLower_[i] != -inf && simplex_info.workUpper_[i] != inf) {
	    // Boxed variable = flip
	    flip_bound(highs_model_object, i);
	  } else {
	    // Other variable = shift
	    simplex_info.costs_perturbed = 1;
	    if (basis.nonbasicMove_[i] == 1) {
	      double random_v = random.fraction();
	      double dual = (1 + random_v) * tau_d;
	      double shift = dual - simplex_info.workDual_[i];
	      simplex_info.workDual_[i] = dual;
	      simplex_info.workCost_[i] = simplex_info.workCost_[i] + shift;
	    } else {
	      double dual = -(1 + random.fraction()) * tau_d;
	      double shift = dual - simplex_info.workDual_[i];
	      simplex_info.workDual_[i] = dual;
	      simplex_info.workCost_[i] = simplex_info.workCost_[i] + shift;
	    }
	  }
	}
      }
    }
    *free_infeasibility_count = workCount;
  }
  
  void compute_dual_infeasible_in_dual(HighsModelObject &highs_model_object, int* dual_infeasibility_count) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    int work_count = 0;
    const double inf = HIGHS_CONST_INF;
    const double tau_d = simplex_info.dual_feasibility_tolerance;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) {
      // Only for non basic variables
      if (!basis.nonbasicFlag_[i]) continue;
      // Free
      if (simplex_info.workLower_[i] == -inf && simplex_info.workUpper_[i] == inf)
	work_count += (fabs(simplex_info.workDual_[i]) >= tau_d);
      // In dual, assuming that boxed variables will be flipped
      if (simplex_info.workLower_[i] == -inf || simplex_info.workUpper_[i] == inf)
	work_count += (basis.nonbasicMove_[i] * simplex_info.workDual_[i] <= -tau_d);
    }
    *dual_infeasibility_count = work_count;
  }
  
  void compute_dual_infeasible_in_primal(HighsModelObject &highs_model_object, int* dual_infeasibility_count) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    int work_count = 0;
    const double inf = HIGHS_CONST_INF;
    const double tau_d = simplex_info.dual_feasibility_tolerance;
    const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
    for (int i = 0; i < numTot; i++) {
      // Only for non basic variables
      if (!basis.nonbasicFlag_[i]) continue;
      // Free
      if (simplex_info.workLower_[i] == -inf && simplex_info.workUpper_[i] == inf)
	work_count += (fabs(simplex_info.workDual_[i]) >= tau_d);
      // In primal don't assume flip
      work_count += (basis.nonbasicMove_[i] * simplex_info.workDual_[i] <= -tau_d);
    }
    *dual_infeasibility_count = work_count;
  }
  
  // Compute the primal values (in baseValue) and set the lower and upper bounds
  // of basic variables
  int set_source_out_from_bound(HighsModelObject &highs_model_object, const int column_out) {
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    int source_out = 0;
    if (simplex_info.workLower_[column_out] != simplex_info.workUpper_[column_out]) {
      if (!highs_isInfinity(-simplex_info.workLower_[column_out])) {
	// Finite LB so source_out = -1 ensures value set to LB if LB < UB
	source_out = -1;
	//      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g]
	//      so setting source_out = -1\n", column_out, simplex_info.workLower_[column_out],
	//      simplex_info.workUpper_[column_out]);
      } else {
	// Infinite LB so source_out = 1 ensures value set to UB
	source_out = 1;
	if (!highs_isInfinity(simplex_info.workUpper_[column_out])) {
	  // Free variable => trouble!
	  printf("TROUBLE: variable %d leaving the basis is free!\n", column_out);
	}
      }
    }
    return source_out;
  }
  
  double compute_primal_objective_function_value(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsScale &scale = highs_model_object.scale_;
    double primal_objective_function_value = 0;
    for (int row = 0; row < solver_lp.numRow_; row++) {
      int var = basis.basicIndex_[row];
      if (var < solver_lp.numCol_)
	primal_objective_function_value += simplex_info.baseValue_[row] * solver_lp.colCost_[var];
    }
    for (int col = 0; col < solver_lp.numCol_; col++) {
      if (basis.nonbasicFlag_[col])
	primal_objective_function_value += simplex_info.workValue_[col] * solver_lp.colCost_[col];
    }
    primal_objective_function_value *= scale.cost_;
    return primal_objective_function_value;
  }
  
  // Record the shift in the cost of a particular column
  double shift_cost(HighsModelObject &highs_model_object, int iCol, double amount) {
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    simplex_info.costs_perturbed = 1;
    assert(simplex_info.workShift_[iCol] == 0);
    simplex_info.workShift_[iCol] = amount;
  }
  
  // Undo the shift in the cost of a particular column
  double shift_back(HighsModelObject &highs_model_object, int iCol) {
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    simplex_info.workDual_[iCol] -= simplex_info.workShift_[iCol];
    simplex_info.workShift_[iCol] = 0;
  }
  
  // The major model updates. Factor calls factor.update; Matrix
  // calls matrix.update; updatePivots does everything---and is
  // called from the likes of HDual::updatePivots
  void update_factor(HighsModelObject &highs_model_object, 
		     HVector *column,
		     HVector *row_ep,
		     int *iRow,
		     int *hint) {
    //    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HFactor &factor = highs_model_object.factor_;
    HighsTimer &timer = highs_model_object.timer_;
    
    timer.start(simplex_info.clock_[UpdateFactorClock]);
    factor.update(column, row_ep, iRow, hint);
    // Now have a representation of B^{-1}, but it is not fresh
    simplex_info.solver_lp_has_invert = true;
    if (simplex_info.update_count >= simplex_info.update_limit) *hint = INVERT_HINT_UPDATE_LIMIT_REACHED;
    timer.stop(simplex_info.clock_[UpdateFactorClock]);
  }
  
  void update_pivots(HighsModelObject &highs_model_object, int columnIn, int rowOut, int sourceOut) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsTimer &timer = highs_model_object.timer_;
    
    timer.start(simplex_info.clock_[UpdatePivotsClock]);
    int columnOut = basis.basicIndex_[rowOut];
    
    // Incoming variable
    basis.basicIndex_[rowOut] = columnIn;
    basis.nonbasicFlag_[columnIn] = 0;
    basis.nonbasicMove_[columnIn] = 0;
    simplex_info.baseLower_[rowOut] = simplex_info.workLower_[columnIn];
    simplex_info.baseUpper_[rowOut] = simplex_info.workUpper_[columnIn];
    
    // Outgoing variable
    basis.nonbasicFlag_[columnOut] = 1;
    //  double dlValue;
    //  double vrLb = simplex_info.workLower_[columnOut];
    //  double vrV = simplex_info.workValue_[columnOut];
    //  double vrUb = simplex_info.workUpper_[columnOut];
    if (simplex_info.workLower_[columnOut] == simplex_info.workUpper_[columnOut]) {
      //    dlValue = simplex_info.workLower_[columnOut]-simplex_info.workValue_[columnOut];
      simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
      basis.nonbasicMove_[columnOut] = 0;
    } else if (sourceOut == -1) {
      //    dlValue = simplex_info.workLower_[columnOut]-simplex_info.workValue_[columnOut];
      simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
      basis.nonbasicMove_[columnOut] = 1;
    } else {
      //    dlValue = simplex_info.workUpper_[columnOut]-simplex_info.workValue_[columnOut];
      simplex_info.workValue_[columnOut] = simplex_info.workUpper_[columnOut];
      basis.nonbasicMove_[columnOut] = -1;
    }
    double nwValue = simplex_info.workValue_[columnOut];
    double vrDual = simplex_info.workDual_[columnOut];
    double dlDualObjectiveValue = nwValue*vrDual;
    //  if (abs(nwValue))
    //    printf("update_pivots columnOut = %6d (%2d): [%11.4g, %11.4g, %11.4g], nwValue = %11.4g, dual = %11.4g, dlObj = %11.4g\n",
    //			   columnOut, basis.nonbasicMove_[columnOut], vrLb, vrV, vrUb, nwValue, vrDual, dlDualObjectiveValue);
    simplex_info.updatedDualObjectiveValue += dlDualObjectiveValue;
    simplex_info.update_count++;
    // Update the number of basic logicals
    if (columnOut < solver_lp.numCol_) simplex_info.num_basic_logicals -= 1;
    if (columnIn < solver_lp.numCol_) simplex_info.num_basic_logicals += 1;
    // No longer have a representation of B^{-1}, and certainly not
    // fresh!
    simplex_info.solver_lp_has_invert = false;
    simplex_info.solver_lp_has_fresh_invert = false;
    // Data are no longer fresh from rebuild
    simplex_info.solver_lp_has_fresh_rebuild = false;
    timer.stop(simplex_info.clock_[UpdatePivotsClock]);
  }
  
  void update_matrix(HighsModelObject &highs_model_object, int columnIn, int columnOut) {
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HMatrix &matrix = highs_model_object.matrix_;
    HighsTimer &timer = highs_model_object.timer_;
    
    timer.start(simplex_info.clock_[UpdateMatrixClock]);
    matrix.update(columnIn, columnOut);
    timer.stop(simplex_info.clock_[UpdateMatrixClock]);
  }
  
#ifdef HiGHSDEV
  void util_analyse_lp_solution(HighsModelObject &highs_model_object) {
    HighsLp &solver_lp = highs_model_object.solver_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    HighsBasis &basis = highs_model_object.basis_;
    HighsScale &scale = highs_model_object.scale_;
    if (simplex_info.solution_status != SimplexSolutionStatus::OPTIMAL) return;
    printf("\nAnalysing the model solution\n");
    fflush(stdout);
    const double inf = HIGHS_CONST_INF;
    const double tlValueEr = 1e-8;
    const double tlPrRsduEr = 1e-8;
    const double tlDuRsduEr = 1e-8;
    const double tlPrIfs = simplex_info.primal_feasibility_tolerance;
    const double tlDuIfs = simplex_info.dual_feasibility_tolerance;
    
    // Copy the values of (nonbasic) primal variables and scatter values of primal
    // variables which are basic
    vector<double> value = simplex_info.workValue_;
    for (int iRow = 0; iRow < solver_lp.numRow_; iRow++)
      value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
    
    // Copy the values of (nonbasic) dual variables and zero values of dual
    // variables which are basic
    vector<double> dual = simplex_info.workDual_;
    for (int iRow = 0; iRow < solver_lp.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;
    
    // Allocate and zero values of row primal activites and column dual activities
    // to check the residuals
    vector<double> sclRowPrAct;
    vector<double> rowPrAct;
    sclRowPrAct.assign(solver_lp.numRow_, 0);
    rowPrAct.assign(solver_lp.numRow_, 0);
    vector<double> sclColDuAct;
    vector<double> colDuAct;
    sclColDuAct.assign(solver_lp.numCol_, 0);
    colDuAct.assign(solver_lp.numCol_, 0);
    
    // Determine row primal activites and column dual activities
    for (int iCol = 0; iCol < solver_lp.numCol_; iCol++) {
      //    printf("\nCol %2d\n", iCol);
      double lcSclColDuAct = -solver_lp.colCost_[iCol];
      double lcColDuAct = -(solver_lp.colCost_[iCol] * scale.cost_) / scale.col_[iCol];
      for (int en = solver_lp.Astart_[iCol]; en < solver_lp.Astart_[iCol + 1]; en++) {
	int iRow = solver_lp.Aindex_[en];
	double Avalue_En = solver_lp.Avalue_[en];
	double unsclAvalue_En = Avalue_En / (scale.col_[iCol] * scale.row_[iRow]);
	sclRowPrAct[iRow] += Avalue_En * value[iCol];
	rowPrAct[iRow] += unsclAvalue_En * value[iCol] * scale.col_[iCol];
	//      double lcSum = lcSclColDuAct - Avalue_En*dual[solver_lp.numCol_+iRow];
	//      printf("Row %2d: %11.4g - (%11.4g*%11.4g=%11.4g) = %11.4g\n",
	//      iRow, lcSclColDuAct, Avalue_En, dual[solver_lp.numCol_+iRow],
	//      Avalue_En*dual[solver_lp.numCol_+iRow], lcSum);
	lcSclColDuAct -= Avalue_En * dual[solver_lp.numCol_ + iRow];
	lcColDuAct -=
          unsclAvalue_En * dual[solver_lp.numCol_ + iRow] * scale.cost_ * scale.row_[iRow];
      }
      sclColDuAct[iCol] = lcSclColDuAct;
      colDuAct[iCol] = lcColDuAct;
    }
    
    // Look for column residual errors and infeasibilities - primal and dual
    if (solver_lp.offset_) printf("Primal objective offset is %11.4g\n", solver_lp.offset_);
    double lcPrObjV = 0;
    double lcValue = 0;
    
    int numRpFreeRowEr = 0;
    int maxRpFreeRowEr = 100;
    int numRpFreeColEr = 0;
    int maxRpFreeColEr = 100;
    
    bool rpAllCol = false;
    int numRpCol = 0;
    int mxRpCol = 100;
    bool rpNoCol = false;
    int numColPrIfs = 0;
    double maxColPrIfs = 0;
    double sumColPrIfs = 0;
    int numSclColPrIfs = 0;
    double maxSclColPrIfs = 0;
    double sumSclColPrIfs = 0;
    int numColDuIfs = 0;
    double maxColDuIfs = 0;
    double sumColDuIfs = 0;
    int numSclColDuIfs = 0;
    double maxSclColDuIfs = 0;
    double sumSclColDuIfs = 0;
    int numColDuRsduEr = 0;
    double sumColDuRsduEr = 0;
    double maxColDuRsduEr = 0;
    int numSclColDuRsduEr = 0;
    double sumSclColDuRsduEr = 0;
    double maxSclColDuRsduEr = 0;
    for (int iCol = 0; iCol < solver_lp.numCol_; iCol++) {
      double sclColValue;
      double sclColDuIfs;
      // Get the unscaled column bounds
      double unsclColLower = solver_lp.colLower_[iCol];
      double unsclColUpper = solver_lp.colUpper_[iCol];
      unsclColLower *= unsclColLower == -inf ? 1 : scale.col_[iCol];
      unsclColUpper *= unsclColUpper == +inf ? 1 : scale.col_[iCol];
      // Determine the column primal values given nonbasicMove and the bounds -
      // and check the dual residual errors and infeasibilities
      if (basis.nonbasicFlag_[iCol]) {
	// Nonbasic variable - check that the value array is correct given
	// nonbasicMove and the bounds
	if (basis.nonbasicMove_[iCol] == NONBASIC_MOVE_UP) {
	  // At lower bound
	  sclColValue = solver_lp.colLower_[iCol];
	  sclColDuIfs = max(-dual[iCol], 0.);
	} else if (basis.nonbasicMove_[iCol] == NONBASIC_MOVE_DN) {
	  // At upper bound
	  sclColValue = solver_lp.colUpper_[iCol];
	  sclColDuIfs = max(dual[iCol], 0.);
	} else {
	  // Fixed or free
	  if (solver_lp.colLower_[iCol] == solver_lp.colUpper_[iCol]) {
	    sclColValue = solver_lp.colUpper_[iCol];
	    sclColDuIfs = 0;
	  } else {
	    // Free
	    //	  bool freeEr = false;
	    if (!highs_isInfinity(-solver_lp.colLower_[iCol])) {
	      // freeEr = true;
	      if (numRpFreeColEr < maxRpFreeColEr) {
		numRpFreeColEr++;
		printf(
		       "Column %7d supposed to be free but has lower bound of %g\n",
		       iCol, solver_lp.colLower_[iCol]);
	      }
	    }
	    if (!highs_isInfinity(solver_lp.colUpper_[iCol])) {
	      // freeEr = true;
	      if (numRpFreeColEr < maxRpFreeColEr) {
		numRpFreeColEr++;
		printf(
		       "Column %7d supposed to be free but has upper bound of %g\n",
		       iCol, solver_lp.colUpper_[iCol]);
	      }
	    }
	    sclColValue = value[iCol];
	    sclColDuIfs = abs(dual[iCol]);
	    //	  if (!freeEr) {printf("Column %7d is free with value %g\n",
	    // iCol ,sclColValue);}
	  }
	}
	double valueEr = abs(sclColValue - value[iCol]);
	if (valueEr > tlValueEr) {
	  printf(
		 "Column %7d has value error of %11.4g for sclColValue = %11.4g and "
		 "value[iCol] = %11.4g\n",
		 iCol, valueEr, sclColValue, value[iCol]);
	  sclColValue = value[iCol];
	}
	
      } else {
	// Basic variable
	sclColValue = value[iCol];
	sclColDuIfs = abs(dual[iCol]);
      }
      
      lcPrObjV += sclColValue * solver_lp.colCost_[iCol];
      
      double unsclColValue = sclColValue * scale.col_[iCol];
      //      assert(highs_isInfinity(-sclColValue));
      //      assert(highs_isInfinity(sclColValue));
      // Assess primal infeasibility
      // For scaled values
      double sclColPrIfs = max(
			       max(solver_lp.colLower_[iCol] - sclColValue, sclColValue - solver_lp.colUpper_[iCol]), 0.0);
      if (sclColPrIfs > tlPrIfs) {
	numSclColPrIfs++;
	sumSclColPrIfs += sclColPrIfs;
      }
      maxSclColPrIfs = max(sclColPrIfs, maxSclColPrIfs);
      // For unscaled values
      double colPrIfs = max(
			    max(unsclColLower - unsclColValue, unsclColValue - unsclColUpper), 0.0);
      if (colPrIfs > tlPrIfs) {
	numColPrIfs++;
	sumColPrIfs += colPrIfs;
      }
      maxColPrIfs = max(colPrIfs, maxColPrIfs);
      
      // Assess dual infeasibility
      // In scaled values
      if (sclColDuIfs > tlDuIfs) {
	numSclColDuIfs++;
	sumSclColDuIfs += sclColDuIfs;
      }
      maxSclColDuIfs = max(sclColDuIfs, maxSclColDuIfs);
      // In unscaled values
      double colDuIfs = sclColDuIfs * scale.cost_ / scale.col_[iCol];
      if (colDuIfs > tlDuIfs) {
	numColDuIfs++;
	sumColDuIfs += colDuIfs;
      }
      maxColDuIfs = max(colDuIfs, maxColDuIfs);
      
      // Check column residual errors
      // Using scaled column activities
      double sclColDual = dual[iCol];
      double sclColDuRsduEr = abs(sclColDuAct[iCol] + sclColDual);
      if (sclColDuRsduEr > tlDuRsduEr) {
	/*
	  bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
	  if (rpCol) {
	  numRpCol++;
	  printf("Col    %7d has a   dual residual error of %11.4g for
	  sclColDuAct[iCol] = %11.4g and -sclColDual = %11.4g\n", iCol,
	  sclColDuRsduEr, sclColDuAct[iCol], -sclColDual);
	  }
	*/
	numSclColDuRsduEr++;
	sumSclColDuRsduEr += sclColDuRsduEr;
      }
      maxSclColDuRsduEr = max(sclColDuRsduEr, maxSclColDuRsduEr);
      // Using unscaled column activities
      double colDual = sclColDual * scale.cost_ / scale.col_[iCol];
      double colDuRsduEr = abs(colDuAct[iCol] + colDual);
      if (colDuRsduEr > tlDuRsduEr) {
	/*
	  bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
	  if (rpCol) {
	  numRpCol++;
	  printf("Col    %7d has a   dual residual error of %11.4g for
	  colDuAct[iCol] = %11.4g and -colDual = %11.4g\n", iCol, colDuRsduEr,
	  colDuAct[iCol], -colDual);
	  }
	*/
	numColDuRsduEr++;
	sumColDuRsduEr += colDuRsduEr;
      }
      maxColDuRsduEr = max(colDuRsduEr, maxColDuRsduEr);
      
      bool erFd = sclColPrIfs > tlPrIfs || colPrIfs > tlPrIfs ||
	sclColDuIfs > tlDuIfs || colDuIfs > tlDuIfs ||
	sclColDuRsduEr > tlDuRsduEr || colDuRsduEr > tlDuRsduEr;
      bool rpCol = (rpAllCol || (numRpCol < mxRpCol && erFd)) && !rpNoCol;
      if (rpCol) {
	numRpCol++;
	printf("\nCol %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iCol,
	       basis.nonbasicFlag_[iCol], basis.nonbasicMove_[iCol], scale.col_[iCol]);
	printf(
	       "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
	       "%11.4g)\n",
	       solver_lp.colLower_[iCol], sclColValue, solver_lp.colUpper_[iCol], sclColPrIfs, sclColDuIfs,
	       sclColDuRsduEr);
	printf(
	       "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: %11.4g) "
	       "\n",
	       unsclColLower, unsclColValue, unsclColUpper, colPrIfs, colDuIfs,
	       colDuRsduEr);
      }
    }
    
    printf(
	   "Found %6d   scaled column primal infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs);
    printf(
	   "Found %6d unscaled column primal infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numColPrIfs, sumColPrIfs, maxColPrIfs);
    printf(
	   "Found %6d   scaled column   dual infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs);
    printf(
	   "Found %6d unscaled column   dual infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numColDuIfs, sumColDuIfs, maxColDuIfs);
    printf(
	   "Found %6d   scaled column   dual residual errors: sum %11.4g; max "
	   "%11.4g\n",
	   numSclColDuRsduEr, sumSclColDuRsduEr, maxSclColDuRsduEr);
    printf(
	   "Found %6d unscaled column   dual residual errors: sum %11.4g; max "
	   "%11.4g\n",
	   numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);
    
    printf(
	   "grep_AnMlSolIfsRsduEr,Col,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
	   "d,%g,%g\n",
	   numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs, numColPrIfs, sumColPrIfs,
	   maxColPrIfs, numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs, numColDuIfs,
	   sumColDuIfs, maxColDuIfs, numSclColDuRsduEr, sumSclColDuRsduEr,
	   maxSclColDuRsduEr, numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);
    
    bool rpAllRow = false;
    int numRpRow = 0;
    int mxRpRow = 100;
    bool rpNoRow = false;
    int numRowPrIfs = 0;
    double sumRowPrIfs = 0;
    double maxRowPrIfs = 0;
    int numSclRowPrIfs = 0;
    double sumSclRowPrIfs = 0;
    double maxSclRowPrIfs = 0;
    int numRowDuIfs = 0;
    double maxRowDuIfs = 0;
    double sumRowDuIfs = 0;
    int numSclRowDuIfs = 0;
    double maxSclRowDuIfs = 0;
    double sumSclRowDuIfs = 0;
    int numRowPrRsduEr = 0;
    double sumRowPrRsduEr = 0;
    double maxRowPrRsduEr = 0;
    int numSclRowPrRsduEr = 0;
    double sumSclRowPrRsduEr = 0;
    double maxSclRowPrRsduEr = 0;
    for (int iRow = 0; iRow < solver_lp.numRow_; iRow++) {
      double sclRowValue;
      double sclRowDuIfs;
      // Get the unscaled row bounds
      double unsclRowLower = solver_lp.rowLower_[iRow];
      double unsclRowUpper = solver_lp.rowUpper_[iRow];
      unsclRowLower *= unsclRowLower == -inf ? 1 : scale.row_[iRow];
      unsclRowUpper *= unsclRowUpper == +inf ? 1 : scale.row_[iRow];
      // Determine the row primal values given nonbasicMove and the bounds - and
      // check the dual residual errors and infeasibilities
      if (basis.nonbasicFlag_[solver_lp.numCol_ + iRow]) {
	// Nonbasic variable
	if (basis.nonbasicMove_[solver_lp.numCol_ + iRow] == NONBASIC_MOVE_DN) {
	  // At lower bound
	  sclRowValue = solver_lp.rowLower_[iRow];
	  sclRowDuIfs = max(dual[solver_lp.numCol_ + iRow], 0.);
	} else if (basis.nonbasicMove_[solver_lp.numCol_ + iRow] == NONBASIC_MOVE_UP) {
	  // At upper bound
	  sclRowValue = solver_lp.rowUpper_[iRow];
	  sclRowDuIfs = max(-dual[solver_lp.numCol_ + iRow], 0.);
	} else {
	  // Fixed or free
	  if (solver_lp.rowLower_[iRow] == solver_lp.rowUpper_[iRow]) {
	    sclRowValue = solver_lp.rowUpper_[iRow];
	    sclRowDuIfs = 0.;
	  } else {
	    // Free
	    //	  bool freeEr = false;
	    if (!highs_isInfinity(-solver_lp.rowLower_[iRow])) {
	      // freeEr = true;
	      if (numRpFreeRowEr < maxRpFreeRowEr) {
		numRpFreeRowEr++;
		printf(
		       "Row    %7d supposed to be free but has lower bound of %g\n",
		       iRow, solver_lp.rowLower_[iRow]);
	      }
	    }
	    if (!highs_isInfinity(solver_lp.rowUpper_[iRow])) {
	      // freeEr = true;
	      if (numRpFreeRowEr < maxRpFreeRowEr) {
		numRpFreeRowEr++;
		printf(
		       "Row    %7d supposed to be free but has upper bound of %g\n",
		       iRow, solver_lp.rowUpper_[iRow]);
	      }
	    }
	    sclRowValue = -value[solver_lp.numCol_ + iRow];
	    sclRowDuIfs = abs(dual[solver_lp.numCol_ + iRow]);
	    //	  if (!freeEr) {printf("Row    %7d is free with value %g\n",
	    // iRow, sclRowValue);}
	  }
	}
	double valueEr = abs(sclRowValue + value[solver_lp.numCol_ + iRow]);
	if (valueEr > tlValueEr) {
	  printf(
		 "Row    %7d has value error of %11.4g for sclRowValue = %11.4g and "
		 "-value[solver_lp.numCol_+iRow] = %11.4g\n",
		 iRow, valueEr, sclRowValue, -value[solver_lp.numCol_ + iRow]);
	  sclRowValue = -value[solver_lp.numCol_ + iRow];
	}
      } else {
	// Basic variable
	sclRowValue = -value[solver_lp.numCol_ + iRow];
	sclRowDuIfs = abs(dual[solver_lp.numCol_ + iRow]);
      }
      //      assert(highs_isInfinity(-sclRowValue));
      //      assert(highs_isInfinity(sclRowValue));
      double unsclRowValue = sclRowValue * scale.row_[iRow];
      
      // Assess primal infeasibility
      // For scaled values
      double sclRowPrIfs = max(
			       max(solver_lp.rowLower_[iRow] - sclRowValue, sclRowValue - solver_lp.rowUpper_[iRow]), 0.0);
      if (sclRowPrIfs > tlPrIfs) {
	numSclRowPrIfs++;
	sumSclRowPrIfs += sclRowPrIfs;
      }
      maxSclRowPrIfs = max(sclRowPrIfs, maxSclRowPrIfs);
      // For unscaled values
      double rowPrIfs = max(
			    max(unsclRowLower - unsclRowValue, unsclRowValue - unsclRowUpper), 0.0);
      if (rowPrIfs > tlPrIfs) {
	numRowPrIfs++;
	sumRowPrIfs += rowPrIfs;
      }
      maxRowPrIfs = max(rowPrIfs, maxRowPrIfs);
      
      // Assess dual infeasibility
      // In scaled values
      if (sclRowDuIfs > tlDuIfs) {
	numSclRowDuIfs++;
	sumSclRowDuIfs += sclRowDuIfs;
      }
      maxSclRowDuIfs = max(sclRowDuIfs, maxSclRowDuIfs);
      // In unscaled values
      double rowDuIfs = sclRowDuIfs * scale.cost_ / scale.row_[iRow];
      if (rowDuIfs > tlDuIfs) {
	numRowDuIfs++;
	sumRowDuIfs += rowDuIfs;
      }
      maxRowDuIfs = max(rowDuIfs, maxRowDuIfs);
      
      // Check row residual errors
      // Using scaled row activities
      double sclRowPrRsduEr = abs(sclRowPrAct[iRow] - sclRowValue);
      if (sclRowPrRsduEr > tlPrRsduEr) {
	/*
	  bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
	  if (rpRow) {
	  numRpRow++;
	  printf("Row    %7d has a primal residual error of %11.4g for
	  sclRowPrAct[iRow] = %11.4g and sclRowValue = %11.4g\n", iRow,
	  sclRowPrRsduEr, sclRowPrAct[iRow], sclRowValue);
	  }
	*/
	numSclRowPrRsduEr++;
	sumSclRowPrRsduEr += sclRowPrRsduEr;
      }
      maxSclRowPrRsduEr = max(sclRowPrRsduEr, maxSclRowPrRsduEr);
      // Using unscaled row activities
      double rowValue = sclRowValue / scale.row_[iRow];
      double rowPrRsduEr = abs(rowPrAct[iRow] - rowValue);
      if (rowPrRsduEr > tlPrRsduEr) {
	/*
	  bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
	  if (rpRow) {
	  numRpRow++;
	  printf("Row    %7d has a primal residual error of %11.4g for
	  rowPrAct[iRow] = %11.4g and rowValue = %11.4g\n", iRow, rowPrRsduEr,
	  rowPrAct[iRow], rowValue);
	  }
	*/
	numRowPrRsduEr++;
	sumRowPrRsduEr += rowPrRsduEr;
      }
      maxRowPrRsduEr = max(rowPrRsduEr, maxRowPrRsduEr);
      
      bool erFd = sclRowPrIfs > tlPrIfs || rowPrIfs > tlPrIfs ||
	sclRowDuIfs > tlDuIfs || rowDuIfs > tlDuIfs ||
	sclRowPrRsduEr > tlPrRsduEr || rowPrRsduEr > tlPrRsduEr;
      bool rpRow = (rpAllRow || (numRpRow < mxRpRow && erFd)) && !rpNoRow;
      if (rpRow) {
	numRpRow++;
	printf("Row %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iRow,
	       basis.nonbasicFlag_[solver_lp.numCol_ + iRow], basis.nonbasicMove_[solver_lp.numCol_ + iRow],
	       scale.row_[iRow]);
	printf(
	       "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
	       "%11.4g)\n",
	       solver_lp.rowLower_[iRow], sclRowValue, solver_lp.rowUpper_[iRow], sclRowPrIfs, sclRowDuIfs,
	       sclRowPrRsduEr);
	printf(
	       "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
	       "%11.4g)\n",
	       unsclRowLower, unsclRowValue, unsclRowUpper, rowPrIfs, rowDuIfs,
	       rowPrRsduEr);
      }
    }
    printf(
	   "Found %6d   scaled    row primal infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs);
    printf(
	   "Found %6d unscaled    row primal infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numRowPrIfs, sumRowPrIfs, maxRowPrIfs);
    printf(
	   "Found %6d   scaled    row   dual infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs);
    printf(
	   "Found %6d unscaled    row   dual infeasibilities: sum %11.4g; max "
	   "%11.4g\n",
	   numRowDuIfs, sumRowDuIfs, maxRowDuIfs);
    printf(
	   "Found %6d   scaled    row primal residual errors: sum %11.4g; max "
	   "%11.4g\n",
	   numSclRowPrRsduEr, sumSclRowPrRsduEr, maxSclRowPrRsduEr);
    printf(
	   "Found %6d unscaled    row primal residual errors: sum %11.4g; max "
	   "%11.4g\n",
	   numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);
    
    printf(
	   "grep_AnMlSolIfsRsduEr,Row,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
	   "d,%g,%g\n",
	   numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs, numRowPrIfs, sumRowPrIfs,
	   maxRowPrIfs, numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs, numRowDuIfs,
	   sumRowDuIfs, maxRowDuIfs, numSclRowPrRsduEr, sumSclRowPrRsduEr,
	   maxSclRowPrRsduEr, numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);
    
    lcPrObjV *= scale.cost_;
    lcPrObjV += solver_lp.offset_;
    double dualObjectiveValue = simplex_info.dualObjectiveValue;
    double ObjEr = abs(dualObjectiveValue - lcPrObjV) / max(1.0, fabs(dualObjectiveValue));
    printf(
	   "Relative objective error of %11.4g: dualObjectiveValue = %g; lcPrObjV = %g\n",
	   ObjEr, dualObjectiveValue, lcPrObjV);
    
  }
#endif
  
  void report_iteration_count_dual_objective_value(HighsModelObject &highs_model_object, int i_v) {
    int iteration_count = highs_model_object.simplex_info_.iteration_count;
    double dual_objective_value = highs_model_object.simplex_info_.dualObjectiveValue;
    HighsPrintMessage(ML_MINIMAL, "%10d  %20.10e  %2d\n", iteration_count, dual_objective_value, i_v);
  }

};
#endif // SIMPLEX_HSIMPLEX_H_
