/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
//#include <cassert>
////#include <iostream>

#include "simplex/HEkk.h"

#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "simplex/HEkkDebug.h"
#include "simplex/HEkkDual.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HFactorDebug.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsRandom.h"

#ifdef OPENMP
#include "omp.h"
#endif

// using std::cout;
// using std::endl;

HighsStatus HEkk::passLp(const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;

  simplex_lp_ = lp;
  // Shouldn't have to check the incoming LP since this is an internal
  // call, but it may be an LP that's set up internally with errors
  // :-) ...
  if (options_.highs_debug_level > HIGHS_DEBUG_LEVEL_NONE) {
    // ... so, if debugging, check the LP.
    call_status = assessLp(simplex_lp_, options_);
    return_status = interpretCallStatus(call_status, return_status, "assessLp");
    if (return_status == HighsStatus::Error) return return_status;
  }
  initialiseForNewLp();
  return HighsStatus::OK;
}

HighsStatus HEkk::solve() {
  initialiseAnalysis();
  if (analysis_.analyse_simplex_time)
    analysis_.simplexTimerStart(SimplexTotalClock);
  iteration_count_ = 0;
  if (initialiseForSolve() == HighsStatus::Error) return HighsStatus::Error;

  assert(simplex_lp_status_.has_basis);
  assert(simplex_lp_status_.has_invert);
  assert(simplex_lp_status_.valid);
  if (scaled_model_status_ == HighsModelStatus::OPTIMAL) return HighsStatus::OK;

  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  std::string algorithm;

  // Indicate that dual and primal rays are not known
  simplex_lp_status_.has_dual_ray = false;
  simplex_lp_status_.has_primal_ray = false;

  chooseSimplexStrategyThreads(options_, simplex_info_);
  int simplex_strategy = simplex_info_.simplex_strategy;

  // Initial solve according to strategy
  if (simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
    algorithm = "primal";
    reportSimplexPhaseIterations(options_.logfile, 
				 options_.output_flag,
				 options_.log_to_console,
				 iteration_count_, simplex_info_, true);
    HighsLogMessage(options_.logfile, HighsMessageType::INFO,
                    "Using EKK primal simplex solver");
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkPrimal::solve");
  } else {
    algorithm = "dual";
    reportSimplexPhaseIterations(options_.logging_file, 
				 options_.output_flag,
				 options_.log_to_console,
				 iteration_count_, simplex_info_, true);
    HEkkDual dual_solver(*this);
    dual_solver.options();
    //
    // Solve, depending on the particular strategy
    if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
      HighsLogMessage(
          options_.logfile, HighsMessageType::INFO,
          "Using EKK parallel dual simplex solver - SIP with %d threads",
          simplex_info_.num_threads);
    } else if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
      HighsLogMessage(
          options_.logfile, HighsMessageType::INFO,
          "Using EKK parallel dual simplex solver - PAMI with %d threads",
          simplex_info_.num_threads);
    } else {
      HighsLogMessage(options_.logfile, HighsMessageType::INFO,
                      "Using EKK dual simplex solver - serial");
    }
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkDual::solve");
  }
  reportSimplexPhaseIterations(options_.logging_file,
			       options_.output_flag,
			       options_.log_to_console,
			       iteration_count_, simplex_info_);
  if (return_status == HighsStatus::Error) return return_status;
  HighsLogMessage(
      options_.logfile, HighsMessageType::INFO,
      "EKK %s simplex solver returns %d primal and %d dual infeasibilities: "
      "Status %s",
      algorithm.c_str(), simplex_info_.num_primal_infeasibility,
      simplex_info_.num_dual_infeasibility,
      utilHighsModelStatusToString(scaled_model_status_).c_str());
  if (scaled_model_status_ == HighsModelStatus::NOTSET) {
    call_status = cleanup();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkDual::solve");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (analysis_.analyse_simplex_time) {
    analysis_.simplexTimerStop(SimplexTotalClock);
    analysis_.reportSimplexTimer();
  }
  if (analysis_.analyse_simplex_data) analysis_.summaryReport();
  if (analysis_.analyse_factor_data) analysis_.reportInvertFormData();
  if (analysis_.analyse_factor_time) analysis_.reportFactorTimer();
  return return_status;
}

HighsStatus HEkk::cleanup() {
  // Clean up from a point with either primal or dual
  // infeasiblilities, but not both
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  if (simplex_info_.num_primal_infeasibility) {
    // Primal infeasibilities, so should be just dual phase 2
    assert(!simplex_info_.num_dual_infeasibility);
    // Use dual simplex (phase 2) with Devex pricing and no perturbation
    simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_DUAL;
    simplex_info_.dual_simplex_cost_perturbation_multiplier = 0;
    simplex_info_.dual_edge_weight_strategy =
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX;
    HEkkDual dual_solver(*this);
    dual_solver.options();
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkDual::solve");
    if (return_status == HighsStatus::Error) return return_status;
  } else {
    // Dual infeasibilities, so should be just primal phase 2
    assert(!simplex_info_.num_primal_infeasibility);
    // Use primal simplex (phase 2) with no perturbation
    simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
    simplex_info_.primal_simplex_bound_perturbation_multiplier = 0;
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkPrimal::solve");
    if (return_status == HighsStatus::Error) return return_status;
  }
  return return_status;
}

HighsStatus HEkk::setBasis() {
  // Set up nonbasicFlag and basicIndex for a logical basis
  const int num_col = simplex_lp_.numCol_;
  const int num_row = simplex_lp_.numRow_;
  const int num_tot = num_col + num_row;
  simplex_basis_.nonbasicFlag_.resize(num_tot);
  simplex_basis_.nonbasicMove_.resize(num_tot);
  simplex_basis_.basicIndex_.resize(num_row);
  for (int iCol = 0; iCol < num_col; iCol++) {
    simplex_basis_.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
    double lower = simplex_lp_.colLower_[iCol];
    double upper = simplex_lp_.colUpper_[iCol];
    int move = illegal_move_value;
    if (lower == upper) {
      // Fixed
      move = NONBASIC_MOVE_ZE;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed. Set to bound of LP that is closer to
        // zero
        if (move == illegal_move_value) {
          if (fabs(lower) < fabs(upper)) {
            move = NONBASIC_MOVE_UP;
          } else {
            move = NONBASIC_MOVE_DN;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = NONBASIC_MOVE_UP;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = NONBASIC_MOVE_DN;
    } else {
      // FREE
      move = NONBASIC_MOVE_ZE;
    }
    assert(move != illegal_move_value);
    simplex_basis_.nonbasicMove_[iCol] = move;
  }
  for (int iRow = 0; iRow < num_row; iRow++) {
    int iVar = num_col + iRow;
    simplex_basis_.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
    simplex_basis_.basicIndex_[iRow] = iVar;
  }
  simplex_info_.num_basic_logicals = num_row;
  simplex_lp_status_.has_basis = true;
  return HighsStatus::OK;
}

HighsStatus HEkk::setBasis(const HighsBasis& basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  if (debugBasisConsistent(options_, simplex_lp_, basis) ==
      HighsDebugStatus::LOGICAL_ERROR) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Supposed to be a Highs basis, but not valid");
    return HighsStatus::Error;
  }
  int num_col = simplex_lp_.numCol_;
  int num_row = simplex_lp_.numRow_;
  int num_tot = num_col + num_row;
  // Resize the basis in case none has yet been defined for this LP
  simplex_basis_.nonbasicFlag_.resize(num_tot);
  simplex_basis_.nonbasicMove_.resize(num_tot);
  simplex_basis_.basicIndex_.resize(num_row);
  int num_basic_variables = 0;
  for (int iCol = 0; iCol < num_col; iCol++) {
    int iVar = iCol;
    const double lower = simplex_lp_.colLower_[iCol];
    const double upper = simplex_lp_.colUpper_[iCol];
    if (basis.col_status[iCol] == HighsBasisStatus::BASIC) {
      simplex_basis_.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
      simplex_basis_.nonbasicMove_[iVar] = 0;
      simplex_basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      simplex_basis_.nonbasicFlag_[iVar] = NONBASIC_FLAG_TRUE;
      if (basis.col_status[iCol] == HighsBasisStatus::LOWER) {
        if (lower == upper) {
          simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
        } else {
          simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_UP;
        }
      } else if (basis.col_status[iCol] == HighsBasisStatus::UPPER) {
        simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_DN;
      } else {
        assert(basis.col_status[iCol] == HighsBasisStatus::ZERO);
        simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      }
    }
  }
  for (int iRow = 0; iRow < num_row; iRow++) {
    int iVar = num_col + iRow;
    const double lower = simplex_lp_.rowLower_[iRow];
    const double upper = simplex_lp_.rowUpper_[iRow];
    if (basis.row_status[iRow] == HighsBasisStatus::BASIC) {
      simplex_basis_.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
      simplex_basis_.nonbasicMove_[iVar] = 0;
      simplex_basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      simplex_basis_.nonbasicFlag_[iVar] = NONBASIC_FLAG_TRUE;
      if (basis.row_status[iRow] == HighsBasisStatus::LOWER) {
        if (lower == upper) {
          simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
        } else {
          simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_DN;
        }
      } else if (basis.row_status[iRow] == HighsBasisStatus::UPPER) {
        simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_UP;
      } else {
        assert(basis.row_status[iRow] == HighsBasisStatus::ZERO);
        simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      }
    }
  }
  simplex_lp_status_.has_basis = true;
  return HighsStatus::OK;
}

HighsStatus HEkk::setBasis(const SimplexBasis& basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  if (debugBasisConsistent(options_, simplex_lp_, basis) ==
      HighsDebugStatus::LOGICAL_ERROR) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Supposed to be a Highs basis, but not valid");
    return HighsStatus::Error;
  }
  simplex_basis_.nonbasicFlag_ = basis.nonbasicFlag_;
  simplex_basis_.nonbasicMove_ = basis.nonbasicMove_;
  simplex_basis_.basicIndex_ = basis.basicIndex_;
  simplex_lp_status_.has_basis = true;
  return HighsStatus::OK;
}

HighsSolution HEkk::getSolution() {
  HighsSolution solution;
  // Scatter the basic primal values
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++)
    simplex_info_.workValue_[simplex_basis_.basicIndex_[iRow]] =
        simplex_info_.baseValue_[iRow];
  // Zero the basic dual values
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++)
    simplex_info_.workDual_[simplex_basis_.basicIndex_[iRow]] = 0;

  // Now we can get the solution
  solution.col_value.resize(simplex_lp_.numCol_);
  solution.col_dual.resize(simplex_lp_.numCol_);
  solution.row_value.resize(simplex_lp_.numRow_);
  solution.row_dual.resize(simplex_lp_.numRow_);

  for (int iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    solution.col_value[iCol] = simplex_info_.workValue_[iCol];
    solution.col_dual[iCol] =
        (int)simplex_lp_.sense_ * simplex_info_.workDual_[iCol];
  }
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    solution.row_value[iRow] =
        -simplex_info_.workValue_[simplex_lp_.numCol_ + iRow];
    solution.row_dual[iRow] =
        (int)simplex_lp_.sense_ *
        simplex_info_.workDual_[simplex_lp_.numCol_ + iRow];
  }
  return solution;
}

HighsBasis HEkk::getHighsBasis() {
  int num_col = simplex_lp_.numCol_;
  int num_row = simplex_lp_.numRow_;
  HighsBasis basis;
  basis.col_status.resize(num_col);
  basis.row_status.resize(num_row);
  assert(simplex_lp_status_.has_basis);
  basis.valid_ = false;
  for (int iCol = 0; iCol < num_col; iCol++) {
    int iVar = iCol;
    const double lower = simplex_lp_.colLower_[iCol];
    const double upper = simplex_lp_.colUpper_[iCol];
    HighsBasisStatus basis_status = HighsBasisStatus::NONBASIC;
    if (!simplex_basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::BASIC;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_UP) {
      basis_status = HighsBasisStatus::LOWER;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_DN) {
      basis_status = HighsBasisStatus::UPPER;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_ZE) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::LOWER;
      } else {
        basis_status = HighsBasisStatus::ZERO;
      }
    }
    basis.col_status[iCol] = basis_status;
  }
  for (int iRow = 0; iRow < num_row; iRow++) {
    int iVar = num_col + iRow;
    const double lower = simplex_lp_.rowLower_[iRow];
    const double upper = simplex_lp_.rowUpper_[iRow];
    HighsBasisStatus basis_status = HighsBasisStatus::NONBASIC;
    if (!simplex_basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::BASIC;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_UP) {
      basis_status = HighsBasisStatus::UPPER;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_DN) {
      basis_status = HighsBasisStatus::LOWER;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_ZE) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::LOWER;
      } else {
        basis_status = HighsBasisStatus::ZERO;
      }
    }
    basis.row_status[iRow] = basis_status;
  }
  basis.valid_ = true;
  return basis;
}

int HEkk::initialiseSimplexLpBasisAndFactor(const bool only_from_known_basis) {
  // If there's no basis, return error if the basis has to be known,
  // otherwise set a logical basis
  if (!simplex_lp_status_.has_basis) {
    if (only_from_known_basis) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "Simplex basis should be known but isn't");
      return -(int)HighsStatus::Error;
    }
    setBasis();
  }
  const int rank_deficiency = computeFactor();
  if (rank_deficiency) {
    // Basis is rank deficient
    if (only_from_known_basis) {
      // If only this basis should be used, then return error
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                      "Supposed to be a full-rank basis, but incorrect");
      return rank_deficiency;
    }
    // Account for rank deficiency by correcing nonbasicFlag
    handleRankDeficiency();
    updateSimplexLpStatus(simplex_lp_status_, LpAction::NEW_BASIS);
    setNonbasicMove();
    simplex_lp_status_.has_basis = true;
    simplex_lp_status_.has_invert = true;
    simplex_lp_status_.has_fresh_invert = true;
  }
  assert(simplex_lp_status_.has_invert);
  return 0;
}

void HEkk::handleRankDeficiency() {
  int rank_deficiency = factor_.rank_deficiency;
  vector<int>& noPvC = factor_.noPvC;
  vector<int>& noPvR = factor_.noPvR;
  for (int k = 0; k < rank_deficiency; k++) {
    int variable_in = simplex_lp_.numCol_ + noPvR[k];
    int variable_out = noPvC[k];
    simplex_basis_.nonbasicFlag_[variable_in] = NONBASIC_FLAG_FALSE;
    simplex_basis_.nonbasicFlag_[variable_out] = NONBASIC_FLAG_TRUE;
  }
  simplex_lp_status_.has_matrix = false;
}

HighsSolutionParams HEkk::getSolutionParams() {
  HighsSolutionParams solution_params;
  solution_params.primal_feasibility_tolerance =
      options_.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options_.dual_feasibility_tolerance;
  if (scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    solution_params.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
    solution_params.dual_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    solution_params.primal_status = PrimalDualStatus::STATUS_NOTSET;
    solution_params.dual_status = PrimalDualStatus::STATUS_NOTSET;
  }
  // Output from solution analysis method
  solution_params.objective_function_value =
      simplex_info_.primal_objective_value;
  solution_params.num_primal_infeasibility =
      simplex_info_.num_primal_infeasibility;
  solution_params.max_primal_infeasibility =
      simplex_info_.max_primal_infeasibility;
  solution_params.sum_primal_infeasibility =
      simplex_info_.sum_primal_infeasibility;
  solution_params.num_dual_infeasibility = simplex_info_.num_dual_infeasibility;
  solution_params.max_dual_infeasibility = simplex_info_.max_dual_infeasibility;
  solution_params.sum_dual_infeasibility = simplex_info_.sum_dual_infeasibility;
  return solution_params;
}

// Private methods

void HEkk::initialiseForNewLp() {
  setSimplexOptions();
  initialiseControl();
  initialiseSimplexLpRandomVectors();
  simplex_lp_status_.initialised = true;
}

HighsStatus HEkk::initialiseForSolve() {
  const int error_return = initialiseSimplexLpBasisAndFactor();
  assert(!error_return);
  if (error_return) return HighsStatus::Error;
  assert(simplex_lp_status_.has_basis);

  initialiseMatrix();  // Timed
  allocateWorkAndBaseArrays();
  initialiseCost(SimplexAlgorithm::PRIMAL, SOLVE_PHASE_UNKNOWN, false);
  initialiseBound(SimplexAlgorithm::PRIMAL, SOLVE_PHASE_UNKNOWN, false);
  initialiseNonbasicValueAndMove();
  computePrimal();                // Timed
  computeDual();                  // Timed
  computeSimplexInfeasible();     // Timed
  computeDualObjectiveValue();    // Timed
  computePrimalObjectiveValue();  // Timed
  simplex_lp_status_.valid = true;

  bool primal_feasible = simplex_info_.num_primal_infeasibility == 0;
  bool dual_feasible = simplex_info_.num_dual_infeasibility == 0;
  scaled_model_status_ = HighsModelStatus::NOTSET;
  if (primal_feasible && dual_feasible)
    scaled_model_status_ = HighsModelStatus::OPTIMAL;
  return HighsStatus::OK;
}

void HEkk::setSimplexOptions() {
  // Copy values of HighsOptions for the simplex solver
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  simplex_info_.simplex_strategy = options_.simplex_strategy;
  simplex_info_.dual_edge_weight_strategy =
      options_.simplex_dual_edge_weight_strategy;
  simplex_info_.price_strategy = options_.simplex_price_strategy;
  simplex_info_.dual_simplex_cost_perturbation_multiplier =
      options_.dual_simplex_cost_perturbation_multiplier;
  simplex_info_.primal_simplex_bound_perturbation_multiplier =
      options_.primal_simplex_bound_perturbation_multiplier;
  simplex_info_.factor_pivot_threshold = options_.factor_pivot_threshold;
  simplex_info_.update_limit = options_.simplex_update_limit;

  // Set values of internal options
  simplex_info_.store_squared_primal_infeasibility = true;
}

void HEkk::initialiseSimplexLpRandomVectors() {
  const int num_col = simplex_lp_.numCol_;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  if (!num_tot) return;
  // Instantiate and (re-)initialise the random number generator
  //  HighsRandom random;
  HighsRandom& random = random_;
  random.initialise();

  if (num_col) {
    // Generate a random permutation of the column indices
    simplex_info_.numColPermutation_.resize(num_col);
    vector<int>& numColPermutation = simplex_info_.numColPermutation_;
    for (int i = 0; i < num_col; i++) numColPermutation[i] = i;
    for (int i = num_col - 1; i >= 1; i--) {
      int j = random.integer() % (i + 1);
      std::swap(numColPermutation[i], numColPermutation[j]);
    }
  }

  // Re-initialise the random number generator and generate the
  // random vectors in the same order as hsol to maintain repeatable
  // performance
  random.initialise();
  //
  // Generate a random permutation of all the indices
  simplex_info_.numTotPermutation_.resize(num_tot);
  vector<int>& numTotPermutation = simplex_info_.numTotPermutation_;
  for (int i = 0; i < num_tot; i++) numTotPermutation[i] = i;
  for (int i = num_tot - 1; i >= 1; i--) {
    int j = random.integer() % (i + 1);
    std::swap(numTotPermutation[i], numTotPermutation[j]);
  }

  // Generate a vector of random reals
  simplex_info_.numTotRandomValue_.resize(num_tot);
  vector<double>& numTotRandomValue = simplex_info_.numTotRandomValue_;
  for (int i = 0; i < num_tot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

void HEkk::chooseSimplexStrategyThreads(const HighsOptions& options,
                                        HighsSimplexInfo& simplex_info) {
  // Given a simplex basis and solution, use the number of primal and
  // dual infeasibilities to determine which simplex variant to use.
  //
  // 1. If it is "CHOOSE", in which case an approapriate stratgy is
  // used
  //
  // 2. If re-solving choose the strategy appropriate to primal or
  // dual feasibility
  //
  int simplex_strategy = options.simplex_strategy;
  if (simplex_info.num_primal_infeasibility > 0) {
    // Not primal feasible, so use dual simplex if choice is permitted
    if (simplex_strategy == SIMPLEX_STRATEGY_CHOOSE)
      simplex_strategy = SIMPLEX_STRATEGY_DUAL;
  } else {
    // Primal feasible - so must be dual infeasible
    assert(simplex_info.num_dual_infeasibility > 0);
    // Use primal simplex if choice is permitted
    if (simplex_strategy == SIMPLEX_STRATEGY_CHOOSE)
      simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
  }
  // Set min/max_threads to correspond to serial code. They will be
  // set to other values if parallel options are used.
  simplex_info.min_threads = 1;
  simplex_info.max_threads = 1;
  // Record the min/max minimum number of HiGHS threads in the options
  const int highs_min_threads = options.highs_min_threads;
  const int highs_max_threads = options.highs_max_threads;
  int omp_max_threads = 0;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  if (options.parallel == on_string &&
      simplex_strategy == SIMPLEX_STRATEGY_DUAL) {
    // The parallel strategy is on and the simplex strategy is dual so use
    // PAMI if there are enough OMP threads
    if (omp_max_threads >= DUAL_MULTI_MIN_THREADS)
      simplex_strategy = SIMPLEX_STRATEGY_DUAL_MULTI;
  }
  //
  // If parallel stratgies are used, the minimum number of HiGHS threads used
  // will be set to be at least the minimum required for the strategy
  //
  // All this is independent of the number of OMP threads available,
  // since code with multiple HiGHS threads can be run in serial.
#ifdef OPENMP
  if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
    simplex_info.min_threads = max(DUAL_TASKS_MIN_THREADS, highs_min_threads);
    simplex_info.max_threads = max(simplex_info.min_threads, highs_max_threads);
  } else if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
    simplex_info.min_threads = max(DUAL_MULTI_MIN_THREADS, highs_min_threads);
    simplex_info.max_threads = max(simplex_info.min_threads, highs_max_threads);
  }
#endif
  // Set the number of HiGHS threads to be used to be the maximum
  // number to be used
  simplex_info.num_threads = simplex_info.max_threads;
  // Give a warning if the number of threads to be used is fewer than
  // the minimum number of HiGHS threads allowed
  if (simplex_info.num_threads < highs_min_threads) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Using %d HiGHS threads for parallel strategy rather than "
                    "minimum number (%d) specified in options",
                    simplex_info.num_threads, highs_min_threads);
  }
  // Give a warning if the number of threads to be used is more than
  // the maximum number of HiGHS threads allowed
  if (simplex_info.num_threads > highs_max_threads) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Using %d HiGHS threads for parallel strategy rather than "
                    "maximum number (%d) specified in options",
                    simplex_info.num_threads, highs_max_threads);
  }
  // Give a warning if the number of threads to be used is fewer than
  // the number of OMP threads available
  if (simplex_info.num_threads > omp_max_threads) {
    HighsLogMessage(
        options.logfile, HighsMessageType::WARNING,
        "Number of OMP threads available = %d < %d = Number of HiGHS threads "
        "to be used: Parallel performance will be less than anticipated",
        omp_max_threads, simplex_info.num_threads);
  }
  // Simplex strategy is now fixed - so set the value to be referred
  // to in the simplex solver
  simplex_info.simplex_strategy = simplex_strategy;
  // Official start of solver Start the solve clock - because
  // setupForSimplexSolve has simplex computations

  if (simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Primal simplex solver unavailable");
    simplex_strategy = SIMPLEX_STRATEGY_DUAL;
  }
}

bool HEkk::getNonsingularInverse(const int solve_phase) {
  assert(simplex_lp_status_.has_basis);
  const vector<int>& basicIndex = simplex_basis_.basicIndex_;
  // Take a copy of basicIndex from before INVERT to be used as the
  // saved ordering of basic variables - so reinvert will run
  // identically.
  const vector<int> basicIndex_before_compute_factor = basicIndex;
  // Save the number of updates performed in case it has to be used to determine
  // a limit
  const int simplex_update_count = simplex_info_.update_count;
  // Dual simplex edge weights are identified with rows, so must be
  // permuted according to INVERT. This must be done if workEdWt_ is
  // not NULL.
  const bool handle_edge_weights = workEdWt_ != NULL;
  // Scatter the edge weights so that, after INVERT, they can be
  // gathered according to the new permutation of basicIndex
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (int i = 0; i < simplex_lp_.numRow_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }

  // Call computeFactor to perform INVERT
  int rank_deficiency = computeFactor();
  const bool artificial_rank_deficiency = false;  //  true;//
  if (artificial_rank_deficiency) {
    if (!simplex_info_.phase1_backtracking_test_done &&
        solve_phase == SOLVE_PHASE_1) {
      // Claim rank deficiency to test backtracking
      printf("Phase1 (Iter %d) Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      simplex_info_.phase1_backtracking_test_done = true;
    } else if (!simplex_info_.phase2_backtracking_test_done &&
               solve_phase == SOLVE_PHASE_2) {
      // Claim rank deficiency to test backtracking
      printf("Phase2 (Iter %d) Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      simplex_info_.phase2_backtracking_test_done = true;
    }
  }
  if (rank_deficiency) {
    // Rank deficient basis, so backtrack to last full rank basis
    //
    // Get the last nonsingular basis - so long as there is one
    if (!getBacktrackingBasis(workEdWtFull_)) return false;
    // Record that backtracking is taking place
    simplex_info_.backtracking_ = true;
    updateSimplexLpStatus(simplex_lp_status_, LpAction::BACKTRACKING);
    int backtrack_rank_deficiency = computeFactor();
    // This basis has previously been inverted successfully, so it shouldn't be
    // singular
    if (backtrack_rank_deficiency) return false;
    // simplex update limit will be half of the number of updates
    // performed, so make sure that at least one update was performed
    if (simplex_update_count <= 1) return false;
    int use_simplex_update_limit = simplex_info_.update_limit;
    int new_simplex_update_limit = simplex_update_count / 2;
    simplex_info_.update_limit = new_simplex_update_limit;
    HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
                    "Rank deficiency of %d after %d simplex updates, so "
                    "backtracking: max updates reduced from %d to %d",
                    rank_deficiency, simplex_update_count,
                    use_simplex_update_limit, new_simplex_update_limit);
  } else {
    // Current basis is full rank so save it
    putBacktrackingBasis(basicIndex_before_compute_factor, workEdWtFull_);
    // Indicate that backtracking is not taking place
    simplex_info_.backtracking_ = false;
    // Reset the update limit in case this is the first successful
    // inversion after backtracking
    simplex_info_.update_limit = options_.simplex_update_limit;
  }
  if (handle_edge_weights) {
    // Gather the edge weights according to the permutation of
    // basicIndex after INVERT
    analysis_.simplexTimerStart(PermWtClock);
    for (int i = 0; i < simplex_lp_.numRow_; i++)
      workEdWt_[i] = workEdWtFull_[basicIndex[i]];
    analysis_.simplexTimerStop(PermWtClock);
  }
  return true;
}

bool HEkk::getBacktrackingBasis(double* scattered_edge_weights) {
  if (!simplex_info_.valid_backtracking_basis_) return false;
  simplex_basis_ = simplex_info_.backtracking_basis_;
  simplex_info_.costs_perturbed =
      simplex_info_.backtracking_basis_costs_perturbed_;
  simplex_info_.workShift_ = simplex_info_.backtracking_basis_workShift_;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (int iVar = 0; iVar < num_tot; iVar++)
      scattered_edge_weights[iVar] =
          simplex_info_.backtracking_basis_edge_weights_[iVar];
  }
  return true;
}

void HEkk::putBacktrackingBasis() {
  const vector<int>& basicIndex = simplex_basis_.basicIndex_;
  const bool handle_edge_weights = workEdWt_ != NULL;
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (int i = 0; i < simplex_lp_.numRow_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }
  putBacktrackingBasis(basicIndex, workEdWtFull_);
}

void HEkk::putBacktrackingBasis(
    const vector<int>& basicIndex_before_compute_factor,
    double* scattered_edge_weights) {
  simplex_info_.valid_backtracking_basis_ = true;
  simplex_info_.backtracking_basis_ = simplex_basis_;
  simplex_info_.backtracking_basis_.basicIndex_ =
      basicIndex_before_compute_factor;
  simplex_info_.backtracking_basis_costs_perturbed_ =
      simplex_info_.costs_perturbed;
  simplex_info_.backtracking_basis_workShift_ = simplex_info_.workShift_;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (int iVar = 0; iVar < num_tot; iVar++)
      simplex_info_.backtracking_basis_edge_weights_[iVar] =
          scattered_edge_weights[iVar];
  }
}

void HEkk::computePrimalObjectiveValue() {
  analysis_.simplexTimerStart(ComputePrObjClock);
  simplex_info_.primal_objective_value = 0;
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    int iVar = simplex_basis_.basicIndex_[iRow];
    if (iVar < simplex_lp_.numCol_) {
      simplex_info_.primal_objective_value +=
          simplex_info_.baseValue_[iRow] * simplex_lp_.colCost_[iVar];
    }
  }
  for (int iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    if (simplex_basis_.nonbasicFlag_[iCol])
      simplex_info_.primal_objective_value +=
          simplex_info_.workValue_[iCol] * simplex_lp_.colCost_[iCol];
  }
  simplex_info_.primal_objective_value *= cost_scale_;
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  simplex_info_.primal_objective_value += simplex_lp_.offset_;
  // Now have primal objective value
  simplex_lp_status_.has_primal_objective_value = true;
  analysis_.simplexTimerStop(ComputePrObjClock);
}

void HEkk::computeDualObjectiveValue(const int phase) {
  analysis_.simplexTimerStart(ComputeDuObjClock);
  simplex_info_.dual_objective_value = 0;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (simplex_basis_.nonbasicFlag_[iCol]) {
      const double term =
          simplex_info_.workValue_[iCol] * simplex_info_.workDual_[iCol];
      if (term) {
        simplex_info_.dual_objective_value +=
            simplex_info_.workValue_[iCol] * simplex_info_.workDual_[iCol];
      }
    }
  }
  simplex_info_.dual_objective_value *= cost_scale_;
  if (phase != 1) {
    // In phase 1 the dual objective has no objective
    // shift. Otherwise, if minimizing the shift is added. If
    // maximizing, workCost (and hence workDual) are negated, so the
    // shift is subtracted. Hence the shift is added according to the
    // sign implied by sense_
    simplex_info_.dual_objective_value +=
        ((int)simplex_lp_.sense_) * simplex_lp_.offset_;
  }
  // Now have dual objective value
  simplex_lp_status_.has_dual_objective_value = true;
  analysis_.simplexTimerStop(ComputeDuObjClock);
}

int HEkk::computeFactor() {
  if (!simplex_lp_status_.has_factor_arrays) {
    assert(simplex_info_.factor_pivot_threshold >=
           options_.factor_pivot_threshold);
    factor_.setup(simplex_lp_.numCol_, simplex_lp_.numRow_,
                  &simplex_lp_.Astart_[0], &simplex_lp_.Aindex_[0],
                  &simplex_lp_.Avalue_[0], &simplex_basis_.basicIndex_[0],
                  options_.highs_debug_level, options_.logfile, options_.output,
                  options_.message_level, simplex_info_.factor_pivot_threshold,
                  options_.factor_pivot_tolerance);
    simplex_lp_status_.has_factor_arrays = true;
  }
  analysis_.simplexTimerStart(InvertClock);
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  if (analysis_.analyse_factor_time) {
    int thread_id = 0;
#ifdef OPENMP
    thread_id = omp_get_thread_num();
#endif
    factor_timer_clock_pointer =
        analysis_.getThreadFactorTimerClockPtr(thread_id);
  }
  const int rank_deficiency = factor_.build(factor_timer_clock_pointer);
  if (analysis_.analyse_factor_data) analysis_.updateInvertFormData(factor_);

  const bool force = rank_deficiency;
  debugCheckInvert(options_, factor_, force);

  if (rank_deficiency) {
    // Have an invertible representation, but of B with column(s)
    // replacements due to singularity. So no (fresh) representation of
    // B^{-1}
    simplex_lp_status_.has_invert = false;
    simplex_lp_status_.has_fresh_invert = false;
  } else {
    // Now have a representation of B^{-1}, and it is fresh!
    simplex_lp_status_.has_invert = true;
    simplex_lp_status_.has_fresh_invert = true;
  }
  // Set the update count to zero since the corrected invertible
  // representation may be used for an initial basis. In any case the
  // number of updates shouldn't be positive
  simplex_info_.update_count = 0;

  analysis_.simplexTimerStop(InvertClock);
  return rank_deficiency;
}

void HEkk::initialiseMatrix() {
  if (!simplex_lp_status_.has_matrix) {
    analysis_.simplexTimerStart(matrixSetupClock);
    matrix_.setup(simplex_lp_.numCol_, simplex_lp_.numRow_,
                  &simplex_lp_.Astart_[0], &simplex_lp_.Aindex_[0],
                  &simplex_lp_.Avalue_[0], &simplex_basis_.nonbasicFlag_[0]);
    simplex_lp_status_.has_matrix = true;
    analysis_.simplexTimerStop(matrixSetupClock);
  }
}

void HEkk::setNonbasicMove() {
  const bool have_solution = false;
  // Don't have a simplex basis since nonbasicMove is not set up.

  // Assign nonbasicMove using as much information as is available
  double lower;
  double upper;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  simplex_basis_.nonbasicMove_.resize(num_tot);

  for (int iVar = 0; iVar < num_tot; iVar++) {
    if (!simplex_basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      continue;
    }
    // Nonbasic variable
    if (iVar < simplex_lp_.numCol_) {
      lower = simplex_lp_.colLower_[iVar];
      upper = simplex_lp_.colUpper_[iVar];
    } else {
      int iRow = iVar - simplex_lp_.numCol_;
      lower = -simplex_lp_.rowUpper_[iRow];
      upper = -simplex_lp_.rowLower_[iRow];
    }
    int move = illegal_move_value;
    if (lower == upper) {
      // Fixed
      move = NONBASIC_MOVE_ZE;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Determine the bound to set the value to according to, in order of
        // priority
        //
        // 1. Any solution value
        if (have_solution) {
          double midpoint = 0.5 * (lower + upper);
          double value = simplex_info_.workValue_[iVar];
          if (value < midpoint) {
            move = NONBASIC_MOVE_UP;
          } else {
            move = NONBASIC_MOVE_DN;
          }
        }
        // 2. Bound of original LP that is closer to zero
        if (move == illegal_move_value) {
          if (fabs(lower) < fabs(upper)) {
            move = NONBASIC_MOVE_UP;
          } else {
            move = NONBASIC_MOVE_DN;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = NONBASIC_MOVE_UP;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = NONBASIC_MOVE_DN;
    } else {
      // FREE
      move = NONBASIC_MOVE_ZE;
    }
    assert(move != illegal_move_value);
    simplex_basis_.nonbasicMove_[iVar] = move;
  }
}

void HEkk::allocateWorkAndBaseArrays() {
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  simplex_info_.workCost_.resize(num_tot);
  simplex_info_.workDual_.resize(num_tot);
  simplex_info_.workShift_.resize(num_tot);

  simplex_info_.workLower_.resize(num_tot);
  simplex_info_.workUpper_.resize(num_tot);
  simplex_info_.workRange_.resize(num_tot);
  simplex_info_.workValue_.resize(num_tot);
  simplex_info_.workLowerShift_.resize(num_tot);
  simplex_info_.workUpperShift_.resize(num_tot);

  // Feel that it should be possible to resize this with in dual
  // solver, and only if Devex is being used, but a pointer to it
  // needs to be set up when constructing HDual
  simplex_info_.devex_index_.resize(num_tot);

  simplex_info_.baseLower_.resize(simplex_lp_.numRow_);
  simplex_info_.baseUpper_.resize(simplex_lp_.numRow_);
  simplex_info_.baseValue_.resize(simplex_lp_.numRow_);
}

void HEkk::initialiseLpColBound() {
  for (int iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    simplex_info_.workLower_[iCol] = simplex_lp_.colLower_[iCol];
    simplex_info_.workUpper_[iCol] = simplex_lp_.colUpper_[iCol];
    simplex_info_.workRange_[iCol] =
        simplex_info_.workUpper_[iCol] - simplex_info_.workLower_[iCol];
    simplex_info_.workLowerShift_[iCol] = 0;
    simplex_info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowBound() {
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    int iCol = simplex_lp_.numCol_ + iRow;
    simplex_info_.workLower_[iCol] = -simplex_lp_.rowUpper_[iRow];
    simplex_info_.workUpper_[iCol] = -simplex_lp_.rowLower_[iRow];
    simplex_info_.workRange_[iCol] =
        simplex_info_.workUpper_[iCol] - simplex_info_.workLower_[iCol];
    simplex_info_.workLowerShift_[iCol] = 0;
    simplex_info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseCost(const SimplexAlgorithm algorithm,
                          const int solvePhase, const bool perturb) {
  // Copy the cost
  initialiseLpColCost();
  initialiseLpRowCost();
  simplex_info_.costs_perturbed = 0;
  // Primal simplex costs are either from the LP or set specially in phase 1
  if (algorithm == SimplexAlgorithm::PRIMAL) return;
  // Dual simplex costs are either from the LP or perturbed
  if (!perturb || simplex_info_.dual_simplex_cost_perturbation_multiplier == 0)
    return;
  // Perturb the original costs, scale down if is too big
  int num_original_nonzero_cost = 0;
  if (analysis_.analyse_simplex_data)
    printf("grep_DuPtrb: Cost perturbation for %s\n",
           simplex_lp_.model_name_.c_str());
  double bigc = 0;
  for (int i = 0; i < simplex_lp_.numCol_; i++) {
    const double abs_cost = fabs(simplex_info_.workCost_[i]);
    bigc = max(bigc, abs_cost);
    if (analysis_.analyse_simplex_data && abs_cost) num_original_nonzero_cost++;
  }
  const int pct0 = (100 * num_original_nonzero_cost) / simplex_lp_.numCol_;
  double average_cost = 0;
  if (analysis_.analyse_simplex_data) {
    if (num_original_nonzero_cost) {
      average_cost = bigc / num_original_nonzero_cost;
    } else {
      printf("grep_DuPtrb:    STRANGE initial workCost has non nonzeros\n");
    }
    printf(
        "grep_DuPtrb:    Initially have %d nonzero costs (%3d%%) with bigc = "
        "%g "
        "and average = %g\n",
        num_original_nonzero_cost, pct0, bigc, average_cost);
  }
  if (bigc > 100) {
    bigc = sqrt(sqrt(bigc));
    if (analysis_.analyse_simplex_data)
      printf("grep_DuPtrb:    Large so set bigc = sqrt(bigc) = %g\n", bigc);
  }

  // If there are few boxed variables, we will just use simple perturbation
  double boxedRate = 0;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int i = 0; i < num_tot; i++)
    boxedRate += (simplex_info_.workRange_[i] < 1e30);
  boxedRate /= num_tot;
  if (boxedRate < 0.01) {
    bigc = min(bigc, 1.0);
    if (analysis_.analyse_simplex_data)
      printf(
          "grep_DuPtrb:    small boxedRate (%g) so set bigc = min(bigc, 1.0) = "
          "%g\n",
          boxedRate, bigc);
  }
  // Determine the perturbation base
  double base = 5e-7 * bigc;
  if (analysis_.analyse_simplex_data)
    printf("grep_DuPtrb:    Perturbation base = %g\n", base);

  // Now do the perturbation
  for (int i = 0; i < simplex_lp_.numCol_; i++) {
    double lower = simplex_lp_.colLower_[i];
    double upper = simplex_lp_.colUpper_[i];
    double xpert = (fabs(simplex_info_.workCost_[i]) + 1) * base *
                   simplex_info_.dual_simplex_cost_perturbation_multiplier *
                   (1 + simplex_info_.numTotRandomValue_[i]);
    const double previous_cost = simplex_info_.workCost_[i];
    if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper >= HIGHS_CONST_INF) {  // Lower
      simplex_info_.workCost_[i] += xpert;
    } else if (lower <= -HIGHS_CONST_INF) {  // Upper
      simplex_info_.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex_info_.workCost_[i] +=
          (simplex_info_.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
    if (analysis_.analyse_simplex_data) {
      const double perturbation1 =
          fabs(simplex_info_.workCost_[i] - previous_cost);
      if (perturbation1)
        updateValueDistribution(perturbation1,
                                analysis_.cost_perturbation1_distribution);
    }
  }
  for (int i = simplex_lp_.numCol_; i < num_tot; i++) {
    double perturbation2 =
        (0.5 - simplex_info_.numTotRandomValue_[i]) *
        simplex_info_.dual_simplex_cost_perturbation_multiplier * 1e-12;
    simplex_info_.workCost_[i] += perturbation2;
    if (analysis_.analyse_simplex_data) {
      perturbation2 = fabs(perturbation2);
      updateValueDistribution(perturbation2,
                              analysis_.cost_perturbation2_distribution);
    }
  }
  simplex_info_.costs_perturbed = 1;
}

void HEkk::initialiseBound(const SimplexAlgorithm algorithm,
                           const int solve_phase, const bool perturb) {
  initialiseLpColBound();
  initialiseLpRowBound();
  simplex_info_.bounds_perturbed = 0;
  // Primal simplex bounds are either from the LP or perturbed
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    if (!perturb ||
        simplex_info_.primal_simplex_bound_perturbation_multiplier == 0)
      return;
    // Perturb the bounds
    // Determine the smallest and largest finite lower/upper bounds
    int num_col = simplex_lp_.numCol_;
    int num_row = simplex_lp_.numRow_;
    int num_tot = num_col + num_row;
    double min_abs_lower = HIGHS_CONST_INF;
    double max_abs_lower = -1;
    double min_abs_upper = HIGHS_CONST_INF;
    double max_abs_upper = -1;
    for (int iVar = 0; iVar < num_tot; iVar++) {
      double abs_lower = fabs(simplex_info_.workLower_[iVar]);
      double abs_upper = fabs(simplex_info_.workUpper_[iVar]);
      if (abs_lower && abs_lower < HIGHS_CONST_INF) {
        min_abs_lower = min(abs_lower, min_abs_lower);
        max_abs_lower = max(abs_lower, max_abs_lower);
      }
      if (abs_upper && abs_upper < HIGHS_CONST_INF) {
        min_abs_upper = min(abs_upper, min_abs_upper);
        max_abs_upper = max(abs_upper, max_abs_upper);
      }
    }
    // printf(
    //     "Nonzero finite lower bounds in [%9.4g, %9.4g]; upper bounds in "
    //     "[%9.4g, %9.4g]\n",
    //     min_abs_lower, max_abs_lower, min_abs_upper, max_abs_upper);

    const double base =
        simplex_info_.primal_simplex_bound_perturbation_multiplier * 5e-7;
    for (int iVar = 0; iVar < num_tot; iVar++) {
      double lower = simplex_info_.workLower_[iVar];
      double upper = simplex_info_.workUpper_[iVar];
      const bool fixed = lower == upper;
      // Don't perturb bounds of nonbasic fixed variables as they stay nonbasic
      if (simplex_basis_.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE && fixed)
        continue;
      double random_value = simplex_info_.numTotRandomValue_[iVar];
      if (lower > -HIGHS_CONST_INF) {
        if (lower < -1) {
          lower -= random_value * base * (-lower);
        } else if (lower < 1) {
          lower -= random_value * base;
        } else {
          lower -= random_value * base * lower;
        }
        simplex_info_.workLower_[iVar] = lower;
      }
      if (upper < HIGHS_CONST_INF) {
        if (upper < -1) {
          upper += random_value * base * (-upper);
        } else if (upper < 1) {
          upper += random_value * base;
        } else {
          upper += random_value * base * upper;
        }
        simplex_info_.workUpper_[iVar] = upper;
      }
      simplex_info_.workRange_[iVar] =
          simplex_info_.workUpper_[iVar] - simplex_info_.workLower_[iVar];
      if (simplex_basis_.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) continue;
      // Set values of nonbasic variables
      if (simplex_basis_.nonbasicMove_[iVar] > 0) {
        simplex_info_.workValue_[iVar] = lower;
      } else if (simplex_basis_.nonbasicMove_[iVar] < 0) {
        simplex_info_.workValue_[iVar] = upper;
      }
    }
    for (int iRow = 0; iRow < num_row; iRow++) {
      int iVar = simplex_basis_.basicIndex_[iRow];
      simplex_info_.baseLower_[iRow] = simplex_info_.workLower_[iVar];
      simplex_info_.baseUpper_[iRow] = simplex_info_.workUpper_[iVar];
    }
    simplex_info_.bounds_perturbed = 1;
    return;
  }
  // Dual simplex costs are either from the LP or set to special values in phase
  // 1
  assert(algorithm == SimplexAlgorithm::DUAL);
  if (solve_phase == SOLVE_PHASE_2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = HIGHS_CONST_INF;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (simplex_info_.workLower_[iCol] == -inf &&
        simplex_info_.workUpper_[iCol] == inf) {
      // Don't change for row variables: they should never become
      // nonbasic when starting from a logical basis, and no crash
      // should make a free row nonbasic, but could an advanced basis
      // make a free row nonbasic.
      // But what it it happened?
      if (iCol >= simplex_lp_.numCol_) continue;
      simplex_info_.workLower_[iCol] = -1000,
      simplex_info_.workUpper_[iCol] = 1000;  // FREE
    } else if (simplex_info_.workLower_[iCol] == -inf) {
      simplex_info_.workLower_[iCol] = -1,
      simplex_info_.workUpper_[iCol] = 0;  // UPPER
    } else if (simplex_info_.workUpper_[iCol] == inf) {
      simplex_info_.workLower_[iCol] = 0,
      simplex_info_.workUpper_[iCol] = 1;  // LOWER
    } else {
      simplex_info_.workLower_[iCol] = 0,
      simplex_info_.workUpper_[iCol] = 0;  // BOXED or FIXED
    }
    simplex_info_.workRange_[iCol] =
        simplex_info_.workUpper_[iCol] - simplex_info_.workLower_[iCol];
  }
}

void HEkk::initialiseLpColCost() {
  for (int iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    simplex_info_.workCost_[iCol] =
        (int)simplex_lp_.sense_ * simplex_lp_.colCost_[iCol];
    simplex_info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowCost() {
  for (int iCol = simplex_lp_.numCol_;
       iCol < simplex_lp_.numCol_ + simplex_lp_.numRow_; iCol++) {
    simplex_info_.workCost_[iCol] = 0;
    simplex_info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseNonbasicValueAndMove() {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int iVar = 0; iVar < num_tot; iVar++) {
    if (!simplex_basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      simplex_basis_.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
      continue;
    }
    // Nonbasic variable
    const double lower = simplex_info_.workLower_[iVar];
    const double upper = simplex_info_.workUpper_[iVar];
    const int original_move = simplex_basis_.nonbasicMove_[iVar];
    double value;
    int move = illegal_move_value;
    if (lower == upper) {
      // Fixed
      value = lower;
      move = NONBASIC_MOVE_ZE;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (original_move == NONBASIC_MOVE_UP) {
          // Set at lower
          value = lower;
          move = NONBASIC_MOVE_UP;
        } else if (original_move == NONBASIC_MOVE_DN) {
          // Set at upper
          value = upper;
          move = NONBASIC_MOVE_DN;
        } else {
          // Invalid nonbasicMove: correct and set value at lower
          value = lower;
          move = NONBASIC_MOVE_UP;
        }
      } else {
        // Lower
        value = lower;
        move = NONBASIC_MOVE_UP;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      value = upper;
      move = NONBASIC_MOVE_DN;
    } else {
      // FREE
      value = 0;
      move = NONBASIC_MOVE_ZE;
    }
    assert(move != illegal_move_value);
    simplex_basis_.nonbasicMove_[iVar] = move;
    simplex_info_.workValue_[iVar] = value;
  }
}

void HEkk::pivotColumnFtran(const int iCol, HVector& col_aq) {
  analysis_.simplexTimerStart(FtranClock);
  col_aq.clear();
  col_aq.packFlag = true;
  matrix_.collect_aj(col_aq, iCol, 1);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis_.col_aq_density);
  factor_.ftran(col_aq, analysis_.col_aq_density,
                analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
  int num_row = simplex_lp_.numRow_;
  const double local_col_aq_density = (double)col_aq.count / num_row;
  analysis_.updateOperationResultDensity(local_col_aq_density,
                                         analysis_.col_aq_density);
  updateOperationResultDensity(local_col_aq_density,
                               simplex_info_.col_aq_density);
  analysis_.simplexTimerStop(FtranClock);
}

void HEkk::unitBtran(const int iRow, HVector& row_ep) {
  analysis_.simplexTimerStart(BtranClock);
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = iRow;
  row_ep.array[iRow] = 1;
  row_ep.packFlag = true;
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep,
                                    analysis_.row_ep_density);
  factor_.btran(row_ep, analysis_.row_ep_density,
                analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
  int num_row = simplex_lp_.numRow_;
  const double local_row_ep_density = (double)row_ep.count / num_row;
  analysis_.updateOperationResultDensity(local_row_ep_density,
                                         analysis_.row_ep_density);
  updateOperationResultDensity(local_row_ep_density,
                               simplex_info_.row_ep_density);
  analysis_.simplexTimerStop(BtranClock);
}

void HEkk::fullBtran(HVector& buffer) {
  // Performs BTRAN on the buffer supplied. Make sure that
  // buffer.count is large (>simplex_lp_.numRow_ to be sure) rather
  // than 0 if the indices of the RHS (and true value of buffer.count)
  // isn't known.
  analysis_.simplexTimerStart(BtranFullClock);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_FULL, buffer,
                                    analysis_.dual_col_density);
  factor_.btran(buffer, analysis_.dual_col_density,
                analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_FULL, buffer);
  const double local_dual_col_density =
      (double)buffer.count / simplex_lp_.numRow_;
  analysis_.updateOperationResultDensity(local_dual_col_density,
                                         analysis_.dual_col_density);
  updateOperationResultDensity(local_dual_col_density,
                               simplex_info_.dual_col_density);
  analysis_.simplexTimerStop(BtranFullClock);
}

void HEkk::choosePriceTechnique(const int price_strategy,
                                const double row_ep_density,
                                bool& use_col_price,
                                bool& use_row_price_w_switch) {
  // By default switch to column PRICE when pi_p has at least this
  // density
  const double density_for_column_price_switch = 0.75;
  use_col_price =
      (price_strategy == SIMPLEX_PRICE_STRATEGY_COL) ||
      (price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH &&
       row_ep_density > density_for_column_price_switch);
  use_row_price_w_switch =
      price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH ||
      price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH;
}

void HEkk::tableauRowPrice(const HVector& row_ep, HVector& row_ap) {
  analysis_.simplexTimerStart(PriceClock);
  const int solver_num_row = simplex_lp_.numRow_;
  const int solver_num_col = simplex_lp_.numCol_;
  const double local_density = 1.0 * row_ep.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  choosePriceTechnique(simplex_info_.price_strategy, local_density,
                       use_col_price, use_row_price_w_switch);
  if (analysis_.analyse_simplex_data) {
    if (use_col_price) {
      const double historical_density_for_non_hypersparse_operation = 1;
      analysis_.operationRecordBefore(
          ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
          historical_density_for_non_hypersparse_operation);
      analysis_.num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis_.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
                                      analysis_.row_ep_density);
      analysis_.num_row_price_with_switch++;
    } else {
      analysis_.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
                                      analysis_.row_ep_density);
      analysis_.num_row_price++;
    }
  }
  row_ap.clear();
  if (use_col_price) {
    // Perform column-wise PRICE
    matrix_.priceByColumn(row_ap, row_ep);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of row_ap
    // becomes extreme
    const double switch_density = matrix_.hyperPRICE;
    matrix_.priceByRowSparseResultWithSwitch(
        row_ap, row_ep, analysis_.row_ap_density, 0, switch_density);
  } else {
    // Perform hyper-sparse row-wise PRICE
    matrix_.priceByRowSparseResult(row_ap, row_ep);
  }
  if (use_col_price) {
    // Column-wise PRICE computes components corresponding to basic
    // variables, so zero these by exploiting the fact that, for basic
    // variables, nonbasicFlag[*]=0
    const int* nonbasicFlag = &simplex_basis_.nonbasicFlag_[0];
    for (int iCol = 0; iCol < solver_num_col; iCol++)
      row_ap.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_ap density
  const double local_row_ap_density = (double)row_ap.count / solver_num_col;
  analysis_.updateOperationResultDensity(local_row_ap_density,
                                         analysis_.row_ap_density);
  updateOperationResultDensity(local_row_ap_density,
                               simplex_info_.row_ap_density);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ap);
  analysis_.simplexTimerStop(PriceClock);
}

void HEkk::fullPrice(const HVector& full_col, HVector& full_row) {
  analysis_.simplexTimerStart(PriceFullClock);
  full_row.clear();
  if (analysis_.analyse_simplex_data) {
    const double historical_density_for_non_hypersparse_operation = 1;
    analysis_.operationRecordBefore(
        ANALYSIS_OPERATION_TYPE_PRICE_FULL, full_col,
        historical_density_for_non_hypersparse_operation);
  }
  matrix_.priceByColumn(full_row, full_col);
  if (analysis_.analyse_simplex_data)
    analysis_.operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_FULL,
                                   full_row);
  analysis_.simplexTimerStop(PriceFullClock);
}

void HEkk::computePrimal() {
  analysis_.simplexTimerStart(ComputePrimalClock);
  const int num_row = simplex_lp_.numRow_;
  const int num_col = simplex_lp_.numCol_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(num_row);
  primal_col.clear();
  for (int i = 0; i < num_col + num_row; i++) {
    if (simplex_basis_.nonbasicFlag_[i] && simplex_info_.workValue_[i] != 0) {
      matrix_.collect_aj(primal_col, i, simplex_info_.workValue_[i]);
    }
  }
  // It's possible that the buffer has no nonzeros, so performing
  // FTRAN is unnecessary. Not much of a saving, but the zero density
  // looks odd in the analysis!
  if (primal_col.count) {
    factor_.ftran(primal_col, analysis_.primal_col_density,
                  analysis_.pointer_serial_factor_clocks);
    const double local_primal_col_density = (double)primal_col.count / num_row;
    analysis_.updateOperationResultDensity(local_primal_col_density,
                                           analysis_.primal_col_density);
    updateOperationResultDensity(local_primal_col_density,
                                 simplex_info_.primal_col_density);
  }
  for (int i = 0; i < num_row; i++) {
    int iCol = simplex_basis_.basicIndex_[i];
    simplex_info_.baseValue_[i] = -primal_col.array[i];
    simplex_info_.baseLower_[i] = simplex_info_.workLower_[iCol];
    simplex_info_.baseUpper_[i] = simplex_info_.workUpper_[iCol];
  }
  // Indicate that the primal infeasiblility information isn't known
  simplex_info_.num_primal_infeasibility = illegal_infeasibility_count;
  simplex_info_.max_primal_infeasibility = illegal_infeasibility_measure;
  simplex_info_.sum_primal_infeasibility = illegal_infeasibility_measure;

  // Now have basic primals
  simplex_lp_status_.has_basic_primal_values = true;
  analysis_.simplexTimerStop(ComputePrimalClock);
}

void HEkk::computeDual() {
  analysis_.simplexTimerStart(ComputeDualClock);
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(simplex_lp_.numRow_);
  dual_col.clear();
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    const double value =
        simplex_info_.workCost_[simplex_basis_.basicIndex_[iRow]] +
        simplex_info_.workShift_[simplex_basis_.basicIndex_[iRow]];
    if (value) {
      dual_col.index[dual_col.count++] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // Copy the costs in case the basic costs are all zero
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int i = 0; i < num_tot; i++)
    simplex_info_.workDual_[i] = simplex_info_.workCost_[i];

  if (dual_col.count) {
    fullBtran(dual_col);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(simplex_lp_.numCol_);
    fullPrice(dual_col, dual_row);
    for (int i = 0; i < simplex_lp_.numCol_; i++)
      simplex_info_.workDual_[i] -= dual_row.array[i];
    for (int i = simplex_lp_.numCol_; i < num_tot; i++)
      simplex_info_.workDual_[i] -= dual_col.array[i - simplex_lp_.numCol_];
  }
  // Indicate that the dual infeasiblility information isn't known
  simplex_info_.num_dual_infeasibility = illegal_infeasibility_count;
  simplex_info_.max_dual_infeasibility = illegal_infeasibility_measure;
  simplex_info_.sum_dual_infeasibility = illegal_infeasibility_measure;

  // Now have nonbasic duals
  simplex_lp_status_.has_nonbasic_dual_values = true;
  analysis_.simplexTimerStop(ComputeDualClock);
}

void HEkk::computeDualInfeasibleWithFlips() {
  // Computes num/max/sum of dual infeasibliities according to
  // nonbasicMove, using the bounds only to identify free variables
  // and non-boxed. Fixed variables are assumed to have nonbasicMove=0
  // so that no dual infeasibility is counted for them. Indeed, when
  // called from cleanup() at the end of dual phase 1, nonbasicMove
  // relates to the phase 1 bounds, but workLower and workUpper will
  // have been set to phase 2 values!
  const double scaled_dual_feasibility_tolerance =
      options_.dual_feasibility_tolerance;
  // Possibly verify that nonbasicMove is correct for fixed variables
  //  debugFixedNonbasicMove(ekk_instance_);

  int num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;

  for (int iVar = 0; iVar < num_tot; iVar++) {
    if (!simplex_basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = simplex_info_.workLower_[iVar];
    const double upper = simplex_info_.workUpper_[iVar];
    const double dual = simplex_info_.workDual_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else if (highs_isInfinity(-lower) || highs_isInfinity(upper)) {
      // Not free or boxed: any dual infeasibility is given by value
      // signed by nonbasicMove.
      //
      // For boxed variables, nonbasicMove may have the wrong sign for
      // dual, but nonbasicMove and the primal value can be flipped to
      // achieve dual feasiblility.
      dual_infeasibility = -simplex_basis_.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  simplex_info_.num_dual_infeasibility = num_dual_infeasibility;
  simplex_info_.max_dual_infeasibility = max_dual_infeasibility;
  simplex_info_.sum_dual_infeasibility = sum_dual_infeasibility;
}

double HEkk::computeDualForTableauColumn(const int iVar,
                                         const HVector& tableau_column) {
  const vector<double>& workCost = simplex_info_.workCost_;
  const vector<int>& basicIndex = simplex_basis_.basicIndex_;

  double dual = simplex_info_.workCost_[iVar];
  for (int i = 0; i < tableau_column.count; i++) {
    int iRow = tableau_column.index[i];
    dual -= tableau_column.array[iRow] * workCost[basicIndex[iRow]];
  }
  return dual;
}

void HEkk::correctDual(int* free_infeasibility_count) {
  const double tau_d = options_.dual_feasibility_tolerance;
  const double inf = HIGHS_CONST_INF;
  int workCount = 0;
  double flip_dual_objective_value_change = 0;
  double shift_dual_objective_value_change = 0;
  int num_flip = 0;
  int num_shift = 0;
  double sum_flip = 0;
  double sum_shift = 0;
  const int num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (int i = 0; i < num_tot; i++) {
    if (simplex_basis_.nonbasicFlag_[i]) {
      if (simplex_info_.workLower_[i] == -inf &&
          simplex_info_.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(simplex_info_.workDual_[i]) >= tau_d);
      } else if (simplex_basis_.nonbasicMove_[i] * simplex_info_.workDual_[i] <=
                 -tau_d) {
        if (simplex_info_.workLower_[i] != -inf &&
            simplex_info_.workUpper_[i] != inf) {
          // Boxed variable = flip
          const int move = simplex_basis_.nonbasicMove_[i];
          flipBound(i);
          double flip =
              simplex_info_.workUpper_[i] - simplex_info_.workLower_[i];
          // Negative dual at lower bound (move=1): flip to upper
          // bound so objective contribution is change in value (flip)
          // times dual, being move*flip*dual
          //
          // Positive dual at upper bound (move=-1): flip to lower
          // bound so objective contribution is change in value
          // (-flip) times dual, being move*flip*dual
          double local_dual_objective_change =
              move * flip * simplex_info_.workDual_[i];
          local_dual_objective_change *= cost_scale_;
          flip_dual_objective_value_change += local_dual_objective_change;
          num_flip++;
          sum_flip += fabs(flip);
        } else if (simplex_info_.allow_cost_perturbation) {
          // Other variable = shift
          //
          // Before 07/01/20, these shifts were always done, but doing
          // it after cost perturbation has been removed can lead to
          // cycling when primal infeasibility has been detecteed in
          // Phase 2, since the shift below removes dual
          // infeasibilities, which are then reinstated after the dual
          // values are recomputed.
          //
          // ToDo: Not shifting leads to dual infeasibilities when an
          // LP is declared to be (primal) infeasible. Should go to
          // phase 1 primal simplex to "prove" infeasibility.
          simplex_info_.costs_perturbed = 1;
          std::string direction;
          double shift;
          if (simplex_basis_.nonbasicMove_[i] == 1) {
            direction = "  up";
            double dual = (1 + random_.fraction()) * tau_d;
            shift = dual - simplex_info_.workDual_[i];
            simplex_info_.workDual_[i] = dual;
            simplex_info_.workCost_[i] = simplex_info_.workCost_[i] + shift;
          } else {
            direction = "down";
            double dual = -(1 + random_.fraction()) * tau_d;
            shift = dual - simplex_info_.workDual_[i];
            simplex_info_.workDual_[i] = dual;
            simplex_info_.workCost_[i] = simplex_info_.workCost_[i] + shift;
          }
          double local_dual_objective_change =
              shift * simplex_info_.workValue_[i];
          local_dual_objective_change *= cost_scale_;
          shift_dual_objective_value_change += local_dual_objective_change;
          num_shift++;
          sum_shift += fabs(shift);
          HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
                            "Move %s: cost shift = %g; objective change = %g\n",
                            direction.c_str(), shift,
                            local_dual_objective_change);
        }
      }
    }
  }
  if (num_flip)
    HighsPrintMessage(
        options_.output, options_.message_level, ML_VERBOSE,
        "Performed %d flip(s): total = %g; objective change = %g\n", num_flip,
        sum_flip, flip_dual_objective_value_change);
  if (num_shift)
    HighsPrintMessage(
        options_.output, options_.message_level, ML_DETAILED,
        "Performed %d cost shift(s): total = %g; objective change = %g\n",
        num_shift, sum_shift, shift_dual_objective_value_change);
  *free_infeasibility_count = workCount;
}

void HEkk::flipBound(const int iCol) {
  int* nonbasicMove = &simplex_basis_.nonbasicMove_[0];
  const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  simplex_info_.workValue_[iCol] = move == 1 ? simplex_info_.workLower_[iCol]
                                             : simplex_info_.workUpper_[iCol];
}

bool HEkk::reinvertOnNumericalTrouble(
    const std::string method_name, double& numerical_trouble_measure,
    const double alpha_from_col, const double alpha_from_row,
    const double numerical_trouble_tolerance) {
  double abs_alpha_from_col = fabs(alpha_from_col);
  double abs_alpha_from_row = fabs(alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  numerical_trouble_measure = abs_alpha_diff / min_abs_alpha;
  const int update_count = simplex_info_.update_count;
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool reinvert = numerical_trouble && update_count > 0;
  ekkDebugReportReinvertOnNumericalTrouble(
      method_name, *this, numerical_trouble_measure, alpha_from_col,
      alpha_from_row, numerical_trouble_tolerance, reinvert);
  if (reinvert) {
    // Consider increasing the Markowitz multiplier
    const double current_pivot_threshold = simplex_info_.factor_pivot_threshold;
    double new_pivot_threshold = 0;
    if (current_pivot_threshold < default_pivot_threshold) {
      // Threshold is below default value, so increase it
      new_pivot_threshold =
          min(current_pivot_threshold * pivot_threshold_change_factor,
              default_pivot_threshold);
    } else if (current_pivot_threshold < max_pivot_threshold) {
      // Threshold is below max value, so increase it if few updates have been
      // performed
      if (update_count < 10)
        new_pivot_threshold =
            min(current_pivot_threshold * pivot_threshold_change_factor,
                max_pivot_threshold);
    }
    if (new_pivot_threshold) {
      HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
                      "   Increasing Markowitz threshold to %g",
                      new_pivot_threshold);
      simplex_info_.factor_pivot_threshold = new_pivot_threshold;
      factor_.setPivotThreshold(new_pivot_threshold);
    }
  }
  return reinvert;
}

// The major model updates. Factor calls factor_.update; Matrix
// calls matrix_.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HEkk::updateFactor(HVector* column, HVector* row_ep, int* iRow,
                        int* hint) {
  analysis_.simplexTimerStart(UpdateFactorClock);
  factor_.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  simplex_lp_status_.has_invert = true;
  if (simplex_info_.update_count >= simplex_info_.update_limit)
    *hint = REBUILD_REASON_UPDATE_LIMIT_REACHED;

  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock = total_syntheticTick_ >= build_syntheticTick_;
  const bool performed_min_updates =
      simplex_info_.update_count >= synthetic_tick_reinversion_min_update_count;
  if (reinvert_syntheticClock && performed_min_updates)
    *hint = REBUILD_REASON_SYNTHETIC_CLOCK_SAYS_INVERT;

  analysis_.simplexTimerStop(UpdateFactorClock);
}

void HEkk::updatePivots(const int variable_in, const int row_out,
                        const int move_out) {
  analysis_.simplexTimerStart(UpdatePivotsClock);
  int variable_out = simplex_basis_.basicIndex_[row_out];

  // Incoming variable
  simplex_basis_.basicIndex_[row_out] = variable_in;
  simplex_basis_.nonbasicFlag_[variable_in] = 0;
  simplex_basis_.nonbasicMove_[variable_in] = 0;
  simplex_info_.baseLower_[row_out] = simplex_info_.workLower_[variable_in];
  simplex_info_.baseUpper_[row_out] = simplex_info_.workUpper_[variable_in];

  // Outgoing variable
  simplex_basis_.nonbasicFlag_[variable_out] = 1;
  if (simplex_info_.workLower_[variable_out] ==
      simplex_info_.workUpper_[variable_out]) {
    simplex_info_.workValue_[variable_out] =
        simplex_info_.workLower_[variable_out];
    simplex_basis_.nonbasicMove_[variable_out] = 0;
  } else if (move_out == -1) {
    simplex_info_.workValue_[variable_out] =
        simplex_info_.workLower_[variable_out];
    simplex_basis_.nonbasicMove_[variable_out] = 1;
  } else {
    simplex_info_.workValue_[variable_out] =
        simplex_info_.workUpper_[variable_out];
    simplex_basis_.nonbasicMove_[variable_out] = -1;
  }
  // Update the dual objective value
  double nwValue = simplex_info_.workValue_[variable_out];
  double vrDual = simplex_info_.workDual_[variable_out];
  double dl_dual_objective_value = nwValue * vrDual;
  simplex_info_.updated_dual_objective_value += dl_dual_objective_value;
  simplex_info_.update_count++;
  // Update the number of basic logicals
  if (variable_out < simplex_lp_.numCol_) simplex_info_.num_basic_logicals++;
  if (variable_in < simplex_lp_.numCol_) simplex_info_.num_basic_logicals--;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  simplex_lp_status_.has_invert = false;
  simplex_lp_status_.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  simplex_lp_status_.has_fresh_rebuild = false;
  analysis_.simplexTimerStop(UpdatePivotsClock);
}

void HEkk::updateMatrix(const int variable_in, const int variable_out) {
  analysis_.simplexTimerStart(UpdateMatrixClock);
  matrix_.update(variable_in, variable_out);
  analysis_.simplexTimerStop(UpdateMatrixClock);
}

void HEkk::computeSimplexInfeasible() {
  computeSimplexPrimalInfeasible();
  computeSimplexDualInfeasible();
}

void HEkk::computeSimplexPrimalInfeasible() {
  // Computes num/max/sum of primal infeasibliities according to the
  // simplex bounds. This is used to determine optimality in dual
  // phase 1 and dual phase 2, albeit using different bounds in
  // workLower/Upper.
  analysis_.simplexTimerStart(ComputePrIfsClock);
  const double scaled_primal_feasibility_tolerance =
      options_.primal_feasibility_tolerance;
  int& num_primal_infeasibility = simplex_info_.num_primal_infeasibility;
  double& max_primal_infeasibility = simplex_info_.max_primal_infeasibility;
  double& sum_primal_infeasibility = simplex_info_.sum_primal_infeasibility;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (int i = 0; i < simplex_lp_.numCol_ + simplex_lp_.numRow_; i++) {
    if (simplex_basis_.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = simplex_info_.workValue_[i];
      double lower = simplex_info_.workLower_[i];
      double upper = simplex_info_.workUpper_[i];
      // @primal_infeasibility calculation
      double primal_infeasibility = 0;
      if (value < lower - scaled_primal_feasibility_tolerance) {
        primal_infeasibility = lower - value;
      } else if (value > upper + scaled_primal_feasibility_tolerance) {
        primal_infeasibility = value - upper;
      }
      if (primal_infeasibility > 0) {
        if (primal_infeasibility > scaled_primal_feasibility_tolerance)
          num_primal_infeasibility++;
        max_primal_infeasibility =
            std::max(primal_infeasibility, max_primal_infeasibility);
        sum_primal_infeasibility += primal_infeasibility;
      }
    }
  }
  for (int i = 0; i < simplex_lp_.numRow_; i++) {
    // Basic variable
    double value = simplex_info_.baseValue_[i];
    double lower = simplex_info_.baseLower_[i];
    double upper = simplex_info_.baseUpper_[i];
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (value < lower - scaled_primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + scaled_primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > scaled_primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputePrIfsClock);
}

void HEkk::computeSimplexDualInfeasible() {
  analysis_.simplexTimerStart(ComputeDuIfsClock);
  // Computes num/max/sum of dual infeasibilities in phase 1 and phase
  // 2 according to nonbasicMove. The bounds are only used to identify
  // free variables. Fixed variables are assumed to have
  // nonbasicMove=0 so that no dual infeasibility is counted for them.
  const double scaled_dual_feasibility_tolerance =
      options_.dual_feasibility_tolerance;
  int& num_dual_infeasibility = simplex_info_.num_dual_infeasibility;
  double& max_dual_infeasibility = simplex_info_.max_dual_infeasibility;
  double& sum_dual_infeasibility = simplex_info_.sum_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (int iCol = 0; iCol < simplex_lp_.numCol_ + simplex_lp_.numRow_; iCol++) {
    if (!simplex_basis_.nonbasicFlag_[iCol]) continue;
    // Nonbasic column
    const double dual = simplex_info_.workDual_[iCol];
    const double lower = simplex_info_.workLower_[iCol];
    const double upper = simplex_info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -simplex_basis_.nonbasicMove_[iCol] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputeDuIfsClock);
}

void HEkk::computeSimplexLpDualInfeasible() {
  // Compute num/max/sum of dual infeasibliities according to the
  // bounds of the simplex LP. Assumes that boxed variables have
  // primal variable at the bound corresponding to the sign of the
  // dual so should only be used in dual phase 1 - where it's only
  // used for reporting after rebuilds.
  const double scaled_dual_feasibility_tolerance =
      options_.dual_feasibility_tolerance;
  int& num_dual_infeasibility =
      analysis_.num_dual_phase_1_lp_dual_infeasibility;
  double& max_dual_infeasibility =
      analysis_.max_dual_phase_1_lp_dual_infeasibility;
  double& sum_dual_infeasibility =
      analysis_.sum_dual_phase_1_lp_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (int iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    int iVar = iCol;
    if (!simplex_basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = simplex_info_.workDual_[iVar];
    const double lower = simplex_lp_.colLower_[iCol];
    const double upper = simplex_lp_.colUpper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  for (int iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    int iVar = simplex_lp_.numCol_ + iRow;
    if (!simplex_basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic row
    const double dual = -simplex_info_.workDual_[iVar];
    const double lower = simplex_lp_.rowLower_[iRow];
    const double upper = simplex_lp_.rowUpper_[iRow];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
}

bool HEkk::sparseLoopStyle(const int count, const int dim, int& to_entry) {
  // Parameter to decide whether to use just the values in a HVector, or
  // use the indices of their nonzeros
  const double density_for_indexing = 0.4;
  const bool use_indices = count >= 0 && count < density_for_indexing * dim;
  if (use_indices) {
    to_entry = count;
  } else {
    to_entry = dim;
  }
  return use_indices;
}

void HEkk::invalidatePrimalMaxSumInfeasibilityRecord() {
  simplex_info_.max_primal_infeasibility = illegal_infeasibility_measure;
  simplex_info_.sum_primal_infeasibility = illegal_infeasibility_measure;
}

void HEkk::invalidatePrimalInfeasibilityRecord() {
  simplex_info_.num_primal_infeasibility = illegal_infeasibility_count;
  invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkk::invalidateDualMaxSumInfeasibilityRecord() {
  simplex_info_.max_dual_infeasibility = illegal_infeasibility_measure;
  simplex_info_.sum_dual_infeasibility = illegal_infeasibility_measure;
}

void HEkk::invalidateDualInfeasibilityRecord() {
  simplex_info_.num_dual_infeasibility = illegal_infeasibility_count;
  invalidateDualMaxSumInfeasibilityRecord();
}

bool HEkk::bailoutReturn() {
  if (solve_bailout_) {
    // If bailout has already been decided: check that it's for one of
    // these reasons
    assert(scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT ||
           scaled_model_status_ == HighsModelStatus::REACHED_ITERATION_LIMIT ||
           scaled_model_status_ ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  }
  return solve_bailout_;
}

bool HEkk::bailoutOnTimeIterations() {
  if (solve_bailout_) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT ||
           scaled_model_status_ == HighsModelStatus::REACHED_ITERATION_LIMIT ||
           scaled_model_status_ ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  } else if (timer_.readRunHighsClock() > options_.time_limit) {
    solve_bailout_ = true;
    scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
  } else if (iteration_count_ >= options_.simplex_iteration_limit) {
    solve_bailout_ = true;
    scaled_model_status_ = HighsModelStatus::REACHED_ITERATION_LIMIT;
  }
  return solve_bailout_;
}

HighsStatus HEkk::returnFromSolve(const HighsStatus return_status) {
  simplex_info_.valid_backtracking_basis_ = false;
  return return_status;
}

double HEkk::computeBasisCondition() {
  int solver_num_row = simplex_lp_.numRow_;
  int solver_num_col = simplex_lp_.numCol_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  const int* Astart = &simplex_lp_.Astart_[0];
  const double* Avalue = &simplex_lp_.Avalue_[0];
  // Compute the Hager condition number estimate for the basis matrix
  const double NoDensity = 1;
  bs_cond_x.resize(solver_num_row);
  bs_cond_y.resize(solver_num_row);
  bs_cond_z.resize(solver_num_row);
  bs_cond_w.resize(solver_num_row);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / solver_num_row;
  double norm_Binv;
  for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    double value = bs_cond_x[r_n];
    if (value) {
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = value;
      row_ep.count++;
    }
  }
  for (int ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor_.ftran(row_ep, NoDensity);
    // zeta = sign(y);
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_y[r_n] = row_ep.array[r_n];
      if (bs_cond_y[r_n] > 0)
        bs_cond_w[r_n] = 1.0;
      else if (bs_cond_y[r_n] < 0)
        bs_cond_w[r_n] = -1.0;
      else
        bs_cond_w[r_n] = 0.0;
    }
    // z=A'\zeta;
    row_ep.clear();
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      double value = bs_cond_w[r_n];
      if (value) {
        row_ep.index[row_ep.count] = r_n;
        row_ep.array[r_n] = value;
        row_ep.count++;
      }
    }
    row_ep.packFlag = false;
    factor_.btran(row_ep, NoDensity);
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    int argmax_z = -1;
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_z[r_n] = row_ep.array[r_n];
      double abs_z_v = fabs(bs_cond_z[r_n]);
      if (abs_z_v > norm_z) {
        norm_z = abs_z_v;
        argmax_z = r_n;
      }
      ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
      norm_Binv += fabs(bs_cond_y[r_n]);
    }
    if (norm_z <= ztx) break;
    // x = zeros(n,1);
    // x(fd_i) = 1;
    for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    int vr_n = simplex_basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (int el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += fabs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  return cond_B;
}

void HEkk::initialiseAnalysis() {
  analysis_.setup(simplex_lp_, options_, iteration_count_);
}
