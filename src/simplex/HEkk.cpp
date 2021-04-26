/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
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
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;

  simplex_lp_ = lp;
  // Shouldn't have to check the incoming LP since this is an internal
  // call, but it may be an LP that's set up internally with errors
  // :-) ...
  if (options_.highs_debug_level > kHighsDebugLevelNone) {
    // ... so, if debugging, check the LP.
    call_status = assessLp(simplex_lp_, options_);
    return_status = interpretCallStatus(call_status, return_status, "assessLp");
    if (return_status == HighsStatus::kError) return return_status;
  }
  initialiseForNewLp();
  return HighsStatus::kOk;
}

HighsStatus HEkk::solve() {
  initialiseAnalysis();
  if (analysis_.analyse_simplex_time)
    analysis_.simplexTimerStart(SimplexTotalClock);
  iteration_count_ = 0;
  if (initialiseForSolve() == HighsStatus::kError) return HighsStatus::kError;

  assert(status_.has_basis);
  assert(status_.has_invert);
  assert(status_.valid);
  if (model_status_ == HighsModelStatus::kOptimal)
    return HighsStatus::kOk;

  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  std::string algorithm;

  // Indicate that dual and primal rays are not known
  status_.has_dual_ray = false;
  status_.has_primal_ray = false;

  // Allow primal and dual perturbations in case a block on them is
  // hanging over from a previous call
  info_.allow_cost_perturbation = true;
  info_.allow_bound_perturbation = true;

  chooseSimplexStrategyThreads(options_, info_);
  HighsInt& simplex_strategy = info_.simplex_strategy;

  // Initial solve according to strategy
  if (simplex_strategy == kSimplexStrategyPrimal) {
    algorithm = "primal";
    reportSimplexPhaseIterations(options_.log_options, iteration_count_, info_,
                                 true);
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Using EKK primal simplex solver\n");
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkPrimal::solve");
  } else {
    algorithm = "dual";
    reportSimplexPhaseIterations(options_.log_options, iteration_count_, info_,
                                 true);
    // Solve, depending on the particular strategy
    if (simplex_strategy == kSimplexStrategyDualTasks) {
      highsLogUser(
          options_.log_options, HighsLogType::kInfo,
          "Using EKK parallel dual simplex solver - SIP with %" HIGHSINT_FORMAT
          " threads\n",
          info_.num_threads);
    } else if (simplex_strategy == kSimplexStrategyDualMulti) {
      highsLogUser(
          options_.log_options, HighsLogType::kInfo,
          "Using EKK parallel dual simplex solver - PAMI with %" HIGHSINT_FORMAT
          " threads\n",
          info_.num_threads);
    } else {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Using EKK dual simplex solver - serial\n");
    }
    HEkkDual dual_solver(*this);
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkDual::solve");

    // Dual simplex solver may set model_status to be
    // kUnboundedOrInfeasible, and Highs::run() may not allow that to
    // be returned, so use primal simplex to distinguish
    if (model_status_ == HighsModelStatus::kUnboundedOrInfeasible &&
        !options_.allow_unbounded_or_infeasible) {
      HEkkPrimal primal_solver(*this);
      call_status = primal_solver.solve();
      assert(called_return_from_solve_);
      return_status =
          interpretCallStatus(call_status, return_status, "HEkkPrimal::solve");
    }
  }
  reportSimplexPhaseIterations(options_.log_options, iteration_count_, info_);
  if (return_status == HighsStatus::kError) return return_status;
  highsLogDev(options_.log_options, HighsLogType::kInfo,
              "EKK %s simplex solver returns %" HIGHSINT_FORMAT
              " primal and %" HIGHSINT_FORMAT
              " dual infeasibilities: "
              "Status %s\n",
              algorithm.c_str(), info_.num_primal_infeasibility,
              info_.num_dual_infeasibility,
              utilModelStatusToString(model_status_).c_str());
  // Can model_status_ = HighsModelStatus::kNotset be returned?
  assert(model_status_ != HighsModelStatus::kNotset);

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
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  if (info_.num_primal_infeasibility) {
    // Primal infeasibilities, so should be just dual phase 2
    assert(!info_.num_dual_infeasibility);
    // Use dual simplex (phase 2) with Devex pricing and no perturbation
    info_.simplex_strategy = kSimplexStrategyDual;
    info_.dual_simplex_cost_perturbation_multiplier = 0;
    info_.dual_edge_weight_strategy = kSimplexDualEdgeWeightStrategyDevex;
    HEkkDual dual_solver(*this);
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkDual::solve");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // Dual infeasibilities, so should be just primal phase 2
    assert(!info_.num_primal_infeasibility);
    // Use primal simplex (phase 2) with no perturbation
    info_.simplex_strategy = kSimplexStrategyPrimal;
    info_.primal_simplex_bound_perturbation_multiplier = 0;
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(call_status, return_status, "HEkkPrimal::solve");
    if (return_status == HighsStatus::kError) return return_status;
  }
  return return_status;
}

HighsStatus HEkk::setBasis() {
  // Set up nonbasicFlag and basicIndex for a logical basis
  const HighsInt num_col = simplex_lp_.numCol_;
  const HighsInt num_row = simplex_lp_.numRow_;
  const HighsInt num_tot = num_col + num_row;
  basis_.nonbasicFlag_.resize(num_tot);
  basis_.nonbasicMove_.resize(num_tot);
  basis_.basicIndex_.resize(num_row);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    basis_.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
    double lower = simplex_lp_.colLower_[iCol];
    double upper = simplex_lp_.colUpper_[iCol];
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed. Set to bound of LP that is closer to
        // zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iCol] = move;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
    basis_.basicIndex_[iRow] = iVar;
  }
  info_.num_basic_logicals = num_row;
  status_.has_basis = true;
  return HighsStatus::kOk;
}

HighsStatus HEkk::setBasis(const HighsBasis& highs_basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  if (debugBasisConsistent(options_, simplex_lp_, highs_basis) ==
      HighsDebugStatus::kLogicalError) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Supposed to be a Highs basis, but not valid\n");
    return HighsStatus::kError;
  }
  HighsInt num_col = simplex_lp_.numCol_;
  HighsInt num_row = simplex_lp_.numRow_;
  HighsInt num_tot = num_col + num_row;
  // Resize the basis in case none has yet been defined for this LP
  basis_.nonbasicFlag_.resize(num_tot);
  basis_.nonbasicMove_.resize(num_tot);
  basis_.basicIndex_.resize(num_row);

  HighsInt num_basic_variables = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;

    const double lower = simplex_lp_.colLower_[iCol];
    const double upper = simplex_lp_.colUpper_[iCol];
    if (highs_basis.col_status[iCol] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (highs_basis.col_status[iCol] == HighsBasisStatus::kLower) {
        if (lower == upper) {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
        } else {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
        }
      } else if (highs_basis.col_status[iCol] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
      } else {
        assert(highs_basis.col_status[iCol] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = simplex_lp_.rowLower_[iRow];
    const double upper = simplex_lp_.rowUpper_[iRow];
    if (highs_basis.row_status[iRow] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (highs_basis.row_status[iRow] == HighsBasisStatus::kLower) {
        if (lower == upper) {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
        } else {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
        }
      } else if (highs_basis.row_status[iRow] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
      } else {
        assert(highs_basis.row_status[iRow] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  status_.has_basis = true;
  return HighsStatus::kOk;
}

HighsStatus HEkk::setBasis(const SimplexBasis& basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  if (debugBasisConsistent(options_, simplex_lp_, basis) ==
      HighsDebugStatus::kLogicalError) {
    highsLogUser(options_.log_options, HighsLogType::kError,
                 "Supposed to be a Highs basis, but not valid\n");
    return HighsStatus::kError;
  }
  basis_.nonbasicFlag_ = basis.nonbasicFlag_;
  basis_.nonbasicMove_ = basis.nonbasicMove_;
  basis_.basicIndex_ = basis.basicIndex_;
  status_.has_basis = true;
  return HighsStatus::kOk;
}

HighsSolution HEkk::getSolution() {
  HighsSolution solution;
  // Scatter the basic primal values
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++)
    info_.workValue_[basis_.basicIndex_[iRow]] = info_.baseValue_[iRow];
  // Zero the basic dual values
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++)
    info_.workDual_[basis_.basicIndex_[iRow]] = 0;

  // Now we can get the solution
  solution.col_value.resize(simplex_lp_.numCol_);
  solution.col_dual.resize(simplex_lp_.numCol_);
  solution.row_value.resize(simplex_lp_.numRow_);
  solution.row_dual.resize(simplex_lp_.numRow_);

  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    solution.col_value[iCol] = info_.workValue_[iCol];
    solution.col_dual[iCol] =
        (HighsInt)simplex_lp_.sense_ * info_.workDual_[iCol];
  }
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    solution.row_value[iRow] = -info_.workValue_[simplex_lp_.numCol_ + iRow];
    solution.row_dual[iRow] = (HighsInt)simplex_lp_.sense_ *
                              info_.workDual_[simplex_lp_.numCol_ + iRow];
  }
  return solution;
}

HighsBasis HEkk::getHighsBasis() {
  HighsInt num_col = simplex_lp_.numCol_;
  HighsInt num_row = simplex_lp_.numRow_;
  HighsBasis highs_basis;
  highs_basis.col_status.resize(num_col);
  highs_basis.row_status.resize(num_row);
  assert(status_.has_basis);
  highs_basis.valid_ = false;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;
    const double lower = simplex_lp_.colLower_[iCol];
    const double upper = simplex_lp_.colUpper_[iCol];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::kLower;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.col_status[iCol] = basis_status;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = simplex_lp_.rowLower_[iRow];
    const double upper = simplex_lp_.rowUpper_[iRow];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::kLower;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.row_status[iRow] = basis_status;
  }
  highs_basis.valid_ = true;
  return highs_basis;
}

HighsInt HEkk::initialiseSimplexLpBasisAndFactor(
    const bool only_from_known_basis) {
  // If there's no basis, return error if the basis has to be known,
  // otherwise set a logical basis
  if (!status_.has_basis) {
    if (only_from_known_basis) {
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Simplex basis should be known but isn't\n");
      return -(HighsInt)HighsStatus::kError;
    }
    setBasis();
  }
  const HighsInt rank_deficiency = computeFactor();
  if (rank_deficiency) {
    // Basis is rank deficient
    if (only_from_known_basis) {
      // If only this basis should be used, then return error
      highsLogUser(options_.log_options, HighsLogType::kError,
                   "Supposed to be a full-rank basis, but incorrect\n");
      return rank_deficiency;
    }
    // Account for rank deficiency by correcing nonbasicFlag
    handleRankDeficiency();
    updateSimplexLpStatus(status_, LpAction::kNewBasis);
    setNonbasicMove();
    status_.has_basis = true;
    status_.has_invert = true;
    status_.has_fresh_invert = true;
  }
  assert(status_.has_invert);
  return 0;
}

void HEkk::handleRankDeficiency() {
  HighsInt rank_deficiency = factor_.rank_deficiency;
  vector<HighsInt>& noPvC = factor_.noPvC;
  vector<HighsInt>& noPvR = factor_.noPvR;
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt variable_in = simplex_lp_.numCol_ + noPvR[k];
    HighsInt variable_out = noPvC[k];
    basis_.nonbasicFlag_[variable_in] = kNonbasicFlagFalse;
    basis_.nonbasicFlag_[variable_out] = kNonbasicFlagTrue;
  }
  status_.has_matrix = false;
}

// Private methods

void HEkk::initialiseForNewLp() {
  setSimplexOptions();
  initialiseControl();
  initialiseSimplexLpRandomVectors();
  status_.initialised = true;
}

bool HEkk::isUnconstrainedLp() {
  bool is_unconstrained_lp = simplex_lp_.numRow_ <= 0;
  if (is_unconstrained_lp)
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "HEkkDual::solve called for LP with non-positive (%" HIGHSINT_FORMAT
        ") number of constraints\n",
        simplex_lp_.numRow_);
  assert(!is_unconstrained_lp);
  return is_unconstrained_lp;
}

HighsStatus HEkk::initialiseForSolve() {
  const HighsInt error_return = initialiseSimplexLpBasisAndFactor();
  assert(!error_return);
  if (error_return) return HighsStatus::kError;
  assert(status_.has_basis);

  updateSimplexOptions();
  initialiseMatrix();  // Timed
  allocateWorkAndBaseArrays();
  initialiseCost(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseBound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseNonbasicValueAndMove();
  computePrimal();                // Timed
  computeDual();                  // Timed
  computeSimplexInfeasible();     // Timed
  computeDualObjectiveValue();    // Timed
  computePrimalObjectiveValue();  // Timed
  status_.valid = true;

  bool primal_feasible = info_.num_primal_infeasibility == 0;
  bool dual_feasible = info_.num_dual_infeasibility == 0;
  model_status_ = HighsModelStatus::kNotset;
  if (primal_feasible && dual_feasible)
    model_status_ = HighsModelStatus::kOptimal;
  return HighsStatus::kOk;
}

void HEkk::setSimplexOptions() {
  // Copy values of HighsOptions for the simplex solver
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_edge_weight_strategy = options_.simplex_dual_edge_weight_strategy;
  info_.price_strategy = options_.simplex_price_strategy;
  info_.dual_simplex_cost_perturbation_multiplier =
      options_.dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_.primal_simplex_bound_perturbation_multiplier;
  info_.factor_pivot_threshold = options_.factor_pivot_threshold;
  info_.update_limit = options_.simplex_update_limit;
  random_.initialise(options_.highs_random_seed);

  // Set values of internal options
  info_.store_squared_primal_infeasibility = true;
}

void HEkk::updateSimplexOptions() {
  // Update some simplex option values from HighsOptions when
  // (re-)solving an LP. Others aren't changed because better values
  // may have been learned due to solving this LP (possibly with some
  // modification) before.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_simplex_cost_perturbation_multiplier =
      options_.dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_.primal_simplex_bound_perturbation_multiplier;
}

void HEkk::initialiseSimplexLpRandomVectors() {
  const HighsInt num_col = simplex_lp_.numCol_;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  if (!num_tot) return;
  // Instantiate and (re-)initialise the random number generator
  //  HighsRandom random;
  HighsRandom& random = random_;
  //  random.initialise();

  if (num_col) {
    // Generate a random permutation of the column indices
    vector<HighsInt>& numColPermutation = info_.numColPermutation_;
    numColPermutation.resize(num_col);
    for (HighsInt i = 0; i < num_col; i++) numColPermutation[i] = i;
    random.shuffle(numColPermutation.data(), num_col);
  }

  // Re-initialise the random number generator and generate the
  // random vectors in the same order as hsol to maintain repeatable
  // performance
  // random.initialise();

  // Generate a random permutation of all the indices
  vector<HighsInt>& numTotPermutation = info_.numTotPermutation_;
  numTotPermutation.resize(num_tot);
  for (HighsInt i = 0; i < num_tot; i++) numTotPermutation[i] = i;
  random.shuffle(numTotPermutation.data(), num_tot);

  // Generate a vector of random reals
  info_.numTotRandomValue_.resize(num_tot);
  vector<double>& numTotRandomValue = info_.numTotRandomValue_;
  for (HighsInt i = 0; i < num_tot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

void HEkk::chooseSimplexStrategyThreads(const HighsOptions& options,
                                        HighsSimplexInfo& info) {
  // Ensure that this is not called with an optimal basis
  assert(info.num_dual_infeasibility > 0 || info.num_primal_infeasibility > 0);
  // Set the internal simplex strategy and number of threads for dual
  // simplex
  HighsInt& simplex_strategy = info.simplex_strategy;
  // By default, use the HighsOptions strategy. If this is
  // kSimplexStrategyChoose, then the strategy used will depend on
  // whether the current basis is primal feasible.
  simplex_strategy = options.simplex_strategy;
  if (simplex_strategy == kSimplexStrategyChoose) {
    // HiGHS is left to choose the simplex strategy
    if (info.num_primal_infeasibility > 0) {
      // Not primal feasible, so use dual simplex
      simplex_strategy = kSimplexStrategyDual;
    } else {
      // Primal feasible. so use primal simplex
      simplex_strategy = kSimplexStrategyPrimal;
    }
  }
  // Set min/max_threads to correspond to serial code. They will be
  // set to other values if parallel options are used.
  info.min_threads = 1;
  info.max_threads = 1;
  // Record the min/max minimum number of HiGHS threads in the options
  const HighsInt highs_min_threads = options.highs_min_threads;
  const HighsInt highs_max_threads = options.highs_max_threads;
  HighsInt omp_max_threads = 0;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  if (options.parallel == kHighsOnString &&
      simplex_strategy == kSimplexStrategyDual) {
    // The parallel strategy is on and the simplex strategy is dual so use
    // PAMI if there are enough OMP threads
    if (omp_max_threads >= kDualMultiMinThreads)
      simplex_strategy = kSimplexStrategyDualMulti;
  }
  //
  // If parallel stratgies are used, the minimum number of HiGHS threads used
  // will be set to be at least the minimum required for the strategy
  //
  // All this is independent of the number of OMP threads available,
  // since code with multiple HiGHS threads can be run in serial.
#ifdef OPENMP
  if (simplex_strategy == kSimplexStrategyDualTasks) {
    info.min_threads = max(kDualTasksMinThreads, highs_min_threads);
    info.max_threads = max(info.min_threads, highs_max_threads);
  } else if (simplex_strategy == kSimplexStrategyDualMulti) {
    info.min_threads = max(kDualMultiMinThreads, highs_min_threads);
    info.max_threads = max(info.min_threads, highs_max_threads);
  }
#endif
  // Set the number of HiGHS threads to be used to be the maximum
  // number to be used
  info.num_threads = info.max_threads;
  // Give a warning if the number of threads to be used is fewer than
  // the minimum number of HiGHS threads allowed
  if (info.num_threads < highs_min_threads) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using %" HIGHSINT_FORMAT
                 " HiGHS threads for parallel strategy rather than "
                 "minimum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_threads, highs_min_threads);
  }
  // Give a warning if the number of threads to be used is more than
  // the maximum number of HiGHS threads allowed
  if (info.num_threads > highs_max_threads) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using %" HIGHSINT_FORMAT
                 " HiGHS threads for parallel strategy rather than "
                 "maximum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_threads, highs_max_threads);
  }
  // Give a warning if the number of threads to be used is fewer than
  // the number of OMP threads available
  if (info.num_threads > omp_max_threads) {
    highsLogUser(
        options.log_options, HighsLogType::kWarning,
        "Number of OMP threads available = %" HIGHSINT_FORMAT
        " < %" HIGHSINT_FORMAT
        " = Number of HiGHS threads "
        "to be used: Parallel performance will be less than anticipated\n",
        omp_max_threads, info.num_threads);
  }
}

bool HEkk::getNonsingularInverse(const HighsInt solve_phase) {
  assert(status_.has_basis);
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  // Take a copy of basicIndex from before INVERT to be used as the
  // saved ordering of basic variables - so reinvert will run
  // identically.
  const vector<HighsInt> basicIndex_before_compute_factor = basicIndex;
  // Save the number of updates performed in case it has to be used to determine
  // a limit
  const HighsInt simplex_update_count = info_.update_count;
  // Dual simplex edge weights are identified with rows, so must be
  // permuted according to INVERT. This must be done if workEdWt_ is
  // not NULL.
  const bool handle_edge_weights = workEdWt_ != NULL;
  // Scatter the edge weights so that, after INVERT, they can be
  // gathered according to the new permutation of basicIndex
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < simplex_lp_.numRow_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }

  // Call computeFactor to perform INVERT
  HighsInt rank_deficiency = computeFactor();
  const bool artificial_rank_deficiency = false;  //  true;//
  if (artificial_rank_deficiency) {
    if (!info_.phase1_backtracking_test_done && solve_phase == kSolvePhase1) {
      // Claim rank deficiency to test backtracking
      printf("Phase1 (Iter %" HIGHSINT_FORMAT
             ") Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      info_.phase1_backtracking_test_done = true;
    } else if (!info_.phase2_backtracking_test_done &&
               solve_phase == kSolvePhase2) {
      // Claim rank deficiency to test backtracking
      printf("Phase2 (Iter %" HIGHSINT_FORMAT
             ") Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      info_.phase2_backtracking_test_done = true;
    }
  }
  if (rank_deficiency) {
    // Rank deficient basis, so backtrack to last full rank basis
    //
    // Get the last nonsingular basis - so long as there is one
    if (!getBacktrackingBasis(workEdWtFull_)) return false;
    // Record that backtracking is taking place
    info_.backtracking_ = true;
    updateSimplexLpStatus(status_, LpAction::kBacktracking);
    HighsInt backtrack_rank_deficiency = computeFactor();
    // This basis has previously been inverted successfully, so it shouldn't be
    // singular
    if (backtrack_rank_deficiency) return false;
    // simplex update limit will be half of the number of updates
    // performed, so make sure that at least one update was performed
    if (simplex_update_count <= 1) return false;
    HighsInt use_simplex_update_limit = info_.update_limit;
    HighsInt new_simplex_update_limit = simplex_update_count / 2;
    info_.update_limit = new_simplex_update_limit;
    highsLogUser(options_.log_options, HighsLogType::kWarning,
                 "Rank deficiency of %" HIGHSINT_FORMAT
                 " after %" HIGHSINT_FORMAT
                 " simplex updates, so "
                 "backtracking: max updates reduced from %" HIGHSINT_FORMAT
                 " to %" HIGHSINT_FORMAT "\n",
                 rank_deficiency, simplex_update_count,
                 use_simplex_update_limit, new_simplex_update_limit);
  } else {
    // Current basis is full rank so save it
    putBacktrackingBasis(basicIndex_before_compute_factor, workEdWtFull_);
    // Indicate that backtracking is not taking place
    info_.backtracking_ = false;
    // Reset the update limit in case this is the first successful
    // inversion after backtracking
    info_.update_limit = options_.simplex_update_limit;
  }
  if (handle_edge_weights) {
    // Gather the edge weights according to the permutation of
    // basicIndex after INVERT
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < simplex_lp_.numRow_; i++)
      workEdWt_[i] = workEdWtFull_[basicIndex[i]];
    analysis_.simplexTimerStop(PermWtClock);
  }
  return true;
}

bool HEkk::getBacktrackingBasis(double* scattered_edge_weights) {
  if (!info_.valid_backtracking_basis_) return false;
  basis_ = info_.backtracking_basis_;
  info_.costs_perturbed = info_.backtracking_basis_costs_perturbed_;
  info_.workShift_ = info_.backtracking_basis_workShift_;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (HighsInt iVar = 0; iVar < num_tot; iVar++)
      scattered_edge_weights[iVar] =
          info_.backtracking_basis_edge_weights_[iVar];
  }
  return true;
}

void HEkk::putBacktrackingBasis() {
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  const bool handle_edge_weights = workEdWt_ != NULL;
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < simplex_lp_.numRow_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }
  putBacktrackingBasis(basicIndex, workEdWtFull_);
}

void HEkk::putBacktrackingBasis(
    const vector<HighsInt>& basicIndex_before_compute_factor,
    double* scattered_edge_weights) {
  info_.valid_backtracking_basis_ = true;
  info_.backtracking_basis_ = basis_;
  info_.backtracking_basis_.basicIndex_ = basicIndex_before_compute_factor;
  info_.backtracking_basis_costs_perturbed_ = info_.costs_perturbed;
  info_.backtracking_basis_workShift_ = info_.workShift_;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (HighsInt iVar = 0; iVar < num_tot; iVar++)
      info_.backtracking_basis_edge_weights_[iVar] =
          scattered_edge_weights[iVar];
  }
}

void HEkk::computePrimalObjectiveValue() {
  analysis_.simplexTimerStart(ComputePrObjClock);
  info_.primal_objective_value = 0;
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    HighsInt iVar = basis_.basicIndex_[iRow];
    if (iVar < simplex_lp_.numCol_) {
      info_.primal_objective_value +=
          info_.baseValue_[iRow] * simplex_lp_.colCost_[iVar];
    }
  }
  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    if (basis_.nonbasicFlag_[iCol])
      info_.primal_objective_value +=
          info_.workValue_[iCol] * simplex_lp_.colCost_[iCol];
  }
  info_.primal_objective_value *= cost_scale_;
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  info_.primal_objective_value += simplex_lp_.offset_;
  // Now have primal objective value
  status_.has_primal_objective_value = true;
  analysis_.simplexTimerStop(ComputePrObjClock);
}

void HEkk::computeDualObjectiveValue(const HighsInt phase) {
  analysis_.simplexTimerStart(ComputeDuObjClock);
  info_.dual_objective_value = 0;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (basis_.nonbasicFlag_[iCol]) {
      const double term = info_.workValue_[iCol] * info_.workDual_[iCol];
      if (term) {
        info_.dual_objective_value +=
            info_.workValue_[iCol] * info_.workDual_[iCol];
      }
    }
  }
  info_.dual_objective_value *= cost_scale_;
  if (phase != 1) {
    // In phase 1 the dual objective has no objective
    // shift. Otherwise, if minimizing the shift is added. If
    // maximizing, workCost (and hence workDual) are negated, so the
    // shift is subtracted. Hence the shift is added according to the
    // sign implied by sense_
    info_.dual_objective_value +=
        ((HighsInt)simplex_lp_.sense_) * simplex_lp_.offset_;
  }
  // Now have dual objective value
  status_.has_dual_objective_value = true;
  analysis_.simplexTimerStop(ComputeDuObjClock);
}

HighsInt HEkk::computeFactor() {
  if (!status_.has_factor_arrays) {
    // todo @ Julian: this fails on glass4
    assert(info_.factor_pivot_threshold >= options_.factor_pivot_threshold);
    factor_.setup(simplex_lp_.numCol_, simplex_lp_.numRow_,
                  &simplex_lp_.Astart_[0], &simplex_lp_.Aindex_[0],
                  &simplex_lp_.Avalue_[0], &basis_.basicIndex_[0],
                  info_.factor_pivot_threshold, options_.factor_pivot_tolerance,
                  options_.highs_debug_level, options_.output_flag,
                  options_.log_file_stream, options_.log_to_console,
                  options_.log_dev_level);
    status_.has_factor_arrays = true;
  }
  analysis_.simplexTimerStart(InvertClock);
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  if (analysis_.analyse_factor_time) {
    HighsInt thread_id = 0;
#ifdef OPENMP
    thread_id = omp_get_thread_num();
#endif
    factor_timer_clock_pointer =
        analysis_.getThreadFactorTimerClockPtr(thread_id);
  }
  const HighsInt rank_deficiency = factor_.build(factor_timer_clock_pointer);
  if (analysis_.analyse_factor_data) analysis_.updateInvertFormData(factor_);

  const bool force = rank_deficiency;
  debugCheckInvert(options_, factor_, force);

  if (rank_deficiency) {
    // Have an invertible representation, but of B with column(s)
    // replacements due to singularity. So no (fresh) representation of
    // B^{-1}
    status_.has_invert = false;
    status_.has_fresh_invert = false;
  } else {
    // Now have a representation of B^{-1}, and it is fresh!
    status_.has_invert = true;
    status_.has_fresh_invert = true;
  }
  // Set the update count to zero since the corrected invertible
  // representation may be used for an initial basis. In any case the
  // number of updates shouldn't be positive
  info_.update_count = 0;

  analysis_.simplexTimerStop(InvertClock);
  return rank_deficiency;
}

void HEkk::initialiseMatrix() {
  if (!status_.has_matrix) {
    analysis_.simplexTimerStart(matrixSetupClock);
    matrix_.setup(simplex_lp_.numCol_, simplex_lp_.numRow_,
                  &simplex_lp_.Astart_[0], &simplex_lp_.Aindex_[0],
                  &simplex_lp_.Avalue_[0], &basis_.nonbasicFlag_[0]);
    status_.has_matrix = true;
    analysis_.simplexTimerStop(matrixSetupClock);
  }
}

void HEkk::setNonbasicMove() {
  const bool have_solution = false;
  // Don't have a simplex basis since nonbasicMove is not set up.

  // Assign nonbasicMove using as much information as is available
  double lower;
  double upper;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  basis_.nonbasicMove_.resize(num_tot);

  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    if (iVar < simplex_lp_.numCol_) {
      lower = simplex_lp_.colLower_[iVar];
      upper = simplex_lp_.colUpper_[iVar];
    } else {
      HighsInt iRow = iVar - simplex_lp_.numCol_;
      lower = -simplex_lp_.rowUpper_[iRow];
      upper = -simplex_lp_.rowLower_[iRow];
    }
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
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
          double value = info_.workValue_[iVar];
          if (value < midpoint) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
        // 2. Bound of original LP that is closer to zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
  }
}

void HEkk::allocateWorkAndBaseArrays() {
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  info_.workCost_.resize(num_tot);
  info_.workDual_.resize(num_tot);
  info_.workShift_.resize(num_tot);

  info_.workLower_.resize(num_tot);
  info_.workUpper_.resize(num_tot);
  info_.workRange_.resize(num_tot);
  info_.workValue_.resize(num_tot);
  info_.workLowerShift_.resize(num_tot);
  info_.workUpperShift_.resize(num_tot);

  // Feel that it should be possible to resize this with in dual
  // solver, and only if Devex is being used, but a pointer to it
  // needs to be set up when constructing HDual
  info_.devex_index_.resize(num_tot);

  info_.baseLower_.resize(simplex_lp_.numRow_);
  info_.baseUpper_.resize(simplex_lp_.numRow_);
  info_.baseValue_.resize(simplex_lp_.numRow_);
}

void HEkk::initialiseLpColBound() {
  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    info_.workLower_[iCol] = simplex_lp_.colLower_[iCol];
    info_.workUpper_[iCol] = simplex_lp_.colUpper_[iCol];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowBound() {
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    HighsInt iCol = simplex_lp_.numCol_ + iRow;
    info_.workLower_[iCol] = -simplex_lp_.rowUpper_[iRow];
    info_.workUpper_[iCol] = -simplex_lp_.rowLower_[iRow];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseCost(const SimplexAlgorithm algorithm,
                          const HighsInt solve_phase, const bool perturb) {
  // Copy the cost
  initialiseLpColCost();
  initialiseLpRowCost();
  info_.costs_perturbed = 0;
  // Primal simplex costs are either from the LP or set specially in phase 1
  if (algorithm == SimplexAlgorithm::kPrimal) return;
  // Dual simplex costs are either from the LP or perturbed
  if (!perturb || info_.dual_simplex_cost_perturbation_multiplier == 0) return;
  // Perturb the original costs, scale down if is too big
  HighsInt num_original_nonzero_cost = 0;
  if (analysis_.analyse_simplex_data)
    printf("grep_DuPtrb: Cost perturbation for %s\n",
           simplex_lp_.model_name_.c_str());
  double bigc = 0;
  for (HighsInt i = 0; i < simplex_lp_.numCol_; i++) {
    const double abs_cost = fabs(info_.workCost_[i]);
    bigc = max(bigc, abs_cost);
    if (analysis_.analyse_simplex_data && abs_cost) num_original_nonzero_cost++;
  }
  const HighsInt pct0 = (100 * num_original_nonzero_cost) / simplex_lp_.numCol_;
  double average_cost = 0;
  if (analysis_.analyse_simplex_data) {
    if (num_original_nonzero_cost) {
      average_cost = bigc / num_original_nonzero_cost;
    } else {
      printf("grep_DuPtrb:    STRANGE initial workCost has non nonzeros\n");
    }
    printf("grep_DuPtrb:    Initially have %" HIGHSINT_FORMAT
           " nonzero costs (%3" HIGHSINT_FORMAT
           "%%) with bigc = "
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
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt i = 0; i < num_tot; i++)
    boxedRate += (info_.workRange_[i] < 1e30);
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
  for (HighsInt i = 0; i < simplex_lp_.numCol_; i++) {
    double lower = simplex_lp_.colLower_[i];
    double upper = simplex_lp_.colUpper_[i];
    double xpert = (fabs(info_.workCost_[i]) + 1) * base *
                   info_.dual_simplex_cost_perturbation_multiplier *
                   (1 + info_.numTotRandomValue_[i]);
    const double previous_cost = info_.workCost_[i];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free - no perturb
    } else if (upper >= kHighsInf) {  // Lower
      info_.workCost_[i] += xpert;
    } else if (lower <= -kHighsInf) {  // Upper
      info_.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      info_.workCost_[i] += (info_.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
    if (analysis_.analyse_simplex_data) {
      const double perturbation1 = fabs(info_.workCost_[i] - previous_cost);
      if (perturbation1)
        updateValueDistribution(perturbation1,
                                analysis_.cost_perturbation1_distribution);
    }
  }
  for (HighsInt i = simplex_lp_.numCol_; i < num_tot; i++) {
    double perturbation2 = (0.5 - info_.numTotRandomValue_[i]) *
                           info_.dual_simplex_cost_perturbation_multiplier *
                           1e-12;
    info_.workCost_[i] += perturbation2;
    if (analysis_.analyse_simplex_data) {
      perturbation2 = fabs(perturbation2);
      updateValueDistribution(perturbation2,
                              analysis_.cost_perturbation2_distribution);
    }
  }
  info_.costs_perturbed = 1;
}

void HEkk::initialiseBound(const SimplexAlgorithm algorithm,
                           const HighsInt solve_phase, const bool perturb) {
  initialiseLpColBound();
  initialiseLpRowBound();
  info_.bounds_perturbed = 0;
  // Primal simplex bounds are either from the LP or perturbed
  if (algorithm == SimplexAlgorithm::kPrimal) {
    if (!perturb || info_.primal_simplex_bound_perturbation_multiplier == 0)
      return;
    // Perturb the bounds
    // Determine the smallest and largest finite lower/upper bounds
    HighsInt num_col = simplex_lp_.numCol_;
    HighsInt num_row = simplex_lp_.numRow_;
    HighsInt num_tot = num_col + num_row;
    double min_abs_lower = kHighsInf;
    double max_abs_lower = -1;
    double min_abs_upper = kHighsInf;
    double max_abs_upper = -1;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double abs_lower = fabs(info_.workLower_[iVar]);
      double abs_upper = fabs(info_.workUpper_[iVar]);
      if (abs_lower && abs_lower < kHighsInf) {
        min_abs_lower = min(abs_lower, min_abs_lower);
        max_abs_lower = max(abs_lower, max_abs_lower);
      }
      if (abs_upper && abs_upper < kHighsInf) {
        min_abs_upper = min(abs_upper, min_abs_upper);
        max_abs_upper = max(abs_upper, max_abs_upper);
      }
    }
    // printf(
    //     "Nonzero finite lower bounds in [%9.4g, %9.4g]; upper bounds in "
    //     "[%9.4g, %9.4g]\n",
    //     min_abs_lower, max_abs_lower, min_abs_upper, max_abs_upper);

    const double base =
        info_.primal_simplex_bound_perturbation_multiplier * 5e-7;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double lower = info_.workLower_[iVar];
      double upper = info_.workUpper_[iVar];
      const bool fixed = lower == upper;
      // Don't perturb bounds of nonbasic fixed variables as they stay nonbasic
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagTrue && fixed) continue;
      double random_value = info_.numTotRandomValue_[iVar];
      if (lower > -kHighsInf) {
        if (lower < -1) {
          lower -= random_value * base * (-lower);
        } else if (lower < 1) {
          lower -= random_value * base;
        } else {
          lower -= random_value * base * lower;
        }
        info_.workLower_[iVar] = lower;
      }
      if (upper < kHighsInf) {
        if (upper < -1) {
          upper += random_value * base * (-upper);
        } else if (upper < 1) {
          upper += random_value * base;
        } else {
          upper += random_value * base * upper;
        }
        info_.workUpper_[iVar] = upper;
      }
      info_.workRange_[iVar] = info_.workUpper_[iVar] - info_.workLower_[iVar];
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
      // Set values of nonbasic variables
      if (basis_.nonbasicMove_[iVar] > 0) {
        info_.workValue_[iVar] = lower;
      } else if (basis_.nonbasicMove_[iVar] < 0) {
        info_.workValue_[iVar] = upper;
      }
    }
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iVar = basis_.basicIndex_[iRow];
      info_.baseLower_[iRow] = info_.workLower_[iVar];
      info_.baseUpper_[iRow] = info_.workUpper_[iVar];
    }
    info_.bounds_perturbed = 1;
    return;
  }
  // Dual simplex bounds are either from the LP or set to special values in
  // phase
  // 1
  assert(algorithm == SimplexAlgorithm::kDual);
  if (solve_phase == kSolvePhase2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = kHighsInf;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (info_.workLower_[iCol] == -inf && info_.workUpper_[iCol] == inf) {
      // Don't change for row variables: they should never become
      // nonbasic when starting from a logical basis, and no crash
      // should make a free row nonbasic, but could an advanced basis
      // make a free row nonbasic.
      // But what it it happened?
      if (iCol >= simplex_lp_.numCol_) continue;
      info_.workLower_[iCol] = -1000,
      info_.workUpper_[iCol] = 1000;  // FREE
    } else if (info_.workLower_[iCol] == -inf) {
      info_.workLower_[iCol] = -1,
      info_.workUpper_[iCol] = 0;  // UPPER
    } else if (info_.workUpper_[iCol] == inf) {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 1;  // LOWER
    } else {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 0;  // BOXED or FIXED
    }
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
  }
}

void HEkk::initialiseLpColCost() {
  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    info_.workCost_[iCol] =
        (HighsInt)simplex_lp_.sense_ * simplex_lp_.colCost_[iCol];
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowCost() {
  for (HighsInt iCol = simplex_lp_.numCol_;
       iCol < simplex_lp_.numCol_ + simplex_lp_.numRow_; iCol++) {
    info_.workCost_[iCol] = 0;
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseNonbasicValueAndMove() {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    const double lower = info_.workLower_[iVar];
    const double upper = info_.workUpper_[iVar];
    const HighsInt original_move = basis_.nonbasicMove_[iVar];
    double value;
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      value = lower;
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (original_move == kNonbasicMoveUp) {
          // Set at lower
          value = lower;
          move = kNonbasicMoveUp;
        } else if (original_move == kNonbasicMoveDn) {
          // Set at upper
          value = upper;
          move = kNonbasicMoveDn;
        } else {
          // Invalid nonbasicMove: correct and set value at lower
          value = lower;
          move = kNonbasicMoveUp;
        }
      } else {
        // Lower
        value = lower;
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      value = upper;
      move = kNonbasicMoveDn;
    } else {
      // FREE
      value = 0;
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
    info_.workValue_[iVar] = value;
  }
}

void HEkk::pivotColumnFtran(const HighsInt iCol, HVector& col_aq) {
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
  HighsInt num_row = simplex_lp_.numRow_;
  const double local_col_aq_density = (double)col_aq.count / num_row;
  analysis_.updateOperationResultDensity(local_col_aq_density,
                                         analysis_.col_aq_density);
  updateOperationResultDensity(local_col_aq_density, info_.col_aq_density);
  analysis_.simplexTimerStop(FtranClock);
}

void HEkk::unitBtran(const HighsInt iRow, HVector& row_ep) {
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
  HighsInt num_row = simplex_lp_.numRow_;
  const double local_row_ep_density = (double)row_ep.count / num_row;
  analysis_.updateOperationResultDensity(local_row_ep_density,
                                         analysis_.row_ep_density);
  updateOperationResultDensity(local_row_ep_density, info_.row_ep_density);
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
  updateOperationResultDensity(local_dual_col_density, info_.dual_col_density);
  analysis_.simplexTimerStop(BtranFullClock);
}

void HEkk::choosePriceTechnique(const HighsInt price_strategy,
                                const double row_ep_density,
                                bool& use_col_price,
                                bool& use_row_price_w_switch) {
  // By default switch to column PRICE when pi_p has at least this
  // density
  const double density_for_column_price_switch = 0.75;
  use_col_price = (price_strategy == kSimplexPriceStrategyCol) ||
                  (price_strategy == kSimplexPriceStrategyRowSwitchColSwitch &&
                   row_ep_density > density_for_column_price_switch);
  use_row_price_w_switch =
      price_strategy == kSimplexPriceStrategyRowSwitch ||
      price_strategy == kSimplexPriceStrategyRowSwitchColSwitch;
}

void HEkk::tableauRowPrice(const HVector& row_ep, HVector& row_ap) {
  analysis_.simplexTimerStart(PriceClock);
  const HighsInt solver_num_row = simplex_lp_.numRow_;
  const HighsInt solver_num_col = simplex_lp_.numCol_;
  const double local_density = 1.0 * row_ep.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  choosePriceTechnique(info_.price_strategy, local_density, use_col_price,
                       use_row_price_w_switch);
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
    const int8_t* nonbasicFlag = &basis_.nonbasicFlag_[0];
    for (HighsInt iCol = 0; iCol < solver_num_col; iCol++)
      row_ap.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_ap density
  const double local_row_ap_density = (double)row_ap.count / solver_num_col;
  analysis_.updateOperationResultDensity(local_row_ap_density,
                                         analysis_.row_ap_density);
  updateOperationResultDensity(local_row_ap_density, info_.row_ap_density);
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
  const HighsInt num_row = simplex_lp_.numRow_;
  const HighsInt num_col = simplex_lp_.numCol_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(num_row);
  primal_col.clear();
  for (HighsInt i = 0; i < num_col + num_row; i++) {
    if (basis_.nonbasicFlag_[i] && info_.workValue_[i] != 0) {
      matrix_.collect_aj(primal_col, i, info_.workValue_[i]);
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
                                 info_.primal_col_density);
  }
  for (HighsInt i = 0; i < num_row; i++) {
    HighsInt iCol = basis_.basicIndex_[i];
    info_.baseValue_[i] = -primal_col.array[i];
    info_.baseLower_[i] = info_.workLower_[iCol];
    info_.baseUpper_[i] = info_.workUpper_[iCol];
  }
  // Indicate that the primal infeasiblility information isn't known
  info_.num_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;

  // Now have basic primals
  status_.has_basic_primal_values = true;
  analysis_.simplexTimerStop(ComputePrimalClock);
}

void HEkk::computeDual() {
  analysis_.simplexTimerStart(ComputeDualClock);
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(simplex_lp_.numRow_);
  dual_col.clear();
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    const double value = info_.workCost_[basis_.basicIndex_[iRow]] +
                         info_.workShift_[basis_.basicIndex_[iRow]];
    if (value) {
      dual_col.index[dual_col.count++] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // Copy the costs in case the basic costs are all zero
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt i = 0; i < num_tot; i++)
    info_.workDual_[i] = info_.workCost_[i];

  if (dual_col.count) {
    fullBtran(dual_col);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(simplex_lp_.numCol_);
    fullPrice(dual_col, dual_row);
    for (HighsInt i = 0; i < simplex_lp_.numCol_; i++)
      info_.workDual_[i] -= dual_row.array[i];
    for (HighsInt i = simplex_lp_.numCol_; i < num_tot; i++)
      info_.workDual_[i] -= dual_col.array[i - simplex_lp_.numCol_];
  }
  // Indicate that the dual infeasiblility information isn't known
  info_.num_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;

  // Now have nonbasic duals
  status_.has_nonbasic_dual_values = true;
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

  HighsInt num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;

  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = info_.workLower_[iVar];
    const double upper = info_.workUpper_[iVar];
    const double dual = info_.workDual_[iVar];
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
      dual_infeasibility = -basis_.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  info_.num_dual_infeasibility = num_dual_infeasibility;
  info_.max_dual_infeasibility = max_dual_infeasibility;
  info_.sum_dual_infeasibility = sum_dual_infeasibility;
}

double HEkk::computeDualForTableauColumn(const HighsInt iVar,
                                         const HVector& tableau_column) {
  const vector<double>& workCost = info_.workCost_;
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;

  double dual = info_.workCost_[iVar];
  for (HighsInt i = 0; i < tableau_column.count; i++) {
    HighsInt iRow = tableau_column.index[i];
    dual -= tableau_column.array[iRow] * workCost[basicIndex[iRow]];
  }
  return dual;
}

bool HEkk::correctDual(HighsInt* free_infeasibility_count) {
  const double tau_d = options_.dual_feasibility_tolerance;
  const double inf = kHighsInf;
  HighsInt workCount = 0;
  double flip_dual_objective_value_change = 0;
  double shift_dual_objective_value_change = 0;
  HighsInt num_flip = 0;
  HighsInt num_shift = 0;
  double sum_flip = 0;
  double sum_shift = 0;
  HighsInt num_shift_skipped = 0;
  const HighsInt num_tot = simplex_lp_.numCol_ + simplex_lp_.numRow_;
  for (HighsInt i = 0; i < num_tot; i++) {
    if (basis_.nonbasicFlag_[i]) {
      if (info_.workLower_[i] == -inf && info_.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(info_.workDual_[i]) >= tau_d);
      } else if (basis_.nonbasicMove_[i] * info_.workDual_[i] <= -tau_d) {
        if (info_.workLower_[i] != -inf && info_.workUpper_[i] != inf) {
          // Boxed variable = flip
          const HighsInt move = basis_.nonbasicMove_[i];
          flipBound(i);
          double flip = info_.workUpper_[i] - info_.workLower_[i];
          // Negative dual at lower bound (move=1): flip to upper
          // bound so objective contribution is change in value (flip)
          // times dual, being move*flip*dual
          //
          // Positive dual at upper bound (move=-1): flip to lower
          // bound so objective contribution is change in value
          // (-flip) times dual, being move*flip*dual
          double local_dual_objective_change = move * flip * info_.workDual_[i];
          local_dual_objective_change *= cost_scale_;
          flip_dual_objective_value_change += local_dual_objective_change;
          num_flip++;
          sum_flip += fabs(flip);
        } else {
          if (info_.allow_cost_perturbation) {
            // Other variable = shift
            info_.costs_perturbed = 1;
            std::string direction;
            double shift;
            if (basis_.nonbasicMove_[i] == 1) {
              direction = "  up";
              double dual = (1 + random_.fraction()) * tau_d;
              shift = dual - info_.workDual_[i];
              info_.workDual_[i] = dual;
              info_.workCost_[i] = info_.workCost_[i] + shift;
            } else {
              direction = "down";
              double dual = -(1 + random_.fraction()) * tau_d;
              shift = dual - info_.workDual_[i];
              info_.workDual_[i] = dual;
              info_.workCost_[i] = info_.workCost_[i] + shift;
            }
            double local_dual_objective_change = shift * info_.workValue_[i];
            local_dual_objective_change *= cost_scale_;
            shift_dual_objective_value_change += local_dual_objective_change;
            num_shift++;
            sum_shift += fabs(shift);
            highsLogDev(options_.log_options, HighsLogType::kVerbose,
                        "Move %s: cost shift = %g; objective change = %g\n",
                        direction.c_str(), shift, local_dual_objective_change);
          } else {
            // Shifting not permitted
            //
            // Before 07/01/20, these shifts were always done, but
            // doing it after cost perturbation has been removed can
            // lead to cycling when dual unboundedness (=> primal
            // infeasibility) has been detecteed in Phase 2, since the
            // shift removes dual infeasibilities, which are then
            // reinstated after the dual values are recomputed.
            //
            // ToDo: Not shifting leads to dual infeasibilities when an
            // LP is declared to be infeasible. Should go to
            // phase 1 primal simplex to "prove" infeasibility.
            //
            num_shift_skipped++;
          }
        }
      }
    }
  }
  if (num_shift_skipped) {
    highsLogDev(options_.log_options, HighsLogType::kError,
                "correctDual: Missed %d cost shifts\n", num_shift_skipped);
    assert(!num_shift_skipped);
    return false;
  }
  if (num_flip)
    highsLogDev(options_.log_options, HighsLogType::kVerbose,
                "Performed %" HIGHSINT_FORMAT
                " flip(s): total = %g; objective change = %g\n",
                num_flip, sum_flip, flip_dual_objective_value_change);
  if (num_shift)
    highsLogDev(options_.log_options, HighsLogType::kDetailed,
                "Performed %" HIGHSINT_FORMAT
                " cost shift(s): total = %g; objective change = %g\n",
                num_shift, sum_shift, shift_dual_objective_value_change);
  *free_infeasibility_count = workCount;
  return true;
}

void HEkk::flipBound(const HighsInt iCol) {
  int8_t* nonbasicMove = &basis_.nonbasicMove_[0];
  const int8_t move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  info_.workValue_[iCol] =
      move == 1 ? info_.workLower_[iCol] : info_.workUpper_[iCol];
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
  const HighsInt update_count = info_.update_count;
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
    const double current_pivot_threshold = info_.factor_pivot_threshold;
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
      highsLogUser(options_.log_options, HighsLogType::kWarning,
                   "   Increasing Markowitz threshold to %g\n",
                   new_pivot_threshold);
      info_.factor_pivot_threshold = new_pivot_threshold;
      factor_.setPivotThreshold(new_pivot_threshold);
    }
  }
  return reinvert;
}

// The major model updates. Factor calls factor_.update; Matrix
// calls matrix_.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HEkk::updateFactor(HVector* column, HVector* row_ep, HighsInt* iRow,
                        HighsInt* hint) {
  analysis_.simplexTimerStart(UpdateFactorClock);
  factor_.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  status_.has_invert = true;
  if (info_.update_count >= info_.update_limit)
    *hint = kRebuildReasonUpdateLimitReached;

  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock = total_syntheticTick_ >= build_syntheticTick_;
  const bool performed_min_updates =
      info_.update_count >= synthetic_tick_reinversion_min_update_count;
  if (reinvert_syntheticClock && performed_min_updates)
    *hint = kRebuildReasonSyntheticClockSaysInvert;

  analysis_.simplexTimerStop(UpdateFactorClock);
}

void HEkk::updatePivots(const HighsInt variable_in, const HighsInt row_out,
                        const HighsInt move_out) {
  analysis_.simplexTimerStart(UpdatePivotsClock);
  HighsInt variable_out = basis_.basicIndex_[row_out];

  // Incoming variable
  basis_.basicIndex_[row_out] = variable_in;
  basis_.nonbasicFlag_[variable_in] = 0;
  basis_.nonbasicMove_[variable_in] = 0;
  info_.baseLower_[row_out] = info_.workLower_[variable_in];
  info_.baseUpper_[row_out] = info_.workUpper_[variable_in];

  // Outgoing variable
  basis_.nonbasicFlag_[variable_out] = 1;
  if (info_.workLower_[variable_out] == info_.workUpper_[variable_out]) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 0;
  } else if (move_out == -1) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 1;
  } else {
    info_.workValue_[variable_out] = info_.workUpper_[variable_out];
    basis_.nonbasicMove_[variable_out] = -1;
  }
  // Update the dual objective value
  double nwValue = info_.workValue_[variable_out];
  double vrDual = info_.workDual_[variable_out];
  double dl_dual_objective_value = nwValue * vrDual;
  info_.updated_dual_objective_value += dl_dual_objective_value;
  info_.update_count++;
  // Update the number of basic logicals
  if (variable_out < simplex_lp_.numCol_) info_.num_basic_logicals++;
  if (variable_in < simplex_lp_.numCol_) info_.num_basic_logicals--;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  status_.has_invert = false;
  status_.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  status_.has_fresh_rebuild = false;
  analysis_.simplexTimerStop(UpdatePivotsClock);
}

void HEkk::updateMatrix(const HighsInt variable_in,
                        const HighsInt variable_out) {
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
  HighsInt& num_primal_infeasibility = info_.num_primal_infeasibility;
  double& max_primal_infeasibility = info_.max_primal_infeasibility;
  double& sum_primal_infeasibility = info_.sum_primal_infeasibility;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (HighsInt i = 0; i < simplex_lp_.numCol_ + simplex_lp_.numRow_; i++) {
    if (basis_.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = info_.workValue_[i];
      double lower = info_.workLower_[i];
      double upper = info_.workUpper_[i];
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
  for (HighsInt i = 0; i < simplex_lp_.numRow_; i++) {
    // Basic variable
    double value = info_.baseValue_[i];
    double lower = info_.baseLower_[i];
    double upper = info_.baseUpper_[i];
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
  HighsInt& num_dual_infeasibility = info_.num_dual_infeasibility;
  double& max_dual_infeasibility = info_.max_dual_infeasibility;
  double& sum_dual_infeasibility = info_.sum_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_ + simplex_lp_.numRow_;
       iCol++) {
    if (!basis_.nonbasicFlag_[iCol]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iCol];
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -basis_.nonbasicMove_[iCol] * dual;
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
  HighsInt& num_dual_infeasibility =
      analysis_.num_dual_phase_1_lp_dual_infeasibility;
  double& max_dual_infeasibility =
      analysis_.max_dual_phase_1_lp_dual_infeasibility;
  double& sum_dual_infeasibility =
      analysis_.sum_dual_phase_1_lp_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < simplex_lp_.numCol_; iCol++) {
    HighsInt iVar = iCol;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iVar];
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
  for (HighsInt iRow = 0; iRow < simplex_lp_.numRow_; iRow++) {
    HighsInt iVar = simplex_lp_.numCol_ + iRow;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic row
    const double dual = -info_.workDual_[iVar];
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

bool HEkk::sparseLoopStyle(const HighsInt count, const HighsInt dim,
                           HighsInt& to_entry) {
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
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidatePrimalInfeasibilityRecord() {
  info_.num_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkk::invalidateDualMaxSumInfeasibilityRecord() {
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidateDualInfeasibilityRecord() {
  info_.num_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  invalidateDualMaxSumInfeasibilityRecord();
}

bool HEkk::bailoutOnTimeIterations() {
  if (solve_bailout_) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveCutoff);
  } else if (timer_.readRunHighsClock() > options_.time_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kTimeLimit;
  } else if (iteration_count_ >= options_.simplex_iteration_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kIterationLimit;
  }
  return solve_bailout_;
}

HighsStatus HEkk::returnFromSolve(const HighsStatus return_status) {
  // Always called before returning from HEkkPrimal/Dual::solve()
  if (solve_bailout_) {
    // If bailout has already been decided: check that it's for one of
    // these reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveCutoff);
  }
  // Check that returnFromSolve has not already been called: it should
  // be called exactly once per solve
  assert(!called_return_from_solve_);
  called_return_from_solve_ = true;
  info_.valid_backtracking_basis_ = false;

  // Initialise the status of the primal and dual solutions
  return_primal_solution_status = kHighsPrimalDualStatusUnknown;
  return_dual_solution_status = kHighsPrimalDualStatusUnknown;
  // Nothing more is known about the solve after an error return
  if (return_status == HighsStatus::kError) return return_status;

  // Check that an invert exists
  assert(status_.has_invert);

  // Determine a primal and possibly a dual solution, removing the
  // effects of perturbations and shifts
  switch (model_status_) {
    case HighsModelStatus::kOptimal: {
      return_primal_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      return_dual_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      break;
    }
    case HighsModelStatus::kInfeasible: {
      // Primal simplex has identified primal infeasibility in phase 1, or
      // dual simplex has identified dual unboundedness in phase 2. In
      // both cases there should be no primal or dual perturbations
      assert(!info_.costs_perturbed && !info_.bounds_perturbed);
      if (exit_algorithm == SimplexAlgorithm::kPrimal) {
        // Reset the simplex costs and recompute duals after primal
        // phase 1
        initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
        computeDual();
      }
      computeSimplexInfeasible();
      // Primal solution is valid, but shouldn't be feasible
      assert(info_.num_primal_infeasibility > 0);
      return_primal_solution_status = kHighsPrimalDualStatusInfeasiblePoint;
      // Dual solution is valid, but is unlikely to be feasible!
      if (info_.num_dual_infeasibility == 0) {
        return_dual_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      } else {
        return_dual_solution_status = kHighsPrimalDualStatusInfeasiblePoint;
      }
      break;
    }
    case HighsModelStatus::kUnboundedOrInfeasible: {
      // Dual simplex has identified dual infeasibility in phase
      // 1. There should be no dual perturbations
      assert(exit_algorithm == SimplexAlgorithm::kDual);
      assert(!info_.costs_perturbed);
      // Reset the simplex bounds and recompute primals
      initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
      computePrimal();
      computeSimplexInfeasible();
      // Primal solution is valid, but is unlikely to be feasible!
      if (info_.num_primal_infeasibility == 0) {
        return_primal_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      } else {
        return_primal_solution_status = kHighsPrimalDualStatusInfeasiblePoint;
      }
      // Dual solution is valid, but shouldn't be feasible
      assert(info_.num_dual_infeasibility > 0);
      return_dual_solution_status = kHighsPrimalDualStatusInfeasiblePoint;
      break;
    }
    case HighsModelStatus::kUnbounded: {
      // Primal simplex has identified unboundedness in phase 2. There
      // should be no primal or dual perturbations
      assert(exit_algorithm == SimplexAlgorithm::kPrimal);
      assert(!info_.costs_perturbed && !info_.bounds_perturbed);
      computeSimplexInfeasible();
      // Primal solution is valid, and should be feasible
      assert(info_.num_primal_infeasibility == 0);
      return_primal_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      // Dual solution is valid, but is unlikely to be feasible!
      if (info_.num_dual_infeasibility == 0) {
        return_dual_solution_status = kHighsPrimalDualStatusFeasiblePoint;
      } else {
        return_dual_solution_status = kHighsPrimalDualStatusInfeasiblePoint;
      }
      break;
    }
    case HighsModelStatus::kObjectiveCutoff:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit: {
      break;
    }
    default: {
      printf("What is default here? Status %s\n",
             utilModelStatusToString(model_status_).c_str());
      break;
    }
  }

  computePrimalObjectiveValue();

  return return_status;
}

double HEkk::computeBasisCondition() {
  HighsInt solver_num_row = simplex_lp_.numRow_;
  HighsInt solver_num_col = simplex_lp_.numCol_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  const HighsInt* Astart = &simplex_lp_.Astart_[0];
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
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    double value = bs_cond_x[r_n];
    if (value) {
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = value;
      row_ep.count++;
    }
  }
  for (HighsInt ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor_.ftran(row_ep, NoDensity);
    // zeta = sign(y);
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
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
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
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
    HighsInt argmax_z = -1;
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
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
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    HighsInt vr_n = basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (HighsInt el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
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
