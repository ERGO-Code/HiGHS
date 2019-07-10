/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplex.h"
#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsStatus.h"
#include "simplex/HCrash.h"
#include "simplex/HVector.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"  // For simplex strategy constants
#include "simplex/SimplexTimer.h"
#include "util/HighsUtils.h"

using std::runtime_error;
#include <cassert>
#include <vector>

void setSimplexOptions(HighsModelObject& highs_model_object) {
  HighsOptions& options = highs_model_object.options_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  //
  // Copy values of HighsOptions for the simplex solver
  //
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  simplex_info.simplex_strategy = options.simplex_strategy;
  simplex_info.dual_edge_weight_strategy =
      options.simplex_dual_edge_weight_strategy;
  simplex_info.price_strategy = options.simplex_price_strategy;
  simplex_info.primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  simplex_info.dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  simplex_info.perturb_costs = options.simplex_perturb_costs;
  simplex_info.update_limit = options.simplex_update_limit;

  // Set values of internal options
  simplex_info.allow_primal_flips_for_dual_feasibility = true;
  if (options.run_as_hsol) simplex_info.allow_primal_flips_for_dual_feasibility = true;
  // Option for analysing the LP solution
  simplex_info.analyseLpSolution = true;
#ifdef HiGHSDEV
  bool useful_analysis = false;
  bool full_timing = false;
  // Options for reporting timing
  simplex_info.report_simplex_inner_clock = useful_analysis;
  simplex_info.report_simplex_outer_clock = full_timing;
  simplex_info.report_simplex_phases_clock = full_timing;
  // Options for analysing the LP and simplex iterations
  simplex_info.analyseLp = useful_analysis;
  simplex_info.analyseSimplexIterations = useful_analysis;
  simplex_info.analyse_invert_time = full_timing;
  simplex_info.analyseRebuildTime = full_timing;
#endif
}

SimplexSolutionStatus transition(HighsModelObject& highs_model_object) {
  // Perform the transition from whatever information is known about
  // the LP to a status where simplex data are set up for the initial
  // rebuild() of the chosen solver - primal, scalar dual or parallel
  // dual.
  //
  // First look at what basis and solution information is known. If a
  // simplex basis is known, then it's used, and there's no need for
  // the solution values. This will usually correspond to hot start
  // when solving MIP problems. Otherwise, generate a simplex basis,
  // thus:
  //
  // If there is a HiGHS basis: use it to determine what's basic and nonbasic
  // (nonbasicFlag).
  //
  // If there's no HiGHS basis: generate nonbasicFlag, possibly by dualising and
  // performing a crash.
  //
  // Use nonbasicFlag to generate basicIndex
  //
  // Use nonbasicFlag and any HiGHS solution to determine nonbasicMove
  //
  HighsTimer& timer = highs_model_object.timer_;
  const HighsOptions& options = highs_model_object.options_;
  const HighsSolution& solution = highs_model_object.solution_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HFactor& factor = highs_model_object.factor_;
  HMatrix& matrix = highs_model_object.matrix_;
  // First determine whether the HiGHS solution space has been
  // allocated, a necessary condition for its values to be used later
  bool have_highs_solution =
    (int)solution.col_value.size() == highs_model_object.lp_.numCol_ &&
    (int)solution.col_dual.size() == highs_model_object.lp_.numCol_ &&
    (int)solution.row_value.size() == highs_model_object.lp_.numRow_ &&
    (int)solution.row_dual.size() == highs_model_object.lp_.numRow_;
  if (!simplex_lp_status.valid) {
    // Simplex LP is not valid, so ensure that it is fully invalidated
    invalidateSimplexLp(simplex_lp_status);
    // Copy the LP to the structure to be used by the solver
    simplex_lp = highs_model_object.lp_;
    // Initialise the real and integer random vectors
    initialiseSimplexLpRandomVectors(highs_model_object);
  }
  if (simplex_lp_status.has_basis) {
    // There is a simplex basis: it should be valid - since it's set internally - but check
    int num_basic_variables = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) {
	num_basic_variables++;
      } else {
	simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_TRUE;
      }
    }
    assert(num_basic_variables == simplex_lp.numRow_);
    if (num_basic_variables != simplex_lp.numRow_) simplex_lp_status.has_basis = false;
  }
  // Now we know whether the simplex basis at least has the right number
  // of basic and nonbasic variables
  if (!simplex_lp_status.has_basis) {
    // There is no simplex basis (or it was found to be invalid) so try to identify one
    if (basis.valid_) {
      // There is is HiGHS basis: use it to construct nonbasicFlag,
      // checking that it has the right number of basic variables
      //
      // Allocate memory for nonbasicFlag
      simplex_basis.nonbasicFlag_.resize(highs_model_object.lp_.numCol_ +
                                         highs_model_object.lp_.numRow_);
      int num_basic_variables = 0;
      for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
        int iVar = iCol;
        if (basis.col_status[iCol] == HighsBasisStatus::BASIC) {
          num_basic_variables++;
          simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
        } else {
          simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_TRUE;
        }
      }
      for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
        int iVar = simplex_lp.numCol_ + iRow;
        if (basis.row_status[iRow] == HighsBasisStatus::BASIC) {
          num_basic_variables++;
          simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_FALSE;
        } else {
          simplex_basis.nonbasicFlag_[iVar] = NONBASIC_FLAG_TRUE;
        }
      }
      assert(num_basic_variables == simplex_lp.numRow_);
      if (num_basic_variables != simplex_lp.numRow_) basis.valid_ = false;
    }
    // nonbasicFlag is valid if the HiGHS basis exists and has the correct
    // number of basic variables
    bool nonbasicFlag_valid = basis.valid_;
    if (!nonbasicFlag_valid) {
      // So, nonbasicFlag is not valid - either because there is no
      // simplex or HiGHS basis, or because what was claimed to be
      // valid has been found to have the wrong number of basic and
      // nonbasic variables
      //
      // This is taken to imply that this is a "new" LP to be solved, so
      //
      // 1. Set simplex options from HiGHS options. This is only done with a new LP so that strategy
      // and knowledge based on run-time experience with the same LP should be preserved.
      //      setSimplexOptions(highs_model_object);
      //
      // 2. Initialise the simplex timing 
      //      SimplexTimer simplex_timer; simplex_timer.initialiseSimplexClocks(highs_model_object);
      //
      // 3. Generate a simplex basis, possibly by performing a crash,
      // and possibly after dualising

      /*
      // Possibly dualise, making sure that no simplex or other data are used to
      initialise if (options.simplex_dualise_strategy !=
      SimplexDualiseStrategy::OFF) { dualiseSimplexLp(highs_model_object);
      have_highs_solution = false;
        // Initialise the real and integer random vectors
        initialiseSimplexLpRandomVectors(highs_model_object);
      }
      */
      // Possibly permute the columns of the LP to be used by the solver.
      if (options.simplex_permute_strategy != SimplexPermuteStrategy::OFF)
        permuteSimplexLp(highs_model_object);

      // Allocate memory for nonbasicFlag
      simplex_basis.nonbasicFlag_.resize(highs_model_object.lp_.numCol_ +
                                         highs_model_object.lp_.numRow_);
      // Set up nonbasicFlag for a logical basis
      for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++)
        simplex_basis.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
      for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
        simplex_basis.nonbasicFlag_[simplex_lp.numCol_ + iRow] =
            NONBASIC_FLAG_FALSE;

      // Possibly find a crash basis
      if (options.simplex_crash_strategy != SimplexCrashStrategy::OFF) {
        HCrash crash(highs_model_object);
        timer.start(simplex_info.clock_[CrashClock]);
        crash.crash(options.simplex_crash_strategy);
        timer.stop(simplex_info.clock_[CrashClock]);
        int num_basic_structurals = 0;
        for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
          if (simplex_basis.nonbasicFlag_[iCol] == NONBASIC_FLAG_FALSE)
            num_basic_structurals++;
        }
        HighsLogMessage(HighsMessageType::INFO,
                        "Crash has created a basis with %d/%d structurals",
                        num_basic_structurals, simplex_lp.numRow_);
      }
    }
    // Now that the dimensions of the LP to be solved by the simplex
    // method are known, make sure that there is a postive number of
    // rows. ToDo: Ensure that LPs with no rows can still be solved
    assert(simplex_lp.numRow_ > 0);
    if (simplex_lp.numRow_ == 0) {
      printf(
          "Cannot currently solve LPs with no rows using the simplex method\n");
      return SimplexSolutionStatus::FAILED;
    }

    // There is now a nonbasicFlag that should be valid - have the
    // right number of basic variables - so check this
    int num_basic_variables = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE)
        num_basic_variables++;
    }
    nonbasicFlag_valid = num_basic_variables == simplex_lp.numRow_;
    assert(nonbasicFlag_valid);
    if (!nonbasicFlag_valid) {
      // Something's gone wrong: any HiGHS basis has been checked and,
      // if there isn't one or it's been found to be invalid, a
      // logical or crash basis has been set up. Both should guarantee
      // the right number of basic variables
      for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++)
        simplex_basis.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
      for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
        simplex_basis.nonbasicFlag_[simplex_lp.numCol_ + iRow] =
            NONBASIC_FLAG_FALSE;
      nonbasicFlag_valid = true;
      // The HiGHS basis shouldn't be valid at this point
      assert(!basis.valid_);
    }
    // Use nonbasicFlag to form basicIndex
    // Allocate memory for basicIndex
    simplex_basis.basicIndex_.resize(highs_model_object.lp_.numRow_);
    num_basic_variables = 0;
    simplex_info.num_basic_logicals = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) {
        simplex_basis.basicIndex_[num_basic_variables] = iVar;
        if (iVar >= simplex_lp.numCol_) simplex_info.num_basic_logicals++;
        num_basic_variables++;
      }
    }
    // Double-check that we have the right number of basic variables
    nonbasicFlag_valid = num_basic_variables == simplex_lp.numRow_;
    assert(nonbasicFlag_valid);
    updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  }
  // Execute from here for all calls
  //
  // Note whether a HiGHS basis can be used to (try to) choose the better bound
  // for boxed variables
  bool have_highs_basis = basis.valid_;
  //
  // Possibly scale the LP to be solved
  //
  // If the LP to be solved isn't scaled then initialise unit scaling
  // factors, to simplify things if no scaling is performed. ToDo This
  // is inefficient if the LP isn't to be scales and is repeatedly
  // hot-started - but is this really going to happen?
  if (!simplex_lp_status.scaling_tried) scaleHighsModelInit(highs_model_object);
  //
  // Scale the LP to be used by the solver if scaling is to be used and the LP
  // is not already scaled
  bool scale_lp = options.simplex_scale_strategy != SimplexScaleStrategy::OFF &&
                  !simplex_lp_status.scaling_tried;
  if (scale_lp) {
    timer.start(simplex_info.clock_[ScaleClock]);    
    scaleSimplexLp(highs_model_object);
    timer.stop(simplex_info.clock_[ScaleClock]);    
#ifdef HiGHSDEV
    // Analyse the scaled LP
    if (simplex_info.analyseLp) {
      analyseLp(highs_model_object.lp_, "Unscaled");
      HighsScale& scale = highs_model_object.scale_;
      if (scale.is_scaled_) {
        analyseVectorValues("Column scaling factors", simplex_lp.numCol_,
                            scale.col_, false);
        analyseVectorValues("Row    scaling factors", simplex_lp.numRow_,
                            scale.row_, false);
        analyseLp(simplex_lp, "Scaled");
      }
    }
#endif
  }
  // Now there is a valid nonbasicFlag and basicIndex, possibly
  // reinvert to check for basis condition/singularity
  //
  // First setup the factor arrays if they don't exist
  if (!simplex_lp_status.has_factor_arrays) {
    factor.setup(simplex_lp.numCol_, simplex_lp.numRow_, &simplex_lp.Astart_[0],
                 &simplex_lp.Aindex_[0], &simplex_lp.Avalue_[0],
                 &simplex_basis.basicIndex_[0]);
    simplex_lp_status.has_factor_arrays = true;
  }
  // Reinvert if there isn't a fresh INVERT. ToDo Override this for MIP hot
  // start
  bool reinvert = !simplex_lp_status.has_fresh_invert;
  if (reinvert) {
    int rankDeficiency = compute_factor(highs_model_object);
    if (rankDeficiency) {
      // ToDo Handle rank deficiency by replacing singular columns with logicals
      throw runtime_error("Transition has singular basis matrix");
    }
    simplex_lp_status.has_fresh_invert = true;
  }
  // Possibly check for basis condition. ToDo Override this for MIP hot start
  bool basis_condition_ok = true;
  if (highs_model_object.options_.simplex_initial_condition_check) {
    timer.start(simplex_info.clock_[BasisConditionClock]);
    double basis_condition = computeBasisCondition(highs_model_object);
    timer.stop(simplex_info.clock_[BasisConditionClock]);
    double basis_condition_tolerance =
        highs_model_object.options_.simplex_initial_condition_tolerance;
    basis_condition_ok = basis_condition < basis_condition_tolerance;
    HighsMessageType message_type = HighsMessageType::INFO;
    std::string condition_comment;
    if (basis_condition_ok) {
      condition_comment = "is within";
    } else {
      message_type = HighsMessageType::WARNING;
      condition_comment = "exceeds";
    }
    HighsLogMessage(
        message_type,
        "Initial basis condition estimate of %g %s the tolerance of %g",
        basis_condition, condition_comment.c_str(), basis_condition_tolerance);
  }
  // ToDo Handle ill-conditioned basis with basis crash, in which case
  // ensure that HiGHS and simplex basis are invalidated and simplex
  // work and base arrays are re-populated
  //  assert(basis_condition_ok);
  if (!basis_condition_ok) {
    HCrash crash(highs_model_object);
    timer.start(simplex_info.clock_[CrashClock]);
    crash.crash(SimplexCrashStrategy::BASIC);
    timer.stop(simplex_info.clock_[CrashClock]);
    HighsLogMessage(HighsMessageType::INFO,
                    "Performed crash to prioritise previously basic variables "
                    "in well-conditioned basis");
    // Use nonbasicFlag to form basicIndex
    // Allocate memory for basicIndex
    simplex_basis.basicIndex_.resize(highs_model_object.lp_.numRow_);
    int num_basic_variables = 0;
    simplex_info.num_basic_logicals = 0;
    for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) {
        simplex_basis.basicIndex_[num_basic_variables] = iVar;
        if (iVar >= simplex_lp.numCol_) simplex_info.num_basic_logicals++;
        num_basic_variables++;
      }
    }
    // Double-check that we have the right number of basic variables
    assert(num_basic_variables == simplex_lp.numRow_);
    updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
    // Report on the outcome of crash
    int num_basic_structurals =
        simplex_lp.numRow_ - simplex_info.num_basic_logicals;
    HighsLogMessage(HighsMessageType::INFO,
                    "Crash has created a basis with %d/%d structurals",
                    num_basic_structurals, simplex_lp.numRow_);
    // Now reinvert
    int rankDeficiency = compute_factor(highs_model_object);
    if (rankDeficiency) {
      // ToDo Handle rank deficiency by replacing singular columns with logicals
      throw runtime_error("Transition has singular basis matrix");
    }
    simplex_lp_status.has_fresh_invert = true;

    // Check the condition after the basis crash
    timer.start(simplex_info.clock_[BasisConditionClock]);
    double basis_condition = computeBasisCondition(highs_model_object);
    timer.stop(simplex_info.clock_[BasisConditionClock]);
    double basis_condition_tolerance =
        highs_model_object.options_.simplex_initial_condition_tolerance;
    basis_condition_ok = basis_condition < basis_condition_tolerance;
    HighsMessageType message_type = HighsMessageType::INFO;
    std::string condition_comment;
    if (basis_condition_ok) {
      condition_comment = "is within";
    } else {
      message_type = HighsMessageType::WARNING;
      condition_comment = "exceeds";
    }
    HighsLogMessage(
        message_type,
        "Initial basis condition estimate of %11.4g %s the tolerance of %g",
        basis_condition, condition_comment.c_str(), basis_condition_tolerance);
  }

  // Now there are nonbasicFlag and basicIndex corresponding to a
  // basis with well-conditioned invertible representation
  //
  // Possibly set up the column-wise and row-wise copies of the matrix
  if (!simplex_lp_status.has_matrix_col_wise ||
      simplex_lp_status.has_matrix_row_wise) {
    matrix.setup(simplex_lp.numCol_, simplex_lp.numRow_, &simplex_lp.Astart_[0],
                 &simplex_lp.Aindex_[0], &simplex_lp.Avalue_[0],
                 &simplex_basis.nonbasicFlag_[0]);
    simplex_lp_status.has_matrix_col_wise = true;
    simplex_lp_status.has_matrix_row_wise = true;
  }
  // Possibly set up the simplex work and base arrays
  // ToDo Stop doing this always
  //  if (!simplex_lp_status.has_basis) {
  // Allocate memory for nonbasicMove
  simplex_basis.nonbasicMove_.resize(simplex_lp.numCol_ + simplex_lp.numRow_);
  allocate_work_and_base_arrays(highs_model_object);
  initialise_cost(highs_model_object);
  initialise_bound(highs_model_object);
  // Don't have a simplex basis since nonbasicMove is not set up.
  const int illegal_move_value = -99;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
      // Nonbasic variable
      double lower = simplex_info.workLower_[iVar];
      double upper = simplex_info.workUpper_[iVar];
      int move = illegal_move_value;
      double value;
      if (lower == upper) {
        // Fixed
        value = lower;
        move = NONBASIC_MOVE_ZE;
      } else if (!highs_isInfinity(-lower)) {
        // Finite lower bound so boxed or lower
        if (!highs_isInfinity(upper)) {
          // Finite upper bound so boxed
          // Determine the bound to set the value to according to, in order of
          // priority
          // 1. Any valid HiGHS basis status
          if (have_highs_basis) {
            if (iVar < simplex_lp.numCol_) {
              // Column
              if (basis.col_status[iVar] == HighsBasisStatus::LOWER) {
                // Set to lower bound
                move = NONBASIC_MOVE_UP;
                value = lower;
              } else if (basis.col_status[iVar] == HighsBasisStatus::UPPER) {
                // Set to upper bound
                move = NONBASIC_MOVE_DN;
                value = upper;
              }
            } else {
              // Row
              int iRow = iVar - simplex_lp.numCol_;
              if (basis.row_status[iRow] == HighsBasisStatus::LOWER) {
                // Set to upper bound
                move = NONBASIC_MOVE_DN;
                value = upper;
              } else if (basis.row_status[iRow] == HighsBasisStatus::UPPER) {
                // Set to lower bound
                move = NONBASIC_MOVE_UP;
                value = lower;
              }
            }
          }
          // 2. Any HiGHS solution value
          if (move == illegal_move_value && have_highs_solution) {
            double midpoint = 0.5 * (lower + upper);
            if (iVar < simplex_lp.numCol_) {
              // Column
              if (solution.col_value[iVar] < midpoint) {
                // Set to lower bound
                move = NONBASIC_MOVE_UP;
                value = lower;
              } else {
                // Set to upper bound
                move = NONBASIC_MOVE_DN;
                value = upper;
              }
            } else {
              // Row
              if (solution.col_value[iVar] < midpoint) {
                // Set to upper bound
                move = NONBASIC_MOVE_DN;
                value = upper;
              } else {
                // Set to lower bound
                move = NONBASIC_MOVE_UP;
                value = lower;
              }
            }
          }
          // 3. Lower bound for original LP
          if (move == illegal_move_value) {
            if (iVar < simplex_lp.numCol_) {
              // Set to lower bound
              move = NONBASIC_MOVE_UP;
              value = lower;
            } else {
              // Row
              // Set to upper bound
              move = NONBASIC_MOVE_DN;
              value = upper;
            }
          }
        } else {
          // Lower (since upper bound is infinite)
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
      simplex_info.workValue_[iVar] = value;
      simplex_basis.nonbasicMove_[iVar] = move;
    } else {
      // Basic variable
      simplex_basis.nonbasicMove_[iVar] = NONBASIC_MOVE_ZE;
    }
  }
  //  } else {}

  // Simplex basis is now valid
  simplex_lp_status.has_basis = true;

  // Possibly solve for the basic primal and nonbasic dual values to determine
  // which simplex solver to use, unless it's forced
  //  if (simplex_lp_status.has_basic_primal_values) {
  compute_primal(highs_model_object);
  simplex_lp_status.has_basic_primal_values = true;
  //}
  //  if (simplex_lp_status.has_basic_dual_values) {
  compute_dual(highs_model_object);
  simplex_lp_status.has_nonbasic_dual_values = true;
  //}
  computeDualObjectiveValue(highs_model_object);
  computePrimalObjectiveValue(highs_model_object);
  simplex_lp_status.valid = true;

  // Store, analyse and possibly report the number of primal and dual
  // infeasiblities and the simplex status
  computePrimalInfeasible(highs_model_object);
  computeDualInfeasible(highs_model_object);
  SimplexSolutionStatus solution_status;
  bool primal_feasible =
      simplex_info.num_primal_infeasibilities ==
      0;  // && max_primal_residual < primal_feasibility_tolerance;
  bool dual_feasible = simplex_info.num_dual_infeasibilities ==
                       0;  // && max_dual_residual < dual_feasibility_tolerance;
  if (primal_feasible) {
    if (dual_feasible) {
      solution_status = SimplexSolutionStatus::OPTIMAL;
    } else {
      solution_status = SimplexSolutionStatus::PRIMAL_FEASIBLE;
    }
  } else {
    if (dual_feasible) {
      solution_status = SimplexSolutionStatus::DUAL_FEASIBLE;
    } else {
      solution_status = SimplexSolutionStatus::UNSET;
    }
  }
  simplex_lp_status.solution_status = solution_status;
  //
#ifdef HiGHSDEV
  // If there is a HiGHS solution then determine the changes in basic
  // and nonbasic values and duals for columns and rows
  if (have_highs_solution) {
    // Go through the columns, finding the differences in nonbasic column values and duals
    int num_nonbasic_col_value_differences = 0;
    double sum_nonbasic_col_value_differences = 0;
    int num_nonbasic_col_dual_differences = 0;
    double sum_nonbasic_col_dual_differences = 0;
    HighsScale& scale = highs_model_object.scale_;
    for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
      int iVar = iCol;
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
	// Consider this nonbasic column
	double local_col_value = simplex_info.workValue_[iVar] * scale.col_[iCol];
	double local_col_dual = simplex_lp.sense_ * simplex_info.workDual_[iVar] / (scale.col_[iCol] / scale.cost_);
	double value_difference = fabs(local_col_value - solution.col_value[iCol]);
	double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
	if (value_difference > simplex_info.primal_feasibility_tolerance) num_nonbasic_col_value_differences++;
	sum_nonbasic_col_value_differences += value_difference;
	if (value_difference > simplex_info.dual_feasibility_tolerance) num_nonbasic_col_dual_differences++;
	sum_nonbasic_col_dual_differences += dual_difference; 
      }
    }
    // Go through the rows, finding the differences in nonbasic and
    // basic row values and duals, as well as differences in basic
    // column values and duals
    int num_nonbasic_row_value_differences = 0;
    double sum_nonbasic_row_value_differences = 0;
    int num_nonbasic_row_dual_differences = 0;
    double sum_nonbasic_row_dual_differences = 0;
    int num_basic_col_value_differences = 0;
    double sum_basic_col_value_differences = 0;
    int num_basic_col_dual_differences = 0;
    double sum_basic_col_dual_differences = 0;
    int num_basic_row_value_differences = 0;
    double sum_basic_row_value_differences = 0;
    int num_basic_row_dual_differences = 0;
    double sum_basic_row_dual_differences = 0;

    for (int ix = 0; ix < simplex_lp.numRow_; ix++) {
      int iRow = ix;
      int iVar = simplex_lp.numCol_ + iRow;
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
	// Consider this nonbasic row
	double local_row_value = - simplex_info.workValue_[iVar] / scale.row_[iRow];
	double local_row_dual = simplex_lp.sense_ * simplex_info.workDual_[iVar] * (scale.row_[iRow] * scale.cost_);
	double value_difference = fabs(local_row_value - solution.row_value[iRow]);
	double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
	if (value_difference > simplex_info.primal_feasibility_tolerance) num_nonbasic_row_value_differences++;
	sum_nonbasic_row_value_differences += value_difference;
	if (value_difference > simplex_info.dual_feasibility_tolerance) num_nonbasic_row_dual_differences++;
	sum_nonbasic_row_dual_differences += dual_difference; 
      }
      // Consider the basic variable associated with this row index
      iVar = simplex_basis.basicIndex_[ix];
      if (iVar < simplex_lp.numCol_) {
	// Consider this basic column
	int iCol = iVar;
	double local_col_value = simplex_info.baseValue_[ix] * scale.col_[iCol];
	double local_col_dual = 0;
	double value_difference = fabs(local_col_value - solution.col_value[iCol]);
	double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
	if (value_difference > simplex_info.primal_feasibility_tolerance) num_basic_col_value_differences++;
	sum_basic_col_value_differences += value_difference;
	if (value_difference > simplex_info.dual_feasibility_tolerance) num_basic_col_dual_differences++;
	sum_basic_col_dual_differences += dual_difference; 
      } else {
	// Consider this basic row
	iRow = iVar - simplex_lp.numCol_;
	double local_row_value = -simplex_info.baseValue_[ix] / scale.row_[iRow];
	double local_row_dual = 0;
	double value_difference = fabs(local_row_value - solution.row_value[iRow]);
	double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
	if (value_difference > simplex_info.primal_feasibility_tolerance) num_basic_row_value_differences++;
	sum_basic_row_value_differences += value_difference;
	if (value_difference > simplex_info.dual_feasibility_tolerance) num_basic_row_dual_differences++;
	sum_basic_row_dual_differences += dual_difference;
      }
    }	
    double acceptable_difference_sum = simplex_info.primal_feasibility_tolerance + simplex_info.dual_feasibility_tolerance;
    bool significant_nonbasic_value_differences = sum_nonbasic_col_value_differences + sum_nonbasic_row_value_differences > 0;
    bool significant_basic_value_differences = sum_basic_col_value_differences + sum_basic_row_value_differences > acceptable_difference_sum;      
    bool significant_nonbasic_col_dual_differences = sum_nonbasic_col_dual_differences > acceptable_difference_sum;
    bool significant_nonbasic_row_dual_differences = sum_nonbasic_row_dual_differences > acceptable_difference_sum;
    bool significant_basic_dual_differences = sum_basic_col_dual_differences + sum_basic_row_dual_differences > 0;
    if (significant_nonbasic_value_differences ||
	significant_basic_value_differences ||
	significant_nonbasic_col_dual_differences ||
	significant_nonbasic_row_dual_differences ||
	significant_basic_dual_differences) {
      printf("In transition(): There are significant value and dual differences\n");
    } else {
      printf("In transition(): There are no significant value and dual differences\n");
    }
    if (significant_nonbasic_value_differences) {
      printf("Nonbasic column value differences: %6d (%11.4g)\n", num_nonbasic_col_value_differences, sum_nonbasic_col_value_differences);
      printf("Nonbasic row    value differences: %6d (%11.4g)\n", num_nonbasic_row_value_differences, sum_nonbasic_row_value_differences);
    }
    if (significant_basic_value_differences) {
      printf("Basic    column value differences: %6d (%11.4g)\n", num_basic_col_value_differences, sum_basic_col_value_differences);
      printf("Basic    row    value differences: %6d (%11.4g)\n", num_basic_row_value_differences, sum_basic_row_value_differences);
    }
    if (significant_nonbasic_col_dual_differences) 
      printf("Nonbasic column  dual differences: %6d (%11.4g)\n", num_nonbasic_col_dual_differences, sum_nonbasic_col_dual_differences);
    if (significant_nonbasic_row_dual_differences)
      printf("Nonbasic row     dual differences: %6d (%11.4g)\n", num_nonbasic_row_dual_differences, sum_nonbasic_row_dual_differences);
    if (significant_basic_dual_differences) {
      printf("Basic    column  dual differences: %6d (%11.4g)\n", num_basic_col_dual_differences, sum_basic_col_dual_differences);
      printf("Basic    row     dual differences: %6d (%11.4g)\n", num_basic_row_dual_differences, sum_basic_row_dual_differences);
    }
    printf("grep_transition,%s,%.15g,%d,%g,%d,%g,%s,%d,%g,%d,%g,%d,%g,%d,%g,%d,%g,%d,%g,%d,%g,%d,%g\n",
	   simplex_lp.model_name_.c_str(),
	   simplex_info.primal_objective_value,
	   simplex_info.num_primal_infeasibilities,
	   simplex_info.sum_primal_infeasibilities,
	   simplex_info.num_dual_infeasibilities,
	   simplex_info.sum_dual_infeasibilities,
	   SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str(),
	   num_nonbasic_col_value_differences, sum_nonbasic_col_value_differences,
	   num_nonbasic_row_value_differences, sum_nonbasic_row_value_differences,
	   num_basic_col_value_differences, sum_basic_col_value_differences,
	   num_basic_row_value_differences, sum_basic_row_value_differences,
	   num_nonbasic_col_dual_differences, sum_nonbasic_col_dual_differences,
	   num_nonbasic_row_dual_differences, sum_nonbasic_row_dual_differences,
	   num_basic_col_dual_differences, sum_basic_col_dual_differences,
	   num_basic_row_dual_differences, sum_basic_row_dual_differences);
  }
#endif  
  HighsLogMessage(
      HighsMessageType::INFO,
      "Initial basic solution: Objective = %.15g; "
      "Infeasibilities Pr %d(%g); Du %d(%g); Status: %s",
      simplex_info.primal_objective_value,
      simplex_info.num_primal_infeasibilities,
      simplex_info.sum_primal_infeasibilities,
      simplex_info.num_dual_infeasibilities,
      simplex_info.sum_dual_infeasibilities,
      SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str());
	 
  return solution_status;
}

bool dual_infeasible(const double value, const double lower, const double upper,
                     const double dual, const double value_tolerance,
                     const double dual_tolerance) {
  double midpoint = (lower + upper) * 0.5;
  double residual = max(lower - value, value - upper);
  bool infeasible = false;
  if (highs_isInfinity(-lower)) {
    // Infinite lower bound
    if (highs_isInfinity(upper)) {
      // Infinite upper bound
      // Free
      infeasible = fabs(dual) >= dual_tolerance;
    } else {
      // Finite upper bound
      // Upper bounded - and assumed to be nonbasic at that bound
      if (fabs(residual) >= value_tolerance) {
        printf("dual_infeasible: %12g %12g %12g %12g %12g\n", value, lower,
               upper, residual, value_tolerance);
      }
      assert(fabs(residual) < value_tolerance);
      infeasible = dual >= dual_tolerance;
    }
  } else {
    // Finite lower bound
    if (highs_isInfinity(upper)) {
      // Infinite upper bound
      // Lower bounded - and assumed to be nonbasic at that bound
      assert(fabs(residual) < value_tolerance);
      infeasible = dual <= -dual_tolerance;
    } else {
      // Finite upper bound
      // Assumed to be nonbasic at that bound
      assert(fabs(residual) < value_tolerance);
      if (lower < upper) {
        // Boxed
        if (value < midpoint) {
          // At lower bound
          infeasible = dual <= -dual_tolerance;
        } else {
          // At upper bound
          infeasible = dual >= dual_tolerance;
        }
      } else {
        // Fixed
        infeasible = false;
      }
    }
  }
  return infeasible;
}

void append_nonbasic_cols_to_basis(HighsLp& lp, HighsBasis& basis,
                                   int XnumNewCol) {
#ifdef HiGHSDEV
  printf("!! Don't do this if basis is invalid! !!\n");
#endif
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  basis.col_status.resize(newNumCol);
  // Make any new columns nonbasic
  for (int col = lp.numCol_; col < newNumCol; col++) {
    if (!highs_isInfinity(-lp.colLower_[col])) {
      // Has finite lower bound so set it there
      basis.col_status[col] = HighsBasisStatus::LOWER;
    } else if (!highs_isInfinity(lp.colUpper_[col])) {
      // Has finite upper bound so set it there
      basis.col_status[col] = HighsBasisStatus::UPPER;
    } else {
      // Free variable so set to zero
      basis.col_status[col] = HighsBasisStatus::ZERO;
    }
  }
}

void append_nonbasic_cols_to_basis(HighsLp& lp, SimplexBasis& simplex_basis,
                                   int XnumNewCol) {
#ifdef HiGHSDEV
  printf("!! Don't do this if basis is invalid! !!\n");
#endif
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  int newNumTot = newNumCol + lp.numRow_;
  simplex_basis.nonbasicFlag_.resize(newNumTot);
  // Shift the row data in basicIndex and nonbasicFlag if necessary
  for (int row = lp.numRow_ - 1; row >= 0; row--) {
    simplex_basis.basicIndex_[row] += XnumNewCol;
    simplex_basis.nonbasicFlag_[newNumCol + row] =
        simplex_basis.nonbasicFlag_[lp.numCol_ + row];
  }
  // Make any new columns nonbasic
  for (int col = lp.numCol_; col < newNumCol; col++) {
    simplex_basis.nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
  }
}

void append_basic_rows_to_basis(HighsLp& lp, HighsBasis& basis,
                                int XnumNewRow) {
#ifdef HiGHSDEV
  printf("!! Don't do this if basis is invalid! !!\n");
#endif
  // Add basic logicals
  if (XnumNewRow == 0) return;
  int newNumRow = lp.numRow_ + XnumNewRow;
  basis.row_status.resize(newNumRow);
  // Make any new rows basic
  for (int row = lp.numRow_; row < newNumRow; row++) {
    basis.row_status[row] = HighsBasisStatus::BASIC;
  }
}

bool highs_basis_ok(
		    // HighsLp& lp, HighsBasis& basis
		    ) {
#ifdef HiGHSDEV
  printf("!! Don't check if basis is invalid! !!\n");
  printf("!! WRITE highs_basis_ok for HighsBasis !!\n");
#endif
  return false;
}

bool nonbasic_flag_basic_index_ok(HighsLp& lp, SimplexBasis& simplex_basis) {
#ifdef HiGHSDEV
  printf("!! Don't check if basis is invalid! !!\n");
#endif
  int numTot = lp.numCol_ + lp.numRow_;
  int num_basic_variables = 0;
  for (int var = 0; var < numTot; var++)
    if (!simplex_basis.nonbasicFlag_[var]) num_basic_variables++;
  assert(num_basic_variables == lp.numRow_);
  if (num_basic_variables != lp.numRow_) return false;
  for (int row = 0; row < lp.numRow_; row++) {
    int flag = simplex_basis.nonbasicFlag_[simplex_basis.basicIndex_[row]];
    assert(!flag);
    if (flag) return false;
  }
  return true;
}

#ifdef HiGHSDEV
void report_basis(HighsLp& lp, HighsBasis& basis) {
#ifdef HiGHSDEV
  printf("!! WRITE report_basis for HighsBasis !!\n");
#endif
  if (lp.numCol_ > 0) printf("   Col          Flag   Move\n");
  for (int col = 0; col < lp.numCol_; col++) {
    printf("%6d         %6d\n", col, (int)basis.col_status[col]);
  }
  if (lp.numRow_ > 0) printf("   Row  Basic   Flag   Move\n");
  for (int row = 0; row < lp.numRow_; row++) {
    printf("%6d         %6d\n", row, (int)basis.row_status[row]);
  }
}

void report_basis(HighsLp& lp, SimplexBasis& simplex_basis) {
  if (lp.numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < lp.numCol_; col++) {
    int var = col;
    if (simplex_basis.nonbasicFlag_[var])
      printf("%6d %6d        %6d\n", var, col,
             simplex_basis.nonbasicFlag_[var]);
    // simplex_basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d\n", var, col, simplex_basis.nonbasicFlag_[var]);
  }
  if (lp.numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    if (simplex_basis.nonbasicFlag_[var])
      printf("%6d %6d %6d %6d\n", var, row, simplex_basis.basicIndex_[row],
             simplex_basis.nonbasicFlag_[var]);
    // simplex_basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, simplex_basis.basicIndex_[row],
             simplex_basis.nonbasicFlag_[var]);
  }
}
#endif

/**
 * @brief Simplex utilities
 */

/*
// Increment iteration count (here!) and (possibly) store the pivots for
// debugging NLA
void record_pivots(int columnIn, int columnOut, double alpha) {
  // NB This is where the iteration count is updated!
  if (columnIn >= 0) simplex_info.iteration_count++;
#ifdef HiGHSDEV
  historyColumnIn.push_back(columnIn);
  historyColumnOut.push_back(columnOut);
  historyAlpha.push_back(alpha);
#endif
}
#ifdef HiGHSDEV
// Store and write out the pivots for debugging NLA
void writePivots(const char* suffix) {
  string filename = "z-" + simplex_lp_->model_name_ + "-" + suffix;
  ofstream output(filename.c_str());
  int count = historyColumnIn.size();
  double current_run_highs_time = timer_->readRunHighsClock();
  output << simplex_lp_->model_name_ << " " << count << "\t" <<
current_run_highs_time << endl; output << setprecision(12); for (int i = 0; i <
count; i++) { output << historyColumnIn[i] << "\t"; output <<
historyColumnOut[i] << "\t"; output << historyAlpha[i] << endl;
  }
  output.close();
}
#endif
*/

void computeDualObjectiveValue(HighsModelObject& highs_model_object,
                               int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;

  simplex_info.dual_objective_value = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (highs_model_object.simplex_basis_.nonbasicFlag_[i]) {
      simplex_info.dual_objective_value +=
          simplex_info.workValue_[i] * simplex_info.workDual_[i];
    }
  }
  if (phase != 1) {
    simplex_info.dual_objective_value *= highs_model_object.scale_.cost_;
    simplex_info.dual_objective_value -= simplex_lp.offset_;
  }
  // Now have dual objective value
  simplex_lp_status.has_dual_objective_value = true;
}

void computePrimalObjectiveValue(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  simplex_info.primal_objective_value = 0;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_basis.basicIndex_[row];
    if (var < simplex_lp.numCol_)
      simplex_info.primal_objective_value +=
          simplex_info.baseValue_[row] * simplex_lp.colCost_[var];
  }
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    if (simplex_basis.nonbasicFlag_[col])
      simplex_info.primal_objective_value +=
          simplex_info.workValue_[col] * simplex_lp.colCost_[col];
  }
  simplex_info.primal_objective_value *= highs_model_object.scale_.cost_;
  simplex_info.primal_objective_value -= simplex_lp.offset_;
  // Now have primal objective value
  simplex_lp_status.has_primal_objective_value = true;
}

void initialiseSimplexLpRandomVectors(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const int numCol = highs_model_object.simplex_lp_.numCol_;
  const int numTot = highs_model_object.simplex_lp_.numCol_ +
                     highs_model_object.simplex_lp_.numRow_;
  // Instantiate and (re-)initialise the random number generator
  HighsRandom& random = highs_model_object.random_;
  random.initialise();
  //
  // Generate a random permutation of the column indices
  simplex_info.numColPermutation_.resize(numCol);
  int* numColPermutation = &simplex_info.numColPermutation_[0];
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
  simplex_info.numTotPermutation_.resize(numTot);
  int* numTotPermutation = &simplex_info.numTotPermutation_[0];
  for (int i = 0; i < numTot; i++) numTotPermutation[i] = i;
  for (int i = numTot - 1; i >= 1; i--) {
    int j = random.integer() % (i + 1);
    std::swap(numTotPermutation[i], numTotPermutation[j]);
  }

  // Generate a vector of random reals
  simplex_info.numTotRandomValue_.resize(numTot);
  double* numTotRandomValue = &simplex_info.numTotRandomValue_[0];
  for (int i = 0; i < numTot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

// SCALING:
#ifdef HiGHSDEV
// Information on large costs
const double tlLargeCo = 1e5;
int numLargeCo;
vector<int> largeCostFlag;
double largeCostScale;
#endif

void scaleHighsModelInit(HighsModelObject& highs_model_object) {
  HighsScale& scale = highs_model_object.scale_;
  scale.is_scaled_ = false;
  scale.col_.assign(highs_model_object.simplex_lp_.numCol_, 1);
  scale.row_.assign(highs_model_object.simplex_lp_.numRow_, 1);
  scale.cost_ = 1;
  scale.extreme_equilibration_improvement_ = 1;
  scale.mean_equilibration_improvement_ = 1;
#ifdef HiGHSDEV
  //  largeCostScale = 1;
#endif
}

void scaleCosts(HighsModelObject& highs_model_object) {
  // Scale the costs by no less than minAlwCostScale
  double max_allowed_cost_scale =
      pow(2.0, highs_model_object.options_.allowed_simplex_scale_factor);
  double cost_scale;
  double max_nonzero_cost = 0;
  for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    if (highs_model_object.simplex_lp_.colCost_[iCol]) {
      max_nonzero_cost =
          max(fabs(highs_model_object.simplex_lp_.colCost_[iCol]),
              max_nonzero_cost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  cost_scale = 1;
  const double ln2 = log(2.0);
  // Scale if the max cost is positive and outside the range [1/16, 16]
  if ((max_nonzero_cost > 0) &&
      ((max_nonzero_cost < (1.0 / 16)) || (max_nonzero_cost > 16))) {
    cost_scale = max_nonzero_cost;
    cost_scale = pow(2.0, floor(log(cost_scale) / ln2 + 0.5));
    cost_scale = min(cost_scale, max_allowed_cost_scale);
  }
  highs_model_object.scale_.cost_ = cost_scale;
  if (cost_scale == 1) return;
  // Scale the costs (and record of max_nonzero_cost) by cost_scale, being at
  // most max_allowed_cost_scale
  for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    highs_model_object.simplex_lp_.colCost_[iCol] /= cost_scale;
  }
  max_nonzero_cost /= cost_scale;

#ifdef HiGHSDEV
  /*
  bool alwLargeCostScaling = false;
    if (alwLargeCostScaling && (numLargeCo > 0)) {
    // Scale any large costs by largeCostScale, being at most (a further)
    // max_allowed_cost_scale
    largeCostScale = max_nonzero_cost;
    largeCostScale = pow(2.0, floor(log(largeCostScale) / ln2 + 0.5));
    largeCostScale = min(largeCostScale, max_allowed_cost_scale);
    printf(
    "   Scaling all |cost| > %11.4g by %11.4g\ngrep_LargeCostScale,%g,%g\n",
    tlLargeCo, largeCostScale, tlLargeCo, largeCostScale);
    for (int iCol = 0; iCol < highs_model_object.simplex_lp_.numCol_; iCol++) {
    if (largeCostFlag[iCol]) {
    highs_model_object.simplex_lp_.colCost_[iCol] /= largeCostScale;
    }
    }
    }
  */
  //  utils.analyseVectorValues("Column costs",
  //  highs_model_object.simplex_lp_.numCol_,
  //  highs_model_object.simplex_lp_.colCost_, false);
#endif
}

void scaleFactorRanges(HighsModelObject& highs_model_object,
                       double& min_col_scale, double& max_col_scale,
                       double& min_row_scale, double& max_row_scale) {
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  double* colScale = &highs_model_object.scale_.col_[0];
  double* rowScale = &highs_model_object.scale_.row_[0];
  // Determine the max and min row and column scaling factors
  min_col_scale = HIGHS_CONST_INF;
  max_col_scale = 1 / HIGHS_CONST_INF;
  min_row_scale = HIGHS_CONST_INF;
  max_row_scale = 1 / HIGHS_CONST_INF;
  for (int iCol = 0; iCol < numCol; iCol++) {
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    min_row_scale = min(rowScale[iRow], min_row_scale);
    max_row_scale = max(rowScale[iRow], max_row_scale);
  }
}

void scaleSimplexLp(HighsModelObject& highs_model_object) {
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.scaling_tried) return;
  // Scale the LP highs_model_object.simplex_lp_, assuming all data are in place
  // Reset all scaling to 1
  HighsScale& scale = highs_model_object.scale_;
  scaleHighsModelInit(highs_model_object);
  int numCol = highs_model_object.simplex_lp_.numCol_;
  int numRow = highs_model_object.simplex_lp_.numRow_;
  double* colScale = &highs_model_object.scale_.col_[0];
  double* rowScale = &highs_model_object.scale_.row_[0];
  int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
  int* Aindex = &highs_model_object.simplex_lp_.Aindex_[0];
  double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
  double* colCost = &highs_model_object.simplex_lp_.colCost_[0];
  double* colLower = &highs_model_object.simplex_lp_.colLower_[0];
  double* colUpper = &highs_model_object.simplex_lp_.colUpper_[0];
  double* rowLower = &highs_model_object.simplex_lp_.rowLower_[0];
  double* rowUpper = &highs_model_object.simplex_lp_.rowUpper_[0];

  // Allow a switch to/from the original scaling rules
  bool original_scaling = highs_model_object.options_.simplex_scale_strategy ==
                          SimplexScaleStrategy::HSOL;
  bool allow_cost_scaling = false;
  if (original_scaling) allow_cost_scaling = false;
  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  double min_matrix_value = HIGHS_CONST_INF, max_matrix_value = 0;
  for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
    double value = fabs(Avalue[k]);
    min_matrix_value = min(min_matrix_value, value);
    max_matrix_value = max(max_matrix_value, value);
  }
  bool no_scaling = min_matrix_value >= 0.2 && max_matrix_value <= 5;
  //   no_scaling = false; printf("!!!! FORCE SCALING !!!!\n");
  if (no_scaling) {
    // No matrix scaling, but possible cost scaling
    HighsLogMessage(HighsMessageType::INFO, "Scaling: Matrix has min(max) values of %g(%g) so none performed", min_matrix_value, max_matrix_value);
    // Possibly scale the costs
    if (allow_cost_scaling) {
      scaleCosts(highs_model_object);
      // Simplex LP is now only scaled if there is a cost scaling factor
      scale.is_scaled_ = scale.cost_ != 1;
    }
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::SCALE);
    return;
  }
  // Include cost in scaling if minimum nonzero cost is less than 0.1
  double min_nonzero_cost = HIGHS_CONST_INF;
  for (int i = 0; i < numCol; i++) {
    if (colCost[i]) min_nonzero_cost = min(fabs(colCost[i]), min_nonzero_cost);
  }
  bool include_cost_in_scaling = false;
  //  if (original_scaling)
  include_cost_in_scaling = min_nonzero_cost < 0.1;

  // Limits on scaling factors
  double max_allow_scale;
  double min_allow_scale;
  if (original_scaling) {
    max_allow_scale = HIGHS_CONST_INF;
  } else {
    max_allow_scale =
        pow(2.0, highs_model_object.options_.allowed_simplex_scale_factor);
  }
  min_allow_scale = 1 / max_allow_scale;

  double min_allow_col_scale = min_allow_scale;
  double max_allow_col_scale = max_allow_scale;
  double min_allow_row_scale = min_allow_scale;
  double max_allow_row_scale = max_allow_scale;

  // Search up to 6 times
  vector<double> row_min_value(numRow, HIGHS_CONST_INF);
  vector<double> row_max_value(numRow, 1 / HIGHS_CONST_INF);
  for (int search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < numCol; iCol++) {
      // For column scale (find)
      double col_min_value = HIGHS_CONST_INF;
      double col_max_value = 1 / HIGHS_CONST_INF;
      double abs_col_cost = fabs(colCost[iCol]);
      if (include_cost_in_scaling && abs_col_cost != 0) {
        col_min_value = min(col_min_value, abs_col_cost);
        col_max_value = max(col_max_value, abs_col_cost);
      }
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
        col_min_value = min(col_min_value, value);
        col_max_value = max(col_max_value, value);
      }
      double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
      // Ensure that column scale factor is not excessively large or small
      colScale[iCol] =
          min(max(min_allow_col_scale, col_equilibration), max_allow_col_scale);
      // For row scale (only collect)
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        int iRow = Aindex[k];
        double value = fabs(Avalue[k]) * colScale[iCol];
        row_min_value[iRow] = min(row_min_value[iRow], value);
        row_max_value[iRow] = max(row_max_value[iRow], value);
      }
    }
    // For row scale (find)
    for (int iRow = 0; iRow < numRow; iRow++) {
      double row_equilibration = 1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
      // Ensure that row scale factor is not excessively large or small
      rowScale[iRow] =
          min(max(min_allow_row_scale, row_equilibration), max_allow_row_scale);
    }
    row_min_value.assign(numRow, HIGHS_CONST_INF);
    row_max_value.assign(numRow, 1 / HIGHS_CONST_INF);
  }
  // Make it numerical better
  // Also determine the max and min row and column scaling factors
  double min_col_scale = HIGHS_CONST_INF;
  double max_col_scale = 1 / HIGHS_CONST_INF;
  double min_row_scale = HIGHS_CONST_INF;
  double max_row_scale = 1 / HIGHS_CONST_INF;
  const double log2 = log(2.0);
  for (int iCol = 0; iCol < numCol; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / log2 + 0.5));
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / log2 + 0.5));
    min_row_scale = min(rowScale[iRow], min_row_scale);
    max_row_scale = max(rowScale[iRow], max_row_scale);
  }
  // Apply scaling to matrix and bounds
  double min_original_col_equilibration = HIGHS_CONST_INF;
  double sum_original_log_col_equilibration = 0;
  double max_original_col_equilibration = 0;
  double min_original_row_equilibration = HIGHS_CONST_INF;
  double sum_original_log_row_equilibration = 0;
  double max_original_row_equilibration = 0;
  double min_col_equilibration = HIGHS_CONST_INF;
  double sum_log_col_equilibration = 0;
  double max_col_equilibration = 0;
  double min_row_equilibration = HIGHS_CONST_INF;
  double sum_log_row_equilibration = 0;
  double max_row_equilibration = 0;
  vector<double> original_row_min_value(numRow, HIGHS_CONST_INF);
  vector<double> original_row_max_value(numRow, 1 / HIGHS_CONST_INF);
  row_min_value.assign(numRow, HIGHS_CONST_INF);
  row_max_value.assign(numRow, 1 / HIGHS_CONST_INF);
  for (int iCol = 0; iCol < numCol; iCol++) {
    double original_col_min_value = HIGHS_CONST_INF;
    double original_col_max_value = 1 / HIGHS_CONST_INF;
    double col_min_value = HIGHS_CONST_INF;
    double col_max_value = 1 / HIGHS_CONST_INF;
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      int iRow = Aindex[k];
      double original_value = fabs(Avalue[k]);
      original_col_min_value = min(original_value, original_col_min_value);
      original_col_max_value = max(original_value, original_col_max_value);
      original_row_min_value[iRow] = min(original_row_min_value[iRow], original_value);
      original_row_max_value[iRow] = max(original_row_max_value[iRow], original_value);
      Avalue[k] *= (colScale[iCol] * rowScale[iRow]);
      double value = fabs(Avalue[k]);
      col_min_value = min(value, col_min_value);
      col_max_value = max(value, col_max_value);
      row_min_value[iRow] = min(row_min_value[iRow], value);
      row_max_value[iRow] = max(row_max_value[iRow], value);
    }
    double original_col_equilibration = 1 / sqrt(original_col_min_value * original_col_max_value);
    min_original_col_equilibration = min(original_col_equilibration, min_original_col_equilibration);
    sum_original_log_col_equilibration += log(original_col_equilibration);
    max_original_col_equilibration = max(original_col_equilibration, max_original_col_equilibration);
    double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
    min_col_equilibration = min(col_equilibration, min_col_equilibration);
    sum_log_col_equilibration += log(col_equilibration);
    max_col_equilibration = max(col_equilibration, max_col_equilibration);
  }

  for (int iRow = 0; iRow < numRow; iRow++) {
    double original_row_equilibration = 1 / sqrt(original_row_min_value[iRow] * original_row_max_value[iRow]);
    min_original_row_equilibration = min(original_row_equilibration, min_original_row_equilibration);
    sum_original_log_row_equilibration += log(original_row_equilibration);
    max_original_row_equilibration = max(original_row_equilibration, max_original_row_equilibration);
    double row_equilibration = 1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
    min_row_equilibration = min(row_equilibration, min_row_equilibration);
    sum_log_row_equilibration += log(row_equilibration);
    max_row_equilibration = max(row_equilibration, max_row_equilibration);
  }
  double geomean_original_col_equilibration = exp(sum_original_log_col_equilibration/numCol);
  double geomean_original_row_equilibration = exp(sum_original_log_row_equilibration/numRow);
#ifdef HiGHSDEV
  HighsLogMessage(HighsMessageType::INFO, "Scaling: Original equilibration: min/mean/max %11.4f/%11.4f/%11.4f (cols); min/mean/max %11.4f/%11.4f/%11.4f (rows)",
	 min_original_col_equilibration,
	 geomean_original_col_equilibration,
	 max_original_col_equilibration,
	 min_original_row_equilibration,
	 geomean_original_row_equilibration,
	 max_original_row_equilibration);
#endif
  double geomean_col_equilibration = exp(sum_log_col_equilibration/numCol);
  double geomean_row_equilibration = exp(sum_log_row_equilibration/numRow);
#ifdef HiGHSDEV
  HighsLogMessage(HighsMessageType::INFO, "Scaling: Final    equilibration: min/mean/max %11.4f/%11.4f/%11.4f (cols); min/mean/max %11.4f/%11.4f/%11.4f (rows)",
	 min_col_equilibration,
	 geomean_col_equilibration,
	 max_col_equilibration,
	 min_row_equilibration,
	 geomean_row_equilibration,
	 max_row_equilibration);
#endif
  scale.extreme_equilibration_improvement_ =
    (max_original_col_equilibration/min_original_col_equilibration +
     max_original_row_equilibration/min_original_row_equilibration)/
    (max_col_equilibration/min_col_equilibration +
     max_row_equilibration/min_row_equilibration);
  scale.mean_equilibration_improvement_ =
    (max(geomean_original_col_equilibration, 1/geomean_original_col_equilibration)*
     max(geomean_original_row_equilibration, 1/geomean_original_row_equilibration))/
    (max(geomean_col_equilibration, 1/geomean_col_equilibration)*
     max(geomean_row_equilibration, 1/geomean_row_equilibration));
  if (!original_scaling) {
    // Abandon scaling if it's not improved equlibration significantly
    // Unscale the matrix
    if (scale.extreme_equilibration_improvement_ < 10 &&
	scale.mean_equilibration_improvement_ < 1.1) {
      for (int iCol = 0; iCol < numCol; iCol++) {
	for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	  int iRow = Aindex[k];
	  Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
	}
      }
      scaleHighsModelInit(highs_model_object);
      HighsLogMessage(HighsMessageType::INFO, "Scaling: Extreme equilibration improved by a factor of only %g and mean equilibration by factor of only %g so no scaling applied",
		      scale.extreme_equilibration_improvement_, scale.mean_equilibration_improvement_);
      // Possibly scale the costs
      if (allow_cost_scaling) {
	scaleCosts(highs_model_object);
	// Simplex LP is now only scaled if there is a cost scaling factor
	scale.is_scaled_ = scale.cost_ != 1;
	  }
      updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
			    LpAction::SCALE);
      return;
    }
  }
  scale.is_scaled_ = true;
  HighsLogMessage(HighsMessageType::INFO, "Scaling: Improved extreme equilibration by factor %g and mean equilibration by factor %g",
		  scale.extreme_equilibration_improvement_, scale.mean_equilibration_improvement_);

  for (int iCol = 0; iCol < numCol; iCol++) {
    colLower[iCol] /= colLower[iCol] == -HIGHS_CONST_INF ? 1 : colScale[iCol];
    colUpper[iCol] /= colUpper[iCol] == +HIGHS_CONST_INF ? 1 : colScale[iCol];
    colCost[iCol] *= colScale[iCol];
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    rowLower[iRow] *= rowLower[iRow] == -HIGHS_CONST_INF ? 1 : rowScale[iRow];
    rowUpper[iRow] *= rowUpper[iRow] == +HIGHS_CONST_INF ? 1 : rowScale[iRow];
  }
  // Deduce the consequences of scaling the LP
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_, LpAction::SCALE);
  // Possibly scale the costs
  if (allow_cost_scaling) scaleCosts(highs_model_object);
}

// PERMUTE:

void permuteSimplexLp(HighsModelObject& highs_model_object) {
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
#ifdef HiGHSDEV
  printf("Called permuteSimplexLp: simplex_lp_status.is_permuted = %d\n",
         simplex_lp_status.is_permuted);
#endif
  if (simplex_lp_status.is_permuted) return;

  int numCol = highs_model_object.simplex_lp_.numCol_;
  vector<int>& numColPermutation =
      highs_model_object.simplex_info_.numColPermutation_;
  vector<int>& Astart = highs_model_object.simplex_lp_.Astart_;
  vector<int>& Aindex = highs_model_object.simplex_lp_.Aindex_;
  vector<double>& Avalue = highs_model_object.simplex_lp_.Avalue_;
  vector<double>& colCost = highs_model_object.simplex_lp_.colCost_;
  vector<double>& colLower = highs_model_object.simplex_lp_.colLower_;
  vector<double>& colUpper = highs_model_object.simplex_lp_.colUpper_;
  vector<double>& colScale = highs_model_object.scale_.col_;

  // 2. Duplicate the original data to copy from
  vector<int> saveAstart = highs_model_object.simplex_lp_.Astart_;
  vector<int> saveAindex = highs_model_object.simplex_lp_.Aindex_;
  vector<double> saveAvalue = highs_model_object.simplex_lp_.Avalue_;
  vector<double> saveColCost = highs_model_object.simplex_lp_.colCost_;
  vector<double> saveColLower = highs_model_object.simplex_lp_.colLower_;
  vector<double> saveColUpper = highs_model_object.simplex_lp_.colUpper_;
  vector<double> saveColScale = highs_model_object.scale_.col_;

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
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::PERMUTE);
}

void initialise_basic_index(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int num_basic_variables = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; var++) {
    if (!simplex_basis.nonbasicFlag_[var]) {
      assert(num_basic_variables < simplex_lp.numRow_);
      simplex_basis.basicIndex_[num_basic_variables] = var;
      num_basic_variables++;
    }
  }
  /*
  if (num_basic_variables != simplex_lp.numRow_) {
    printf("STRANGE: %d = num_basic_variables != simplex_lp.numRow_ = %d\n",
  num_basic_variables, simplex_lp.numRow_); fflush(stdout);
  }
  */
  assert(num_basic_variables == simplex_lp.numRow_);
}

void allocate_work_and_base_arrays(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Allocate bounds and solution spaces
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  simplex_info.workCost_.resize(numTot);
  simplex_info.workDual_.resize(numTot);
  simplex_info.workShift_.resize(numTot);

  simplex_info.workLower_.resize(numTot);
  simplex_info.workUpper_.resize(numTot);
  simplex_info.workRange_.resize(numTot);
  simplex_info.workValue_.resize(numTot);

  simplex_info.baseLower_.resize(simplex_lp.numRow_);
  simplex_info.baseUpper_.resize(simplex_lp.numRow_);
  simplex_info.baseValue_.resize(simplex_lp.numRow_);
}

void initialise_from_nonbasic(HighsModelObject& highs_model_object) {
  // Initialise basicIndex from nonbasic* then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  initialise_basic_index(highs_model_object);
  allocate_work_and_base_arrays(highs_model_object);
  populate_work_arrays(highs_model_object);

  // Deduce the consequences of a new basis
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_BASIS);
}

void replace_from_nonbasic(HighsModelObject& highs_model_object) {
  // Initialise basicIndex using nonbasic* then populate (where possible)
  // work* arrays
  initialise_basic_index(highs_model_object);
  populate_work_arrays(highs_model_object);

  // Deduce the consequences of a new basis
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_BASIS);
}

void initialise_with_logical_basis(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays

  for (int row = 0; row < simplex_lp.numRow_; row++)
    simplex_basis.basicIndex_[row] = simplex_lp.numCol_ + row;
  for (int col = 0; col < simplex_lp.numCol_; col++)
    simplex_basis.nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
  simplex_lp_status.has_basis = true;
  simplex_info.num_basic_logicals = simplex_lp.numRow_;

  allocate_work_and_base_arrays(highs_model_object);
  populate_work_arrays(highs_model_object);

  // Deduce the consequences of a new basis
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_BASIS);
}

void initialise_value_from_nonbasic(HighsModelObject& highs_model_object,
                                    int firstvar, int lastvar) {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  assert(firstvar >= 0);
  assert(lastvar < highs_model_object.simplex_lp_.numCol_ + highs_model_object.simplex_lp_.numRow_);
  // double dl_pr_act, norm_dl_pr_act;
  // norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++) {
    if (simplex_basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      // double prev_pr_act = simplex_info.workValue_[var];
      if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
        // Fixed
        simplex_info.workValue_[var] = simplex_info.workLower_[var];
        simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      } else if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        // Finite lower bound so boxed or lower
        if (!highs_isInfinity(simplex_info.workUpper_[var])) {
          // Finite upper bound so boxed
          if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
            // Set at lower
            simplex_info.workValue_[var] = simplex_info.workLower_[var];
          } else if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
            // Set at upper
            simplex_info.workValue_[var] = simplex_info.workUpper_[var];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
            simplex_info.workValue_[var] = simplex_info.workLower_[var];
          }
        } else {
          // Lower
          simplex_info.workValue_[var] = simplex_info.workLower_[var];
          simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
        }
      } else if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        // Upper
        simplex_info.workValue_[var] = simplex_info.workUpper_[var];
        simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        simplex_info.workValue_[var] = 0;
        simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      }
      // dl_pr_act = simplex_info.workValue_[var] - prev_pr_act;
      // norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (fabs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
      //      %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var,
      // simplex_info.workLower_[var],
      // simplex_info.workValue_[var], simplex_info.workUpper_[var],
      // simplex_info.workDual_[var], dl_pr_act);
    } else {
      // Basic variable
      simplex_basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
    }
  }
  //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
  //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is
  //  %g\n", norm_dl_pr_act);
}

void initialise_value(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  initialise_value_from_nonbasic(highs_model_object, 0, numTot - 1);
}

void initialise_phase2_col_bound(HighsModelObject& highs_model_object,
                                 int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  assert(firstcol >= 0);
  assert(lastcol < simplex_lp.numCol_);
  for (int col = firstcol; col <= lastcol; col++) {
    simplex_info.workLower_[col] = simplex_lp.colLower_[col];
    simplex_info.workUpper_[col] = simplex_lp.colUpper_[col];
    simplex_info.workRange_[col] =
        simplex_info.workUpper_[col] - simplex_info.workLower_[col];
  }
}

void initialise_phase2_row_bound(HighsModelObject& highs_model_object,
                                 int firstrow, int lastrow) {
  // Copy bounds and compute ranges
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  assert(firstrow >= 0);
  assert(lastrow < simplex_lp.numRow_);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = simplex_lp.numCol_ + row;
    simplex_info.workLower_[var] = -simplex_lp.rowUpper_[row];
    simplex_info.workUpper_[var] = -simplex_lp.rowLower_[row];
    simplex_info.workRange_[var] =
        simplex_info.workUpper_[var] - simplex_info.workLower_[var];
  }
}

void initialise_bound(HighsModelObject& highs_model_object, int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initialise_phase2_col_bound(highs_model_object, 0, simplex_lp.numCol_ - 1);
  initialise_phase2_row_bound(highs_model_object, 0, simplex_lp.numRow_ - 1);
  if (phase == 2) return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HIGHS_CONST_INF;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_info.workLower_[i] == -inf &&
        simplex_info.workUpper_[i] == inf) {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= simplex_lp.numCol_) continue;
      simplex_info.workLower_[i] = -1000,
      simplex_info.workUpper_[i] = 1000;  // FREE
    } else if (simplex_info.workLower_[i] == -inf) {
      simplex_info.workLower_[i] = -1, simplex_info.workUpper_[i] = 0;  // UPPER
    } else if (simplex_info.workUpper_[i] == inf) {
      simplex_info.workLower_[i] = 0, simplex_info.workUpper_[i] = 1;  // LOWER
    } else {
      simplex_info.workLower_[i] = 0,
      simplex_info.workUpper_[i] = 0;  // BOXED or FIXED
    }
    simplex_info.workRange_[i] =
        simplex_info.workUpper_[i] - simplex_info.workLower_[i];
  }
}

void initialise_phase2_col_cost(HighsModelObject& highs_model_object,
                                int firstcol, int lastcol) {
  // Copy the Phase 2 cost and zero the shift
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    simplex_info.workCost_[var] = simplex_lp.sense_ * simplex_lp.colCost_[col];
    simplex_info.workShift_[var] = 0.;
  }
}

void initialise_phase2_row_cost(HighsModelObject& highs_model_object,
                                int firstrow, int lastrow) {
  // Zero the cost and shift
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  for (int row = firstrow; row <= lastrow; row++) {
    int var = simplex_lp.numCol_ + row;
    simplex_info.workCost_[var] = 0;
    simplex_info.workShift_[var] = 0.;
  }
}

void initialise_cost(HighsModelObject& highs_model_object, int perturb) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Copy the cost
  initialise_phase2_col_cost(highs_model_object, 0, simplex_lp.numCol_ - 1);
  initialise_phase2_row_cost(highs_model_object, 0, simplex_lp.numRow_ - 1);
  // See if we want to skip perturbation
  simplex_info.costs_perturbed = 0;
  if (perturb == 0 || simplex_info.perturb_costs == 0) return;
  simplex_info.costs_perturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < simplex_lp.numCol_; i++)
    bigc = max(bigc, fabs(simplex_info.workCost_[i]));
  if (bigc > 100) bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++)
    boxedRate += (simplex_info.workRange_[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) bigc = min(bigc, 1.0);
  if (bigc < 1) {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < simplex_lp.numCol_; i++) {
    double lower = simplex_lp.colLower_[i];
    double upper = simplex_lp.colUpper_[i];
    double xpert = (fabs(simplex_info.workCost_[i]) + 1) * base *
                   (1 + simplex_info.numTotRandomValue_[i]);
    if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper == HIGHS_CONST_INF) {  // Lower
      simplex_info.workCost_[i] += xpert;
    } else if (lower == -HIGHS_CONST_INF) {  // Upper
      simplex_info.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex_info.workCost_[i] +=
          (simplex_info.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }

  for (int i = simplex_lp.numCol_; i < numTot; i++) {
    simplex_info.workCost_[i] +=
        (0.5 - simplex_info.numTotRandomValue_[i]) * 1e-12;
  }
}

int get_nonbasicMove(HighsModelObject& highs_model_object, int var) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  assert(var >= 0);
  assert(var < highs_model_object.simplex_lp_.numCol_ + highs_model_object.simplex_lp_.numRow_);
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

void populate_work_arrays(HighsModelObject& highs_model_object) {
  // Initialize the values
  initialise_cost(highs_model_object);
  initialise_bound(highs_model_object);
  initialise_value(highs_model_object);
}

void replace_with_logical_basis(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Replace basis with a logical basis then populate (where possible)
  // work* arrays
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_lp.numCol_ + row;
    simplex_basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    simplex_basis.basicIndex_[row] = var;
  }
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    simplex_basis.nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
  }
  simplex_info.num_basic_logicals = simplex_lp.numRow_;

  populate_work_arrays(highs_model_object);

  // Deduce the consequences of a new basis
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_BASIS);
}

void replace_with_new_basis(HighsModelObject& highs_model_object,
                            const int* XbasicIndex) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  // Replace basis with a new basis then populate (where possible)
  // work* arrays
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; var++) {
    simplex_basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
  }
  simplex_info.num_basic_logicals = 0;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = XbasicIndex[row];
    if (var >= simplex_lp.numCol_) simplex_info.num_basic_logicals++;
    simplex_basis.basicIndex_[row] = var;
    simplex_basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
  }

  populate_work_arrays(highs_model_object);

  // Deduce the consequences of a new basis
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_BASIS);
}

void setup_num_basic_logicals(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  simplex_info.num_basic_logicals = 0;
  for (int i = 0; i < simplex_lp.numRow_; i++)
    if (simplex_basis.basicIndex_[i] >= simplex_lp.numCol_)
      simplex_info.num_basic_logicals += 1;
#ifdef HiGHSDEV
  printf("Determined num_basic_logicals = %d of %d\n",
         simplex_info.num_basic_logicals, simplex_lp.numRow_);
#endif
}

#ifdef HiGHSDEV
void reportSimplexProfiling(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexTimer simplex_timer;
  HighsTimer& timer = highs_model_object.timer_;

  if (simplex_info.simplex_strategy == SimplexStrategy::PRIMAL) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportSimplexInnerClock(highs_model_object);
    }
  } else if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_PLAIN) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportSimplexInnerClock(highs_model_object);
    }
    if (simplex_info.report_simplex_outer_clock) {
      simplex_timer.reportDualSimplexIterateClock(highs_model_object);
      simplex_timer.reportDualSimplexOuterClock(highs_model_object);
    }
  }

  //  if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
  //    int reportList[] = {
  //        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
  //        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
  //        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
  //        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
  //        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
  //        HTICK_GROUP1};
  //    int reportCount = sizeof(reportList) / sizeof(int);
  //    timer.report(reportCount, reportList, 0.0);
  //  }

  if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
    //    int reportList[] = {
    //        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
    //        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
    //        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
    //        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
    //        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
    //        HTICK_UPDATE_ROW_EP};
    //    int reportCount = sizeof(reportList) / sizeof(int);
    //    timer.report(reportCount, reportList, 0.0);
    printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
           highs_model_object.lp_.model_name_.c_str(), simplex_info.pami_cutoff,
           simplex_info.iteration_count / (1.0 + simplex_info.multi_iteration));
  }

  if (simplex_info.report_simplex_phases_clock) {
    simplex_timer.reportSimplexTotalClock(highs_model_object);
    simplex_timer.reportSimplexPhasesClock(highs_model_object);
  }

  if (simplex_info.analyse_invert_time) {
    double current_run_highs_time = timer.readRunHighsClock();
    int iClock = simplex_info.clock_[InvertClock];
    simplex_info.total_inverts = timer.clock_num_call[iClock];
    simplex_info.total_invert_time = timer.clock_time[iClock];

    printf(
        "Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time "
        "= %11.4g",
        simplex_info.total_inverts, simplex_info.total_invert_time,
        current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n",
             (100 * simplex_info.total_invert_time) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
  if (simplex_info.analyseRebuildTime) {
    double current_run_highs_time = timer.readRunHighsClock();
    HighsClockRecord totalRebuildClock;
    timer.clockInit(totalRebuildClock);
    timer.clockAdd(totalRebuildClock,
                   simplex_info.clock_[IterateDualRebuildClock]);
    timer.clockAdd(totalRebuildClock,
                   simplex_info.clock_[IteratePrimalRebuildClock]);
    int totalRebuilds = 0;
    double totalRebuildTime = 0;
    printf("Time: Total rebuild time = %11.4g (%4d) of Total time = %11.4g",
           totalRebuildTime, totalRebuilds, current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n", (100 * totalRebuildTime) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
}
#endif

double computeBasisCondition(HighsModelObject& highs_model_object) {
  int solver_num_row = highs_model_object.simplex_lp_.numRow_;
  int solver_num_col = highs_model_object.simplex_lp_.numCol_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  HFactor& factor = highs_model_object.factor_;
  const int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
  const double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
  // Compute the Hager condition number estimate for the basis matrix
  double NoDensity = 1;
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
  row_ep.count = solver_num_row;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    row_ep.index[r_n] = r_n;
    row_ep.array[r_n] = bs_cond_x[r_n];
  }
  for (int ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor.ftran(row_ep, NoDensity);
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
    row_ep.count = solver_num_row;
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      row_ep.index[r_n] = r_n;
      row_ep.array[r_n] = bs_cond_w[r_n];
    }
    row_ep.packFlag = false;
    factor.btran(row_ep, NoDensity);
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
    int vr_n = highs_model_object.simplex_basis_.basicIndex_[r_n];
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

bool work_arrays_ok(HighsModelObject& highs_model_object, int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  //  printf("Called work_arrays_ok(%d)\n", phase);cout << flush;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == simplex_lp.colLower_[col];
        if (!ok) {
          printf("For col %d, simplex_info.workLower_ should be %g but is %g\n",
                 col, simplex_lp.colLower_[col], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == simplex_lp.colUpper_[col];
        if (!ok) {
          printf("For col %d, simplex_info.workUpper_ should be %g but is %g\n",
                 col, simplex_lp.colUpper_[col], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == -simplex_lp.rowUpper_[row];
        if (!ok) {
          printf("For row %d, simplex_info.workLower_ should be %g but is %g\n",
                 row, -simplex_lp.rowUpper_[row], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == -simplex_lp.rowLower_[row];
        if (!ok) {
          printf("For row %d, simplex_info.workUpper_ should be %g but is %g\n",
                 row, -simplex_lp.rowLower_[row], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
  }
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    ok = simplex_info.workRange_[var] ==
         (simplex_info.workUpper_[var] - simplex_info.workLower_[var]);
    if (!ok) {
      printf(
          "For variable %d, simplex_info.workRange_ should be %g = %g - %g "
          "but is %g\n",
          var, simplex_info.workUpper_[var] - simplex_info.workLower_[var],
          simplex_info.workUpper_[var], simplex_info.workLower_[var],
          simplex_info.workRange_[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!simplex_info.costs_perturbed) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      ok = simplex_info.workCost_[var] ==
           simplex_lp.sense_ * simplex_lp.colCost_[col];
      if (!ok) {
        printf("For col %d, simplex_info.workLower_ should be %g but is %g\n",
               col, simplex_lp.colLower_[col], simplex_info.workCost_[var]);
        return ok;
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      ok = simplex_info.workCost_[var] == 0.;
      if (!ok) {
        printf("For row %d, simplex_info.workCost_ should be zero but is %g\n",
               row, simplex_info.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool one_nonbasic_move_vs_work_arrays_ok(HighsModelObject& highs_model_object,
                                         int var) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  assert(var >= 0);
  assert(var < simplex_lp.numCol_ + simplex_lp.numRow_);
  // Make sure we're not checking a basic variable
  if (!simplex_basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-simplex_info.workLower_[var])) {
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
        // Fixed variable
        ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          printf(
              "Fixed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] so nonbasic "
              "move should be zero but is %d\n",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
        if (!ok) {
          printf(
              "Fixed variable %d (simplex_lp.numCol_ = %d) so "
              "simplex_info.work value should be %g but "
              "is %g\n",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          printf(
              "Boxed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %d\n",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_info.workUpper_[var] - simplex_info.workLower_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                "NONBASIC_MOVE_UP so work "
                "value should be %g but is %g\n",
                var, simplex_lp.numCol_, simplex_info.workLower_[var],
                simplex_info.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                "NONBASIC_MOVE_DN so work "
                "value should be %g but is %g\n",
                var, simplex_lp.numCol_, simplex_info.workUpper_[var],
                simplex_info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d\n",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            NONBASIC_MOVE_UP, simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d\n",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, simplex_lp.numCol_, simplex_info.workUpper_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        printf(
            "Free variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, %11g] "
            "so nonbasic "
            "move should be zero but is  %d\n",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == 0.0;
      if (!ok) {
        printf(
            "Free variable %d (simplex_lp.numCol_ = %d) so work value should "
            "be zero but "
            "is %g\n",
            var, simplex_lp.numCol_, simplex_info.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool all_nonbasic_move_vs_work_arrays_ok(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  //    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  bool ok;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    printf(
        "NonbasicMoveVsWorkArrays: var = %2d; simplex_basis.nonbasicFlag_[var] "
        "= %2d\n",
        var, simplex_basis.nonbasicFlag_[var]);
    if (!simplex_basis.nonbasicFlag_[var]) continue;
    ok = one_nonbasic_move_vs_work_arrays_ok(highs_model_object, var);
    if (!ok) {
      printf("Error in NonbasicMoveVsWorkArrays for nonbasic variable %d\n",
             var);
      assert(ok);
      return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool ok_to_solve(HighsModelObject& highs_model_object, int level, int phase) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  //  printf("Called ok_to_solve(%1d, %1d)\n", level, phase);
  bool ok;
  // Level 0: Minimal check - just look at flags. This means we trust them!
  ok = simplex_lp_status.has_basis && simplex_lp_status.has_matrix_col_wise &&
       simplex_lp_status.has_matrix_row_wise &&
       simplex_lp_status.has_factor_arrays &&
       simplex_lp_status.has_dual_steepest_edge_weights &&
       simplex_lp_status.has_invert;
  // TODO: Eliminate the following line ASAP!!!
  ok = true;
  if (!ok) {
    if (!simplex_lp_status.has_basis)
      printf("Not OK to solve since simplex_lp_status.has_basis = %d\n",
             simplex_lp_status.has_basis);
    if (!simplex_lp_status.has_matrix_col_wise)
      printf(
          "Not OK to solve since simplex_lp_status.has_matrix_col_wise "
          "= %d\n",
          simplex_lp_status.has_matrix_col_wise);
    if (!simplex_lp_status.has_matrix_row_wise)
      printf(
          "Not OK to solve since simplex_lp_status.has_matrix_row_wise "
          "= %d\n",
          simplex_lp_status.has_matrix_row_wise);
    //    if (!simplex_lp_status.has_factor_arrays)
    //      printf("Not OK to solve since
    //      simplex_lp_status.has_factor_arrays = %d\n",
    //             simplex_lp_status.has_factor_arrays);
    if (!simplex_lp_status.has_dual_steepest_edge_weights)
      printf(
          "Not OK to solve since "
          "simplex_lp_status.has_dual_steepest_edge_weights = %d\n",
          simplex_lp_status.has_dual_steepest_edge_weights);
    if (!simplex_lp_status.has_invert)
      printf("Not OK to solve since simplex_lp_status.has_invert = %d\n",
             simplex_lp_status.has_invert);
  }
  assert(ok);
  if (level <= 0) return ok;
  // Level 1: Basis and data check
  ok = nonbasic_flag_basic_index_ok(simplex_lp,
                                    highs_model_object.simplex_basis_);
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
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    if (simplex_basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      ok = one_nonbasic_move_vs_work_arrays_ok(highs_model_object, var);
      if (!ok) {
        printf("Error in nonbasicMoveVsWorkArrays for variable %d of %d\n", var,
               numTot);
        assert(ok);
        return ok;
      }
    }
  }
  return ok;
}

void flip_bound(HighsModelObject& highs_model_object, int iCol) {
  int* nonbasicMove = &highs_model_object.simplex_basis_.nonbasicMove_[0];
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  simplex_info.workValue_[iCol] =
      move == 1 ? simplex_info.workLower_[iCol] : simplex_info.workUpper_[iCol];
}
/*
int handle_rank_deficiency(HighsModelObject &highs_model_object) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HFactor &factor = highs_model_object.factor_;
  SimplexBasis &simplex_basis = highs_model_object.simplex_basis_;
  int rankDeficiency = factor.rankDeficiency;
  const int *noPvC = factor.getNoPvC();
  printf("Returned %d = factor.build();\n", rankDeficiency);
  fflush(stdout);
  vector<int> basicRows;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  basicRows.resize(numTot);
  //    printf("Before - simplex_basis.basicIndex_:"); for (int iRow=0;
iRow<simplex_lp.numRow_; iRow++)
  //    printf(" %2d", simplex_basis.basicIndex_[iRow]); printf("\n");
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
basicRows[simplex_basis.basicIndex_[iRow]] = iRow; for (int k = 0; k <
rankDeficiency; k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = simplex_lp.numCol_ + factor.noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut,
simplex_info.workLower_[columnOut],
    //      simplex_info.workUpper_[columnOut]);
    if (simplex_basis.basicIndex_[rowOut] != columnOut) {
      printf("%d = simplex_basis.basicIndex_[rowOut] != noPvC[k] = %d\n",
simplex_basis.basicIndex_[rowOut], columnOut); fflush(stdout);
    }
    int sourceOut = setSourceOutFmBd(columnOut);
    updatePivots(columnIn, rowOut, sourceOut);
    updateMatrix(columnIn, columnOut);
  }
  //    printf("After  - simplex_basis.basicIndex_:"); for (int iRow=0;
iRow<simplex_lp.numRow_; iRow++)
  //    printf(" %2d", simplex_basis.basicIndex_[iRow]); printf("\n");
#ifdef HiGHSDEV
  factor.checkInvert();
#endif
  return 0;
}
*/
int compute_factor(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HFactor& factor = highs_model_object.factor_;
#ifdef HiGHSDEV
  HighsTimer& timer = highs_model_object.timer_;
  double tt0 = 0;
  int iClock = simplex_info.clock_[InvertClock];
  if (simplex_info.analyse_invert_time) tt0 = timer.clock_time[iClock];
#endif
  // TODO Understand why handling noPvC and noPvR in what seem to be
  // different ways ends up equivalent.
  int rankDeficiency = factor.build();
  if (rankDeficiency) {
    //    handle_rank_deficiency();
    //    simplex_lp_status.solution_status = SimplexSolutionStatus::SINGULAR;
#ifdef HiGHSDEV
    //    writePivots("failed");
#endif
    //      return rankDeficiency;
  }
  //    printf("INVERT: After %d iterations and %d updates\n",
  //    simplex_info.iteration_count, simplex_info.update_count);
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
        simplex_info.total_inverts, simplex_info.iteration_count, invertTime,
        simplex_info.total_invert_time);
  }
#endif

  // Now have a representation of B^{-1}, and it is fresh!
  simplex_lp_status.has_invert = true;
  simplex_lp_status.has_fresh_invert = true;
  return 0;
}

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void compute_primal(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HMatrix& matrix = highs_model_object.matrix_;
  HFactor& factor = highs_model_object.factor_;
  // Setup a local buffer for the values of basic variables
  HVector buffer;
  buffer.setup(simplex_lp.numRow_);
  buffer.clear();
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    //    if (simplex_basis.nonbasicFlag_[i] && simplex_info.workValue_[i] == 0)
    //    printf("\nAfter adding in %12g*a[%2d]\n", simplex_info.workValue_[i],
    //    i);
    if (simplex_basis.nonbasicFlag_[i] && simplex_info.workValue_[i] != 0) {
      matrix.collect_aj(buffer, i, simplex_info.workValue_[i]);
      //      printf("\nAfter adding in %12g*a[%2d]\nRow     value\n",
      //      simplex_info.workValue_[i], i); for (int iRow = 0; iRow <
      //      simplex_lp.numRow_; iRow++) printf("%3d %12g\n", iRow,
      //      buffer.array[iRow]);
    }
  }
  //  for (int i = 0; i < simplex_lp.numRow_; i++) printf("Bf FTRAN: Row %2d has
  //  value %12g\n", i, buffer.array[i]);
  factor.ftran(buffer, 1);

  for (int i = 0; i < simplex_lp.numRow_; i++) {
    int iCol = simplex_basis.basicIndex_[i];
    simplex_info.baseValue_[i] = -buffer.array[i];
    simplex_info.baseLower_[i] = simplex_info.workLower_[iCol];
    simplex_info.baseUpper_[i] = simplex_info.workUpper_[iCol];
  }
  // Now have basic primals
  simplex_lp_status.has_basic_primal_values = true;
}

void computePrimalInfeasible(HighsModelObject& highs_model_object,
                             const bool report) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int num_nonbasic_primal_infeasibilities = 0;
  int num_basic_primal_infeasibilities = 0;
  double max_nonbasic_primal_infeasibility = 0;
  double max_basic_primal_infeasibility = 0;
  double sum_nonbasic_primal_infeasibilities = 0;
  double sum_basic_primal_infeasibilities = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  //  int nonbasic_ix = 0;
  for (int i = 0; i < numTot; i++) {
    if (simplex_basis.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = simplex_info.workValue_[i];
      double lower = simplex_info.workLower_[i];
      double upper = simplex_info.workUpper_[i];
      double primal_infeasibility = max(lower - value, value - upper);
      //      printf("Nonbasic column %2d is %2d, [%12g, %12g, %12g]
      //      primal_infeasibility = %12g\n", nonbasic_ix, i, lower, value,
      //      upper, primal_infeasibility); nonbasic_ix++;
      if (primal_infeasibility > 0) {
	if (primal_infeasibility > simplex_info.primal_feasibility_tolerance) num_nonbasic_primal_infeasibilities++;
	max_nonbasic_primal_infeasibility =
	  std::max(primal_infeasibility, max_nonbasic_primal_infeasibility);
	sum_nonbasic_primal_infeasibilities += primal_infeasibility;
      }
    }
  }
  for (int i = 0; i < simplex_lp.numRow_; i++) {
    // Basic variable
    double value = simplex_info.baseValue_[i];
    double lower = simplex_info.baseLower_[i];
    double upper = simplex_info.baseUpper_[i];
    double primal_infeasibility = max(lower - value, value - upper);
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > simplex_info.primal_feasibility_tolerance) num_basic_primal_infeasibilities++;
      max_basic_primal_infeasibility =
	std::max(primal_infeasibility, max_basic_primal_infeasibility);
      sum_basic_primal_infeasibilities += primal_infeasibility;
    }
  }
  int num_primal_infeasibilities =
      num_nonbasic_primal_infeasibilities + num_basic_primal_infeasibilities;
  double max_primal_infeasibility = std::max(max_nonbasic_primal_infeasibility,
                                             max_basic_primal_infeasibility);
  double sum_primal_infeasibilities =
      sum_nonbasic_primal_infeasibilities + sum_basic_primal_infeasibilities;
  if (report) {
#ifdef HiGHSDEV
    if (num_primal_infeasibilities) {
      int num_iter = simplex_info.iteration_count;
      printf(
	     "Iter %d has %d (%d+%d) primal infeasibilities (max = %g = max[%g, "
	     "%g]) summing to %g (%g+%g)\n",
	     num_iter, num_primal_infeasibilities,
	     num_nonbasic_primal_infeasibilities, num_basic_primal_infeasibilities,
	     max_primal_infeasibility, max_nonbasic_primal_infeasibility,
	     max_basic_primal_infeasibility, sum_primal_infeasibilities,
	     sum_nonbasic_primal_infeasibilities, sum_basic_primal_infeasibilities);
    }
#endif
  }
  simplex_info.num_primal_infeasibilities = num_primal_infeasibilities;
  simplex_info.max_primal_infeasibility = max_primal_infeasibility;
  simplex_info.sum_primal_infeasibilities = sum_primal_infeasibilities;
}

void computeDualInfeasible(HighsModelObject& highs_model_object,
                           const bool report) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int num_fixed_variable_move_errors = 0;
  int num_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibilities = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    double lower = simplex_info.workLower_[iVar];
    double upper = simplex_info.workUpper_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(simplex_info.workDual_[iVar]);
    } else {
      // Not fixed: any dual infeasibility is given by value signed by
      // nonbasicMove. This assumes that nonbasicMove=0 for fixed
      // variables
      dual_infeasibility =
          -simplex_basis.nonbasicMove_[iVar] * simplex_info.workDual_[iVar];
      if (lower == upper && simplex_basis.nonbasicMove_[iVar]) num_fixed_variable_move_errors++;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= simplex_info.dual_feasibility_tolerance) num_dual_infeasibilities++;
      max_dual_infeasibility =
	std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  // Check that there are no fixed variables with nonzero nonbasicMove
  if (num_fixed_variable_move_errors) {
    HighsLogMessage(HighsMessageType::ERROR,
		    "In computeDualInfeasible there are %d fixed variables with nonzero nonbasicMove",
		    num_fixed_variable_move_errors);
  }
  assert(num_fixed_variable_move_errors==0);

  if (report) {
#ifdef HiGHSDEV
    if (num_dual_infeasibilities) {
      int num_iter = simplex_info.iteration_count;
      printf("Iter %d has %d dual infeasibilities (max = %g) summing to %g\n",
	     num_iter, num_dual_infeasibilities, max_dual_infeasibility,
	     sum_dual_infeasibilities);
    }
#endif
  }
  simplex_info.num_dual_infeasibilities = num_dual_infeasibilities;
  simplex_info.max_dual_infeasibility = max_dual_infeasibility;
  simplex_info.sum_dual_infeasibilities = sum_dual_infeasibilities;
}

void computeDualInfeasibleWithFlips(HighsModelObject& highs_model_object,
				    const bool report) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int num_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibilities = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    double lower = simplex_info.workLower_[iVar];
    double upper = simplex_info.workUpper_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(simplex_info.workDual_[iVar]);
    } else if (highs_isInfinity(-lower) || highs_isInfinity(upper)) {
      // Not boxed: any dual infeasibility is given by value signed by
      // nonbasicMove
      dual_infeasibility =
          -simplex_basis.nonbasicMove_[iVar] * simplex_info.workDual_[iVar];
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= simplex_info.dual_feasibility_tolerance) num_dual_infeasibilities++;
      max_dual_infeasibility =
	std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  if (report) {
#ifdef HiGHSDEV
    if (num_dual_infeasibilities) {
      int num_iter = simplex_info.iteration_count;
      printf("Iter %d has %d dual infeasibilities (max = %g) summing to %g\n",
	     num_iter, num_dual_infeasibilities, max_dual_infeasibility,
	     sum_dual_infeasibilities);
    }
#endif
  }
  simplex_info.num_dual_infeasibilities = num_dual_infeasibilities;
  simplex_info.max_dual_infeasibility = max_dual_infeasibility;
  simplex_info.sum_dual_infeasibilities = sum_dual_infeasibilities;
}

void compute_dual(HighsModelObject& highs_model_object) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HMatrix& matrix = highs_model_object.matrix_;
  HFactor& factor = highs_model_object.factor_;
  bool an_compute_dual_norm2 = false;
  double btran_rhs_norm2;
  double btran_sol_norm2;
  double work_dual_norm2;

  // Create a local buffer for the pi vector
  HVector buffer;
  buffer.setup(simplex_lp.numRow_);
  buffer.clear();
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    buffer.index[iRow] = iRow;
    buffer.array[iRow] =
        simplex_info.workCost_[simplex_basis.basicIndex_[iRow]] +
        simplex_info.workShift_[simplex_basis.basicIndex_[iRow]];
  }
  buffer.count = simplex_lp.numRow_;
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
  bufferLong.setup(simplex_lp.numCol_);
  bufferLong.clear();
  matrix.price_by_col(bufferLong, buffer);
  for (int i = 0; i < simplex_lp.numCol_; i++) {
    simplex_info.workDual_[i] = simplex_info.workCost_[i] - bufferLong.array[i];
  }
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = simplex_lp.numCol_; i < numTot; i++) {
    simplex_info.workDual_[i] =
        simplex_info.workCost_[i] - buffer.array[i - simplex_lp.numCol_];
  }

  if (an_compute_dual_norm2) {
    work_dual_norm2 = 0;
    for (int i = 0; i < numTot; i++)
      work_dual_norm2 += simplex_info.workDual_[i] * simplex_info.workDual_[i];
    work_dual_norm2 = sqrt(work_dual_norm2);
    //  printf("compute_dual: B.pi=c_B has ||c_B||=%11.4g; ||pi||=%11.4g;
    //  ||pi^TA-c||=%11.4g\n", btran_rhs_norm2, btran_sol_norm2,
    //  work_dual_norm2);
    double current_dual_feasibility_tolerance =
        simplex_info.dual_feasibility_tolerance;
    double new_dual_feasibility_tolerance = work_dual_norm2 / 1e16;
    if (new_dual_feasibility_tolerance > 1e-1) {
      printf(
          "Seriously: do you expect to solve an LP with ||pi^TA-c||=%11.4g?\n",
          work_dual_norm2);
    } else if (new_dual_feasibility_tolerance >
               10 * current_dual_feasibility_tolerance) {
      printf(
          "||pi^TA-c|| = %12g so solving with dual_feasibility_tolerance = "
          "%12g\n",
          work_dual_norm2, new_dual_feasibility_tolerance);
      simplex_info.dual_feasibility_tolerance = new_dual_feasibility_tolerance;
    }
  }

  // Now have nonbasic duals
  simplex_lp_status.has_nonbasic_dual_values = true;
}

void correct_dual(HighsModelObject& highs_model_object,
                  int* free_infeasibility_count) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsRandom& random = highs_model_object.random_;
  const double tau_d = simplex_info.dual_feasibility_tolerance;
  const double inf = HIGHS_CONST_INF;
  int workCount = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_basis.nonbasicFlag_[i]) {
      if (simplex_info.workLower_[i] == -inf &&
          simplex_info.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(simplex_info.workDual_[i]) >= tau_d);
      } else if (simplex_basis.nonbasicMove_[i] * simplex_info.workDual_[i] <=
                 -tau_d) {
        if (simplex_info.workLower_[i] != -inf &&
            simplex_info.workUpper_[i] != inf) {
          // Boxed variable = flip
          flip_bound(highs_model_object, i);
        } else {
          // Other variable = shift
          simplex_info.costs_perturbed = 1;
          if (simplex_basis.nonbasicMove_[i] == 1) {
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

// Record the shift in the cost of a particular column
void shift_cost(HighsModelObject& highs_model_object, int iCol, double amount) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  simplex_info.costs_perturbed = 1;
  assert(simplex_info.workShift_[iCol] == 0);
  simplex_info.workShift_[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void shift_back(HighsModelObject& highs_model_object, int iCol) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  simplex_info.workDual_[iCol] -= simplex_info.workShift_[iCol];
  simplex_info.workShift_[iCol] = 0;
}

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void update_factor(HighsModelObject& highs_model_object, HVector* column,
                   HVector* row_ep, int* iRow, int* hint) {
  //    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HFactor& factor = highs_model_object.factor_;
  HighsTimer& timer = highs_model_object.timer_;

  timer.start(simplex_info.clock_[UpdateFactorClock]);
  factor.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  simplex_lp_status.has_invert = true;
  if (simplex_info.update_count >= simplex_info.update_limit)
    *hint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  timer.stop(simplex_info.clock_[UpdateFactorClock]);
}

void update_pivots(HighsModelObject& highs_model_object, int columnIn,
                   int rowOut, int sourceOut) {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsTimer& timer = highs_model_object.timer_;

  timer.start(simplex_info.clock_[UpdatePivotsClock]);
  int columnOut = simplex_basis.basicIndex_[rowOut];

  // Incoming variable
  simplex_basis.basicIndex_[rowOut] = columnIn;
  simplex_basis.nonbasicFlag_[columnIn] = 0;
  simplex_basis.nonbasicMove_[columnIn] = 0;
  simplex_info.baseLower_[rowOut] = simplex_info.workLower_[columnIn];
  simplex_info.baseUpper_[rowOut] = simplex_info.workUpper_[columnIn];

  // Outgoing variable
  simplex_basis.nonbasicFlag_[columnOut] = 1;
  //  double dlValue;
  //  double vrLb = simplex_info.workLower_[columnOut];
  //  double vrV = simplex_info.workValue_[columnOut];
  //  double vrUb = simplex_info.workUpper_[columnOut];
  if (simplex_info.workLower_[columnOut] ==
      simplex_info.workUpper_[columnOut]) {
    //    dlValue =
    //    simplex_info.workLower_[columnOut]-simplex_info.workValue_[columnOut];
    simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = 0;
  } else if (sourceOut == -1) {
    //    dlValue =
    //    simplex_info.workLower_[columnOut]-simplex_info.workValue_[columnOut];
    simplex_info.workValue_[columnOut] = simplex_info.workLower_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = 1;
  } else {
    //    dlValue =
    //    simplex_info.workUpper_[columnOut]-simplex_info.workValue_[columnOut];
    simplex_info.workValue_[columnOut] = simplex_info.workUpper_[columnOut];
    simplex_basis.nonbasicMove_[columnOut] = -1;
  }
  double nwValue = simplex_info.workValue_[columnOut];
  double vrDual = simplex_info.workDual_[columnOut];
  double dl_dual_objective_value = nwValue * vrDual;
  //  if (fabs(nwValue))
  //    printf("update_pivots columnOut = %6d (%2d): [%11.4g, %11.4g, %11.4g],
  //    nwValue = %11.4g, dual = %11.4g, dlObj = %11.4g\n",
  //			   columnOut, simplex_basis.nonbasicMove_[columnOut], vrLb, vrV,
  //vrUb, nwValue, vrDual, dl_dual_objective_value);
  simplex_info.updated_dual_objective_value += dl_dual_objective_value;
  simplex_info.update_count++;
  // Update the number of basic logicals
  if (columnOut < simplex_lp.numCol_) simplex_info.num_basic_logicals -= 1;
  if (columnIn < simplex_lp.numCol_) simplex_info.num_basic_logicals += 1;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  simplex_lp_status.has_invert = false;
  simplex_lp_status.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = false;
  timer.stop(simplex_info.clock_[UpdatePivotsClock]);
}

void update_matrix(HighsModelObject& highs_model_object, int columnIn,
                   int columnOut) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HMatrix& matrix = highs_model_object.matrix_;
  HighsTimer& timer = highs_model_object.timer_;

  timer.start(simplex_info.clock_[UpdateMatrixClock]);
  matrix.update(columnIn, columnOut);
  timer.stop(simplex_info.clock_[UpdateMatrixClock]);
}

void logRebuild(HighsModelObject& highs_model_object, const bool primal,
                const int solve_phase) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  double objective_value;
  string simplex_variant;
  if (primal) {
    simplex_variant = "Pr";
    objective_value = simplex_info.primal_objective_value;
  } else {
    simplex_variant = "Du";
    objective_value = simplex_info.dual_objective_value;
  }
  if (solve_phase < 2) {
    HighsLogMessage(HighsMessageType::INFO, "Iter %10d: %20.10e %sPh%1d",
                    simplex_info.iteration_count, objective_value,
                    simplex_variant.c_str(), solve_phase);
  } else if (!primal && simplex_info.sum_dual_infeasibilities == 0) {
    HighsLogMessage(
        HighsMessageType::INFO, "Iter %10d: %20.10e %sPh%1d Pr: %d(%g)",
        simplex_info.iteration_count, objective_value, simplex_variant.c_str(),
        solve_phase, simplex_info.num_primal_infeasibilities,
        simplex_info.sum_primal_infeasibilities);
  } else {
    HighsLogMessage(HighsMessageType::INFO,
                    "Iter %10d: %20.10e %sPh%1d Pr: %d(%g); Du: %d(%g)",
                    simplex_info.iteration_count, objective_value,
                    simplex_variant.c_str(), solve_phase,
                    simplex_info.num_primal_infeasibilities,
                    simplex_info.sum_primal_infeasibilities,
                    simplex_info.num_dual_infeasibilities,
                    simplex_info.sum_dual_infeasibilities);
  }
}

// Return a string representation of SimplexSolutionStatus.
std::string SimplexSolutionStatusToString(SimplexSolutionStatus status) {
  switch (status) {
    case SimplexSolutionStatus::UNSET:
      return "Unset";
      break;
    case SimplexSolutionStatus::OPTIMAL:
      return "Optimal";
      break;
    case SimplexSolutionStatus::PRIMAL_FEASIBLE:
      return "Primal feasible";
      break;
    case SimplexSolutionStatus::DUAL_FEASIBLE:
      return "Dual feasible";
      break;
    case SimplexSolutionStatus::INFEASIBLE:
      return "Infeasible";
      break;
    case SimplexSolutionStatus::UNBOUNDED:
      return "Primal unbounded";
      break;
    case SimplexSolutionStatus::SINGULAR:
      return "Singular basis";
      break;
    case SimplexSolutionStatus::FAILED:
      return "Failed";
      break;
    case SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return "Reached dual objective value upper bound";
      break;
    case SimplexSolutionStatus::OUT_OF_TIME:
      return "Time limit exceeded";
      break;
  }
  return "";
}

void reportSimplexLpStatus(HighsSimplexLpStatus& simplex_lp_status,
                           const char* message) {
  printf("\nReporting solver status and flags: %s\n\n", message);
  printf("  valid =                          %d\n", simplex_lp_status.valid);
  printf("  is_dualised =                    %d\n",
         simplex_lp_status.is_dualised);
  printf("  is_permuted =                    %d\n",
         simplex_lp_status.is_permuted);
  printf("  is_scaled =                      %d\n",
         simplex_lp_status.scaling_tried);
  printf("  has_basis =                      %d\n",
         simplex_lp_status.has_basis);
  printf("  has_matrix_col_wise =            %d\n",
         simplex_lp_status.has_matrix_col_wise);
  printf("  has_matrix_row_wise =            %d\n",
         simplex_lp_status.has_matrix_row_wise);
  printf("  has_factor_arrays =              %d\n",
         simplex_lp_status.has_factor_arrays);
  printf("  has_dual_steepest_edge_weights = %d\n",
         simplex_lp_status.has_dual_steepest_edge_weights);
  printf("  has_nonbasic_dual_values =       %d\n",
         simplex_lp_status.has_nonbasic_dual_values);
  printf("  has_basic_primal_values =        %d\n",
         simplex_lp_status.has_basic_primal_values);
  printf("  has_invert =                     %d\n",
         simplex_lp_status.has_invert);
  printf("  has_fresh_invert =               %d\n",
         simplex_lp_status.has_fresh_invert);
  printf("  has_fresh_rebuild =              %d\n",
         simplex_lp_status.has_fresh_rebuild);
  printf("  has_dual_objective_value =       %d\n",
         simplex_lp_status.has_dual_objective_value);
  printf("  has_primal_objective_value =     %d\n",
         simplex_lp_status.has_primal_objective_value);
}

void invalidateSimplexLpData(HighsSimplexLpStatus& simplex_lp_status) {
  simplex_lp_status.has_basis = false;
  simplex_lp_status.has_matrix_col_wise = false;
  simplex_lp_status.has_matrix_row_wise = false;
  simplex_lp_status.has_factor_arrays = false;
  simplex_lp_status.has_dual_steepest_edge_weights = false;
  simplex_lp_status.has_nonbasic_dual_values = false;
  simplex_lp_status.has_basic_primal_values = false;
  simplex_lp_status.has_invert = false;
  simplex_lp_status.has_fresh_invert = false;
  simplex_lp_status.has_fresh_rebuild = false;
  simplex_lp_status.has_dual_objective_value = false;
  simplex_lp_status.has_primal_objective_value = false;
}

void invalidateSimplexLp(HighsSimplexLpStatus& simplex_lp_status) {
  simplex_lp_status.valid = false;
  simplex_lp_status.is_dualised = false;
  simplex_lp_status.is_permuted = false;
  simplex_lp_status.scaling_tried = false;
  invalidateSimplexLpData(simplex_lp_status);
}

void updateSimplexLpStatus(HighsSimplexLpStatus& simplex_lp_status,
                           LpAction action) {
  switch (action) {
    case LpAction::DUALISE:
#ifdef HIGHSDEV
      printf(" LpAction::DUALISE\n");
#endif
      simplex_lp_status.is_dualised = true;
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::PERMUTE:
#ifdef HIGHSDEV
      printf(" LpAction::PERMUTE\n");
#endif
      simplex_lp_status.is_permuted = true;
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::SCALE:
#ifdef HIGHSDEV
      printf(" LpAction::SCALE\n");
#endif
      simplex_lp_status.scaling_tried = true;
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::NEW_COSTS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COSTS\n");
#endif
      //      initCost();
      simplex_lp_status.has_nonbasic_dual_values = false;
      simplex_lp_status.has_fresh_rebuild = false;
      simplex_lp_status.has_dual_objective_value = false;
      simplex_lp_status.has_primal_objective_value = false;
      break;
    case LpAction::NEW_BOUNDS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BOUNDS\n");
#endif
      //      simplex_info.simplex_lp_ = true;
      //     initBound();
      //     initValue();
      simplex_lp_status.has_basic_primal_values = false;
      simplex_lp_status.has_fresh_rebuild = false;
      simplex_lp_status.has_dual_objective_value = false;
      simplex_lp_status.has_primal_objective_value = false;
      break;
    case LpAction::NEW_BASIS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_BASIS\n");
#endif
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::NEW_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_COLS\n");
#endif
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::NEW_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::NEW_ROWS\n");
#endif
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::DEL_COLS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_COLS\n");
#endif
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::DEL_ROWS:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS\n");
#endif
      invalidateSimplexLpData(simplex_lp_status);
      break;
    case LpAction::DEL_ROWS_BASIS_OK:
#ifdef HIGHSDEV
      printf(" LpAction::DEL_ROWS_BASIS_OK\n");
#endif
      //      simplex_info.simplex_lp_ = true;
      break;
    default:
#ifdef HIGHSDEV
      printf(" Unrecognised LpAction::%d\n", (int)action);
#endif
      break;
  }
}
