/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.cpp
 * @brief
 * @author Julian Hall
 */
#include "pdlp/HiPdlpWrapper.h"

#include "pdlp/hipdlp/logger.hpp"
#include "pdlp/hipdlp/pdhg.hpp"
#include "pdlp/hipdlp/restart.hpp"

void getHiPpdlpParamsFromOptions(const HighsOptions& options, HighsTimer& timer,
                                 PrimalDualParams& params);

HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object) {
  return solveLpHiPdlp(solver_object.options_, solver_object.timer_,
                       solver_object.lp_, solver_object.basis_,
                       solver_object.solution_, solver_object.model_status_,
                       solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  std::string log_filename = "";
  LogLevel log_level;
  if (options.log_dev_level == kHighsLogDevLevelInfo) {
    log_level = LogLevel::kInfo;
  } else if (options.log_dev_level == kHighsLogDevLevelDetailed) {
    log_level = LogLevel::kVerbose;
  } else if (options.log_dev_level == kHighsLogDevLevelVerbose) {
    log_level = LogLevel::kDebug;
  } else {
    log_level = LogLevel::kNone;
  }
  // --- Initialize Logger ---
  Logger logger(log_level);
  if (!log_filename.empty()) logger.set_log_file(log_filename);
  logger.print_header();

  Timer total_timer;
  /*** Order of operations 
   * Preprocess with HiPdlp
   * Scale with HiPdlp
   * Solve with HiPdlp
   * Unscale with HiPdlp
   * Postprocess with HiPDLP
   * ***/
  // 2. Preprocess with HiPdlp
  PDLPSolver pdlp(logger);
  pdlp.setParams(options, timer);
  HighsLp preprocessed_lp;
  pdlp.passLp(&lp);
  //logger_.info("Preprocessing LP to handle ranged constraints...");
  pdlp.PreprocessLp();

  // 3. Scale with HiPdlp
  pdlp.scaling_.ScaleProblem();

  // 4. Solve with HiPdlp
  std::vector<double> x, y;
  //pdlp.Solve(preprocessed_lp, params, x, y);

  // 5. Unscale with HiPdlp
  pdlp.scaling_.UnscaleSolution(x, y);

  // 6. Postprocess with HiPDLP
  HighsSolution pdlp_solution;
  //pdlp.Postsolve(presolved_lp, preprocessed_lp, x, y, pdlp_solution);


  // --- Print Summary ---
  logger.print_summary(pdlp.GetResults(), pdlp.GetIterationCount(),
                       total_timer.read());

  highs_info.pdlp_iteration_count = pdlp.GetIterationCount();

  model_status = HighsModelStatus::kUnknown;
  highs_solution.clear();
  highs_basis.valid = false;

  const TerminationStatus termination_status = pdlp.GetResults().term_code;
  switch (termination_status) {
    case TerminationStatus::OPTIMAL: {
      model_status = HighsModelStatus::kOptimal;
      break;
    }
    case TerminationStatus::INFEASIBLE: {
      assert(111 == 222);
      model_status = HighsModelStatus::kInfeasible;
      return HighsStatus::kOk;
      break;
    }
    case TerminationStatus::UNBOUNDED: {
      assert(111 == 333);
      model_status = HighsModelStatus::kUnbounded;
      return HighsStatus::kOk;
      break;
    }
    case TerminationStatus::TIMEOUT: {
      // ToDo IterationLimit termination needs to be handled separately
      model_status =
          highs_info.pdlp_iteration_count >= options.pdlp_iteration_limit
              ? HighsModelStatus::kIterationLimit
              : HighsModelStatus::kTimeLimit;
      break;
    }
    case TerminationStatus::WARNING:
    case TerminationStatus::FEASIBLE: {
      assert(111 == 555);
      model_status = HighsModelStatus::kUnknown;
      return HighsStatus::kError;
    }
    default:
      assert(termination_status == TerminationStatus::ERROR);
      return HighsStatus::kError;
  }
  assert(termination_status == TerminationStatus::OPTIMAL ||
         termination_status == TerminationStatus::TIMEOUT);
  highs_solution.col_value = x;
  highs_solution.col_value.resize(lp.num_col_);
  highs_solution.row_dual = y;
  lp.a_matrix_.product(highs_solution.row_value, highs_solution.col_value);
  lp.a_matrix_.productTranspose(highs_solution.col_dual,
                                highs_solution.row_dual);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    highs_solution.col_dual[iCol] =
        lp.col_cost_[iCol] - highs_solution.col_dual[iCol];
  highs_solution.value_valid = true;
  highs_solution.dual_valid = true;
  return HighsStatus::kOk;
}



void PrimalDualParams::initialise() {
  this->eta = 0;
  this->omega = 0;
  this->tolerance = 0;
  this->max_iterations = 0;
  this->device_type = Device::CPU;
  this->time_limit = 3600.0;
  this->restart_strategy = RestartStrategy::NO_RESTART;
  this->fixed_restart_interval = 0;
  this->use_halpern_restart = false;
  this->scaling_method = ScalingMethod::NONE;
  this->use_ruiz_scaling = false;
  this->use_pc_scaling = false;
  this->use_l2_scaling = false;
  this->ruiz_iterations = 10;
  this->ruiz_norm = INFINITY;
  this->pc_alpha = 1.0;
  this->step_size_strategy = StepSizeStrategy::FIXED;
  this->malitsky_pock_params.initialise();
  this->adaptive_linesearch_params.initialise();
}

void MalitskyPockParams::initialise() {
  this->step_size_interpolation = 0.5;  // Between 0 and 1
  this->step_size_downscaling_factor = 0.7;
  this->linesearch_contraction_factor = 0.99;
}

void AdaptiveLinesearchParams::initialise() {
  this->step_size_reduction_exponent = 0.3;
  this->step_size_growth_exponent = 0.6;
}
