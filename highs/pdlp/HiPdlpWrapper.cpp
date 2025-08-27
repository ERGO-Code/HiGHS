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

void getHiPpdlpParamsFromOptions(const HighsOptions& options, HighsTimer& timer, PrimalDualParams& params);
			      
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
  LogLevel log_level = LogLevel::kInfo;

  // --- Initialize Logger ---
  Logger logger(log_level);
  if (!log_filename.empty()) logger.set_log_file(log_filename);
  logger.print_header();
    
  // --- Initialize Parameters with Defaults ---
  PrimalDualParams params{};

  getHiPpdlpParamsFromOptions(options, timer, params);
  std::vector<double> x, y;
  PDLPSolver pdlp(logger);
    
  Timer total_timer;
  pdlp.Solve((HighsLp&)lp, params, x, y);
    
  // --- Print Summary ---
  logger.print_summary(pdlp.GetResults(), pdlp.GetIterationCount(), total_timer.read());
  return HighsStatus::kError;
}

void getHiPpdlpParamsFromOptions(const HighsOptions& options, HighsTimer& timer, PrimalDualParams& params) {
  params.initialise();
  //  params.eta = 0; Not set in parse_options_file
  //  params.omega = 0; Not set in parse_opions_file
  //
  // Is params.tolerance for primal and dual feasibility, and
  // optimality?
  params.tolerance = options.pdlp_optimality_tolerance;
  if (options.kkt_tolerance != kDefaultKktTolerance) {
    params.tolerance = options.kkt_tolerance;
  }
  params.max_iterations = options.pdlp_iteration_limit;
  params.device_type = Device::CPU;
  // HiPDLP has its own timer, so set its time limit according to
  // the time remaining with respect to the HiGHS time limit (if
  // finite)
  double time_limit = options.time_limit;
  if (time_limit < kHighsInf) {
    time_limit -= timer.read();
    time_limit = std::max(0.0, time_limit);
  }
  params.time_limit = time_limit;

  params.scaling_method = ScalingMethod::NONE;
  params.use_ruiz_scaling = false;
  params.use_pc_scaling = false;
  params.use_l2_scaling = false;
  if ((options.pdlp_features_off & kPdlpScalingOff) == 0) {
    // Use scaling: now see which
    params.use_ruiz_scaling = options.pdlp_scaling_mode & kPdlpScalingRuiz;
    params.use_pc_scaling = options.pdlp_scaling_mode & kPdlpScalingPC;
    params.use_l2_scaling = options.pdlp_scaling_mode & kPdlpScalingL2;
  }
  params.ruiz_iterations = options.pdlp_ruiz_iterations;
  //  params.ruiz_norm = INFINITY; Not set in parse_opions_file
  //  params.pc_alpha = 1.0; Not set in parse_opions_file

  // Restart strategy maps 0/1/2 to RestartStrategy
  params.restart_strategy = RestartStrategy::NO_RESTART;
  if ((options.pdlp_features_off & kPdlpRestartOff) == 0) {
    // Use restart: now see which
    if (options.pdlp_restart_strategy == kPdlpRestartStrategyFixed) {
      params.restart_strategy = RestartStrategy::FIXED_RESTART;
    } else if (options.pdlp_restart_strategy == kPdlpRestartStrategyAdaptive) {
      params.restart_strategy = RestartStrategy::ADAPTIVE_RESTART;
    }
  }
  //  params.fixed_restart_interval = 0; Not set in parse_opions_file
  //  params.use_halpern_restart = false; Not set in parse_opions_file
  
  params.step_size_strategy = StepSizeStrategy::FIXED;
  if ((options.pdlp_features_off & kPdlpAdaptiveStepSizeOff) == 0) {
    // Use adaptive step size: now see which
    if (options.pdlp_step_size_strategy == kPdlpStepSizeStrategyAdaptive) {
      params.step_size_strategy = StepSizeStrategy::ADAPTIVE;
    } else if (options.pdlp_step_size_strategy == kPdlpStepSizeStrategyMalitskyPock) {
      params.step_size_strategy = StepSizeStrategy::MALITSKY_POCK;
    }
  }
  //  params.malitsky_pock_params.initialise(); Not set in parse_opions_file
  //  params.adaptive_linesearch_params.initialise(); Not set in parse_opions_file

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

  
