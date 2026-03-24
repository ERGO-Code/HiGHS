/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/logger.cc
 * @brief
 */
#include "logger.hpp"

#include <iostream>
#include <map>
#include <sstream>

// Helper to convert enums to strings for pretty printing
template <typename T>
std::string EnumToString(T e, const std::map<T, std::string>& map) {
  auto it = map.find(e);
  return it != map.end() ? it->second : "UNKNOWN";
}

// Logger implementation
void Logger::SetLevel(const HighsInt log_dev_level) {
  if (log_dev_level == kHighsLogDevLevelVerbose) {
    console_level_ = LogLevel::kDebug;
  } else if (log_dev_level == kHighsLogDevLevelDetailed) {
    console_level_ = LogLevel::kVerbose;
  } else if (log_dev_level == kHighsLogDevLevelInfo) {
    console_level_ = LogLevel::kDetailed;
  } else {
    console_level_ = LogLevel::kInfo;  // None;
  }
}

void Logger::Log(LogLevel level, const std::string& message) const {
  // Now using HiGHS IO
  if (level <= console_level_)
    highsLogUser(log_options_, HighsLogType::kInfo, "%s\n", message.c_str());
}

void Logger::Info(const std::string& message) const {
  Log(LogLevel::kInfo, message);
}
void Logger::Detailed(const std::string& message) const {
  Log(LogLevel::kDetailed, message);
}
void Logger::Verbose(const std::string& message) const {
  Log(LogLevel::kVerbose, message);
}
void Logger::Debug(const std::string& message) const {
  Log(LogLevel::kDebug, message);
}

void Logger::PrintHeader() const {
  Detailed("Using HiPDLP with developer logging");
}

void Logger::PrintParams(const PrimalDualParams& params) const {
  Detailed("\nSolver Parameters:");
  std::stringstream ss;

  std::map<RestartStrategy, std::string> restart_map = {
      {RestartStrategy::NO_RESTART, "None"},
      {RestartStrategy::FIXED_RESTART, "Fixed"},
      {RestartStrategy::ADAPTIVE_RESTART, "Adaptive"}};
  std::map<StepSizeStrategy, std::string> step_size_map = {
      {StepSizeStrategy::FIXED, "Fixed"},
      {StepSizeStrategy::ADAPTIVE, "Adaptive"},
      {StepSizeStrategy::MALITSKY_POCK, "Malitsky-Pock"},
      {StepSizeStrategy::PID, "PID"}};

  ss << "  - Max Iterations: " << params.max_iterations;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Tolerance: " << params.tolerance;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Restart Strategy: "
     << EnumToString(params.restart_strategy, restart_map);
  Detailed(ss.str());
  ss.str("");

  std::string scaling_method = "None";
  if (params.use_ruiz_scaling && params.use_pc_scaling &&
      params.use_l2_scaling) {
    scaling_method = "Combined (Ruiz + PC + L2)";
  } else if (params.use_ruiz_scaling && params.use_pc_scaling) {
    scaling_method = "Ruiz + Pock-Chambolle";
  } else if (params.use_ruiz_scaling && params.use_l2_scaling) {
    scaling_method = "Ruiz + L2-Norm";
  } else if (params.use_pc_scaling && params.use_l2_scaling) {
    scaling_method = "Pock-Chambolle + L2-Norm";
  } else if (params.use_ruiz_scaling) {
    scaling_method = "Ruiz";
  } else if (params.use_pc_scaling) {
    scaling_method = "Pock-Chambolle";
  } else if (params.use_l2_scaling) {
    scaling_method = "L2-Norm";
  }

  ss << "  - Scaling Method: " << scaling_method;
  Detailed(ss.str());
  ss.str("");

  // Optionally, add more details about scaling parameters if any scaling is
  // used
  if (params.use_ruiz_scaling) {
    ss.str("");
    ss << "   * Ruiz iterations: " << params.ruiz_iterations;
    Detailed(ss.str());
  }
  if (params.use_pc_scaling) {
    ss.str("");
    ss << "   * PC alpha: " << params.pc_alpha;
    Detailed(ss.str());
  }
  ss << "  - Step Size Strategy: "
     << EnumToString(params.step_size_strategy, step_size_map);
  Detailed(ss.str());
  Detailed("------------------------------------------------------------");
}

void Logger::PrintIterationHeader() const {
  Info(
      "     Iter    Primal Feas   Dual Feas     P-D Gap   Step Size       "
      "Time");
}

void Logger::PrintIterationStats(HighsInt iter, const SolverResults& results,
                                 const double step_size,
                                 const double time) const {
  std::stringstream ss;
  // clang-format off
  ss << std::fixed << std::setprecision(9) << std::setw(9) << iter << "    "
     << std::scientific << std::setprecision(2)
     << std::setw(11) << results.primal_feasibility << " "
     << std::setw(11) << results.dual_feasibility << " "
     << std::setw(11) << results.duality_gap << "   "
     << std::fixed << std::setprecision(4)
     << std::setw(9) << step_size << "   "
     << std::setprecision(1)
     << std::setw(8) << time;
  // clang-format on
  Info(ss.str());
}

void Logger::PrintSummary(const SolverResults& results, HighsInt total_iter,
                          double total_time) const {
  Detailed("\n-------------------- Solver Summary --------------------");
  std::stringstream ss;

  std::map<TerminationStatus, std::string> term_map = {
      {TerminationStatus::NOTSET, "Not set"},
      {TerminationStatus::OPTIMAL, "Optimal"},
      {TerminationStatus::INFEASIBLE, "Infeasible"},
      {TerminationStatus::UNBOUNDED, "Unbounded"},
      {TerminationStatus::MAXITER, "Maxiter"},
      {TerminationStatus::TIMEOUT, "Timeout"},
      {TerminationStatus::ERROR, "Error"}};

  ss << "  - Termination Status: " << EnumToString(results.term_code, term_map);
  Detailed(ss.str());
  ss.str("");
  ss << "  - Total Iterations: " << total_iter;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Total Time: " << std::fixed << std::setprecision(3) << total_time
     << " seconds";
  Detailed(ss.str());
  ss.str("");
  ss << "  - Final Primal Objective: " << std::scientific
     << std::setprecision(6) << results.primal_obj;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Final Duality Gap: " << std::scientific << std::setprecision(6)
     << results.duality_gap;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Final Primal Feasibility: " << std::scientific
     << std::setprecision(6) << results.primal_feasibility;
  Detailed(ss.str());
  ss.str("");
  ss << "  - Final Dual Feasibility: " << std::scientific
     << std::setprecision(6) << results.dual_feasibility;
  Detailed(ss.str());
  Detailed("------------------------------------------------------------");
}
