/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/logger.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_LOGGER_HPP
#define PDLP_HIPDLP_LOGGER_HPP

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "restart.hpp"  // For PrimalDualParams
#include "solver_results.hpp"

// Log verbosity level
enum class LogLevel {
  kNone,      // No output
  kInfo,      // Standard output: summary, termination, major events
  kDetailed,  // Developer output
  kVerbose,   // Detailed output: iteration-level info
  kDebug      // Verbose + debug info for developers
};

class Logger {
 public:
  void SetLevel(const HighsInt log_dev_level);
  void SetHighsLogOptions(const HighsLogOptions log_options) {
    log_options_ = log_options;
  }
  LogLevel GetLogLevel() const { return console_level_; }
  // Logging methods for different levels
  void Info(const std::string& message) const;
  void Detailed(const std::string& message) const;
  void Verbose(const std::string& message) const;
  void Debug(const std::string& message) const;

  // Formatted printing functions
  void PrintHeader() const;
  void PrintParams(const PrimalDualParams& params) const;
  void PrintIterationHeader() const;
  void PrintIterationStats(HighsInt iter, const SolverResults& current_results,
                           const double current_eta, const double time) const;
  void PrintSummary(const SolverResults& results, HighsInt total_iter,
                    double total_time) const;
  void SetLogOptions(HighsLogOptions log_options) {
    log_options_ = log_options;
  }
  LogLevel GetConsoleLevel() const { return console_level_; }

 private:
  void Log(LogLevel level, const std::string& message) const;
  LogLevel console_level_;
  HighsLogOptions log_options_;
};

#endif  // PDLP_HIPDLP_LOGGER_HPP
