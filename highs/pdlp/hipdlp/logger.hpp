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
  void setLevel(const HighsInt log_dev_level);
  void setHighsLogOptions(const HighsLogOptions log_options) {
    log_options_ = log_options;
  }
  LogLevel getLogLevel() const { return console_level_; }
  // Logging methods for different levels
  void info(const std::string& message) const;
  void detailed(const std::string& message) const;
  void verbose(const std::string& message) const;
  void debug(const std::string& message) const;

  // Formatted printing functions
  void printHeader() const;
  void printParams(const PrimalDualParams& params) const;
  void printIterationHeader() const;
  void printIterationStats(HighsInt iter, const SolverResults& current_results,
                           const double current_eta, const double time) const;
  void printSummary(const SolverResults& results, HighsInt total_iter,
                    double total_time) const;
  void setLogOptions(HighsLogOptions log_options) {
    log_options_ = log_options;
  }
  LogLevel getConsoleLevel() const { return console_level_; }

 private:
  void log(LogLevel level, const std::string& message) const;
  LogLevel console_level_;
  HighsLogOptions log_options_;
};

#endif  // PDLP_HIPDLP_LOGGER_HPP
