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
  kNone,     // No output
  kInfo,     // Standard output: summary, termination, major events
  kVerbose,  // Detailed output: iteration-level info
  kDebug     // Verbose + debug info for developers
};

class Timer {
 public:
  Timer();
  void reset();
  double read();

 private:
  std::chrono::high_resolution_clock::time_point start_time_;
};

class Logger {
 public:
  void setLevel(const HighsInt log_dev_level);
  void passHighsLogOptions(const HighsLogOptions log_options) { log_options_ = log_options;}
  void set_log_file(const std::string& filename);
  LogLevel getLogLevel() const { return console_level_; }
  // Logging methods for different levels
  void info(const std::string& message);
  void verbose(const std::string& message);
  void debug(const std::string& message);

  // Formatted printing functions
  void print_header();
  void print_params(const PrimalDualParams& params);
  void print_iteration_header();
  void print_iteration_stats(int iter, const SolverResults& current_results,
                             double current_eta);
  void print_summary(const SolverResults& results, int total_iter,
                     double total_time);

 private:
  void log(LogLevel level, const std::string& message);
  LogLevel console_level_;
  HighsLogOptions log_options_;
  
};

#endif  // PDLP_HIPDLP_LOGGER_HPP
