/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-08-11 10:52:55
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-11 10:53:02
 * @FilePath: /cupdlp-CPP/include/logger.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <memory>
#include "solver_results.hpp"
#include "restart.hpp" // For PrimalDualParams

// Log verbosity level
enum class LogLevel {
    kNone,    // No output
    kInfo,    // Standard output: summary, termination, major events
    kVerbose, // Detailed output: iteration-level info
    kDebug    // Verbose + debug info for developers
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
    Logger(LogLevel level = LogLevel::kInfo);
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
    void print_iteration_stats(int iter, const SolverResults& current_results, double current_eta);
    void print_summary(const SolverResults& results, int total_iter, double total_time);

private:
    void log(LogLevel level, const std::string& message);
    LogLevel console_level_;
    std::ofstream log_file_;
};

#endif // LOGGER_HPP
