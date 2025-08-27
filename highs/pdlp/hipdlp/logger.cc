#include "logger.hpp"
#include <iostream>
#include <sstream>
#include <map>

// Helper to convert enums to strings for pretty printing
template<typename T>
std::string enum_to_string(T e, const std::map<T, std::string>& map) {
    auto it = map.find(e);
    return it != map.end() ? it->second : "UNKNOWN";
}

// Timer implementation
Timer::Timer() { reset(); }
void Timer::reset() { start_time_ = std::chrono::high_resolution_clock::now(); }
double Timer::read() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time_;
    return elapsed.count();
}

// Logger implementation
Logger::Logger(LogLevel level) : console_level_(level) {}

void Logger::set_log_file(const std::string& filename) {
    log_file_.open(filename, std::ios::out | std::ios::trunc);
    if (!log_file_.is_open()) {
        std::cerr << "Error: Could not open log file: " << filename << std::endl;
    }
}

void Logger::log(LogLevel level, const std::string& message) {
    if (level <= console_level_) {
        std::cout << message << std::endl;
    }
    if (log_file_.is_open()) {
        log_file_ << message << std::endl;
    }
}

void Logger::info(const std::string& message) { log(LogLevel::kInfo, message); }
void Logger::verbose(const std::string& message) { log(LogLevel::kVerbose, message); }
void Logger::debug(const std::string& message) { log(LogLevel::kDebug, message); }

void Logger::print_header() {
    info("------------------------------------------------------------");
    info("          PDLP Solver - C++ Implementation          ");
    info("------------------------------------------------------------");
}

void Logger::print_params(const PrimalDualParams& params) {
    info("\nSolver Parameters:");
    std::stringstream ss;
    
    std::map<RestartStrategy, std::string> restart_map = {{RestartStrategy::NO_RESTART, "None"}, {RestartStrategy::FIXED_RESTART, "Fixed"}, {RestartStrategy::ADAPTIVE_RESTART, "Adaptive"}};
    std::map<ScalingMethod, std::string> scaling_map = {{ScalingMethod::NONE, "None"}, {ScalingMethod::RUIZ, "Ruiz"}, {ScalingMethod::POCK_CHAMBOLLE, "Pock-Chambolle"}, {ScalingMethod::L2_NORM, "L2-Norm"}};
    std::map<StepSizeStrategy, std::string> step_size_map = {{StepSizeStrategy::FIXED, "Fixed"}, {StepSizeStrategy::ADAPTIVE, "Adaptive"}, {StepSizeStrategy::MALITSKY_POCK, "Malitsky-Pock"}};

    ss << "  - Max Iterations: " << params.max_iterations;
    info(ss.str()); ss.str("");
    ss << "  - Tolerance: " << params.tolerance;
    info(ss.str()); ss.str("");
    ss << "  - Restart Strategy: " << enum_to_string(params.restart_strategy, restart_map);
    info(ss.str()); ss.str("");
    ss << "  - Scaling Method: " << enum_to_string(params.scaling_method, scaling_map);
    info(ss.str()); ss.str("");
    ss << "  - Step Size Strategy: " << enum_to_string(params.step_size_strategy, step_size_map);
    info(ss.str());
    info("------------------------------------------------------------");
}

void Logger::print_iteration_header() {
    verbose("\n-------------------------------------------------------------------------------------------------");
    verbose(" Iter   | Primal Feas | Dual Feas  | Duality Gap  | Step Size");
    verbose("-------------------------------------------------------------------------------------------------");
}

void Logger::print_iteration_stats(int iter, const SolverResults& results, double step_size) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4)
       << " " << std::setw(6) << iter << " | "
       << std::scientific << std::setprecision(2)
       << std::setw(11) << results.primal_feasibility << " | "
       << std::setw(10) << results.dual_feasibility << " | "
       << std::setw(11) << results.duality_gap << " | "
       << std::fixed << std::setprecision(4)
       << std::setw(9) << step_size;
    verbose(ss.str());
}

void Logger::print_summary(const SolverResults& results, int total_iter, double total_time) {
    info("\n-------------------- Solver Summary --------------------");
    std::stringstream ss;
    
    std::map<TerminationStatus, std::string> term_map = {{TerminationStatus::OPTIMAL, "Optimal"}, {TerminationStatus::TIMEOUT, "Timeout"}, {TerminationStatus::FEASIBLE, "Feasible"}};
    
    ss << "  - Termination Status: " << enum_to_string(results.term_code, term_map);
    info(ss.str()); ss.str("");
    ss << "  - Total Iterations: " << total_iter;
    info(ss.str()); ss.str("");
    ss << "  - Total Time: " << std::fixed << std::setprecision(3) << total_time << " seconds";
    info(ss.str()); ss.str("");
    ss << "  - Final Primal Objective: " << std::scientific << std::setprecision(6) << results.primal_obj;
    info(ss.str()); ss.str("");
    ss << "  - Final Duality Gap: " << std::scientific << std::setprecision(6) << results.duality_gap;
    info(ss.str()); ss.str("");
    ss << "  - Final Primal Feasibility: " << std::scientific << std::setprecision(6) << results.primal_feasibility;
    info(ss.str()); ss.str("");
    ss << "  - Final Dual Feasibility: " << std::scientific << std::setprecision(6) << results.dual_feasibility;
    info(ss.str());
    info("------------------------------------------------------------");
}
