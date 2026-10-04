// Controlled legacy/rewrite GPU comparison with real HiGHS presolve/postsolve.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "Highs.h"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/algorithm/solver.hpp"

int main(int argc, char** argv) {
  if (argc != 5) return 2;  // model, legacy|rewrite, export prefix, limit
#ifndef CUPDLP_GPU
  std::fprintf(stderr, "Requires a legacy GPU build\n");
  return 2;
#endif
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
  };
  const std::string implementation = argv[2], prefix = argv[3];
  if (implementation != "legacy" && implementation != "rewrite") return 2;
  Highs highs;
  highs.setOptionValue("output_flag", true);
  highs.setOptionValue("threads", (HighsInt)1);
  const auto read_start = Clock::now();
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const double read_seconds = elapsed(read_start);
  if (highs.getModel().isMip() || highs.getModel().isQp()) {
    std::fprintf(stderr, "This benchmark requires continuous LP input\n");
    return 3;
  }
  highs.setOptionValue("solver", "hipdlp");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("hipdlp_backend", "cuda");
  highs.setOptionValue("hipdlp_scaling", "legacy");  // Ruiz+PC, no CR
  highs.setOptionValue("hipdlp_scaling_selector", false);
  highs.setOptionValue("pdlp_scaling_mode", (HighsInt)5);
  highs.setOptionValue("pdlp_ruiz_iterations", (HighsInt)10);
  highs.setOptionValue("pdlp_step_size_strategy", (HighsInt)3);
  highs.setOptionValue("kkt_tolerance", 1e-4);
  highs.setOptionValue("pdlp_optimality_tolerance", 1e-4);
  highs.setOptionValue("time_limit", std::atof(argv[4]));
  const auto start = Clock::now();
  const auto ps = highs.presolve();
  const double presolve_seconds = elapsed(start);
  const auto presolve_status = highs.getModelPresolveStatus();
  if (ps == HighsStatus::kError ||
      (presolve_status != HighsPresolveStatus::kReduced &&
       presolve_status != HighsPresolveStatus::kNotReduced &&
       presolve_status != HighsPresolveStatus::kReducedToEmpty)) return 4;
  const auto& lp = highs.getPresolvedLp();
  HighsSolution solution;
  HighsInfo info;
  HighsBasis basis;
  HighsTimer timer;
  HighsCallback callback(&highs);
  HighsModelStatus native_status = HighsModelStatus::kNotset;
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  HighsStatus status = HighsStatus::kOk;
  const auto algorithm_start = Clock::now();
  timer.start();  // Legacy reads this clock to enforce its native time limit.
  if (presolve_status == HighsPresolveStatus::kReducedToEmpty) {
    solution.value_valid = solution.dual_valid = true;
    native_status = HighsModelStatus::kOptimal;
    info.pdlp_iteration_count = 0;
  } else if (implementation == "legacy") {
    status = hipdlp::api::solve(highs.getOptions(), timer, lp, basis,
        solution, native_status, info, callback,
        hipdlp::api::MigrationOptions(hipdlp::api::SolverImplementation::kLegacy,
                                     hipdlp::backend::Kind::kAuto));
  } else {
    const auto translated = hipdlp::api::translateOptions(highs.getOptions(),
        hipdlp::api::MigrationOptions(hipdlp::api::SolverImplementation::kRewrite,
                                     hipdlp::backend::Kind::kCuda));
    if (!translated.ok()) return 5;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &diagnostics;
    const auto result = hipdlp::algorithm::solve(translated.options,
        hipdlp::model::OriginalModel(lp), context);
    if (!result.ok() || !result.has_solution) {
      std::fprintf(stderr, "Solve failed: %s\n", result.error.message.c_str());
      return 6;
    }
    solution = result.solution;
    info.pdlp_iteration_count = result.iteration_count;
    native_status = result.status == hipdlp::api::SolverStatus::kOptimal ? HighsModelStatus::kOptimal :
        result.status == hipdlp::api::SolverStatus::kTimeLimit ? HighsModelStatus::kTimeLimit :
        result.status == hipdlp::api::SolverStatus::kIterationLimit ? HighsModelStatus::kIterationLimit : HighsModelStatus::kUnknown;
  }
  const double algorithm_seconds = elapsed(algorithm_start);
  timer.stop();
  if (status == HighsStatus::kError) return 6;
  const auto postsolve_start = Clock::now();
  // No basis is passed: postsolve only recovers vectors, without simplex cleanup.
  if (presolve_status != HighsPresolveStatus::kNotReduced) {
    if (highs.postsolve(solution) == HighsStatus::kError) return 7;
    solution = highs.getSolution();
  }
  const double postsolve_seconds = elapsed(postsolve_start);
  const double seconds = elapsed(start);
  std::printf("BENCH_TIMING {\"read_wall_sec\":%.17g,\"solver_wall_sec\":%.17g}\n", read_seconds, seconds);
  const auto dump = [&](const char* suffix, const std::vector<double>& values) {
    FILE* file = std::fopen((prefix+suffix).c_str(), "w");
    if (!file) return false;
    for (double v : values) std::fprintf(file, "%.17g\n", v);
    return std::fclose(file) == 0;
  };
  if (!dump(".primal.txt", solution.col_value) || !dump(".dual.txt", solution.row_dual) ||
      !dump(".reduced_cost.txt", solution.col_dual)) return 8;
  FILE* file = std::fopen((prefix+".json").c_str(), "w");
  if (!file) return 8;
  const char* status_name = native_status == HighsModelStatus::kOptimal ? "OPTIMAL" :
      native_status == HighsModelStatus::kTimeLimit ? "TIME_LIMIT" :
      native_status == HighsModelStatus::kIterationLimit ? "ITERATION_LIMIT" : "OTHER";
  std::fprintf(file, "{\"status\":\"%s\",\"iterations\":%lld,\"time\":%.17g,"
      "\"presolve_seconds\":%.17g,\"algorithm_seconds\":%.17g,\"postsolve_seconds\":%.17g,"
      "\"presolve_status\":\"%s\",\"reduced_rows\":%lld,\"reduced_cols\":%lld,\"reduced_nnz\":%zu,"
      "\"rewrite_preprocessing_seconds\":%.17g,\"backend_setup_seconds\":%.17g,"
      "\"operator_norm_seconds\":%.17g,\"iterations_seconds\":%.17g,"
      "\"host_kkt_checks\":%zu,\"kkt_screens\":%zu}\n",
      status_name, (long long)info.pdlp_iteration_count, seconds, presolve_seconds,
      algorithm_seconds, postsolve_seconds, highs.presolveStatusToString(presolve_status).c_str(),
      (long long)lp.num_row_, (long long)lp.num_col_, lp.a_matrix_.value_.size(),
      diagnostics.preprocessing_seconds, diagnostics.backend_setup_seconds,
      diagnostics.operator_norm_seconds, diagnostics.iterations_seconds,
      diagnostics.original_kkt_verifications, diagnostics.kkt_screens);
  return std::fclose(file) == 0 ? 0 : 9;
}
