#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "Highs.h"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/algorithm/solver.hpp"

int main(int argc, char** argv) {
  if (argc != 4) return 2;
  using Clock = std::chrono::steady_clock;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  const auto read_start = Clock::now();
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const double read_seconds = std::chrono::duration<double>(Clock::now() - read_start).count();
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("hipdlp_scaling", "auto");
  highs.setOptionValue("pdlp_step_size_strategy", (HighsInt)3);
  highs.setOptionValue("kkt_tolerance", 1e-4);
  highs.setOptionValue("pdlp_optimality_tolerance", 1e-4);
  highs.setOptionValue("time_limit", std::atof(argv[3]));
  const auto translated = hipdlp::api::translateOptions(highs.getOptions(),
      hipdlp::api::MigrationOptions(hipdlp::api::SolverImplementation::kRewrite,
                                  hipdlp::backend::Kind::kCuda));
  if (!translated.ok()) return 4;
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const auto start = Clock::now();
  const auto result = hipdlp::algorithm::solve(translated.options,
      hipdlp::model::OriginalModel(highs.getLp()), context);
  const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
  std::printf("BENCH_TIMING {\"read_wall_sec\":%.17g,\"solver_wall_sec\":%.17g}\n",
              read_seconds, seconds);
  const char* status = result.status == hipdlp::api::SolverStatus::kOptimal ? "OPTIMAL" :
      result.status == hipdlp::api::SolverStatus::kTimeLimit ? "TIME_LIMIT" :
      result.status == hipdlp::api::SolverStatus::kIterationLimit ? "ITERATION_LIMIT" : "OTHER";
  if (!result.ok() || !result.has_solution) {
    std::fprintf(stderr, "%s\n", result.error.message.c_str());
    return 5;
  }
  const auto dump = [&](const char* suffix, const std::vector<double>& values) {
    FILE* file = std::fopen((std::string(argv[2])+suffix).c_str(), "w");
    if (!file) return false;
    for (double v : values) std::fprintf(file, "%.17g\n", v);
    return std::fclose(file) == 0;
  };
  if (!dump(".primal.txt", result.solution.col_value) ||
      !dump(".dual.txt", result.solution.row_dual) ||
      !dump(".reduced_cost.txt", result.solution.col_dual)) return 6;
  FILE* file = std::fopen((std::string(argv[2])+".json").c_str(), "w");
  if (!file) return 7;
  std::fprintf(file,
      "{\"status\":\"%s\",\"iterations\":%lld,\"time\":%.17g,"
      "\"preprocessing_seconds\":%.17g,\"host_kkt_checks\":%zu,\"kkt_screens\":%zu}\n",
      status, (long long)result.iteration_count, seconds, diagnostics.preprocessing_seconds,
      diagnostics.original_kkt_verifications, diagnostics.kkt_screens);
  return std::fclose(file) == 0 ? 0 : 8;
}
