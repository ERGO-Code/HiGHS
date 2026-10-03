// Time several fresh rewrite solves of one scratch LP without rereading it.
// Input I/O is shared; solver_wall_sec includes preprocessing, pilots and KKT.
// Usage: large_compare SCRATCH.mps NAME cpu|cuda CONFIGS SECONDS ITERATIONS
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "lp_data/HighsCallback.h"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
std::string quoted(const std::string& value) {
  std::string result = "\"";
  for (char c : value) {
    if (c == '\\' || c == '"') result += '\\';
    if (c == '\n')
      result += "\\n";
    else
      result += c;
  }
  return result + "\"";
}
const char* statusName(HighsModelStatus status) {
  switch (status) {
    case HighsModelStatus::kOptimal:
      return "OPTIMAL";
    case HighsModelStatus::kInfeasible:
      return "INFEASIBLE";
    case HighsModelStatus::kUnbounded:
      return "UNBOUNDED";
    case HighsModelStatus::kUnboundedOrInfeasible:
      return "UNBOUNDED_OR_INFEASIBLE";
    case HighsModelStatus::kIterationLimit:
      return "ITERATION_LIMIT";
    case HighsModelStatus::kTimeLimit:
      return "TIME_LIMIT";
    default:
      return "ERROR";
  }
}
void number(const char* name, double value) {
  if (std::isfinite(value))
    std::printf(",\"%s\":%.17g", name, value);
  else
    std::printf(",\"%s\":null", name);
}
bool configure(Highs& highs, const std::string& config,
               const std::string& backend, double seconds,
               HighsInt iterations) {
  const bool valid = config == "A" || config == "L" || config == "S" ||
                     config == "P" || config == "R" || config == "C" ||
                     config == "O" || config == "AP" || config == "AR";
  if (!valid) return false;
  std::string scaling = "auto";
  if (config == "L" || config == "P" || config == "R") scaling = "legacy";
  if (config == "C") scaling = "curtis_reid";
  if (config == "O") scaling = "off";
  const HighsInt mask = (config == "P" || config == "AP")   ? 4
                        : (config == "R" || config == "AR") ? 1
                                                            : 5;
  return highs.setOptionValue("solver", "hipdlp") == HighsStatus::kOk &&
         highs.setOptionValue("presolve", "off") == HighsStatus::kOk &&
         highs.setOptionValue("hipdlp_backend", backend) == HighsStatus::kOk &&
         highs.setOptionValue("hipdlp_scaling", scaling) == HighsStatus::kOk &&
         highs.setOptionValue("hipdlp_scaling_selector", config == "S") ==
             HighsStatus::kOk &&
         highs.setOptionValue("pdlp_scaling_mode", mask) == HighsStatus::kOk &&
         highs.setOptionValue("pdlp_step_size_strategy",
                              kPdlpStepSizeStrategyFixed) == HighsStatus::kOk &&
         highs.setOptionValue("kkt_tolerance", 1e-4) == HighsStatus::kOk &&
         highs.setOptionValue("pdlp_optimality_tolerance", 1e-4) ==
             HighsStatus::kOk &&
         highs.setOptionValue("time_limit", seconds) == HighsStatus::kOk &&
         highs.setOptionValue("pdlp_iteration_limit", iterations) ==
             HighsStatus::kOk;
}
hipdlp::api::SolveResult run(Highs& highs, hipdlp::backend::Kind backend,
                             HighsModelStatus& model_status, HighsInfo& info,
                             HighsStatus& call_status) {
  HighsBasis basis;
  HighsSolution solution;
  HighsCallback callback(&highs);
  HighsTimer timer;
  hipdlp::api::SolveResult result;
  call_status = hipdlp::api::solve(
      highs.getOptions(), timer, highs.getLp(), basis, solution, model_status,
      info, callback,
      hipdlp::api::MigrationOptions(hipdlp::api::SolverImplementation::kRewrite,
                                    backend),
      &result);
  return result;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 7) return 2;
  hipdlp::backend::Kind backend;
  if (!hipdlp::api::backendKindFromString(argv[3], backend)) return 2;
  const double seconds = std::atof(argv[5]);
  const HighsInt iterations = static_cast<HighsInt>(std::atoll(argv[6]));
  Highs highs;
  highs.setOptionValue("output_flag", false);
  const auto read_start = Clock::now();
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const double read_seconds = elapsed(read_start);

  // Warm the selected backend on a separate one-row LP, without changing the
  // benchmark LP or carrying any iterate into a timed solve.
  Highs warm;
  warm.setOptionValue("output_flag", false);
  HighsLp tiny;
  tiny.num_row_ = tiny.num_col_ = 1;
  tiny.col_cost_ = {1};
  tiny.col_lower_ = {0};
  tiny.col_upper_ = {kHighsInf};
  tiny.row_lower_ = tiny.row_upper_ = {1};
  tiny.a_matrix_.format_ = MatrixFormat::kColwise;
  tiny.a_matrix_.num_row_ = tiny.a_matrix_.num_col_ = 1;
  tiny.a_matrix_.start_ = {0, 1};
  tiny.a_matrix_.index_ = {0};
  tiny.a_matrix_.value_ = {1};
  if (warm.passModel(tiny) == HighsStatus::kError ||
      !configure(warm, "L", argv[3], 30, 40))
    return 4;
  HighsModelStatus warm_status = HighsModelStatus::kNotset;
  HighsInfo warm_info;
  HighsStatus warm_call;
  if (!run(warm, backend, warm_status, warm_info, warm_call).ok()) return 4;

  std::printf(
      "LARGE_MODEL {\"name\":%s,\"rows\":%lld,\"cols\":%lld,\"nnz\":%lld",
      quoted(argv[2]).c_str(), (long long)highs.getLp().num_row_,
      (long long)highs.getLp().num_col_,
      (long long)highs.getLp().a_matrix_.numNz());
  number("read_wall_sec", read_seconds);
  std::printf("}\n");
  std::fflush(stdout);
  highs.setOptionValue("output_flag", true);
  std::stringstream configs(argv[4]);
  std::string config;
  int order = 0;
  while (std::getline(configs, config, ',')) {
    if (!configure(highs, config, argv[3], seconds, iterations)) return 2;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsStatus call_status;
    std::printf("LARGE_BEGIN %s\n", config.c_str());
    std::fflush(stdout);
    const auto start = Clock::now();
    const auto result = run(highs, backend, model_status, info, call_status);
    const double solver_seconds = elapsed(start);
    const auto verify_start = Clock::now();
    const auto kkt =
        hipdlp::verification::verifyKkt(highs.getLp(), result.solution, 1e-4);
    bool certificate_pass = false;
    if (result.certificate.kind ==
        hipdlp::api::CertificateKind::kPrimalUnboundedness)
      certificate_pass =
          hipdlp::verification::verifyPrimalRay(
              highs.getLp(), result.certificate.column_direction, 1e-4)
              .pass;
    else if (result.certificate.kind ==
                 hipdlp::api::CertificateKind::kPrimalInfeasibility ||
             result.certificate.kind ==
                 hipdlp::api::CertificateKind::kAmbiguous)
      certificate_pass = hipdlp::verification::verifyFarkasCertificate(
                             highs.getLp(), result.certificate.row_direction,
                             result.certificate.column_direction, 1e-4)
                             .pass;
    std::printf(
        "LARGE_RESULT {\"name\":%s,\"config\":%s,\"backend\":%s,\"status\":%s,"
        "\"call_status\":%d,\"iterations\":%lld,\"order\":%d,\"kkt_pass\":%s,"
        "\"certificate_pass\":%s,\"kkt_tolerance\":0.0001,\"pdlp_optimality_"
        "tolerance\":0.0001",
        quoted(argv[2]).c_str(), quoted(config).c_str(),
        quoted(argv[3]).c_str(), quoted(statusName(model_status)).c_str(),
        (int)call_status, (long long)result.iteration_count, order++,
        kkt.pass ? "true" : "false", certificate_pass ? "true" : "false");
    number("solver_wall_sec", solver_seconds);
    number("verification_wall_sec", elapsed(verify_start));
    number("time_limit", seconds);
    std::printf(",\"iteration_limit\":%lld", (long long)iterations);
    number("relative_primal_residual", kkt.metrics.relative_primal_residual);
    number("relative_dual_residual", kkt.metrics.relative_dual_residual);
    number("relative_objective_gap", kkt.metrics.relative_objective_gap);
    std::printf("}\n");
    std::fflush(stdout);
  }
}
