// Migration/benchmark adapter for the HiPDLP rewrite.
//
// This executable is deliberately NOT installed and NOT a public interface. It
// exists so that the legacy and rewrite HiPDLP implementations can be compared
// through exactly the same option translation, the same solve entry point and
// the same timer boundary, without adding any public implementation option.
//
// The solver timer starts after the model has been read and stops when the
// solve call returns; end_to_end_sec additionally covers reading and reporting.
// This calls the solver directly, without HiGHS presolve or postsolve. Use the
// highs executable for end-to-end runs with presolve.

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "Highs.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsCallback.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {

struct AdapterOptions {
  std::string instance;
  std::string implementation = "rewrite";
  std::string backend = "auto";
  std::string step_size = "fixed";
  std::string scaling = "auto";
  double time_limit = 3600.0;
  double kkt_tolerance = 1e-4;
  HighsInt iteration_limit = kHighsIInf;
  double verification_tolerance = 1e-4;
  std::string output_dir;
  bool log = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: hipdlp_migration_adapter --instance FILE.mps\n"
               "         [--implementation legacy|rewrite]\n"
               "         [--backend auto|cpu|cuda]\n"
               "         [--step-size fixed|pid]\n"
               "         [--scaling auto|off|curtis_reid|legacy]\n"
               "         [--time-limit SECONDS]\n"
               "         [--kkt-tolerance TOL]\n"
               "         [--iteration-limit N]\n"
               "         [--output-dir DIR]\n"
               "         [--log]\n");
}

bool nextValue(const int argc, char** argv, int& index, std::string& value) {
  if (index + 1 >= argc) return false;
  value = argv[++index];
  return true;
}

bool parseArguments(const int argc, char** argv, AdapterOptions& options) {
  for (int index = 1; index < argc; ++index) {
    const std::string flag = argv[index];
    std::string value;
    if (flag == "--instance") {
      if (!nextValue(argc, argv, index, options.instance)) return false;
    } else if (flag == "--implementation") {
      if (!nextValue(argc, argv, index, options.implementation)) return false;
    } else if (flag == "--backend") {
      if (!nextValue(argc, argv, index, options.backend)) return false;
    } else if (flag == "--step-size") {
      if (!nextValue(argc, argv, index, options.step_size)) return false;
    } else if (flag == "--scaling") {
      if (!nextValue(argc, argv, index, options.scaling)) return false;
    } else if (flag == "--time-limit") {
      if (!nextValue(argc, argv, index, value)) return false;
      options.time_limit = std::atof(value.c_str());
    } else if (flag == "--kkt-tolerance") {
      if (!nextValue(argc, argv, index, value)) return false;
      options.kkt_tolerance = std::atof(value.c_str());
      options.verification_tolerance = options.kkt_tolerance;
    } else if (flag == "--iteration-limit") {
      if (!nextValue(argc, argv, index, value)) return false;
      options.iteration_limit = std::atol(value.c_str());
    } else if (flag == "--output-dir") {
      if (!nextValue(argc, argv, index, options.output_dir)) return false;
    } else if (flag == "--log") {
      options.log = true;
    } else {
      std::fprintf(stderr, "unknown flag %s\n", flag.c_str());
      return false;
    }
  }
  return !options.instance.empty();
}

const char* modelStatusName(const HighsModelStatus status) {
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
    case HighsModelStatus::kSolveError:
      return "SOLVE_ERROR";
    case HighsModelStatus::kNotset:
      return "NOT_SET";
    default:
      return "OTHER";
  }
}

const char* certificateKindName(const hipdlp::api::CertificateKind kind) {
  switch (kind) {
    case hipdlp::api::CertificateKind::kNone:
      return "none";
    case hipdlp::api::CertificateKind::kPrimalInfeasibility:
      return "primal_infeasibility";
    case hipdlp::api::CertificateKind::kPrimalUnboundedness:
      return "primal_unboundedness";
    case hipdlp::api::CertificateKind::kAmbiguous:
      return "ambiguous";
  }
  return "unknown";
}

std::string formatLine(const char* format, ...) {
  char buffer[1024];
  va_list arguments;
  va_start(arguments, format);
  vsnprintf(buffer, sizeof(buffer), format, arguments);
  va_end(arguments);
  std::string line(buffer);
  line += "\n";
  return line;
}

std::string quotedField(const char* key, const std::string& value) {
  return formatLine("  \"%s\": \"%s\",", key, value.c_str());
}

std::string intField(const char* key, const long long value,
                     const bool available = true) {
  if (!available) return formatLine("  \"%s\": null,", key);
  return formatLine("  \"%s\": %lld,", key, value);
}

std::string boolField(const char* key, const bool value) {
  return formatLine("  \"%s\": %s,", key, value ? "true" : "false");
}

std::string doubleField(const char* key, const double value) {
  if (!std::isfinite(value)) return formatLine("  \"%s\": null,", key);
  return formatLine("  \"%s\": %.17g,", key, value);
}

std::string nullableDoubleField(const char* key, const double value,
                                const bool available) {
  if (!available || !std::isfinite(value))
    return formatLine("  \"%s\": null,", key);
  return formatLine("  \"%s\": %.17g,", key, value);
}
}  // namespace

int main(int argc, char** argv) {
  AdapterOptions adapter;
  if (!parseArguments(argc, argv, adapter)) {
    usage();
    return 2;
  }
  hipdlp::api::SolverImplementation implementation =
      hipdlp::api::SolverImplementation::kRewrite;
  if (adapter.implementation == "legacy") {
    implementation = hipdlp::api::SolverImplementation::kLegacy;
  } else if (adapter.implementation != "rewrite") {
    std::fprintf(stderr, "unknown implementation %s\n",
                 adapter.implementation.c_str());
    return 2;
  }
  hipdlp::backend::Kind backend = hipdlp::backend::Kind::kAuto;
  if (!hipdlp::api::backendKindFromString(adapter.backend, backend)) {
    std::fprintf(stderr, "unknown backend %s\n", adapter.backend.c_str());
    return 2;
  }
  HighsInt step_size_strategy = kPdlpStepSizeStrategyFixed;
  if (adapter.step_size == "pid") {
    step_size_strategy = kPdlpStepSizeStrategyPid;
  } else if (adapter.step_size != "fixed") {
    std::fprintf(stderr, "unknown step size %s\n", adapter.step_size.c_str());
    return 2;
  }
  // Legacy has no backend selector of its own, so it only ever runs as auto.
  const bool is_legacy =
      implementation == hipdlp::api::SolverImplementation::kLegacy;
  const std::string effective_backend = is_legacy ? "auto" : adapter.backend;
  const hipdlp::backend::Kind migration_backend =
      is_legacy ? hipdlp::backend::Kind::kAuto : backend;

  const std::chrono::steady_clock::time_point end_to_end_start =
      std::chrono::steady_clock::now();

  Highs highs;
  highs.setOptionValue("output_flag", adapter.log);
  // A warning (for example dropped tiny matrix entries) still yields a model.
  if (highs.readModel(adapter.instance) == HighsStatus::kError) {
    std::fprintf(stderr, "failed to read %s\n", adapter.instance.c_str());
    return 3;
  }
  highs.setOptionValue("solver", kHiPdlpString);
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("hipdlp_backend", effective_backend);
  highs.setOptionValue("hipdlp_scaling", adapter.scaling);
  highs.setOptionValue("pdlp_step_size_strategy", step_size_strategy);
  highs.setOptionValue("kkt_tolerance", adapter.kkt_tolerance);
  // Legacy stops on kkt_tolerance; the rewrite stops on this option, whose
  // default is tighter, so a fair comparison needs both set.
  highs.setOptionValue("pdlp_optimality_tolerance", adapter.kkt_tolerance);
  highs.setOptionValue("time_limit", adapter.time_limit);
  highs.setOptionValue("pdlp_iteration_limit", adapter.iteration_limit);

  const HighsLp& lp = highs.getLp();
  const HighsOptions& options = highs.getOptions();
  HighsBasis basis;
  HighsSolution solution;
  HighsInfo info;
  HighsModelStatus model_status = HighsModelStatus::kNotset;
  HighsCallback callback(&highs);
  HighsTimer timer;
  hipdlp::api::SolveResult rewrite_result;
  const hipdlp::api::MigrationOptions migration(implementation,
                                                migration_backend);

  // Solver timer boundary: identical for legacy and rewrite.
  const std::chrono::steady_clock::time_point solve_start =
      std::chrono::steady_clock::now();
  const HighsStatus status =
      hipdlp::api::solve(options, timer, lp, basis, solution, model_status,
                         info, callback, migration, &rewrite_result);
  const double solver_wall_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    solve_start)
          .count();

  const hipdlp::verification::KktVerificationResult kkt =
      hipdlp::verification::verifyKkt(lp, solution,
                                      adapter.verification_tolerance);

  bool certificate_pass = false;
  const char* certificate_kind =
      certificateKindName(rewrite_result.certificate.kind);
  if (rewrite_result.certificate.kind ==
      hipdlp::api::CertificateKind::kPrimalUnboundedness) {
    certificate_pass = hipdlp::verification::verifyPrimalRay(
                           lp, rewrite_result.certificate.column_direction,
                           adapter.verification_tolerance)
                           .pass;
  } else if (rewrite_result.certificate.kind ==
                 hipdlp::api::CertificateKind::kPrimalInfeasibility ||
             rewrite_result.certificate.kind ==
                 hipdlp::api::CertificateKind::kAmbiguous) {
    certificate_pass = hipdlp::verification::verifyFarkasCertificate(
                           lp, rewrite_result.certificate.row_direction,
                           rewrite_result.certificate.column_direction,
                           adapter.verification_tolerance)
                           .pass;
  }

  double objective = 0.0;
  bool objective_available =
      solution.value_valid && std::isfinite(lp.offset_) &&
      solution.col_value.size() == static_cast<std::size_t>(lp.num_col_);
  if (objective_available) {
    objective = lp.offset_;
    for (HighsInt col = 0; col < lp.num_col_; ++col) {
      if (!std::isfinite(lp.col_cost_[col]) ||
          !std::isfinite(solution.col_value[col])) {
        objective_available = false;
        break;
      }
      objective += lp.col_cost_[col] * solution.col_value[col];
      if (!std::isfinite(objective)) {
        objective_available = false;
        break;
      }
    }
  }
  const bool iterations_available =
      status != HighsStatus::kError && info.pdlp_iteration_count >= 0;

  const double end_to_end_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    end_to_end_start)
          .count();
  const char* return_status_name =
      status == HighsStatus::kOk
          ? "ok"
          : (status == HighsStatus::kWarning ? "warning" : "error");

  std::string json;
  json += "{\n";
  json += formatLine("  \"schema\": \"hipdlp-migration-adapter/1\",");
  json += quotedField("instance", adapter.instance);
  json += quotedField("implementation", adapter.implementation);
  json += quotedField("backend", effective_backend);
  json += quotedField("requested_backend", adapter.backend);
  json += quotedField("step_size", adapter.step_size);
  json += quotedField("status", modelStatusName(model_status));
  json += quotedField("call_status", return_status_name);
  json +=
      intField("iterations", static_cast<long long>(info.pdlp_iteration_count),
               iterations_available);
  json += nullableDoubleField("objective", objective, objective_available);
  json += doubleField("solver_wall_sec", solver_wall_sec);
  json += doubleField("end_to_end_sec", end_to_end_sec);
  json += doubleField("relative_primal_residual",
                      kkt.metrics.relative_primal_residual);
  json +=
      doubleField("relative_dual_residual", kkt.metrics.relative_dual_residual);
  json +=
      doubleField("relative_objective_gap", kkt.metrics.relative_objective_gap);
  json += boolField("kkt_input_valid", kkt.input_valid);
  json += boolField("kkt_pass", kkt.pass);
  json += quotedField("certificate_kind", certificate_kind);
  json += boolField("certificate_pass", certificate_pass);
  json += "  \"options\": {\n";
  json += "    \"presolve\": \"off\",\n";
  json +=
      formatLine("    \"hipdlp_scaling\": \"%s\",", adapter.scaling.c_str());
  json += formatLine("    \"pdlp_step_size_strategy\": %lld,",
                     static_cast<long long>(step_size_strategy));
  json += formatLine("    \"kkt_tolerance\": %.17g,", adapter.kkt_tolerance);
  json += formatLine("    \"time_limit\": %.17g,", adapter.time_limit);
  json += formatLine("    \"pdlp_iteration_limit\": %lld",
                     static_cast<long long>(adapter.iteration_limit));
  json += "  },\n";
  json += quotedField("executable", argv[0]);
  json += formatLine("  \"highs_version\": \"%d.%d.%d\"", HIGHS_VERSION_MAJOR,
                     HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH);
  json += "}\n";

  const std::string solver_name =
      is_legacy ? "hipdlp-legacy"
                : std::string("hipdlp-rewrite-") + effective_backend;
  std::string common_json;
  common_json += "{\n";
  common_json += formatLine("  \"schema\": \"lp-cross-solver-result/1\",");
  common_json += quotedField("solver", solver_name);
  common_json += quotedField("instance", adapter.instance);
  common_json += quotedField("status", modelStatusName(model_status));
  common_json += quotedField("return_status", return_status_name);
  common_json += nullableDoubleField("solver_time_sec", solver_wall_sec, true);
  common_json +=
      intField("iterations", static_cast<long long>(info.pdlp_iteration_count),
               iterations_available);
  common_json +=
      nullableDoubleField("primal_objective", objective, objective_available);
  common_json += nullableDoubleField("relative_primal_residual",
                                     kkt.metrics.relative_primal_residual,
                                     kkt.metrics_finite);
  common_json += nullableDoubleField("relative_dual_residual",
                                     kkt.metrics.relative_dual_residual,
                                     kkt.metrics_finite);
  common_json += nullableDoubleField("relative_objective_gap",
                                     kkt.metrics.relative_objective_gap,
                                     kkt.metrics_finite);
  common_json += boolField("kkt_input_valid", kkt.input_valid);
  common_json += boolField("kkt_pass", kkt.pass);
  common_json += quotedField("certificate_kind", certificate_kind);
  common_json += formatLine("  \"certificate_pass\": %s",
                            certificate_pass ? "true" : "false");
  common_json += "}\n";

  std::fputs(json.c_str(), stdout);
  if (!adapter.output_dir.empty()) {
    const std::string summary_path =
        adapter.output_dir + "/hipdlp_summary.json";
    FILE* summary_file = std::fopen(summary_path.c_str(), "w");
    if (summary_file == nullptr) {
      std::fprintf(stderr, "failed to write %s\n", summary_path.c_str());
      return 4;
    }
    std::fputs(json.c_str(), summary_file);
    std::fclose(summary_file);

    const std::string result_path = adapter.output_dir + "/result.json";
    FILE* result_file = std::fopen(result_path.c_str(), "w");
    if (result_file == nullptr) {
      std::fprintf(stderr, "failed to write %s\n", result_path.c_str());
      return 4;
    }
    std::fputs(common_json.c_str(), result_file);
    std::fclose(result_file);
  }

  return status == HighsStatus::kError ? 1 : 0;
}
