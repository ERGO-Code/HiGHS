#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {

struct AdapterOptions {
  std::string instance;
  std::string solver;
  std::string presolve;
  double time_limit = 0.0;
  double kkt_tolerance = 0.0;
  std::string output_dir;
};

struct CertificateResult {
  std::string kind = "none";
  bool pass = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: lp_solver_validation_adapter --instance FILE\n"
               "         --solver simplex|ipm|pdlp\n"
               "         --presolve on|off\n"
               "         --time-limit SECONDS\n"
               "         --kkt-tolerance TOLERANCE\n"
               "         --output-dir DIRECTORY\n");
}

bool nextValue(const int argc, char** argv, int& index, std::string& value) {
  if (index + 1 >= argc) return false;
  value = argv[++index];
  return true;
}

bool parseDouble(const std::string& text, double& value) {
  char* end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return end != text.c_str() && *end == '\0' && std::isfinite(value);
}

bool parseArguments(const int argc, char** argv, AdapterOptions& options) {
  for (int index = 1; index < argc; ++index) {
    const std::string flag = argv[index];
    std::string value;
    if (flag == "--instance") {
      if (!nextValue(argc, argv, index, options.instance)) return false;
    } else if (flag == "--solver") {
      if (!nextValue(argc, argv, index, options.solver)) return false;
    } else if (flag == "--presolve") {
      if (!nextValue(argc, argv, index, options.presolve)) return false;
    } else if (flag == "--time-limit") {
      if (!nextValue(argc, argv, index, value) ||
          !parseDouble(value, options.time_limit))
        return false;
    } else if (flag == "--kkt-tolerance") {
      if (!nextValue(argc, argv, index, value) ||
          !parseDouble(value, options.kkt_tolerance))
        return false;
    } else if (flag == "--output-dir") {
      if (!nextValue(argc, argv, index, options.output_dir)) return false;
    } else {
      std::fprintf(stderr, "unknown flag %s\n", flag.c_str());
      return false;
    }
  }
  const bool solver_ok = options.solver == "simplex" ||
                         options.solver == "ipm" || options.solver == "pdlp";
  const bool presolve_ok =
      options.presolve == "on" || options.presolve == "off";
  return !options.instance.empty() && solver_ok && presolve_ok &&
         options.time_limit >= 0.0 && options.kkt_tolerance > 0.0 &&
         !options.output_dir.empty();
}

const char* modelStatusName(const HighsModelStatus status) {
  switch (status) {
    case HighsModelStatus::kNotset:
      return "NOT_SET";
    case HighsModelStatus::kLoadError:
      return "LOAD_ERROR";
    case HighsModelStatus::kModelError:
      return "MODEL_ERROR";
    case HighsModelStatus::kPresolveError:
      return "PRESOLVE_ERROR";
    case HighsModelStatus::kSolveError:
      return "SOLVE_ERROR";
    case HighsModelStatus::kPostsolveError:
      return "POSTSOLVE_ERROR";
    case HighsModelStatus::kModelEmpty:
      return "MODEL_EMPTY";
    case HighsModelStatus::kOptimal:
      return "OPTIMAL";
    case HighsModelStatus::kInfeasible:
      return "INFEASIBLE";
    case HighsModelStatus::kUnboundedOrInfeasible:
      return "UNBOUNDED_OR_INFEASIBLE";
    case HighsModelStatus::kUnbounded:
      return "UNBOUNDED";
    case HighsModelStatus::kObjectiveBound:
      return "OBJECTIVE_BOUND";
    case HighsModelStatus::kObjectiveTarget:
      return "OBJECTIVE_TARGET";
    case HighsModelStatus::kTimeLimit:
      return "TIME_LIMIT";
    case HighsModelStatus::kIterationLimit:
      return "ITERATION_LIMIT";
    case HighsModelStatus::kUnknown:
      return "UNKNOWN";
    default:
      return "OTHER";
  }
}

const char* returnStatusName(const HighsStatus status) {
  if (status == HighsStatus::kOk) return "ok";
  if (status == HighsStatus::kWarning) return "warning";
  return "error";
}

std::string formatLine(const char* format, ...) {
  char buffer[1024];
  va_list arguments;
  va_start(arguments, format);
  std::vsnprintf(buffer, sizeof(buffer), format, arguments);
  va_end(arguments);
  return std::string(buffer) + "\n";
}

std::string escapeJson(const std::string& value) {
  std::string escaped;
  for (std::size_t index = 0; index < value.size(); ++index) {
    const unsigned char character =
        static_cast<unsigned char>(value[index]);
    switch (character) {
      case '\"':
        escaped += "\\\"";
        break;
      case '\\':
        escaped += "\\\\";
        break;
      case '\b':
        escaped += "\\b";
        break;
      case '\f':
        escaped += "\\f";
        break;
      case '\n':
        escaped += "\\n";
        break;
      case '\r':
        escaped += "\\r";
        break;
      case '\t':
        escaped += "\\t";
        break;
      default:
        if (character < 0x20) {
          char encoded[7];
          std::snprintf(encoded, sizeof(encoded), "\\u%04x", character);
          escaped += encoded;
        } else {
          escaped += static_cast<char>(character);
        }
    }
  }
  return escaped;
}

std::string quotedField(const char* key, const std::string& value) {
  return formatLine("  \"%s\": \"%s\",", key, escapeJson(value).c_str());
}

std::string boolField(const char* key, const bool value) {
  return formatLine("  \"%s\": %s,", key, value ? "true" : "false");
}

std::string doubleField(const char* key, const double value,
                        const bool available) {
  if (!available || !std::isfinite(value))
    return formatLine("  \"%s\": null,", key);
  return formatLine("  \"%s\": %.17g,", key, value);
}

std::string integerField(const char* key, const HighsInt value,
                         const bool available) {
  if (!available) return formatLine("  \"%s\": null,", key);
  return formatLine("  \"%s\": %lld,", key,
                    static_cast<long long>(value));
}

HighsInt solverIterations(const std::string& solver, const HighsInfo& info) {
  if (solver == "simplex") return info.simplex_iteration_count;
  if (solver == "ipm") return info.ipm_iteration_count;
  return info.pdlp_iteration_count;
}

bool computeObjective(const HighsLp& lp, const HighsSolution& solution,
                      double& objective) {
  if (!solution.value_valid ||
      solution.col_value.size() != static_cast<std::size_t>(lp.num_col_))
    return false;
  objective = lp.offset_;
  for (HighsInt column = 0; column < lp.num_col_; ++column) {
    if (!std::isfinite(solution.col_value[column])) return false;
    objective += lp.col_cost_[column] * solution.col_value[column];
    if (!std::isfinite(objective)) return false;
  }
  return true;
}

CertificateResult verifyPrimalCertificate(Highs& highs, const HighsLp& lp,
                                          const double tolerance) {
  CertificateResult result;
  if (lp.num_col_ <= 0) return result;
  std::vector<double> direction(lp.num_col_, 0.0);
  bool has_ray = false;
  if (highs.getPrimalRay(has_ray, direction.data()) != HighsStatus::kOk ||
      !has_ray)
    return result;
  result.kind = "primal_unboundedness";
  result.pass =
      hipdlp::verification::verifyPrimalRay(lp, direction, tolerance).pass;
  return result;
}

CertificateResult verifyDualCertificate(Highs& highs, const HighsLp& lp,
                                        const double tolerance) {
  CertificateResult result;
  if (lp.num_row_ <= 0) return result;
  std::vector<double> row_direction(lp.num_row_, 0.0);
  bool has_ray = false;
  if (highs.getDualRay(has_ray, row_direction.data()) != HighsStatus::kOk ||
      !has_ray)
    return result;
  std::vector<double> column_direction;
  lp.a_matrix_.productTranspose(column_direction, row_direction);
  for (std::size_t column = 0; column < column_direction.size(); ++column)
    column_direction[column] = -column_direction[column];
  result.kind = "primal_infeasibility";
  result.pass = hipdlp::verification::verifyFarkasCertificate(
                    lp, row_direction, column_direction, tolerance)
                    .pass;
  return result;
}

CertificateResult verifyCertificate(Highs& highs, const HighsLp& lp,
                                    const HighsModelStatus status,
                                    const double tolerance) {
  if (status == HighsModelStatus::kUnbounded)
    return verifyPrimalCertificate(highs, lp, tolerance);
  if (status == HighsModelStatus::kInfeasible)
    return verifyDualCertificate(highs, lp, tolerance);
  if (status == HighsModelStatus::kUnboundedOrInfeasible) {
    CertificateResult result = verifyPrimalCertificate(highs, lp, tolerance);
    if (result.kind != "none") return result;
    return verifyDualCertificate(highs, lp, tolerance);
  }
  return CertificateResult();
}

}  // namespace

int main(int argc, char** argv) {
  AdapterOptions options;
  if (!parseArguments(argc, argv, options)) {
    usage();
    return 2;
  }

  Highs highs;
  if (highs.setOptionValue("output_flag", false) != HighsStatus::kOk ||
      highs.readModel(options.instance) != HighsStatus::kOk) {
    std::fprintf(stderr, "failed to read %s\n", options.instance.c_str());
    return 3;
  }
  if (highs.setOptionValue("solver", options.solver) != HighsStatus::kOk ||
      highs.setOptionValue("presolve", options.presolve) != HighsStatus::kOk ||
      highs.setOptionValue("time_limit", options.time_limit) !=
          HighsStatus::kOk ||
      highs.setOptionValue("kkt_tolerance", options.kkt_tolerance) !=
          HighsStatus::kOk) {
    std::fprintf(stderr, "failed to set solver options\n");
    return 3;
  }

  const std::chrono::steady_clock::time_point solve_start =
      std::chrono::steady_clock::now();
  const HighsStatus return_status = highs.run();
  const double solver_time =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    solve_start)
          .count();

  const HighsModelStatus model_status = highs.getModelStatus();
  const HighsLp lp = highs.getLp();
  const HighsSolution solution = highs.getSolution();
  const HighsInfo info = highs.getInfo();
  const hipdlp::verification::KktVerificationResult kkt =
      hipdlp::verification::verifyKkt(lp, solution, options.kkt_tolerance);
  const CertificateResult certificate =
      verifyCertificate(highs, lp, model_status, options.kkt_tolerance);

  double primal_objective = 0.0;
  const bool objective_available =
      computeObjective(lp, solution, primal_objective);
  const HighsInt iterations = solverIterations(options.solver, info);

  std::string json;
  json += "{\n";
  json += formatLine("  \"schema\": \"lp-cross-solver-result/1\",");
  json += quotedField("solver", options.solver);
  json += quotedField("instance", options.instance);
  json += quotedField("status", modelStatusName(model_status));
  json += quotedField("return_status", returnStatusName(return_status));
  json += doubleField("solver_time_sec", solver_time, true);
  json += integerField("iterations", iterations, info.valid && iterations >= 0);
  json += doubleField("primal_objective", primal_objective,
                      objective_available);
  json += doubleField("relative_primal_residual",
                      kkt.metrics.relative_primal_residual,
                      kkt.metrics_finite);
  json += doubleField("relative_dual_residual",
                      kkt.metrics.relative_dual_residual, kkt.metrics_finite);
  json += doubleField("relative_objective_gap",
                      kkt.metrics.relative_objective_gap, kkt.metrics_finite);
  json += boolField("kkt_input_valid", kkt.input_valid);
  json += boolField("kkt_pass", kkt.pass);
  json += quotedField("certificate_kind", certificate.kind);
  json += formatLine("  \"certificate_pass\": %s",
                     certificate.pass ? "true" : "false");
  json += "}\n";

  const std::string output_path = options.output_dir + "/result.json";
  std::ofstream output(output_path.c_str());
  if (!output) {
    std::fprintf(stderr, "failed to write %s\n", output_path.c_str());
    return 4;
  }
  output << json;
  output.close();
  std::fputs(json.c_str(), stdout);
  return return_status == HighsStatus::kError ? 1 : 0;
}