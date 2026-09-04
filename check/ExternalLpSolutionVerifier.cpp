#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {

struct VerifierOptions {
  std::string instance;
  std::string solution;
  std::string primal_values;
  std::string row_dual_values;
  std::string solver;
  std::string status;
  double solver_time = 0.0;
  HighsInt iterations = 0;
  double tolerance = 0.0;
  std::string output;
  bool derive_reduced_costs = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: external_lp_solution_verifier --instance FILE\n"
               "         (--solution HIGHS_RAW_SOLUTION |\n"
               "          --primal-values FILE --row-dual-values FILE)\n"
               "         --solver NAME --status STATUS\n"
               "         --solver-time SECONDS --iterations COUNT\n"
               "         --tolerance TOLERANCE --output FILE\n"
               "         [--derive-reduced-costs]\n");
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

bool parseInteger(const std::string& text, HighsInt& value) {
  char* end = nullptr;
  const long long parsed = std::strtoll(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || parsed < 0) return false;
  value = static_cast<HighsInt>(parsed);
  return static_cast<long long>(value) == parsed;
}

bool parseArguments(const int argc, char** argv, VerifierOptions& options) {
  for (int index = 1; index < argc; ++index) {
    const std::string flag = argv[index];
    std::string value;
    if (flag == "--instance") {
      if (!nextValue(argc, argv, index, options.instance)) return false;
    } else if (flag == "--solution") {
      if (!nextValue(argc, argv, index, options.solution)) return false;
    } else if (flag == "--primal-values") {
      if (!nextValue(argc, argv, index, options.primal_values)) return false;
    } else if (flag == "--row-dual-values") {
      if (!nextValue(argc, argv, index, options.row_dual_values)) return false;
    } else if (flag == "--solver") {
      if (!nextValue(argc, argv, index, options.solver)) return false;
    } else if (flag == "--status") {
      if (!nextValue(argc, argv, index, options.status)) return false;
    } else if (flag == "--solver-time") {
      if (!nextValue(argc, argv, index, value) ||
          !parseDouble(value, options.solver_time))
        return false;
    } else if (flag == "--iterations") {
      if (!nextValue(argc, argv, index, value) ||
          !parseInteger(value, options.iterations))
        return false;
    } else if (flag == "--tolerance") {
      if (!nextValue(argc, argv, index, value) ||
          !parseDouble(value, options.tolerance))
        return false;
    } else if (flag == "--output") {
      if (!nextValue(argc, argv, index, options.output)) return false;
    } else if (flag == "--derive-reduced-costs") {
      options.derive_reduced_costs = true;
    } else {
      std::fprintf(stderr, "unknown flag %s\n", flag.c_str());
      return false;
    }
  }
    const bool raw_solution = !options.solution.empty();
    const bool vector_solution = !options.primal_values.empty() &&
             !options.row_dual_values.empty();
    return !options.instance.empty() && raw_solution != vector_solution &&
      !options.solver.empty() && !options.status.empty() &&
         options.solver_time >= 0.0 && options.tolerance > 0.0 &&
         !options.output.empty();
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

bool finiteVector(const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

bool readNumericVector(const std::string& path, const HighsInt expected_size,
                       std::vector<double>& values) {
  std::ifstream input(path.c_str());
  if (!input) return false;
  values.clear();
  double value = 0.0;
  while (input >> value) {
    if (!std::isfinite(value)) return false;
    values.push_back(value);
  }
  return input.eof() &&
         values.size() == static_cast<std::size_t>(expected_size);
}

bool countHeaderMatches(const std::string& line, const char* label,
                        const HighsInt expected_count) {
  std::istringstream stream(line);
  std::string hash;
  std::string parsed_label;
  long long parsed_count = -1;
  stream >> hash >> parsed_label >> parsed_count;
  return hash == "#" && parsed_label == label &&
         parsed_count == static_cast<long long>(expected_count);
}

bool hasCompleteRawDualSection(const std::string& path,
                               const HighsInt num_col,
                               const HighsInt num_row) {
  std::ifstream input(path.c_str());
  std::string line;
  while (std::getline(input, line) && line != "# Dual solution values") {
  }
  if (!input || !std::getline(input, line) || line == "None") return false;
  if (!std::getline(input, line) ||
      !countHeaderMatches(line, "Columns", num_col))
    return false;
  for (HighsInt column = 0; column < num_col; ++column) {
    if (!std::getline(input, line)) return false;
  }
  if (!std::getline(input, line) ||
      !countHeaderMatches(line, "Rows", num_row))
    return false;
  for (HighsInt row = 0; row < num_row; ++row) {
    if (!std::getline(input, line)) return false;
  }
  return true;
}

bool computeObjective(const HighsLp& lp, const HighsSolution& solution,
                      double& objective) {
  if (!solution.value_valid ||
      solution.col_value.size() != static_cast<std::size_t>(lp.num_col_) ||
      !finiteVector(solution.col_value))
    return false;
  objective = lp.offset_;
  for (HighsInt column = 0; column < lp.num_col_; ++column) {
    objective += lp.col_cost_[column] * solution.col_value[column];
    if (!std::isfinite(objective)) return false;
  }
  return true;
}

bool deriveReducedCosts(const HighsLp& lp, const HighsSolution& solution,
                        std::vector<double>& reduced_costs) {
  if (!solution.dual_valid ||
      solution.row_dual.size() != static_cast<std::size_t>(lp.num_row_) ||
      solution.col_dual.size() != static_cast<std::size_t>(lp.num_col_) ||
      !finiteVector(solution.row_dual) || !finiteVector(solution.col_dual))
    return false;
  std::vector<double> transpose_product;
  lp.a_matrix_.productTranspose(transpose_product, solution.row_dual);
  if (transpose_product.size() != static_cast<std::size_t>(lp.num_col_) ||
      !finiteVector(transpose_product))
    return false;
  reduced_costs.resize(lp.num_col_);
  for (HighsInt column = 0; column < lp.num_col_; ++column) {
    reduced_costs[column] = lp.col_cost_[column] - transpose_product[column];
    if (!std::isfinite(reduced_costs[column])) return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  VerifierOptions options;
  if (!parseArguments(argc, argv, options)) {
    usage();
    return 2;
  }

  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(options.instance) != HighsStatus::kOk) {
    std::fprintf(stderr, "failed to read %s\n", options.instance.c_str());
    return 3;
  }

  const HighsLp& lp = highs.getLp();
  HighsSolution solution;
  bool reported_reduced_cost_consistency_available = false;
  if (!options.solution.empty()) {
    if (highs.readSolution(options.solution, kSolutionStyleRaw) !=
        HighsStatus::kOk) {
      std::fprintf(stderr, "failed to read %s\n", options.solution.c_str());
      return 3;
    }
    solution = highs.getSolution();
    solution.dual_valid = hasCompleteRawDualSection(
        options.solution, lp.num_col_, lp.num_row_);
    reported_reduced_cost_consistency_available = solution.dual_valid;
  } else {
    if (!readNumericVector(options.primal_values, lp.num_col_,
                           solution.col_value) ||
        !readNumericVector(options.row_dual_values, lp.num_row_,
                           solution.row_dual)) {
      std::fprintf(stderr, "failed to read complete finite solution vectors\n");
      return 3;
    }
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_dual.assign(lp.num_col_, 0.0);
    options.derive_reduced_costs = true;
  }
  solution.row_value.clear();
  if (solution.value_valid &&
      solution.col_value.size() == static_cast<std::size_t>(lp.num_col_))
    lp.a_matrix_.product(solution.row_value, solution.col_value);

  std::vector<double> derived_reduced_costs;
  const bool reduced_cost_consistency_available =
      deriveReducedCosts(lp, solution, derived_reduced_costs);
  double reported_reduced_cost_max_abs_error = 0.0;
  if (reduced_cost_consistency_available) {
    for (HighsInt column = 0; column < lp.num_col_; ++column) {
      reported_reduced_cost_max_abs_error = std::max(
          reported_reduced_cost_max_abs_error,
          std::fabs(solution.col_dual[column] -
                    derived_reduced_costs[column]));
    }
  }
  const bool reported_reduced_cost_consistency_pass =
      reported_reduced_cost_consistency_available &&
      reduced_cost_consistency_available &&
      reported_reduced_cost_max_abs_error <= options.tolerance;
  if (options.derive_reduced_costs && !reduced_cost_consistency_available) {
    std::fprintf(stderr, "cannot derive reduced costs from invalid row duals\n");
    return 3;
  }
  if (options.derive_reduced_costs)
    solution.col_dual = derived_reduced_costs;

  const hipdlp::verification::KktVerificationResult kkt =
      hipdlp::verification::verifyKkt(lp, solution, options.tolerance);
  double primal_objective = 0.0;
  const bool objective_available =
      computeObjective(lp, solution, primal_objective);

  std::string json;
  json += "{\n";
  json += formatLine("  \"schema\": \"lp-cross-solver-result/1\",");
  json += quotedField("solver", options.solver);
  json += quotedField("instance", options.instance);
  json += quotedField("status", options.status);
  json += formatLine("  \"return_status\": \"ok\",");
  json += doubleField("solver_time_sec", options.solver_time, true);
  json += formatLine("  \"iterations\": %lld,",
                     static_cast<long long>(options.iterations));
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
  json += quotedField("reduced_cost_source",
                      options.derive_reduced_costs ? "derived_from_row_dual"
                                                   : "reported");
  json += boolField("reported_reduced_cost_consistency_available",
                    reported_reduced_cost_consistency_available);
  json += boolField("reported_reduced_cost_consistency_pass",
                    reported_reduced_cost_consistency_pass);
  json += doubleField("reported_reduced_cost_max_abs_error",
                      reported_reduced_cost_max_abs_error,
                      reported_reduced_cost_consistency_available &&
                          reduced_cost_consistency_available);
  json += formatLine("  \"certificate_kind\": \"none\",");
  json += formatLine("  \"certificate_pass\": false");
  json += "}\n";

  std::ofstream output(options.output.c_str());
  if (!output) {
    std::fprintf(stderr, "failed to write %s\n", options.output.c_str());
    return 4;
  }
  output << json;
  output.close();
  std::fputs(json.c_str(), stdout);
  return 0;
}