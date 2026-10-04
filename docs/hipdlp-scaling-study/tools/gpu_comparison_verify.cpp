#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Highs.h"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

static bool readValues(const std::string& path, HighsInt count,
                       std::vector<double>& values) {
  std::ifstream file(path);
  double v;
  while (file >> v) {
    if (!std::isfinite(v)) return false;
    values.push_back(v);
  }
  return file.eof() && values.size() == static_cast<std::size_t>(count);
}

// Keep the read-only LP in memory while verifying serial solver runs. Input
// lines contain primal, row-dual, reduced-cost and output paths, tab separated.
int main(int argc, char** argv) {
  if (argc != 2) return 2;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const auto& lp = highs.getLp();
  std::cout << "{\"rows\":" << lp.num_row_ << ",\"cols\":" << lp.num_col_
            << ",\"nnz\":" << lp.a_matrix_.value_.size() << "}" << std::endl;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream request(line);
    std::string x, y, z, output;
    if (!std::getline(request, x, '\t') || !std::getline(request, y, '\t') ||
        !std::getline(request, z, '\t') || !std::getline(request, output)) return 4;
    std::ofstream file(output);
    if (!file) return 5;
    HighsSolution solution;
    if (!readValues(x, lp.num_col_, solution.col_value) ||
        !readValues(y, lp.num_row_, solution.row_dual) ||
        !readValues(z, lp.num_col_, solution.col_dual)) {
      file << "{\"kkt_pass\":false,\"unverified\":true,\"reason\":\"Invalid export\"}\n";
    } else {
      solution.value_valid = solution.dual_valid = true;
      lp.a_matrix_.product(solution.row_value, solution.col_value);
      const auto checked = hipdlp::verification::verifyKkt(lp, solution, 1e-4);
      if (!checked.input_valid || !checked.metrics_finite) {
        file << "{\"kkt_pass\":false,\"unverified\":true,\"reason\":\"Invalid KKT metrics\"}\n";
      } else {
        const auto& m = checked.metrics;
        file.precision(17);
        file << "{\"kkt_contract\":\"hipdlp-original-kkt/2\",\"tolerance\":0.0001,\"kkt_pass\":"
             << (checked.pass ? "true" : "false")
             << ",\"objective\":" << m.primal_objective
             << ",\"primal\":" << m.relative_primal_residual
             << ",\"dual\":" << m.relative_dual_residual
             << ",\"gap\":" << m.relative_objective_gap
             << ",\"complementarity\":" << m.relative_complementarity
             << ",\"stationarity\":" << m.max_stationarity_residual << "}\n";
      }
    }
    file.close();
    std::cout << "DONE" << std::endl;
  }
}
