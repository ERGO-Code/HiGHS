#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include "Highs.h"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

bool read(const char* path, int count, std::vector<double>& values) {
  std::ifstream file(path);
  double value;
  while (file >> value) {
    if (!std::isfinite(value)) return false;
    values.push_back(value);
  }
  return values.size() == static_cast<std::size_t>(count);
}

// Passive verification with the solver's reported bound duals. Reconstructing
// bound duals from rounded row dual exports can magnify roundoff at huge bounds.
int main(int argc, char** argv) {
  if (argc != 6) return 2;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const auto& lp = highs.getLp();
  HighsSolution solution;
  if (!read(argv[2], lp.num_col_, solution.col_value) ||
      !read(argv[3], lp.num_row_, solution.row_dual) ||
      !read(argv[4], lp.num_col_, solution.col_dual)) return 4;
  solution.value_valid = solution.dual_valid = true;
  lp.a_matrix_.product(solution.row_value, solution.col_value);
  const auto kkt = hipdlp::verification::verifyKkt(lp, solution, 1e-4);
  const auto& m = kkt.metrics;
  std::ofstream output(argv[5]);
  output.precision(17);
  output << "{\"kkt_contract\":\"hipdlp-original-kkt/2\",\"kkt_tolerance\":0.0001,\"kkt_pass\":"
         << (kkt.pass ? "true" : "false")
         << ",\"primal_objective\":" << m.primal_objective
         << ",\"relative_primal_residual\":" << m.relative_primal_residual
         << ",\"relative_dual_residual\":" << m.relative_dual_residual
         << ",\"relative_objective_gap\":" << m.relative_objective_gap
         << ",\"relative_complementarity\":" << m.relative_complementarity
         << ",\"max_stationarity_residual\":" << m.max_stationarity_residual
         << ",\"reduced_cost_source\":\"reported\"}\n";
  return 0;
}
