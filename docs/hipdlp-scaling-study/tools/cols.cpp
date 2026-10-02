#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
using namespace hipdlp::model;
int main(int argc, char** argv) {
  Highs highs; highs.setOptionValue("output_flag", false); highs.readModel(argv[1]);
  HighsLp lp = highs.getLp(); lp.a_matrix_.ensureColwise();
  ProblemBuildResult built = buildProblem(lp); const Problem& p = built.problem;
  ScalingOptions o; o.mode = ScalingMode::kAuto; o.legacy_scaling_mode = kPdlpScalingRuiz + kPdlpScalingPC;
  ScalingResult r = computeScalingPipeline(p.matrix(), o);
  const auto& C = r.scaling.col_scale; const auto& Rr = r.scaling.row_scale;
  FILE* f = std::fopen(argv[2], "r"); int nc, nr; if (fscanf(f, "%d %d", &nc, &nr) != 2) return 1;
  std::vector<double> x(nc), d(nc), y(nr);
  for (auto& v : x) if (fscanf(f, "%lf", &v) != 1) return 1;
  for (auto& v : d) if (fscanf(f, "%lf", &v) != 1) return 1;
  for (auto& v : y) if (fscanf(f, "%lf", &v) != 1) return 1;
  std::vector<int> idx(nc); for (int j = 0; j < nc; ++j) idx[j] = j;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return C[a] < C[b]; });
  // transformed row index for each original row
  const auto& tr = p.transformation().original_to_transformed_row;
  for (int k = 0; k < 6; ++k) {
    int j = idx[k];
    std::printf("col %d  C=%.3g  cost=%.4g  bounds=[%.4g, %.4g]  x=%.4g  reduced cost(d)=%.4g\n", j, C[j], lp.col_cost_[j], lp.col_lower_[j], lp.col_upper_[j], x[j], d[j]);
    int n = lp.a_matrix_.start_[j + 1] - lp.a_matrix_.start_[j];
    std::printf("   %d entries:", n);
    for (int e = lp.a_matrix_.start_[j]; e < std::min(lp.a_matrix_.start_[j + 1], lp.a_matrix_.start_[j] + 6); ++e) {
      int row = lp.a_matrix_.index_[e];
      std::printf(" (row %d: a=%.3g, R=%.3g)", row, lp.a_matrix_.value_[e], Rr[tr[row]]);
    }
    std::printf("\n");
  }
  // rows with the largest R: their largest entry magnitude and size
  std::vector<int> ridx(p.numRows()); for (int i = 0; i < p.numRows(); ++i) ridx[i] = i;
  std::sort(ridx.begin(), ridx.end(), [&](int a, int b) { return Rr[a] > Rr[b]; });
  std::printf("\nrows with largest R:\n");
  const SparseMatrix& A = p.matrix();
  for (int k = 0; k < 6; ++k) {
    int i = ridx[k]; double mx = 0; int cnt = A.starts[i + 1] - A.starts[i];
    for (int e = A.starts[i]; e < A.starts[i + 1]; ++e) mx = std::max(mx, std::fabs(A.values[e]));
    std::printf("  transformed row %d: R=%.3g  nnz=%d  max|a|=%.3g  rhs=%.4g\n", i, Rr[i], cnt, mx, p.rowLower()[i] / Rr[i]);
  }
  return 0;
}
