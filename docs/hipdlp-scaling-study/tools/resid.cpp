// Where does the unscaled dual residual live? Join per-column dual violation with the column scale.
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
  // resid MPS SOL CONFIG(A|L)
  Highs highs; highs.setOptionValue("output_flag", false); highs.readModel(argv[1]);
  const HighsLp& lp = highs.getLp();
  FILE* f = std::fopen(argv[2], "r"); int nc, nr; std::fscanf(f, "%d %d", &nc, &nr);
  std::vector<double> x(nc), d(nc), y(nr);
  for (auto& v : x) std::fscanf(f, "%lf", &v);
  for (auto& v : d) std::fscanf(f, "%lf", &v);
  for (auto& v : y) std::fscanf(f, "%lf", &v);
  std::fclose(f);
  ProblemBuildResult built = buildProblem(lp);
  const Problem& p = built.problem;
  ScalingOptions o; o.mode = std::string(argv[3]) == "A" ? ScalingMode::kAuto : ScalingMode::kLegacy;
  o.legacy_scaling_mode = kPdlpScalingRuiz + kPdlpScalingPC;
  ScalingResult r = computeScalingPipeline(p.matrix(), o);
  const auto& C = r.scaling.col_scale;
  const auto& R = r.scaling.row_scale;
  // reduced cost from y: d = c - A^T y (original LP, colwise)
  HighsLp l2 = lp; l2.a_matrix_.ensureColwise();
  std::vector<double> rc(nc);
  double sgn = lp.sense_ == ObjSense::kMinimize ? 1 : -1;
  for (int j = 0; j < nc; ++j) {
    double s = 0;
    for (int e = l2.a_matrix_.start_[j]; e < l2.a_matrix_.start_[j + 1]; ++e) s += l2.a_matrix_.value_[e] * y[l2.a_matrix_.index_[e]];
    rc[j] = lp.col_cost_[j] - s;
  }
  // violation of dual feasibility by bound type, using the solver's reported x for active bound detection
  std::vector<double> v(nc);
  for (int j = 0; j < nc; ++j) {
    double lo = lp.col_lower_[j], up = lp.col_upper_[j]; bool hl = lo > -1e20, hu = up < 1e20;
    double dj = rc[j];
    if (!hl && !hu) v[j] = std::fabs(dj);
    else if (hl && !hu) v[j] = std::max(0.0, -dj);
    else if (!hl && hu) v[j] = std::max(0.0, dj);
    else v[j] = 0;  // boxed: sign handled by complementarity, ignore here
  }
  // buckets by log10 C
  const int NB = 12; double sum_unsc[NB] = {0}, sum_sc[NB] = {0}; int cnt[NB] = {0};
  double tot_u = 0, tot_s = 0;
  for (int j = 0; j < nc; ++j) {
    int b = (int)std::floor(std::log10(C[j])) + 8; b = std::max(0, std::min(NB - 1, b));
    sum_unsc[b] += v[j] * v[j]; sum_sc[b] += (v[j] * C[j]) * (v[j] * C[j]); cnt[b]++;
    tot_u += v[j] * v[j]; tot_s += (v[j] * C[j]) * (v[j] * C[j]);
  }
  std::printf("total dual violation L2: unscaled %.4g  scaled %.4g  (ratio %.3g)\n", std::sqrt(tot_u), std::sqrt(tot_s), std::sqrt(tot_u / std::max(tot_s, 1e-300)));
  std::printf("by decade of column scale C (col counts, share of unscaled sq-violation, share of scaled sq-violation)\n");
  for (int b = 0; b < NB; ++b) if (cnt[b])
    std::printf("  C~1e%-3d cols %7d  unscaled share %6.2f%%  scaled share %6.2f%%\n", b - 8, cnt[b], 100 * sum_unsc[b] / std::max(tot_u, 1e-300), 100 * sum_sc[b] / std::max(tot_s, 1e-300));
  return 0;
}
