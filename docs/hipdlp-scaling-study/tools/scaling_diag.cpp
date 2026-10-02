// Scaling diagnostic: how do the HiPDLP rewrite scaling pipelines change the
// numerical shape of an LP? Uses the rewrite's own buildProblem and scaling
// functions, so the numbers are what the solver actually sees.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

using namespace hipdlp::model;

namespace {

struct Dist {
  double min = 0, p01 = 0, p05 = 0, p50 = 0, p95 = 0, p99 = 0, max = 0;
  size_t n = 0;
};

// Distribution of log10(|v|) over nonzero finite entries.
Dist logDist(std::vector<double> v) {
  std::vector<double> l;
  l.reserve(v.size());
  for (double x : v) {
    double a = std::fabs(x);
    if (a > 0 && std::isfinite(a) && a < 1e20) l.push_back(std::log10(a));
  }
  Dist d;
  d.n = l.size();
  if (l.empty()) return d;
  std::sort(l.begin(), l.end());
  auto q = [&](double f) { return l[std::min(l.size() - 1, (size_t)(f * (l.size() - 1)))]; };
  d.min = l.front(); d.p01 = q(0.01); d.p05 = q(0.05); d.p50 = q(0.5);
  d.p95 = q(0.95); d.p99 = q(0.99); d.max = l.back();
  return d;
}

void printDist(const char* name, const Dist& d) {
  std::printf("  %-14s n=%-9zu log10|.|  min %7.2f p1 %7.2f p5 %7.2f med %7.2f p95 %7.2f p99 %7.2f max %7.2f\n",
              name, d.n, d.min, d.p01, d.p05, d.p50, d.p95, d.p99, d.max);
}

double l2(const std::vector<double>& v, bool finite_only = true) {
  long double s = 0;
  for (double x : v)
    if (!finite_only || (std::isfinite(x) && std::fabs(x) < 1e20)) s += (long double)x * x;
  return std::sqrt((double)s);
}

// sigma_max(A) by power iteration on A^T A.
double sigmaMax(const SparseMatrix& A, const SparseMatrix& AT, int iters = 60) {
  std::vector<double> x(A.num_cols), y(A.num_rows);
  for (int i = 0; i < A.num_cols; ++i) x[i] = 1.0 + 0.37 * std::sin(1.0 * i);
  double nrm = l2(x), est = 0;
  for (auto& v : x) v /= nrm;
  for (int it = 0; it < iters; ++it) {
    for (int r = 0; r < A.num_rows; ++r) {
      double s = 0;
      for (int e = A.starts[r]; e < A.starts[r + 1]; ++e) s += A.values[e] * x[A.indices[e]];
      y[r] = s;
    }
    for (int c = 0; c < A.num_cols; ++c) {
      double s = 0;
      for (int e = AT.starts[c]; e < AT.starts[c + 1]; ++e) s += AT.values[e] * y[AT.indices[e]];
      x[c] = s;
    }
    nrm = l2(x);
    est = std::sqrt(nrm);
    if (nrm == 0) break;
    for (auto& v : x) v /= nrm;
  }
  return est;
}

struct Scaled {
  SparseMatrix A, AT;
  std::vector<double> c, b, lo, up, R, C;
};

SparseMatrix transposeOf(const SparseMatrix& A) {
  SparseMatrix T;
  T.num_rows = A.num_cols; T.num_cols = A.num_rows;
  T.starts.assign(T.num_rows + 1, 0);
  for (int e = 0; e < A.numNonzeros(); ++e) T.starts[A.indices[e] + 1]++;
  for (int i = 0; i < T.num_rows; ++i) T.starts[i + 1] += T.starts[i];
  T.indices.resize(A.numNonzeros()); T.values.resize(A.numNonzeros());
  std::vector<int> pos(T.starts.begin(), T.starts.end() - 1);
  for (int r = 0; r < A.num_rows; ++r)
    for (int e = A.starts[r]; e < A.starts[r + 1]; ++e) {
      int p = pos[A.indices[e]]++;
      T.indices[p] = r; T.values[p] = A.values[e];
    }
  return T;
}

Scaled applyScale(const Problem& p, const std::vector<double>& R, const std::vector<double>& C) {
  Scaled s;
  s.R = R; s.C = C;
  s.A = p.matrix();
  for (int r = 0; r < s.A.num_rows; ++r)
    for (int e = s.A.starts[r]; e < s.A.starts[r + 1]; ++e)
      s.A.values[e] *= R[r] * C[s.A.indices[e]];
  s.AT = transposeOf(s.A);
  s.c = p.costs(); s.lo = p.colLower(); s.up = p.colUpper(); s.b = p.rowLower();
  for (int j = 0; j < (int)s.c.size(); ++j) {
    s.c[j] *= C[j];
    if (std::fabs(s.lo[j]) < 1e20) s.lo[j] /= C[j];
    if (std::fabs(s.up[j]) < 1e20) s.up[j] /= C[j];
  }
  for (int i = 0; i < (int)s.b.size(); ++i)
    if (std::fabs(s.b[i]) < 1e20) s.b[i] *= R[i];
  return s;
}

void report(const char* label, const Scaled& s) {
  std::printf("\n=== %s ===\n", label);
  printDist("A", logDist(s.A.values));
  printDist("rowscale R", logDist(s.R));
  printDist("colscale C", logDist(s.C));
  printDist("c", logDist(s.c));
  printDist("b (rowlower)", logDist(s.b));
  std::vector<double> ub;
  for (double v : s.up) if (std::fabs(v) < 1e20) ub.push_back(v);
  for (double v : s.lo) if (std::fabs(v) < 1e20) ub.push_back(v);
  printDist("finite col bnds", logDist(ub));
  double nA = sigmaMax(s.A, s.AT);
  std::vector<double> bf;
  for (double v : s.b) if (std::fabs(v) < 1e20) bf.push_back(v);
  std::printf("  sigma_max(A)=%.4g   ||c||2=%.4g   ||b||2=%.4g   ||bnd||2=%.4g\n", nA, l2(s.c), l2(bf), l2(ub));
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return 2;
  double drop = argc > 2 ? std::atof(argv[2]) : 0.0;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.readModel(argv[1]);
  HighsLp lp = highs.getLp();
  if (drop > 0) {
    lp.a_matrix_.ensureColwise();
    HighsInt out = 0; long dropped = 0;
    std::vector<HighsInt> ns(lp.num_col_ + 1, 0);
    for (HighsInt c = 0; c < lp.num_col_; ++c) {
      for (HighsInt e = lp.a_matrix_.start_[c]; e < lp.a_matrix_.start_[c + 1]; ++e) {
        if (std::fabs(lp.a_matrix_.value_[e]) >= drop) { lp.a_matrix_.index_[out] = lp.a_matrix_.index_[e]; lp.a_matrix_.value_[out] = lp.a_matrix_.value_[e]; ++out; }
        else ++dropped;
      }
      ns[c + 1] = out;
    }
    lp.a_matrix_.start_ = ns; lp.a_matrix_.index_.resize(out); lp.a_matrix_.value_.resize(out);
    std::printf("dropped %ld entries below %.1e\n", dropped, drop);
  }
  std::printf("LP: %d rows %d cols %d nnz\n", (int)lp.num_row_, (int)lp.num_col_,
              (int)lp.a_matrix_.numNz());
  ProblemBuildResult built = buildProblem(lp);
  if (!built.ok()) { std::printf("buildProblem failed\n"); return 1; }
  const Problem& p = built.problem;
  std::printf("transformed: %d rows %d cols %d nnz\n", p.numRows(), p.numCols(), p.matrix().numNonzeros());

  {
    const SparseMatrix& M = p.matrix();
    long cnt[16] = {0};
    for (double v : M.values) { int d = (int)std::floor(std::log10(std::fabs(v))); d = std::max(-9, std::min(d, 1)); cnt[d + 9]++; }
    std::printf("decade histogram of |a| (fraction of nnz):\n");
    for (int d = -9; d <= 1; ++d) std::printf("  1e%-3d: %6.2f%%\n", d, 100.0 * cnt[d + 9] / M.numNonzeros());
  }
  std::vector<double> R1(p.numRows(), 1.0), C1(p.numCols(), 1.0);
  report("off (transformed problem)", applyScale(p, R1, C1));

  struct Cfg { const char* name; ScalingMode mode; int legacy_mask; int cr; int ruiz; };
  std::vector<Cfg> cfgs = {
      {"curtis_reid only (20 sweeps)", ScalingMode::kCurtisReid, 0, 20, 10},
      {"auto = CR + Ruiz10 + PC", ScalingMode::kAuto, kPdlpScalingRuiz + kPdlpScalingPC, 20, 10},
      {"legacy = Ruiz10 + PC", ScalingMode::kLegacy, kPdlpScalingRuiz + kPdlpScalingPC, 20, 10},
  };
  for (const Cfg& c : cfgs) {
    ScalingOptions o;
    o.mode = c.mode; o.legacy_scaling_mode = c.legacy_mask;
    o.curtis_reid_sweeps = c.cr; o.ruiz_iterations = c.ruiz;
    ScalingResult r = computeScalingPipeline(p.matrix(), o);
    if (!r.ok()) { std::printf("%s failed\n", c.name); continue; }
    report(c.name, applyScale(p, r.scaling.row_scale, r.scaling.col_scale));
  }
  return 0;
}
