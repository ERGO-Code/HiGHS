// Per-(instance, scaling config) feature extraction using the rewrite's own scaling code.
// usage: feat INSTANCE.mps[.gz] NAME  -> CSV rows on stdout (no header; see kHeader)
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

static const char* kHeader =
    "instance,config,rows,cols,nnz,a_lmin,a_lp1,a_lmed,a_lp99,a_lmax,"
    "sig,amax,amin_l,"
    "Rmin_l,Rp1_l,Rp99_l,Rmax_l,Cmin_l,Cp1_l,Cp99_l,Cmax_l,"
    "ampR_mean,ampR_max,ampC_mean,ampC_max,"
    "cn,bn,bndn,cn0,bn0,bndn0,sig0";

static double l2(const std::vector<double>& v) {
  long double s = 0;
  for (double x : v) if (std::isfinite(x) && std::fabs(x) < 1e20) s += (long double)x * x;
  return std::sqrt((double)s);
}
static double pct(std::vector<double> v, double f) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[std::min(v.size() - 1, (size_t)(f * (v.size() - 1)))];
}
static SparseMatrix transposeOf(const SparseMatrix& A) {
  SparseMatrix T; T.num_rows = A.num_cols; T.num_cols = A.num_rows;
  T.starts.assign(T.num_rows + 1, 0);
  for (int e = 0; e < A.numNonzeros(); ++e) T.starts[A.indices[e] + 1]++;
  for (int i = 0; i < T.num_rows; ++i) T.starts[i + 1] += T.starts[i];
  T.indices.resize(A.numNonzeros()); T.values.resize(A.numNonzeros());
  std::vector<int> pos(T.starts.begin(), T.starts.end() - 1);
  for (int r = 0; r < A.num_rows; ++r)
    for (int e = A.starts[r]; e < A.starts[r + 1]; ++e) {
      int p = pos[A.indices[e]]++; T.indices[p] = r; T.values[p] = A.values[e];
    }
  return T;
}
static double sigmaMax(const SparseMatrix& A, const SparseMatrix& AT, int iters = 80) {
  if (A.num_cols == 0 || A.num_rows == 0 || A.numNonzeros() == 0) return 0;
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
    nrm = l2(x); est = std::sqrt(nrm);
    if (nrm == 0) break;
    for (auto& v : x) v /= nrm;
  }
  return est;
}
static double lg(double x) { return x > 0 ? std::log10(x) : -99; }

int main(int argc, char** argv) {
  if (argc < 3) return 2;
  const std::string name = argv[2];
  if (argc > 3 && std::string(argv[3]) == "--header") { std::puts(kHeader); return 0; }
  Highs highs; highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) { std::fprintf(stderr, "read fail\n"); return 3; }
  ProblemBuildResult built = buildProblem(highs.getLp());
  if (!built.ok()) { std::fprintf(stderr, "build fail\n"); return 4; }
  const Problem& p = built.problem;
  const int m = p.numRows(), n = p.numCols();

  struct Cfg { const char* name; ScalingMode mode; int mask; };
  const Cfg cfgs[] = {
      {"A", ScalingMode::kAuto, kPdlpScalingRuiz + kPdlpScalingPC},
      {"L", ScalingMode::kLegacy, kPdlpScalingRuiz + kPdlpScalingPC},
      {"P", ScalingMode::kLegacy, kPdlpScalingPC},
      {"R", ScalingMode::kLegacy, kPdlpScalingRuiz},
      {"C", ScalingMode::kCurtisReid, 0},
      {"O", ScalingMode::kOff, 0},
      {"AP", ScalingMode::kAuto, kPdlpScalingPC},
      {"AR", ScalingMode::kAuto, kPdlpScalingRuiz},
  };

  // unscaled reference norms
  double cn0 = l2(p.costs());
  std::vector<double> b0, bnd0;
  for (double v : p.rowLower()) if (std::fabs(v) < 1e20) b0.push_back(v);
  for (double v : p.colLower()) if (std::fabs(v) < 1e20) bnd0.push_back(v);
  for (double v : p.colUpper()) if (std::fabs(v) < 1e20) bnd0.push_back(v);
  double bn0 = l2(b0), bndn0 = l2(bnd0);
  double sig0 = sigmaMax(p.matrix(), transposeOf(p.matrix()));

  for (const Cfg& c : cfgs) {
    std::vector<double> R(m, 1.0), C(n, 1.0);
    if (c.mode != ScalingMode::kOff) {
      ScalingOptions o; o.mode = c.mode; o.legacy_scaling_mode = c.mask;
      o.curtis_reid_sweeps = 20; o.ruiz_iterations = 10;
      ScalingResult r = computeScalingPipeline(p.matrix(), o);
      if (!r.ok()) { std::fprintf(stderr, "scaling fail %s\n", c.name); continue; }
      R = r.scaling.row_scale; C = r.scaling.col_scale;
    }
    SparseMatrix A = p.matrix();
    std::vector<double> absv;
    absv.reserve(A.values.size());
    double amax = 0, amin = 1e300;
    for (int r = 0; r < m; ++r)
      for (int e = A.starts[r]; e < A.starts[r + 1]; ++e) {
        A.values[e] *= R[r] * C[A.indices[e]];
        double a = std::fabs(A.values[e]); amax = std::max(amax, a); amin = std::min(amin, a);
        absv.push_back(std::log10(a));
      }
    SparseMatrix AT = transposeOf(A);
    double sig = sigmaMax(A, AT);
    std::vector<double> lR(m), lC(n);
    double ampR_mean = 0, ampR_max = 0, ampC_mean = 0, ampC_max = 0;
    for (int i = 0; i < m; ++i) { lR[i] = lg(R[i]); ampR_mean += 1.0 / R[i]; ampR_max = std::max(ampR_max, 1.0 / R[i]); }
    for (int j = 0; j < n; ++j) { lC[j] = lg(C[j]); ampC_mean += 1.0 / C[j]; ampC_max = std::max(ampC_max, 1.0 / C[j]); }
    if (m) ampR_mean /= m;
    if (n) ampC_mean /= n;
    std::vector<double> cs(n), bs, bnds;
    for (int j = 0; j < n; ++j) {
      cs[j] = p.costs()[j] * C[j];
      double lo = p.colLower()[j], up = p.colUpper()[j];
      if (std::fabs(lo) < 1e20) bnds.push_back(lo / C[j]);
      if (std::fabs(up) < 1e20) bnds.push_back(up / C[j]);
    }
    for (int i = 0; i < m; ++i) { double v = p.rowLower()[i]; if (std::fabs(v) < 1e20) bs.push_back(v * R[i]); }
    std::printf("%s,%s,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.6g,%.6g,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
                name.c_str(), c.name, m, n, A.numNonzeros(),
                pct(absv, 0), pct(absv, 0.01), pct(absv, 0.5), pct(absv, 0.99), pct(absv, 1),
                sig, amax, lg(amin),
                pct(lR, 0), pct(lR, 0.01), pct(lR, 0.99), pct(lR, 1),
                pct(lC, 0), pct(lC, 0.01), pct(lC, 0.99), pct(lC, 1),
                ampR_mean, ampR_max, ampC_mean, ampC_max,
                l2(cs), l2(bs), l2(bnds), cn0, bn0, bndn0, sig0);
  }
  return 0;
}
