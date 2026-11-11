// Absolute correctness tests of HiPDLP: LPs whose optimal primal-dual
// solution is known in advance, checked against SolveCpu (every build) and
// the GPU Solve (CUDA builds with a device).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

#include "HConfig.h"
#include "catch.hpp"
#include "pdlp/hipdlp/solver.h"

namespace {
const double kInf = std::numeric_limits<double>::infinity();

// xorshift64, the same generator as TestHiPdlpReference.cpp, so a seed means
// the same LP on every platform (std::uniform_real_distribution does not).
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed)
      : s(seed * 2685821657736338717ULL + 1442695040888963407ULL) {}
  uint64_t next() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  }
  double u() { return double(next() >> 11) * (1.0 / 9007199254740992.0); }
  double range(double a, double b) { return a + (b - a) * u(); }
  bool coin() { return u() < 0.5; }
};

// An LP together with its optimal solution, in the HiGHS sign convention of
// lp.sense (for a maximization the duals and reduced costs flip sign).
struct Planted {
  hipdlp::HostLp lp;
  std::vector<double> x, y, z;
  double objective = 0;
};

// min -x1 - x2  s.t.  x1 + 2 x2 <= 4,  3 x1 + x2 <= 6,  x >= 0.
// Both rows active, both columns interior: x* = (1.6, 1.2), y* = (-0.4, -0.2),
// z* = 0, objective -2.8.
hipdlp::HostLp handLp() {
  hipdlp::HostLp lp;
  lp.rows = 2;
  lp.columns = 2;
  lp.row_start = {0, 2, 4};
  lp.column_index = {0, 1, 0, 1};
  lp.value = {1, 2, 3, 1};
  lp.cost = {-1, -1};
  lp.row_lower = {-kInf, -kInf};
  lp.row_upper = {4, 6};
  lp.column_lower = {0, 0};
  lp.column_upper = {kInf, kInf};
  return lp;
}

// Builds min c'x s.t. row_lower <= Ax <= row_upper, column_lower <= x <=
// column_upper around a chosen optimal (x*, y*, z*), then sets c = A'y* + z*.
// Rows 0..k-1 are active with y != 0, columns 0..k-1 are strictly interior with
// z = 0, every other column sits on a bound with z != 0 and every other row is
// strictly inactive with y = 0. The k x k block A[0:k, 0:k] is strictly
// diagonally dominant, hence nonsingular, so x* and (y*, z*) are unique.
Planted plant(uint64_t seed, int m, int n, int k, int sense) {
  Rng r(seed);
  Planted p;
  hipdlp::HostLp& lp = p.lp;
  lp.rows = m;
  lp.columns = n;
  lp.sense = sense;
  lp.offset = r.range(-5, 5);

  // 1. Columns: bounds, x*, z* (minimization signs).
  std::vector<double> x(n), z(n, 0.0);
  lp.column_lower.resize(n);
  lp.column_upper.resize(n);
  for (int j = 0; j < n; ++j) {
    double l, u, t = r.u();
    if (j < k) {       // interior, z = 0
      if (t < 0.25) {  // boxed
        l = r.range(-3, 0);
        u = l + r.range(1, 4);
        x[j] = l + (u - l) * r.range(0.2, 0.8);
      } else if (t < 0.5) {  // lower bound only
        l = r.range(-2, 2);
        u = kInf;
        x[j] = l + r.range(0.5, 3);
      } else if (t < 0.75) {  // upper bound only
        l = -kInf;
        u = r.range(-2, 2);
        x[j] = u - r.range(0.5, 3);
      } else {  // free
        l = -kInf;
        u = kInf;
        x[j] = r.range(-3, 3);
      }
    } else if (t < 0.4) {  // at the lower bound, z > 0
      l = r.range(-2, 2);
      u = r.coin() ? l + r.range(1, 4) : kInf;
      x[j] = l;
      z[j] = r.range(0.5, 2);
    } else if (t < 0.8) {  // at the upper bound, z < 0
      u = r.range(-2, 2);
      l = r.coin() ? u - r.range(1, 4) : -kInf;
      x[j] = u;
      z[j] = -r.range(0.5, 2);
    } else {  // fixed, z of either sign
      l = r.range(-2, 2);
      u = l;
      x[j] = l;
      z[j] = r.coin() ? r.range(0.5, 2) : -r.range(0.5, 2);
    }
    lp.column_lower[j] = l;
    lp.column_upper[j] = u;
  }

  // 2. Matrix, row by row; the block rows get a dominant diagonal.
  const double density = 0.3;
  lp.row_start.push_back(0);
  for (int i = 0; i < m; ++i) {
    std::vector<double> row(n, 0.0);
    for (int j = 0; j < n; ++j)
      if (r.u() < density) row[j] = (r.coin() ? 1 : -1) * r.range(0.5, 2);
    if (i < k) {
      double off = 0;
      for (int j = 0; j < k; ++j) {
        if (j != i) off += std::abs(row[j]);
      }
      row[i] = (r.coin() ? 1 : -1) * (off + r.range(1, 2));
    }
    for (int j = 0; j < n; ++j) {
      if (row[j] != 0.0) {
        lp.column_index.push_back(j);
        lp.value.push_back(row[j]);
      }
    }
    lp.row_start.push_back(int(lp.value.size()));
  }

  // 3. Row activities A x*, then row bounds and y* (minimization signs).
  std::vector<double> y(m, 0.0);
  lp.row_lower.resize(m);
  lp.row_upper.resize(m);
  for (int i = 0; i < m; ++i) {
    double ax = 0;
    for (int e = lp.row_start[i]; e < lp.row_start[i + 1]; ++e) {
      ax += lp.value[e] * x[lp.column_index[e]];
    }
    double t = r.u();
    if (i < k) {      // active, y != 0
      if (t < 0.4) {  // at the lower bound
        lp.row_lower[i] = ax;
        lp.row_upper[i] = r.coin() ? ax + r.range(1, 3) : kInf;
        y[i] = r.range(0.5, 2);
      } else if (t < 0.8) {  // at the upper bound
        lp.row_lower[i] = r.coin() ? ax - r.range(1, 3) : -kInf;
        lp.row_upper[i] = ax;
        y[i] = -r.range(0.5, 2);
      } else {  // equality
        lp.row_lower[i] = ax;
        lp.row_upper[i] = ax;
        y[i] = r.coin() ? r.range(0.5, 2) : -r.range(0.5, 2);
      }
    } else if (t < 0.3) {  // strictly inactive from here on, y = 0: boxed
      y[i] = 0;
      lp.row_lower[i] = ax - r.range(0.5, 2);
      lp.row_upper[i] = ax + r.range(0.5, 2);
    } else if (t < 0.6) {  // lower bound only
      y[i] = 0;
      lp.row_lower[i] = ax - r.range(0.5, 2);
      lp.row_upper[i] = kInf;
    } else if (t < 0.9) {  // upper bound only
      y[i] = 0;
      lp.row_lower[i] = -kInf;
      lp.row_upper[i] = ax + r.range(0.5, 2);
    } else {  // free
      y[i] = 0;
      lp.row_lower[i] = -kInf;
      lp.row_upper[i] = kInf;
    }
  }

  // 4. c = A'y* + z*, then the sense: a maximization of -c has the same x*,
  // and HiGHS reports its duals with the opposite sign.
  lp.cost = z;
  for (int i = 0; i < m; ++i) {
    for (int e = lp.row_start[i]; e < lp.row_start[i + 1]; ++e) {
      int j = lp.column_index[e];
      double a = lp.value[e];
      lp.cost[j] += a * y[i];
    }
  }

  for (int j = 0; j < n; ++j) {
    lp.cost[j] *= sense;
    z[j] *= sense;
  }
  for (int i = 0; i < m; ++i) y[i] *= sense;
  p.x = x;
  p.y = y;
  p.z = z;
  p.objective = lp.offset + std::inner_product(lp.cost.begin(), lp.cost.end(),
                                               x.begin(), 0.0);
  return p;
}

// Every solver this build has: SolveCpu always, the GPU Solve when CUDA was
// detected at configure time and a usable device answers the hand LP.
typedef hipdlp::SolveResult (*SolveFunction)(const hipdlp::HostLp&,
                                             const hipdlp::SolveOptions&);
struct Backend {
  const char* name;
  SolveFunction solve;
};

std::vector<Backend> backends() {
  std::vector<Backend> list;
  list.push_back({"CPU", hipdlp::SolveCpu});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (hipdlp::Solve(handLp()).status == hipdlp::Status::kOptimal) {
    list.push_back({"GPU", hipdlp::Solve});
  }
#endif
  return list;
}

// The solver's answer against the planted one: every component of x, y and z
// within `tol`, the objective within `tol_objective` relative.
void checkSolution(const Planted& p, const hipdlp::SolveResult& r, double tol,
                   double tol_objective) {
  REQUIRE(r.status == hipdlp::Status::kOptimal);
  for (int j = 0; j < p.lp.columns; ++j) {
    CHECK(r.primal[j] == Approx(p.x[j]).margin(tol));
    CHECK(r.reduced_cost[j] == Approx(p.z[j]).margin(tol));
  }
  for (int i = 0; i < p.lp.rows; ++i) {
    CHECK(r.dual[i] == Approx(p.y[i]).margin(tol));
  }

  CHECK(
      r.metrics.primal_objective ==
      Approx(p.objective).margin(tol_objective * (1 + std::fabs(p.objective))));
}

// The three PDLP termination measures of the point a solver returned,
// recomputed from their definition in solver.h: minimization sense, no offset,
// z nonzero only on a bound x is within sqrt(tolerance) of.
struct L2 {
  double primal, dual, gap;
};

L2 recomputeL2(const hipdlp::HostLp& lp, const hipdlp::SolveResult& r,
               double tolerance) {
  // y in minimization signs, A x and A'y.
  const double sense = lp.sense;
  std::vector<double> y(lp.rows), ax(lp.rows, 0.0), aty(lp.columns, 0.0);
  for (int i = 0; i < lp.rows; ++i) {
    y[i] = sense * r.dual[i];
    for (int e = lp.row_start[i]; e < lp.row_start[i + 1]; ++e) {
      int j = lp.column_index[e];
      double a = lp.value[e];
      ax[i] += a * r.primal[j];
      aty[j] += a * y[i];
    }
  }

  const double near = std::sqrt(tolerance);
  double violation_sq = 0, b_sq = 0, residual_sq = 0, c_sq = 0, p = 0, d = 0;

  // Columns: bound violation, the dual residual left after the part of
  // c - A'y a nearby bound absorbs, and z'x for the dual objective.
  for (int j = 0; j < lp.columns; ++j) {
    const double c = sense * lp.cost[j], x = r.primal[j];
    const double l = lp.column_lower[j], u = lp.column_upper[j];
    const double v = std::max(std::max(l - x, x - u), 0.0);

    violation_sq += v * v;
    const double raw = c - aty[j];
    const bool at_lower = std::isfinite(l) && std::fabs(x - l) <= near;
    const bool at_upper = std::isfinite(u) && std::fabs(x - u) <= near;
    double z = 0;
    if (at_lower && at_upper)
      z = raw;
    else if (at_lower)
      z = std::max(raw, 0.0);
    else if (at_upper)
      z = std::min(raw, 0.0);
    residual_sq += (raw - z) * (raw - z);
    c_sq += c * c;
    p += c * x;
    d += z * x;
  }

  // Rows: bound violation, ||b|| (an equality bound counted once) and each
  // row's bound by the sign of its multiplier for the dual objective.
  for (int i = 0; i < lp.rows; ++i) {
    const double l = lp.row_lower[i], u = lp.row_upper[i];
    violation_sq += std::max(std::max(l - ax[i], ax[i] - u), 0.0) *
                    std::max(std::max(l - ax[i], ax[i] - u), 0.0);
    if (std::isfinite(l) && l != u) b_sq += l * l;
    if (std::isfinite(u)) b_sq += u * u;
    if (y[i] > 0) d += y[i] * l;
    if (y[i] < 0) d += y[i] * u;
  }

  L2 result;
  result.primal = std::sqrt(violation_sq) / (1 + std::sqrt(b_sq));
  result.dual = std::sqrt(residual_sq) / (1 + std::sqrt(c_sq));
  result.gap = std::fabs(p - d) / (1 + std::fabs(p) + std::fabs(d));
  return result;
}

// Each reported measure against the recomputed one (the GPU sums in another
// order: measured differences <= 3e-15 absolute), and OPTIMAL must mean all
// three are within the tolerance.
void checkMetrics(const Planted& p, const hipdlp::SolveResult& r,
                  double tolerance) {
  const L2 mine = recomputeL2(p.lp, r, tolerance);
  CHECK(r.metrics.l2_primal == Approx(mine.primal).epsilon(1e-9).margin(1e-13));
  CHECK(r.metrics.l2_dual == Approx(mine.dual).epsilon(1e-9).margin(1e-13));
  CHECK(r.metrics.l2_gap == Approx(mine.gap).epsilon(1e-9).margin(1e-13));
  if (r.status == hipdlp::Status::kOptimal) {
    CHECK(mine.primal <= tolerance);
    CHECK(mine.dual <= tolerance);
    CHECK(mine.gap <= tolerance);
  }
}

}  // namespace

// The hand-solved LP: signs of duals and objective for min, max and an offset.
TEST_CASE("hipdlp-planted-hand", "[hipdlp][hipdlp-planted]") {
  hipdlp::SolveOptions options;
  options.tolerance = 1e-8;

  SECTION("minimize") {
    const hipdlp::SolveResult r = hipdlp::SolveCpu(handLp(), options);
    REQUIRE(r.status == hipdlp::Status::kOptimal);
    CHECK(r.primal[0] == Approx(1.6).margin(1e-6));
    CHECK(r.primal[1] == Approx(1.2).margin(1e-6));
    CHECK(r.dual[0] == Approx(-0.4).margin(1e-6));
    CHECK(r.dual[1] == Approx(-0.2).margin(1e-6));
    CHECK(r.reduced_cost[0] == Approx(0).margin(1e-6));
    CHECK(r.reduced_cost[1] == Approx(0).margin(1e-6));
    CHECK(r.metrics.primal_objective == Approx(-2.8).margin(1e-6));
  }

  // max x1 + x2 has the same x*; the duals and the objective flip sign.
  SECTION("maximize") {
    hipdlp::HostLp lp = handLp();
    lp.sense = -1;
    lp.cost = {1, 1};
    const hipdlp::SolveResult r = hipdlp::SolveCpu(lp, options);
    REQUIRE(r.status == hipdlp::Status::kOptimal);
    CHECK(r.primal[0] == Approx(1.6).margin(1e-6));
    CHECK(r.primal[1] == Approx(1.2).margin(1e-6));
    CHECK(r.dual[0] == Approx(0.4).margin(1e-6));
    CHECK(r.dual[1] == Approx(0.2).margin(1e-6));
    CHECK(r.reduced_cost[0] == Approx(0).margin(1e-6));
    CHECK(r.reduced_cost[1] == Approx(0).margin(1e-6));
    CHECK(r.metrics.primal_objective == Approx(2.8).margin(1e-6));
  }
  SECTION("offset") {
    hipdlp::HostLp lp = handLp();
    lp.offset = 10;
    const hipdlp::SolveResult r = hipdlp::SolveCpu(lp, options);
    REQUIRE(r.status == hipdlp::Status::kOptimal);
    CHECK(r.metrics.primal_objective == Approx(-2.8 + 10).margin(1e-6));
  }
}

// The generator itself: the planted point satisfies the KKT conditions.
TEST_CASE("hipdlp-planted-generator", "[hipdlp][hipdlp-planted]") {
  for (uint64_t seed = 1; seed <= 20; ++seed)
    for (int sense : {1, -1}) {
      const Planted p = plant(seed, 30, 40, 12, sense);
      const hipdlp::HostLp& lp = p.lp;
      // residual = c - A'y - z and ax = A x.
      std::vector<double> residual = lp.cost;
      for (int j = 0; j < lp.columns; ++j) residual[j] -= p.z[j];
      std::vector<double> ax(lp.rows, 0.0);
      for (int i = 0; i < lp.rows; ++i)
        for (int e = lp.row_start[i]; e < lp.row_start[i + 1]; ++e) {
          const int j = lp.column_index[e];
          residual[j] -= lp.value[e] * p.y[i];
          ax[i] += lp.value[e] * p.x[j];
        }

      // Stationarity and primal feasibility (rows up to rounding).
      for (int j = 0; j < lp.columns; ++j) {
        CHECK(p.x[j] >= lp.column_lower[j] - 1e-12);
        CHECK(p.x[j] <= lp.column_upper[j] + 1e-12);
        CHECK(std::fabs(residual[j]) <= 1e-12);
      }
      for (int i = 0; i < lp.rows; ++i) {
        CHECK(ax[i] >= lp.row_lower[i] - 1e-12);
        CHECK(ax[i] <= lp.row_upper[i] + 1e-12);
      }
      // Complementarity in minimization signs: a positive multiplier sits on
      // the lower bound, a negative one on the upper bound.
      for (int j = 0; j < lp.columns; ++j) {
        const double zmin = p.z[j] * sense;
        if (zmin > 0) CHECK(p.x[j] == lp.column_lower[j]);
        if (zmin < 0) CHECK(p.x[j] == lp.column_upper[j]);
      }
      for (int i = 0; i < lp.rows; ++i) {
        const double ymin = p.y[i] * sense;
        if (ymin > 0) CHECK(ax[i] == Approx(lp.row_lower[i]).margin(1e-12));
        if (ymin < 0) CHECK(ax[i] == Approx(lp.row_upper[i]).margin(1e-12));
      }
    }
}

TEST_CASE("hipdlp-planted-first-solve", "[hipdlp][hipdlp-planted]") {
  const Planted p = plant(1, 30, 40, 12, 1);
  hipdlp::SolveOptions options;
  options.tolerance = 1e-8;
  const hipdlp::SolveResult r = hipdlp::SolveCpu(p.lp, options);
  REQUIRE(r.status == hipdlp::Status::kOptimal);
  for (int j = 0; j < p.lp.columns; ++j) {
    CHECK(r.primal[j] == Approx(p.x[j]).margin(1e-5));
  }
  CHECK(r.metrics.primal_objective == Approx(p.objective).margin(1e-5));
}

// Shapes: wide, tall, square with every row active, no active row at all,
// tall with every column interior. Tolerances measured on SolveCpu over these
// LPs (2026-11-09): stop 1e-8 -> |error| <= 2.2e-7, objective <= 2.5e-8
// relative; stop 1e-4 -> |error| <= 2.6e-3, objective <= 3.9e-4; asserted with
// about 10x margin.
TEST_CASE("hipdlp-planted-solve", "[hipdlp][hipdlp-planted]") {
  const int shapes[][3] = {
      {30, 40, 12}, {40, 30, 12}, {10, 10, 10}, {20, 30, 0}, {50, 20, 20}};
  for (const Backend& backend : backends())
    for (const auto& shape : shapes)
      for (uint64_t seed = 1; seed <= 20; ++seed)
        for (int sense : {1, -1}) {
          const Planted p = plant(seed, shape[0], shape[1], shape[2], sense);
          INFO(backend.name << " m " << shape[0] << " n " << shape[1] << " k "
                            << shape[2] << " seed " << seed << " sense "
                            << sense);
          hipdlp::SolveOptions options;
          options.iteration_limit = 100000;  // a bug must fail, not hang
          options.tolerance = 1e-8;
          checkSolution(p, backend.solve(p.lp, options), 1e-5, 1e-6);
          options.tolerance = 1e-4;
          checkSolution(p, backend.solve(p.lp, options), 2e-2, 3e-3);
        }
}

// The measures are a formula of the returned point, so check them on
// unfinished iterates too (iteration limits 1, 7, 64) where many columns sit
// near, but not at, a bound: that is where a wrong near-bound test shows.
TEST_CASE("hipdlp-planted-metrics", "[hipdlp][hipdlp-planted]") {
  const int shapes[][3] = {
      {30, 40, 12}, {40, 30, 12}, {10, 10, 10}, {20, 30, 0}, {50, 20, 20}};
  for (const Backend& backend : backends())
    for (const auto& shape : shapes)
      for (uint64_t seed = 1; seed <= 5; ++seed)
        for (int sense : {1, -1})
          for (double tolerance : {1e-4, 1e-8})
            for (int64_t limit : {1, 7, 64, 100000}) {
              const Planted p =
                  plant(seed, shape[0], shape[1], shape[2], sense);
              INFO(backend.name << " m " << shape[0] << " n " << shape[1]
                                << " k " << shape[2] << " seed " << seed
                                << " sense " << sense << " tolerance "
                                << tolerance << " limit " << limit);
              hipdlp::SolveOptions options;
              options.tolerance = tolerance;
              options.iteration_limit = limit;
              checkMetrics(p, backend.solve(p.lp, options), tolerance);
            }
}
