// Differential test of the GPU HiPDLP solver against its CPU reference.
//
// Both solvers run the same reflected Halpern PDHG with the same restart rule,
// so after a fixed number of iterations they must be in the same state: the
// discrete quantities (iterations, restarts, restart decisions) equal, the
// vectors equal up to the rounding of a different summation order. A wrong
// formula, a wrong constant in the restart rule or a wrong sign in the primal
// weight update shows up many orders of magnitude above that tolerance.
//
// The CPU reference lives in highs/pdlp/hipdlp/solver.cc, the test hooks of
// both solvers in highs/pdlp/hipdlp/solver_test.h.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "catch.hpp"

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/solver_test.h"

namespace {

const bool dev_run = false;

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
};

// min c'x s.t. lr <= Ax <= ur, lx <= x <= ux, built around a chosen
// primal-dual optimal pair (x*, y*, z*) with c = A'y* + z*, so a finite optimum
// exists. sense = -1 stores -c and maximizes.
hipdlp::HostLp makeLp(uint64_t seed, int m, int n, double density, int sense) {
  Rng r(seed);
  const double inf = std::numeric_limits<double>::infinity();
  hipdlp::HostLp lp;
  lp.rows = m;
  lp.columns = n;
  lp.sense = sense;
  lp.offset = 3.0;
  lp.row_start.push_back(0);
  for (int i = 0; i < m; ++i) {
    size_t before = lp.column_index.size();
    for (int j = 0; j < n; ++j)
      if (r.u() < density) {
        lp.column_index.push_back(j);
        lp.value.push_back((r.u() < 0.5 ? -1 : 1) *
                           std::pow(10.0, r.range(-1.5, 1.5)));
      }
    if (lp.column_index.size() == before) {
      lp.column_index.push_back(int(r.u() * n) % n);
      lp.value.push_back(r.range(0.5, 2));
    }
    lp.row_start.push_back(int(lp.column_index.size()));
  }
  std::vector<double> xs(n), zs(n, 0.0);
  lp.column_lower.resize(n);
  lp.column_upper.resize(n);
  for (int j = 0; j < n; ++j) {
    double t = r.u(), l, u;
    if (t < 0.30) {
      l = r.range(-3, 1);
      u = l + r.range(0.5, 5);
    } else if (t < 0.65) {
      l = r.range(-2, 2);
      u = inf;
    } else if (t < 0.75) {
      l = -inf;
      u = r.range(-2, 2);
    } else if (t < 0.90) {
      l = -inf;
      u = inf;
    } else {
      l = u = r.range(-2, 2);
    }
    double p = r.u(), x, z = 0;
    if (l == u) {
      x = l;
      z = r.range(-2, 2);
    } else if (std::isfinite(l) && p < 0.4) {
      x = l;
      z = r.u() < 0.7 ? r.range(0, 2) : 0;
    } else if (std::isfinite(u) && p < 0.6) {
      x = u;
      z = r.u() < 0.7 ? -r.range(0, 2) : 0;
    } else if (std::isfinite(l) && std::isfinite(u)) {
      x = r.range(l, u);
    } else if (std::isfinite(l)) {
      x = l + r.range(0, 3);
    } else if (std::isfinite(u)) {
      x = u - r.range(0, 3);
    } else {
      x = r.range(-2, 2);
    }
    lp.column_lower[j] = l;
    lp.column_upper[j] = u;
    xs[j] = x;
    zs[j] = z;
  }
  std::vector<double> ax(m, 0.0), ys(m, 0.0);
  for (int i = 0; i < m; ++i)
    for (int k = lp.row_start[i]; k < lp.row_start[i + 1]; ++k)
      ax[i] += lp.value[k] * xs[lp.column_index[k]];
  lp.row_lower.resize(m);
  lp.row_upper.resize(m);
  for (int i = 0; i < m; ++i) {
    double t = r.u(), y = 0;
    if (t < 0.25) {
      lp.row_lower[i] = lp.row_upper[i] = ax[i];
      y = r.range(-2, 2);
    } else if (t < 0.40) {
      lp.row_lower[i] = ax[i] - r.range(0.1, 2);
      lp.row_upper[i] = ax[i] + r.range(0.1, 2);
    } else if (t < 0.65) {
      double s = r.u() < 0.5 ? 0 : r.range(0.1, 2);
      lp.row_lower[i] = ax[i] - s;
      lp.row_upper[i] = inf;
      if (s == 0 && r.u() < 0.7) y = r.range(0, 2);
    } else if (t < 0.90) {
      double s = r.u() < 0.5 ? 0 : r.range(0.1, 2);
      lp.row_lower[i] = -inf;
      lp.row_upper[i] = ax[i] + s;
      if (s == 0 && r.u() < 0.7) y = -r.range(0, 2);
    } else {
      lp.row_lower[i] = -inf;
      lp.row_upper[i] = inf;
    }
    ys[i] = y;
  }
  std::vector<double> c = zs;
  for (int i = 0; i < m; ++i)
    for (int k = lp.row_start[i]; k < lp.row_start[i + 1]; ++k)
      c[lp.column_index[k]] += lp.value[k] * ys[i];
  lp.cost.resize(n);
  for (int j = 0; j < n; ++j) lp.cost[j] = c[j] * sense;
  return lp;
}

// The model of a HiGHS instance as the solvers take it: row-wise CSR, plain
// ints, and IEEE infinities for absent bounds.
hipdlp::HostLp readHostLp(const std::string& instance) {
  Highs h;
  h.setOptionValue("output_flag", false);
  const std::string file =
      std::string(HIGHS_DIR) + "/check/instances/" + instance + ".mps";
  REQUIRE(h.readModel(file) == HighsStatus::kOk);
  const HighsLp& model = h.getLp();
  HighsSparseMatrix matrix = model.a_matrix_;
  matrix.ensureRowwise();
  const double inf = std::numeric_limits<double>::infinity();
  auto bound = [&](double v) {
    return v >= kHighsInf ? inf : (v <= -kHighsInf ? -inf : v);
  };
  hipdlp::HostLp lp;
  lp.rows = int(model.num_row_);
  lp.columns = int(model.num_col_);
  lp.sense = int(model.sense_);
  lp.offset = model.offset_;
  for (HighsInt i = 0; i <= model.num_row_; ++i)
    lp.row_start.push_back(int(matrix.start_[i]));
  const HighsInt nnz = matrix.start_[model.num_row_];
  for (HighsInt k = 0; k < nnz; ++k) {
    lp.column_index.push_back(int(matrix.index_[k]));
    lp.value.push_back(matrix.value_[k]);
  }
  lp.cost = model.col_cost_;
  for (double v : model.row_lower_) lp.row_lower.push_back(bound(v));
  for (double v : model.row_upper_) lp.row_upper.push_back(bound(v));
  for (double v : model.col_lower_) lp.column_lower.push_back(bound(v));
  for (double v : model.col_upper_) lp.column_upper.push_back(bound(v));
  return lp;
}

// max_i |a_i - b_i| / (1 + |b_i|); NaN propagates, a size mismatch is infinite.
double diff(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
  double worst = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double d = std::fabs(a[i] - b[i]) / (1 + std::fabs(b[i]));
    if (!(d <= worst)) worst = d;
  }
  return worst;
}
double diff(double a, double b) {
  return std::fabs(a - b) / (1 + std::fabs(b));
}

// How far the CPU may drift from the GPU after `limit` iterations: the two sum
// in a different order (SpMV, norms). Measured on random LPs and netlib LPs:
// N=1 3e-15, N=40 7e-13, N=200 2e-12, N=1000 2e-11.
double stateTolerance(int64_t limit) {
  return limit <= 3      ? 1e-11
         : limit <= 200  ? 1e-10
         : limit <= 1000 ? 1e-9
                         : 1e-8;
}

// "" when the two decision logs agree in every discrete field and in omega and
// the FPEs to 1e-6, else what differs first.
std::string compareDecisions(const std::vector<hipdlp::Decision>& gpu,
                             const std::vector<hipdlp::Decision>& cpu) {
  std::ostringstream out;
  if (gpu.size() != cpu.size()) {
    out << "decision count " << gpu.size() << " (GPU) vs " << cpu.size()
        << " (CPU)";
    return out.str();
  }
  for (size_t i = 0; i < gpu.size(); ++i) {
    const hipdlp::Decision &g = gpu[i], &c = cpu[i];
    if (g.iteration != c.iteration || g.epoch != c.epoch ||
        g.status != c.status || g.restart != c.restart ||
        g.balance_epoch_limit != c.balance_epoch_limit) {
      out << "decision " << i << " at iteration " << g.iteration << ": restart "
          << g.restart << " (GPU) vs " << c.restart << " (CPU), epoch "
          << g.epoch << " vs " << c.epoch << ", balance_epoch_limit "
          << g.balance_epoch_limit << " vs " << c.balance_epoch_limit;
      return out.str();
    }
    const double worst = std::max({diff(g.omega, c.omega), diff(g.fpe, c.fpe),
                                   diff(g.reference_fpe, c.reference_fpe)});
    if (!(worst <= 1e-6)) {
      out << "decision " << i << " at iteration " << g.iteration
          << ": omega/fpe differ by " << worst;
      return out.str();
    }
  }
  return "";
}

struct Tally {
  int checked = 0;
  int skipped = 0;
};

// Runs both solvers for each iteration limit and both weight modes, and checks
// that their final states agree.
void compareSolvers(const std::string& name, const hipdlp::HostLp& lp,
                    Tally& tally) {
  const int64_t limits[] = {0, 1, 2, 3, 39, 40, 41, 80, 81, 200, 400, 1000};
  for (const bool fixed_weight : {true, false}) {
    for (const int64_t limit : limits) {
      std::ostringstream label;
      label << name << ", primal weight "
            << (fixed_weight ? "fixed" : "updated") << ", " << limit
            << " iterations";
      INFO(label.str());

      hipdlp::SolveOptions options;
      options.tolerance = 1e-15;  // keep the GPU from stopping on optimality
      options.iteration_limit = limit;
      hipdlp::Snapshot gpu, cpu;
      const hipdlp::SolveResult rg =
          hipdlp::SolveForTesting(lp, options, &gpu, fixed_weight);
      if (limit > 0 && rg.status != hipdlp::Status::kIterationLimit) {
        // The GPU converged to machine precision before the limit; the CPU
        // reference has no termination test to match.
        ++tally.skipped;
        continue;
      }
      if (limit > 0 && gpu.fpe < 1e-12) {
        // z - T(z) is rounding noise: the weight update divides two
        // noise-sized movements, so omega (and every restart decision taken on
        // it) is not comparable any more.
        ++tally.skipped;
        continue;
      }
      const hipdlp::SolveResult rc =
          hipdlp::SolveCpuForTesting(lp, options, &cpu, fixed_weight);
      ++tally.checked;

      const double tol = stateTolerance(limit);
      const double scaling =
          std::max({diff(gpu.row_multiplier, cpu.row_multiplier),
                    diff(gpu.column_multiplier, cpu.column_multiplier),
                    diff(gpu.scaled_value, cpu.scaled_value)});
      const double state = std::max(
          {diff(gpu.x, cpu.x), diff(gpu.y, cpu.y), diff(gpu.xp, cpu.xp),
           diff(gpu.yp, cpu.yp), diff(gpu.xa, cpu.xa), diff(gpu.ya, cpu.ya)});
      const double epoch_move = std::max(diff(gpu.epoch_dx, cpu.epoch_dx),
                                         diff(gpu.epoch_dy, cpu.epoch_dy));
      const std::string decisions =
          compareDecisions(gpu.decisions, cpu.decisions);
      if (dev_run)
        printf("%s: scaling %.1e state %.1e epoch %.1e omega %.1e fpe %.1e\n",
               label.str().c_str(), scaling, state, epoch_move,
               diff(gpu.omega, cpu.omega), diff(gpu.fpe, cpu.fpe));

      CHECK(scaling <= 1e-14);
      CHECK(state <= tol);
      CHECK(epoch_move <= 10 * tol);
      CHECK(rc.status == rg.status);
      CHECK(rc.iterations == rg.iterations);
      CHECK(rc.restarts == rg.restarts);
      CHECK(cpu.epoch_iteration == gpu.epoch_iteration);
      CHECK(diff(gpu.omega, cpu.omega) <= 1e-6);
      CHECK(diff(gpu.fpe, cpu.fpe) <= 1e-9);
      CHECK(decisions == "");
    }
  }
}

// Runs Solve (GPU) and SolveCpu to termination at two tolerances and checks
// that they stop at the same iteration with the same solution and measures.
void compareToTermination(const std::string& name, const hipdlp::HostLp& lp,
                          Tally& tally) {
  for (const double tolerance : {1e-4, 1e-6}) {
    std::ostringstream label;
    label << name << ", tolerance " << tolerance;
    INFO(label.str());

    hipdlp::SolveOptions options;
    options.tolerance = tolerance;
    options.iteration_limit = 200000;
    const hipdlp::SolveResult gpu = hipdlp::Solve(lp, options);
    const hipdlp::SolveResult cpu = hipdlp::SolveCpu(lp, options);
    ++tally.checked;

    CHECK(cpu.status == gpu.status);
    CHECK(cpu.iterations == gpu.iterations);
    CHECK(cpu.restarts == gpu.restarts);
    CHECK(cpu.kkt_checks == gpu.kkt_checks);
    CHECK(cpu.certificate_checks == gpu.certificate_checks);
    CHECK(cpu.certificate_valid == gpu.certificate_valid);
    CHECK(cpu.ray_valid == gpu.ray_valid);
    CHECK(cpu.primal_valid == gpu.primal_valid);
    CHECK(cpu.dual_valid == gpu.dual_valid);
    CHECK(cpu.metrics_valid == gpu.metrics_valid);
    CHECK(cpu.rows == gpu.rows);
    CHECK(cpu.columns == gpu.columns);
    CHECK(cpu.nonzeros == gpu.nonzeros);
    // The solutions agree to the rounding of a different summation order,
    // amplified over up to tens of thousands of iterations.
    CHECK(diff(gpu.primal, cpu.primal) <= 1e-6);
    CHECK(diff(gpu.dual, cpu.dual) <= 1e-6);
    CHECK(diff(gpu.reduced_cost, cpu.reduced_cost) <= 1e-6);
    CHECK(diff(gpu.row_activity, cpu.row_activity) <= 1e-6);
    CHECK(diff(gpu.metrics.l2_primal, cpu.metrics.l2_primal) <= 1e-6);
    CHECK(diff(gpu.metrics.l2_dual, cpu.metrics.l2_dual) <= 1e-6);
    CHECK(diff(gpu.metrics.l2_gap, cpu.metrics.l2_gap) <= 1e-6);
    CHECK(diff(gpu.metrics.primal_objective, cpu.metrics.primal_objective) <=
          1e-8);
    CHECK(diff(gpu.metrics.dual_objective, cpu.metrics.dual_objective) <= 1e-8);
  }
}

const double kInf = std::numeric_limits<double>::infinity();

hipdlp::HostLp smallLp(int rows, int columns, std::vector<int> start,
                       std::vector<int> index, std::vector<double> value,
                       std::vector<double> cost, std::vector<double> row_lower,
                       std::vector<double> row_upper,
                       std::vector<double> column_lower,
                       std::vector<double> column_upper, int sense = 1) {
  hipdlp::HostLp lp;
  lp.rows = rows;
  lp.columns = columns;
  lp.sense = sense;
  lp.row_start = start;
  lp.column_index = index;
  lp.value = value;
  lp.cost = cost;
  lp.row_lower = row_lower;
  lp.row_upper = row_upper;
  lp.column_lower = column_lower;
  lp.column_upper = column_upper;
  return lp;
}

// Appends a copy of the first row with finite bound whose range cannot be met
// together with the original row: no point is feasible.
void makeInfeasible(hipdlp::HostLp& lp) {
  for (int i = 0; i < lp.rows; ++i) {
    double lower, upper;
    if (std::isfinite(lp.row_upper[i])) {
      lower = lp.row_upper[i] + 1;
      upper = kInf;
    } else if (std::isfinite(lp.row_lower[i])) {
      lower = -kInf;
      upper = lp.row_lower[i] - 1;
    } else {
      continue;
    }
    const std::vector<int> index(lp.column_index.begin() + lp.row_start[i],
                                 lp.column_index.begin() + lp.row_start[i + 1]);
    const std::vector<double> value(lp.value.begin() + lp.row_start[i],
                                    lp.value.begin() + lp.row_start[i + 1]);
    lp.column_index.insert(lp.column_index.end(), index.begin(), index.end());
    lp.value.insert(lp.value.end(), value.begin(), value.end());
    lp.row_start.push_back(static_cast<int>(lp.column_index.size()));
    lp.row_lower.push_back(lower);
    lp.row_upper.push_back(upper);
    ++lp.rows;
    return;
  }
}

// A machine can have the CUDA toolkit without a usable GPU.
bool deviceAvailable() {
  const hipdlp::HostLp lp = makeLp(1, 3, 4, 0.8, 1);
  hipdlp::SolveOptions options;
  options.iteration_limit = 1;
  hipdlp::Snapshot snapshot;
  const hipdlp::SolveResult result =
      hipdlp::SolveForTesting(lp, options, &snapshot, true);
  return result.status == hipdlp::Status::kIterationLimit;
}

}  // namespace

TEST_CASE("hipdlp-cpu-reference-random", "[hipdlp][hipdlp-reference]") {
  if (!deviceAvailable()) {
    WARN("No usable CUDA/HIP device: skipping the GPU-vs-CPU comparison");
    return;
  }
  struct Case {
    const char* name;
    uint64_t seed;
    int rows, columns;
    double density;
    int sense;
  };
  const Case cases[] = {{"tiny", 1, 3, 4, 0.8, 1},
                        {"small", 2, 30, 50, 0.15, 1},
                        {"small-max", 3, 30, 50, 0.15, -1},
                        {"medium", 4, 300, 500, 0.02, 1}};
  Tally tally;
  for (const Case& c : cases)
    compareSolvers(
        c.name, makeLp(c.seed, c.rows, c.columns, c.density, c.sense), tally);
  // A skip rule that fires everywhere would make the test pass vacuously.
  REQUIRE(tally.checked >= 70);
}

TEST_CASE("hipdlp-cpu-reference-netlib", "[hipdlp][hipdlp-reference]") {
  if (!deviceAvailable()) {
    WARN("No usable CUDA/HIP device: skipping the GPU-vs-CPU comparison");
    return;
  }
  Tally tally;
  for (const char* instance :
       {"afiro", "adlittle", "israel", "stair", "standata"})
    compareSolvers(instance, readHostLp(instance), tally);
  REQUIRE(tally.checked >= 100);
}

TEST_CASE("hipdlp-cpu-solve-random", "[hipdlp][hipdlp-reference]") {
  if (!deviceAvailable()) {
    WARN("No usable CUDA/HIP device: skipping the GPU-vs-CPU comparison");
    return;
  }
  struct Case {
    const char* name;
    uint64_t seed;
    int rows, columns;
    double density;
    int sense;
  };
  const Case cases[] = {{"tiny", 1, 3, 4, 0.8, 1},
                        {"small", 2, 30, 50, 0.15, 1},
                        {"small-max", 3, 30, 50, 0.15, -1},
                        {"medium", 4, 300, 500, 0.02, 1}};
  Tally tally;
  for (const Case& c : cases)
    compareToTermination(
        c.name, makeLp(c.seed, c.rows, c.columns, c.density, c.sense), tally);

  // The start point (x = 0, y = 0) is already optimal: the first check, before
  // any iteration, must stop the solve.
  hipdlp::HostLp origin;
  origin.rows = origin.columns = 1;
  origin.row_start = {0, 1};
  origin.column_index = {0};
  origin.value = {1.0};
  origin.cost = {1.0};
  origin.row_lower = {-std::numeric_limits<double>::infinity()};
  origin.row_upper = {5.0};
  origin.column_lower = {0.0};
  origin.column_upper = {10.0};
  compareToTermination("origin-optimal", origin, tally);
  hipdlp::SolveOptions options;
  CHECK(hipdlp::SolveCpu(origin, options).iterations == 0);
  CHECK(hipdlp::SolveCpu(origin, options).status == hipdlp::Status::kOptimal);
  REQUIRE(tally.checked == 10);
}

TEST_CASE("hipdlp-cpu-solve-netlib", "[hipdlp][hipdlp-reference]") {
  if (!deviceAvailable()) {
    WARN("No usable CUDA/HIP device: skipping the GPU-vs-CPU comparison");
    return;
  }
  Tally tally;
  for (const char* instance :
       {"afiro", "adlittle", "israel", "stair", "standata"})
    compareToTermination(instance, readHostLp(instance), tally);
  REQUIRE(tally.checked == 10);
}

TEST_CASE("hipdlp-cpu-certificates", "[hipdlp][hipdlp-reference]") {
  if (!deviceAvailable()) {
    WARN("No usable CUDA/HIP device: skipping the GPU-vs-CPU comparison");
    return;
  }
  using hipdlp::Status;
  struct Case {
    const char* name;
    hipdlp::HostLp lp;
    Status expected;
  };
  std::vector<Case> cases;
  // A row [1, 2] without entries: the start point contradicts it, so the
  // Farkas certificate exists before the first iteration.
  cases.push_back(
      {"empty-row",
       smallLp(1, 1, {0, 0}, {}, {}, {1.0}, {1.0}, {2.0}, {0.0}, {10.0}),
       Status::kInfeasible});
  // x >= 20 with x <= 10
  cases.push_back(
      {"bound-clash",
       smallLp(1, 1, {0, 1}, {0}, {1.0}, {1.0}, {20.0}, {kInf}, {0.0}, {10.0}),
       Status::kInfeasible});
  // x + y >= 10 and x + y <= 5
  cases.push_back(
      {"two-rows",
       smallLp(2, 2, {0, 2, 4}, {0, 1, 0, 1}, {1, 1, 1, 1}, {1, 1},
               {10.0, -kInf}, {kInf, 5.0}, {0.0, 0.0}, {kInf, kInf}),
       Status::kInfeasible});
  // min -x - y with x - y <= 1: the ray (1, 1) decreases the cost for ever.
  cases.push_back({"ray",
                   smallLp(1, 2, {0, 2}, {0, 1}, {1.0, -1.0}, {-1.0, -1.0},
                           {-kInf}, {1.0}, {0.0, 0.0}, {kInf, kInf}),
                   Status::kUnbounded});
  // The same problem maximized: the sign of the cost must follow the sense.
  cases.push_back({"ray-max",
                   smallLp(1, 2, {0, 2}, {0, 1}, {1.0, -1.0}, {1.0, 1.0},
                           {-kInf}, {1.0}, {0.0, 0.0}, {kInf, kInf}, -1),
                   Status::kUnbounded});
  // min x with x free and x + y >= -1, y >= 0
  cases.push_back({"free-ray",
                   smallLp(1, 2, {0, 2}, {0, 1}, {1.0, 1.0}, {1.0, 0.0}, {-1.0},
                           {kInf}, {-kInf, 0.0}, {kInf, kInf}),
                   Status::kUnbounded});
  // Badly scaled rows and columns: the iterates live in scaled coordinates and
  // only give a certificate after unscaling.
  cases.push_back(
      {"ray-scaled",
       smallLp(2, 2, {0, 2, 4}, {0, 1, 0, 1}, {1000.0, -1.0, 1.0, 0.001},
               {-1.0, -1.0}, {-kInf, -10.0}, {1.0, kInf}, {0.0, 0.0},
               {kInf, kInf}),
       Status::kUnbounded});
  cases.push_back(
      {"farkas-scaled",
       smallLp(2, 2, {0, 2, 4}, {0, 1, 0, 1}, {1000.0, 1000.0, 0.001, 0.001},
               {1, 1}, {10000.0, -kInf}, {kInf, 0.005}, {0.0, 0.0},
               {kInf, kInf}),
       Status::kInfeasible});
  // Feasible only far from the start point: the ray is found while the start
  // point is still infeasible.
  cases.push_back({"ray-late",
                   smallLp(2, 2, {0, 2, 4}, {0, 1, 0, 1}, {1.0, -1.0, 1.0, 1.0},
                           {-1.0, -1.0}, {-kInf, 50.0}, {1.0, kInf}, {0.0, 0.0},
                           {kInf, kInf}),
                   Status::kUnbounded});
  // x <= 0 and y in [0, 1]: the ray decreases x, against an upper bound.
  cases.push_back({"ray-down",
                   smallLp(1, 2, {0, 2}, {0, 1}, {1.0, 1.0}, {1.0, 0.0},
                           {-kInf}, {10.0}, {-kInf, 0.0}, {0.0, 1.0}),
                   Status::kUnbounded});
  // x <= 1 but a row asks x >= 3: the Farkas multiplier of the column sits on
  // an upper bound.
  cases.push_back(
      {"farkas-upper",
       smallLp(1, 1, {0, 1}, {0}, {1.0}, {1.0}, {3.0}, {kInf}, {-kInf}, {1.0}),
       Status::kInfeasible});
  // Random LPs made infeasible by a contradicting copy of a row.
  struct Random {
    const char* name;
    uint64_t seed;
    int rows, columns;
    double density;
    int sense;
  };
  const Random random[] = {{"tiny-infeasible", 1, 3, 4, 0.8, 1},
                           {"small-infeasible", 2, 30, 50, 0.15, 1},
                           {"small-max-infeasible", 3, 30, 50, 0.15, -1},
                           {"medium-infeasible", 4, 300, 500, 0.02, 1}};
  for (const Random& r : random) {
    hipdlp::HostLp lp = makeLp(r.seed, r.rows, r.columns, r.density, r.sense);
    makeInfeasible(lp);
    cases.push_back({r.name, lp, Status::kInfeasible});
  }
  cases.push_back({"woodinfe", readHostLp("woodinfe"), Status::kInfeasible});

  Tally tally;
  for (const Case& c : cases) {
    INFO(c.name);
    hipdlp::SolveOptions options;
    options.iteration_limit = 200000;
    // The CPU solver must reach the verdict on its own, not only agree with
    // the GPU about it.
    CHECK(hipdlp::SolveCpu(c.lp, options).status == c.expected);
    compareToTermination(c.name, c.lp, tally);
  }
  REQUIRE(tally.checked == 2 * static_cast<int>(cases.size()));
}

#endif  // HIGHS_HIPDLP_CUDA_DETECTED
