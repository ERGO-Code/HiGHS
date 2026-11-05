#include "solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "solver_test.h"

namespace hipdlp {
namespace {
constexpr double kEta = 0.99;
constexpr int64_t kCheckInterval = 40;
constexpr int64_t kMinBalanceEpoch = 400;
const double kInfinity = std::numeric_limits<double>::infinity();

struct Csr {
  int rows, cols;
  std::vector<int> start, index;
  std::vector<double> value;
};

// The problem in scaled coordinates: x = D_c x_scaled, y = D_r y_scaled,
// A_scaled = D_r A D_c.
struct Model {
  Csr a;
  Csr original;
  std::vector<double> row_scale, column_scale;
  std::vector<double> cost, column_lower, column_upper, row_lower, row_upper;
};

struct Iterate {
  std::vector<double> x, y, xp, yp, xr, yr, dx, dy, xa, ya, epoch_dx, epoch_dy;
};

struct Weight {
  double omega, best, integral;
  int direction;
  int64_t epoch_limit;
};

struct Report {
  std::vector<double> x, y, z, row_activity;
  double l2_primal = 0, l2_dual = 0, l2_gap = 0, primal_objective = 0,
         dual_objective = 0;
  bool bad = false;
};

Model Prepare(const HostLp& lp) {
  Model m;
  // TODO 1: copy the matrix into m.a; row_scale and column_scale start at 1
  m.a.rows = lp.rows;
  m.a.cols = lp.columns;
  m.a.start = lp.row_start;
  m.a.index = lp.column_index;
  m.a.value = lp.value;
  m.row_scale.assign(lp.rows, 1.0);
  m.column_scale.assign(lp.columns, 1.0);
  m.original = m.a;
  for (int sweep = 0; sweep < 11;
       ++sweep) {  // 10 max-norm sweeps, 1 sum-norm sweep
    const bool sum_norm = sweep == 10;
    std::vector<double> row_norm(lp.rows, 0.0), column_norm(lp.columns, 0.0);
    for (int r = 0; r < lp.rows; ++r) {
      for (int k = m.a.start[r]; k < m.a.start[r + 1]; ++k) {
        const int c = m.a.index[k];
        const double magnitude = std::abs(m.a.value[k]);

        if (sum_norm) {
          row_norm[r] += magnitude;
          column_norm[c] += magnitude;
        } else {
          row_norm[r] = std::max(row_norm[r], magnitude);
          column_norm[c] = std::max(column_norm[c], magnitude);
        }
      }
    }

    std::vector<double> row_factor(lp.rows, 1.0);
    std::vector<double> column_factor(lp.columns, 1.0);

    for (int i = 0; i < lp.rows; ++i) {
      if (row_norm[i] > 0) {
        row_factor[i] = 1.0 / std::sqrt(std::max(1e-30, row_norm[i]));
        m.row_scale[i] *= row_factor[i];
      }
    }

    for (int j = 0; j < lp.columns; ++j) {
      if (column_norm[j] > 0) {
        column_factor[j] = 1.0 / std::sqrt(std::max(1e-30, column_norm[j]));
        m.column_scale[j] *= column_factor[j];
      }
    }

    for (int i = 0; i < lp.rows; ++i) {
      for (int k = m.a.start[i]; k < m.a.start[i + 1]; ++k) {
        const int j = m.a.index[k];
        m.a.value[k] *= row_factor[i] * column_factor[j];
      }
    }
  }

  m.cost = lp.cost;
  // c_s = D_c c
  for (int j = 0; j < lp.columns; ++j) {
    m.cost[j] *= lp.sense * m.column_scale[j];
  }

  // x bounds: l_s^x = D_c^{-1} l^x
  m.column_lower = lp.column_lower;
  m.column_upper = lp.column_upper;
  for (int j = 0; j < lp.columns; ++j) {
    m.column_lower[j] /= m.column_scale[j];
    m.column_upper[j] /= m.column_scale[j];
  }

  // row bounds: l_s^r = D_r l^r
  m.row_lower = lp.row_lower;
  m.row_upper = lp.row_upper;
  for (int i = 0; i < lp.rows; ++i) {
    m.row_lower[i] *= m.row_scale[i];
    m.row_upper[i] *= m.row_scale[i];
  }

  return m;
}

double Project(double v, double lower, double upper) {
  return std::min(std::max(v, lower), upper);
}

void Multiply(const Csr& a, const std::vector<double>& x,
              std::vector<double>& y) {
  y.assign(a.rows, 0.0);
  for (int i = 0; i < a.rows; ++i) {
    double sum = 0.0;
    for (int k = a.start[i]; k < a.start[i + 1]; ++k)
      sum += a.value[k] * x[a.index[k]];
    y[i] = sum;
  }
}

void MultiplyTranspose(const Csr& a, const std::vector<double>& y,
                       std::vector<double>& x) {
  x.assign(a.cols, 0.0);
  for (int i = 0; i < a.rows; ++i) {
    for (int k = a.start[i]; k < a.start[i + 1]; ++k)
      x[a.index[k]] += a.value[k] * y[i];
  }
}

Iterate Initialize(const Model& m) {
  Iterate it;
  it.x.resize(m.a.cols);
  for (int j = 0; j < m.a.cols; ++j)
    it.x[j] = Project(0.0, m.column_lower[j], m.column_upper[j]);
  it.y.assign(m.a.rows, 0.0);
  it.xp = it.xr = it.xa = it.x;
  it.yp.assign(m.a.rows, 0.0);
  it.yr.assign(m.a.rows, 0.0);
  it.dx.assign(m.a.cols, 0.0);
  it.dy.assign(m.a.rows, 0.0);
  it.ya.assign(m.a.rows, 0.0);
  it.epoch_dx.assign(m.a.cols, 0.0);
  it.epoch_dy.assign(m.a.rows, 0.0);

  return it;
}

void UpdateWeight(const Iterate& it, int64_t epoch, Weight& w) {
  double dx2 = 0, dy2 = 0;
  for (size_t j = 0; j < it.x.size(); ++j) {
    dx2 += (it.xp[j] - it.xa[j]) * (it.xp[j] - it.xa[j]);
  }
  for (size_t i = 0; i < it.y.size(); ++i) {
    dy2 += (it.yp[i] - it.ya[i]) * (it.yp[i] - it.ya[i]);
  }
  w.epoch_limit = 0;
  if (dx2 > 0 && dy2 > 0 && std::isfinite(dx2) && std::isfinite(dy2)) {
    double e = 0.5 * (std::log(dy2) - std::log(dx2)) - std::log(w.omega);
    w.integral = 0.3 * w.integral + e;
    double requested = 0.99 * e + 0.01 * w.integral;
    double change = Project(requested, -std::log(10.0), std::log(10.0));
    double next_log =
        Project(std::log(w.omega) + change, std::log(1e-8), std::log(1e8));
    int direction =
        std::abs(requested) >= std::log(10.0) &&
                std::abs(next_log - std::log(w.omega)) >= 0.99 * std::log(10.0)
            ? (requested > 0 ? 1 : -1)
            : 0;
    if (direction) {
      const int64_t length = direction == w.direction ? epoch / 2 : epoch;
      w.epoch_limit = length > kMinBalanceEpoch ? length : kMinBalanceEpoch;
    }
    w.direction = direction;
    w.omega = std::exp(next_log);
    w.best = w.omega;
  } else {
    w.direction = 0;
    w.omega = w.best;
    w.integral = 0;
  }
}

void Evaluate(const Model& m, Iterate& it, double omega) {
  const int n = m.a.cols, rows = m.a.rows;
  std::vector<double> aty, axr;
  MultiplyTranspose(m.a, it.y, aty);
  for (int j = 0; j < n; ++j) {
    double raw = it.x[j] - (kEta / omega) * (m.cost[j] - aty[j]);
    it.xp[j] = Project(raw, m.column_lower[j], m.column_upper[j]);
    it.dx[j] = it.x[j] - it.xp[j];
    it.xr[j] = 2 * it.xp[j] - it.x[j];
  }
  Multiply(m.a, it.xr, axr);
  const double sigma = kEta * omega;
  for (int i = 0; i < rows; ++i) {
    double t = it.y[i] / sigma - axr[i];
    it.yp[i] = sigma * (t - Project(t, -m.row_upper[i], -m.row_lower[i]));
    it.dy[i] = it.y[i] - it.yp[i];
    it.yr[i] = 2 * it.yp[i] - it.y[i];
  }
}

void Blend(const Model& m, Iterate& it, int64_t epoch) {
  double alpha = 1.0 / (epoch + 2);
  for (int i = 0; i < m.a.rows; ++i) {
    it.y[i] = alpha * it.ya[i] + (1 - alpha) * it.yr[i];
  }

  for (int j = 0; j < m.a.cols; ++j) {
    it.x[j] = alpha * it.xa[j] + (1 - alpha) * it.xr[j];
  }
}

double FixedPointError(const Model& m, const Iterate& it, double omega) {
  std::vector<double> atdy;
  MultiplyTranspose(m.a, it.dy, atdy);
  double dx2 = 0, dy2 = 0, cross = 0;
  for (int j = 0; j < m.a.cols; ++j) {
    dx2 += it.dx[j] * it.dx[j];
    cross += it.dx[j] * atdy[j];
  }
  for (int i = 0; i < m.a.rows; ++i) dy2 += it.dy[i] * it.dy[i];
  const double energy = omega * dx2 + dy2 / omega + 2 * kEta * cross;
  return std::sqrt(std::max(0.0, energy));
}

double Multiplier(double v, double lower, double upper) {
  if (!std::isfinite(lower)) {
    v = std::min(v, 0.0);
  }
  if (!std::isfinite(upper)) {
    v = std::max(v, 0.0);
  }
  return v;
}

double BoundProduct(double v, double lower, double upper) {
  if (v == 0) {
    return 0;
  } else if (v > 0) {
    return lower * v;
  } else {
    return upper * v;
  }
}

Report Measure(const HostLp& lp, const Model& m, const Iterate& it,
               double tolerance) {
  const int n = m.a.cols, rows = m.a.rows;
  Report r;
  r.x.resize(n);
  r.y.resize(rows);
  r.z.resize(n);

  for (int j = 0; j < m.a.cols; ++j) {
    r.x[j] = m.column_scale[j] * it.xp[j];
  }
  for (int i = 0; i < m.a.rows; ++i) {
    r.y[i] =
        Multiplier(m.row_scale[i] * it.yp[i], lp.row_lower[i], lp.row_upper[i]);
  }

  std::vector<double> aty;
  Multiply(m.original, r.x, r.row_activity);
  MultiplyTranspose(m.original, r.y, aty);

  double p = 0, d = 0, cross = 0, primal_sq = 0, dual_sq = 0, c_sq = 0,
         b_sq = 0;
  for (int j = 0; j < n; ++j) {
    const double c = lp.sense * lp.cost[j];
    const double l = lp.column_lower[j], u = lp.column_upper[j], x = r.x[j];
    const double raw_z = c - aty[j];
    r.z[j] = Multiplier(raw_z, l, u);
    p += c * x;
    d += BoundProduct(r.z[j], l, u);
    const double violation = std::max(std::max(l - x, x - u), 0.0);
    primal_sq += violation * violation;
    const bool at_lower = std::isfinite(l) && (x - l) * (x - l) <= tolerance;
    const bool at_upper = std::isfinite(u) && (x - u) * (x - u) <= tolerance;
    const double z = at_lower && at_upper ? raw_z
                     : at_lower           ? std::max(raw_z, 0.0)
                     : at_upper           ? std::min(raw_z, 0.0)
                                          : 0.0;
    dual_sq += (raw_z - z) * (raw_z - z);
    c_sq += c * c;
    cross += z * x;
  }
  for (int i = 0; i < rows; ++i) {
    const double l = lp.row_lower[i], u = lp.row_upper[i],
                 ax = r.row_activity[i];
    const double violation = std::max(std::max(l - ax, ax - u), 0.0);
    primal_sq += violation * violation;
    if (std::isfinite(l) && l != u) b_sq += l * l;
    if (std::isfinite(u)) b_sq += u * u;
    const double product = BoundProduct(r.y[i], l, u);
    cross += product;
    d += product;
  }

  r.l2_primal = std::sqrt(primal_sq) / (1 + std::sqrt(b_sq));
  r.l2_dual = std::sqrt(dual_sq) / (1 + std::sqrt(c_sq));
  r.l2_gap = std::fabs(p - cross) / (1 + std::fabs(p) + std::fabs(cross));
  r.primal_objective = p;
  r.dual_objective = d;
  r.bad = !std::isfinite(primal_sq + dual_sq + b_sq + c_sq + p + d + cross);
  return r;
}

double Seconds(std::chrono::steady_clock::time_point since) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - since)
      .count();
}

Status NextStatus(const Report& report, const SolveOptions& options,
                  int64_t iterations, double seconds) {
  const bool optimal = !report.bad && report.l2_primal <= options.tolerance &&
                       report.l2_dual <= options.tolerance &&
                       report.l2_gap <= options.tolerance;
  if (report.bad) return Status::kNumericalError;
  if (optimal) return Status::kOptimal;
  if (seconds >= options.time_limit_seconds) return Status::kTimeLimit;
  if (iterations >= options.iteration_limit) return Status::kIterationLimit;
  return Status::kRunning;
}

SolveResult RunCpu(const HostLp& lp, const SolveOptions& options,
                   Snapshot* snapshot, bool fixed_weight) {
  const std::chrono::steady_clock::time_point start =
      std::chrono::steady_clock::now();
  SolveResult result;
  result.rows = lp.rows;
  result.columns = lp.columns;
  result.nonzeros = lp.value.size();
  const Model m = Prepare(lp);
  Iterate it = Initialize(m);
  Weight w;
  w.omega = w.best = 1.0;
  w.integral = 0;
  w.direction = 0;
  w.epoch_limit = 0;

  Status status = Status::kRunning;
  int64_t iterations = 0, epoch = 0, restarts = 0, epoch_iteration = 0;
  double fpe = 0, reference_fpe = 0, previous_fpe = kInfinity;

  if (snapshot) snapshot->decisions.clear();
  auto log_decision = [&](bool restart) {
    if (!snapshot) return;
    Decision d;
    d.iteration = iterations;
    d.epoch = epoch;
    d.status = status;
    d.restart = restart;
    d.omega = w.omega;
    d.fpe = fpe;
    d.reference_fpe = reference_fpe;
    d.balance_epoch_limit = w.epoch_limit;
    snapshot->decisions.push_back(d);
  };

  Report report = Measure(lp, m, it, options.tolerance);
  int64_t kkt_checks = 1;
  status = NextStatus(report, options, iterations, Seconds(start));
  log_decision(false);
  if (status == Status::kRunning) {
    Evaluate(m, it, w.omega);
    fpe = reference_fpe = FixedPointError(m, it, w.omega);
  }

  while (status == Status::kRunning) {
    Evaluate(m, it, w.omega);
    Blend(m, it, epoch);
    ++iterations;
    ++epoch;
    const bool check = iterations % kCheckInterval == 0 ||
                       iterations >= options.iteration_limit ||
                       Seconds(start) >= options.time_limit_seconds;
    if (!check) continue;

    fpe = FixedPointError(m, it, w.omega);
    report = Measure(lp, m, it, options.tolerance);
    ++kkt_checks;
    status = NextStatus(report, options, iterations, Seconds(start));
    // A balance deadline only interrupts an epoch whose FPE has not fallen by
    // half.
    const bool restart =
        status == Status::kRunning && epoch > 0 &&
        iterations % kCheckInterval == 0 && reference_fpe > 0 &&
        (iterations == kCheckInterval || fpe <= 0.2 * reference_fpe ||
         (fpe <= 0.5 * reference_fpe && fpe > previous_fpe) ||
         double(epoch) >= 0.36 * double(iterations) ||
         (w.epoch_limit > 0 && epoch >= w.epoch_limit &&
          fpe > 0.5 * reference_fpe));

    previous_fpe = fpe;
    log_decision(restart);
    if (!restart) continue;
    if (!fixed_weight) UpdateWeight(it, epoch, w);
    for (int j = 0; j < m.a.cols; ++j) {
      it.epoch_dx[j] = it.x[j] - it.xa[j];
      it.x[j] = it.xa[j] = it.xp[j];
    }
    for (int i = 0; i < m.a.rows; ++i) {
      it.epoch_dy[i] = it.y[i] - it.ya[i];
      it.y[i] = it.ya[i] = it.yp[i];
    }
    ++restarts;
    epoch_iteration = iterations;
    epoch = 0;
    Evaluate(m, it, w.omega);
    fpe = reference_fpe = FixedPointError(m, it, w.omega);
    previous_fpe = kInfinity;
  }

  result.status = status;
  result.iterations = iterations;
  result.restarts = restarts;
  result.kkt_checks = kkt_checks;
  result.metrics_valid = !report.bad;
  result.primal_valid = !report.bad && report.l2_primal <= options.tolerance;
  result.dual_valid = status == Status::kOptimal;
  result.primal = report.x;
  result.dual = report.y;
  result.reduced_cost = report.z;
  for (double& y : result.dual) y *= lp.sense;
  for (double& z : result.reduced_cost) z *= lp.sense;
  result.row_activity = report.row_activity;
  result.metrics.l2_primal = report.l2_primal;
  result.metrics.l2_dual = report.l2_dual;
  result.metrics.l2_gap = report.l2_gap;
  result.metrics.primal_objective =
      report.primal_objective * lp.sense + lp.offset;
  result.metrics.dual_objective = report.dual_objective * lp.sense + lp.offset;
  result.times.total = Seconds(start);
  if (snapshot) {
    snapshot->row_multiplier = m.row_scale;
    snapshot->column_multiplier = m.column_scale;
    snapshot->scaled_value = m.a.value;
    snapshot->x = it.x;
    snapshot->y = it.y;
    snapshot->xp = it.xp;
    snapshot->yp = it.yp;
    snapshot->xa = it.xa;
    snapshot->ya = it.ya;
    snapshot->epoch_dx = it.epoch_dx;
    snapshot->epoch_dy = it.epoch_dy;
    snapshot->epoch_iteration = epoch_iteration;
    snapshot->omega = w.omega;
    snapshot->fpe = fpe;
  }
  return result;
}

}  // unnamed namespace

SolveResult SolveCpu(const HostLp& lp, const SolveOptions& options) {
  return RunCpu(lp, options, nullptr, false);
}

SolveResult SolveCpuForTesting(const HostLp& lp, const SolveOptions& options,
                               Snapshot* snapshot, bool fixed_weight) {
  return RunCpu(lp, options, snapshot, fixed_weight);
}

const char* StatusName(Status status) {
  switch (status) {
    case Status::kRunning:
      return "RUNNING";
    case Status::kOptimal:
      return "OPTIMAL";
    case Status::kInfeasible:
      return "INFEASIBLE";
    case Status::kUnbounded:
      return "UNBOUNDED";
    case Status::kTimeLimit:
      return "TIME_LIMIT";
    case Status::kIterationLimit:
      return "ITERATION_LIMIT";
    case Status::kInvalidModel:
      return "INVALID_MODEL";
    case Status::kUnsupportedModel:
      return "UNSUPPORTED_MODEL";
    case Status::kOutOfMemory:
      return "OUT_OF_MEMORY";
    case Status::kNumericalError:
      return "NUMERICAL_ERROR";
    case Status::kDeviceError:
      return "DEVICE_ERROR";
  }
  return "DEVICE_ERROR";
}

}  // namespace hipdlp
