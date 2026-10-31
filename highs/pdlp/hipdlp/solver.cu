#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "device.h"
#include "kernel_types.h"
#include "solver.h"

namespace hipdlp {
constexpr int kBlock = 256;
constexpr double kEta = 0.99;
constexpr double kCertificateTolerance = 1e-8;
constexpr int64_t kMinBalanceEpoch = 400; 

constexpr double kDeviationLevel = 1e-5;
constexpr int kCheckInterval = 40;  // iterations between two checks

// HIP has no conditional graph nodes, so there the host drives the control flow: every IF /
// WHILE body is a graph of its own, the decision kernels write their branch into a flag in
// device memory, and the host reads the flags between launches. A chunk graph runs
// kCheckInterval iterations and the check, so the host waits once per check. The same kernels
// serve both builds, only SetBranch differs. -DHIPDLP_HOST_LOOP builds this path with CUDA, to
// test it against the conditional-node graph.
#if !defined(HIPDLP_HOST_LOOP) && defined(USE_HIP)
#define HIPDLP_HOST_LOOP
#endif
#ifdef HIPDLP_HOST_LOOP
using Branch = unsigned*;
static __device__ void SetBranch(Branch b, bool value) { *b = value; }
// True when branch a or b (null: none) is taken.
static __device__ bool Taken(const unsigned* a, const unsigned* b = nullptr) {
  return (a && *a) || (b && *b);
}
// A check fell due inside a chunk (time or iteration limit, nonfinite values): the iterate stays
// put until the check that ends the chunk. The products in between recompute the same values.
static __device__ bool Held(const DeviceState* s) { return s->hold; }
#else
using Branch = cudaGraphConditionalHandle;
static __device__ void SetBranch(Branch b, bool value) { cudaGraphSetConditional(b, value); }
static __device__ bool Held(const DeviceState*) { return false; }
#endif

using Clock = std::chrono::steady_clock;
static double Seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
static int Blocks(int size) { return std::max(1, int((int64_t(size) + kBlock - 1) / kBlock)); }

// Blocks for one lane group per row: kBlock / lanes groups per block, at most 4096 blocks.
static int WarpBlocks(int rows, int lanes) {
  int64_t groups = kBlock / lanes;
  return std::max(1, int(std::min(int64_t(4096), (int64_t(rows) + groups - 1) / groups)));
}

#ifdef USE_HIP
// Nanoseconds per wall_clock64() tick, set by CheckDevice from the device's wall clock rate.
__device__ double wall_clock_ns;
static __device__ unsigned long long DeviceTime() {
  return (unsigned long long)(double(wall_clock64()) * wall_clock_ns);
}
#else
static __device__ unsigned long long DeviceTime() {
  unsigned long long value;
  asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(value));
  return value;
}
#endif

static __device__ bool Expired(const DeviceState* s) {
  return double(DeviceTime() - s->begin) * 1e-9 >= s->remaining_seconds;
}

static __device__ double Project(double v, double l, double u) { return fmin(fmax(v, l), u); }

static __device__ double Multiplier(double v, double l, double u) {
  if (!isfinite(l)) v = fmin(v, 0.0);
  if (!isfinite(u)) v = fmax(v, 0.0);
  return v;
}
static __device__ double Recession(double v, double l, double u) {
  if (isfinite(l)) v = fmax(v, 0.0);
  if (isfinite(u)) v = fmin(v, 0.0);
  return v;
}
static __device__ double BoundProduct(double v, double l, double u) {
  return v > 0 ? v * l : (v < 0 ? v * u : 0.0);
}
static __device__ double ConeViolation(double v, double l, double u) {
  double error = 0;
  if (isfinite(l)) error = fmax(error, -v);
  if (isfinite(u)) error = fmax(error, v);
  return error;
}

static __device__ void BoundMetrics(Stats& s, double v, double multiplier, double l, double u,
                                    int index) {
  if (!isfinite(v) || !isfinite(multiplier)) {
    s.bad = 1;
    return;
  }
  double raw = 0, relative = 0;
  if (isfinite(l)) {
    raw = fmax(raw, l - v);
    relative = fmax(relative, (l - v) / (1 + fabs(l)));
  }
  if (isfinite(u)) {
    raw = fmax(raw, v - u);
    relative = fmax(relative, (v - u) / (1 + fabs(u)));
  }
  if (relative > s.primal) {
    s.primal = relative;
    s.primal_index = index;
  }
  if (raw > s.raw) {
    s.raw = raw;
    s.raw_index = index;
  }
  double contribution = BoundProduct(multiplier, l, u);
  double comp = multiplier == 0 ? 0 : fabs(multiplier * v - contribution);
  s.d += contribution;
  if (comp > s.comp) {
    s.comp = comp;
    s.comp_index = index;
  }
  if (!isfinite(contribution) || !isfinite(comp)) s.bad = 1;
}

// lpKktCheck pass 1 for one row activity or column value v with multiplier `mult` and bounds
// [l, u]: a value within sqrt(tol) of a bound is active and puts that bound into the active
// bound norm; the dual objective takes the bound on the side of the midpoint v is on (not the
// bound the sign of the multiplier selects), and a free variable contributes its multiplier.
static __device__ void KktBounds(Stats& s, double v, double mult, double l, double u, double tol) {
  if (!isfinite(v) || !isfinite(mult)) return;
  if (isfinite(l) && (l - v) * (l - v) <= tol)
    s.active_bound = fmax(s.active_bound, fabs(l));
  else if (isfinite(u) && (v - u) * (v - u) <= tol)
    s.active_bound = fmax(s.active_bound, fabs(u));
  double bound = 1;
  if (isfinite(l) || isfinite(u)) bound = v < 0.5 * (l + u) ? l : u;
  s.dual_pos += bound * mult;
}
// lpKktCheck pass 2 (getKktFailures) for the same value: primal violation over 1 + the active
// bound norm (or the violated bound itself when the value is not near a bound and that bound is
// larger), dual sign violation over 1 + the active cost norm (or the cost when the multiplier
// is not near zero). Only a boxed variable of meaningful length has a two-sided sign rule: it
// uses the multiplier sign that the side of the midpoint requires. `cost` is 0 for rows.
static __device__ void KktRelative(Stats& s, double v, double mult, double l, double u,
                                   double cost, double tol, double active_bound,
                                   double active_cost) {
  if (!isfinite(v) || !isfinite(mult)) return;
  double violation = fmax(fmax(l - v, v - u), 0.0);
  bool near_bound = (l - v) * (l - v) <= tol || (v - u) * (v - u) <= tol;
  int side = 0;  // -1: below the midpoint (lower rule), +1: above (upper rule), 0: none
  double sign_violation = 0;
  if (l < u) {
    if (!isfinite(l) && !isfinite(u))
      sign_violation = fabs(mult);
    else if ((u - l) * (u - l) > tol) {
      side = v < 0.5 * (l + u) ? -1 : 1;
      sign_violation = side < 0 ? fmax(-mult, 0.0) : fmax(mult, 0.0);
    }
  }
  if (violation > 0) {
    double measure = active_bound;
    if (!near_bound) measure = fmax(side > 0 ? fabs(u) : fabs(l), measure);
    s.rel_primal = fmax(s.rel_primal, violation / (1 + measure));
  }
  if (sign_violation > 0) {
    double measure = active_cost;
    if (cost != 0 && mult * mult >= tol) measure = fmax(fabs(cost), measure);
    s.rel_dual = fmax(s.rel_dual, sign_violation / (1 + measure));
  }
}
static __device__ void StoreStats(Stats value, Stats* partial) {
  __shared__ typename cub::BlockReduce<Stats, kBlock>::TempStorage storage;
  Stats total = cub::BlockReduce<Stats, kBlock>(storage).Reduce(value, Combine{});
  if (threadIdx.x == 0) partial[blockIdx.x] = total;
}

__global__ void Begin(DeviceState* s) {
  s->decision_count = 0;
  s->farkas_info = CertificateInfo{};
  s->ray_info = CertificateInfo{};
  s->feasible_iteration = 0;
  s->direction_iteration = s->epoch_iteration = 0;
  s->direction_source = -1;
  s->direction_sign = 0;
  s->status = s->remaining_seconds <= 0 ? Status::kTimeLimit : Status::kRunning;
  s->bad = 0;
  s->iterations = s->epoch = s->restarts = s->kkt_checks = s->certificate_checks = s->spmv_count =
      0;
  s->feasible_valid = s->optimal = s->farkas_valid = s->ray_valid = s->restart = 0;
  s->omega = s->best_omega = 1;
  s->balance_epoch_limit = 0;
  s->balance_direction = 0;
  s->deviation = s->deviation_next = 0;
  s->hold = 0;
  s->integral = 0;
  s->previous_fpe = INFINITY;
  s->reference_fpe = s->fpe = 0;
  s->margin = s->stationarity = s->cert_scale = 0;
  s->ray_margin = s->ray_row_violation = s->ray_column_violation = 0;
  s->validation_ns = 0;
  s->metrics = Metrics{};
  s->begin = DeviceTime();
}

// A bounded grid of warps handles both sparse and very long rows/columns.
__global__ void ScalingFactors(int size, const int* start, const double* a, double* factor,
                               double* cumulative, bool pc, int lanes, DeviceState* s) {
  if (s->status != Status::kRunning) return;
  int lane = threadIdx.x % lanes;
  int stride = gridDim.x * blockDim.x / lanes;
  for (int i = (blockIdx.x * blockDim.x + threadIdx.x) / lanes; i < size; i += stride) {
    int active = 0;
    if (lane == 0) active = s->status == Status::kRunning && !Expired(s);
    active = GroupBroadcast(active, lanes);
    if (!active) {
      if (lane == 0) factor[i] = 1;
      continue;
    }
    double norm = 0;
    for (int k = start[i] + lane; k < start[i + 1]; k += lanes) {
      double v = fabs(a[k]);
      if (!isfinite(v)) atomicOr(&s->bad, 1);
      norm = pc ? norm + v : fmax(norm, v);
    }
    for (int offset = lanes / 2; offset > 0; offset /= 2) {
      double other = GroupShflDown(norm, offset, lanes);
      norm = pc ? norm + other : fmax(norm, other);
    }
    if (lane == 0) {
      double f = norm == 0 ? 1 : 1 / sqrt(fmax(norm, 1e-30));
      factor[i] = f;
      cumulative[i] *= f;
      if (!isfinite(norm) || !isfinite(cumulative[i]) || cumulative[i] <= 0) atomicOr(&s->bad, 1);
    }
  }
}

__global__ void ScaleMatrix(int size, const int* start, const int* index, double* a,
                            const double* left, const double* right, int lanes, DeviceState* s) {
  if (s->status != Status::kRunning) return;
  int stride = gridDim.x * blockDim.x / lanes;
  for (int i = (blockIdx.x * blockDim.x + threadIdx.x) / lanes; i < size; i += stride) {
    for (int k = start[i] + threadIdx.x % lanes; k < start[i + 1]; k += lanes) {
      a[k] *= left[i] * right[index[k]];
      if (!isfinite(a[k])) atomicOr(&s->bad, 1);
    }
  }
}

__global__ void ScalingGuard (DeviceState* s) {
  if (s-> bad)
    s->status = Status::kNumericalError;
  else if (Expired(s))
    s->status = Status::kTimeLimit;
}

__global__ void Initialize(DeviceLp lp, Vectors v, DeviceState* s)  {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < lp.n) {
    lp.cs[i] = lp.c[i] * lp.dc[i];
    lp.lxs[i] = lp.lx[i] / lp.dc[i];
    lp.uxs[i] = lp.ux[i] / lp.dc[i];
    double start = s->warm_start ? v.out_x[i] / lp.dc[i] : 0; //out_x holds the start point
    v.x[i] = v.xp[i] = v.xr[i] = v.xa[i] = Project(start, lp.lxs[i], lp.uxs[i]);
    v.dx[i] = v.epoch_dx[i] = 0;
    if (!isfinite(lp.cs[i]) || (isfinite(lp.lx[i]) && !isfinite(lp.lxs[i])) ||
        (isfinite(lp.ux[i]) && !isfinite(lp.uxs[i])) || !isfinite(v.x[i]) || !isfinite(start))
      atomicOr(&s->bad, 1);
    }
    if (i < lp.m) {
      lp.lrs[i] = lp.lr[i] * lp.dr[i];
      lp.urs[i] = lp.ur[i] * lp.dr[i];
      double start = s->warm_start ? v.out_y[i] / lp.dr[i] : 0.0; // already of a legal sign
      v.y[i] = v.yp[i] = v.yr[i] = v.ya[i] = start;
      v.dy[i] = v.epoch_dy[i] = 0;
      if ((isfinite(lp.lr[i]) && !isfinite(lp.lrs[i])) ||
          (isfinite(lp.ur[i]) && !isfinite(lp.urs[i])) || !isfinite(start))
        atomicOr(&s->bad, 1);
  }
}

__global__ void PrimalStep(DeviceLp lp, Vectors v, DeviceState* s) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= lp.n || s->status != Status::kRunning || Held(s)) return;
  double raw;
  if (s->deviation) {
    raw = v.x[j] - (kEta / s->omega) * (v.ga[j] - v.aty[j]);
    v.xp[j] = Project(raw, lp.lxs[j] - v.xa[j], lp.uxs[j] - v.xa[j]);
  } else { // absolute deviation
    raw = v.x[j] - (kEta / s->omega) * (lp.cs[j] - v.aty[j]);
    v.xp[j] = Project(raw, lp.lxs[j], lp.uxs[j]);
  }
  v.dx[j] = v.x[j] - v.xp[j];
  v.xr[j] = 2 * v.xp[j] - v.x[j];
  if (!isfinite(raw) || !isfinite(v.xr[j])) atomicOr(&s->bad, 1);
}

__global__ void DualStep(DeviceLp lp, Vectors v, DeviceState* s, bool blend) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (s->status != Status::kRunning || Held(s)) return;
  double alpha = 1 / (double(s->epoch) + 2);
  if (i < lp.m) {
    double sigma = kEta * s->omega;
    double t = v.y[i] / sigma - v.axr[i];
    if (s->deviation) {
      // t is the deviation of the projection argument. The row binds on its lower (upper) side
      // iff the absolute multiplier ya + sigma (t - lower) is positive (ya + sigma (t - upper) is
      // negative); otherwise the multiplier is 0. An infinite side gives -inf (+inf), never binds.
      double low = sigma * (t - v.lower[i]), up = sigma * (t - v.upper[i]);
      v.yp[i] = v.ya[i] + low > 0 ? low : (v.ya[i] + up < 0 ? up : -v.ya[i]);
    } else
      v.yp[i] = sigma * (t - Project(t, -lp.urs[i], -lp.lrs[i]));
    v.dy[i] = v.y[i] - v.yp[i];
    v.yr[i] = 2 * v.yp[i] - v.y[i];
    if (!isfinite(t) || !isfinite(v.yr[i])) atomicOr(&s->bad, 1);
    if (blend && s->deviation)
      v.y[i] = (1 - alpha) * v.yr[i];  // the anchor's deviation is 0
    else if (blend)
      v.y[i] = alpha * v.ya[i] + (1 - alpha) * v.yr[i];
    if (blend && !isfinite(v.y[i])) atomicOr(&s->bad, 1);
  }
  if (blend && i < lp.n && s->deviation)
    v.x[i] = (1 - alpha) * v.xr[i];
  else if (blend && i < lp.n)
    v.x[i] = alpha * v.xa[i] + (1 - alpha) * v.xr[i];
  if (blend && i < lp.n && !isfinite(v.x[i])) atomicOr(&s->bad, 1);
}

__global__ void StepControl(DeviceState* s, Branch check) {
  if (s->status == Status::kRunning && !Held(s)) {
    ++s->iterations;
    ++s->epoch;
    s->spmv_count += 2;
  }
  bool due = s->bad || s->iterations % kCheckInterval == 0 ||
             s->iterations >= s->iteration_limit || Expired(s);
  SetBranch(check, due);
#ifdef HIPDLP_HOST_LOOP
  s->hold = due;
#endif
}

__global__ void FpeStats(DeviceLp lp, Vectors v, Stats* partial) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  if (i < lp.n) {
    a.dx2 = v.dx[i] * v.dx[i];
    a.cross = v.dx[i] * v.atdy[i];
  }
  if (i < lp.m) a.dy2 = v.dy[i] * v.dy[i];
  a.absolute = fabs(a.cross);
  if (!isfinite(a.dx2) || !isfinite(a.dy2) || !isfinite(a.cross)) a.bad = 1;
  StoreStats(a, partial);
}

__global__ void Unscale(DeviceLp lp, Vectors v, DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < lp.n && s->deviation)
    v.out_x[i] = lp.dc[i] * (v.xa[i] + v.xp[i]);
  else if (i < lp.n)
    v.out_x[i] = lp.dc[i] * v.xp[i];
  if (i < lp.m) {
    double y;
    if (s->deviation)
      y = lp.dr[i] * (v.ya[i] + v.yp[i]);
    else
      y = lp.dr[i] * v.yp[i];
    if (!isfinite(y)) atomicOr(&s->bad, 1);
    v.out_y[i] = Multiplier(y, lp.lr[i], lp.ur[i]);
  }
}

__global__ void KktStats(DeviceLp lp, Vectors v, Stats* partial, double tol) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  if (i < lp.n) {
    double raw_z = lp.c[i] - v.out_aty[i];
    double z = Multiplier(raw_z, lp.lx[i], lp.ux[i]);
    v.out_z[i] = z;
    a.dual = fabs(raw_z - z);
    a.dual_index = i;
    a.p = lp.c[i] * v.out_x[i];
    if (!isfinite(raw_z) || !isfinite(a.p)) a.bad = 1;
    BoundMetrics(a, v.out_x[i], z, lp.lx[i], lp.ux[i], lp.m + i);
    KktBounds(a, v.out_x[i], z, lp.lx[i], lp.ux[i], tol);
    if (isfinite(z) && z * z < tol) a.active_cost = fabs(lp.c[i]);
    double violation = fmax(fmax(lp.lx[i] - v.out_x[i], v.out_x[i] - lp.ux[i]), 0.0);
    a.l2_primal_sq = violation * violation;
    // L2 measures as cuPDLPx: its reduced cost is the primal projection's slack, i.e. nonzero only
    // on a bound the value sits at (here: within sqrt(tol), as KktBounds) and only with the sign
    // that bound allows; any other part of c - A'y is dual residual.
    double x = v.out_x[i];
    bool at_lower = isfinite(lp.lx[i]) && (x - lp.lx[i]) * (x - lp.lx[i]) <= tol;
    bool at_upper = isfinite(lp.ux[i]) && (x - lp.ux[i]) * (x - lp.ux[i]) <= tol;
    double zl2 = at_lower && at_upper ? raw_z
                 : at_lower           ? fmax(raw_z, 0.0)
                 : at_upper           ? fmin(raw_z, 0.0)
                                      : 0.0;
    a.l2_dual_sq = (raw_z - zl2) * (raw_z - zl2);
    a.c_sq = lp.c[i] * lp.c[i];
    a.cross = zl2 * x;  // KKT pass: `cross` carries cuPDLPx's dual objective (l2_gap)
  }
  if (i < lp.m) {
    BoundMetrics(a, v.out_ax[i], v.out_y[i], lp.lr[i], lp.ur[i], i);
    KktBounds(a, v.out_ax[i], v.out_y[i], lp.lr[i], lp.ur[i], tol);
    double violation = fmax(fmax(lp.lr[i] - v.out_ax[i], v.out_ax[i] - lp.ur[i]), 0.0);
    a.l2_primal_sq += violation * violation;
    if (isfinite(lp.lr[i]) && lp.lr[i] != lp.ur[i]) a.b_sq += lp.lr[i] * lp.lr[i];
    if (isfinite(lp.ur[i])) a.b_sq += lp.ur[i] * lp.ur[i];
    a.cross += BoundProduct(v.out_y[i], lp.lr[i], lp.ur[i]);
  }
  StoreStats(a, partial);
}

// Pass 2 of the lpKktCheck measures; reads the norms ReduceStats(mode 0) stored in `s`. out_aty
// still holds A^T y of this check (nothing runs between the two passes).
__global__ void KktRelativeStats(DeviceLp lp, Vectors v, const DeviceState* s, Stats* partial,
                                 double tol) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  double active_bound = s->active_bound_norm, active_cost = s->active_cost_norm;
  if (i < lp.n) {
    KktRelative(a, v.out_x[i], v.out_z[i], lp.lx[i], lp.ux[i], lp.c[i], tol, active_bound,
                active_cost);
    double residual = fabs(lp.c[i] - v.out_aty[i] - v.out_z[i]);
    if (isfinite(residual)) a.rel_residual = residual / (1 + active_cost);
  }
  if (i < lp.m)
    KktRelative(a, v.out_ax[i], v.out_y[i], lp.lr[i], lp.ur[i], 0.0, tol, active_bound,
                active_cost);
  StoreStats(a, partial);
}
__global__ void SaveFeasible(DeviceLp lp, Vectors v, DeviceState* s) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < lp.n && s->metrics.primal <= s->tolerance && !s->bad) v.feasible[j] = v.out_x[j];
}

__global__ void CertificateDirection(DeviceLp lp, Vectors v, DeviceState* s, int source, int sign) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i == 0) {
    s->direction_source = s->iterations == 0 && source != 2 ? 3 : source;
    s->direction_sign = sign;
    s->direction_iteration = source == 2 ? s->epoch_iteration : s->iterations;
  }
  if (i < lp.n) {
    double d = source == 2 ? v.epoch_dx[i]
                           : (source == 0 ? (s->deviation ? v.x[i] : v.x[i] - v.xa[i])
                                          : v.x[i] - v.xr[i]);
    d *= sign * lp.dc[i];
    if (!isfinite(d)) atomicOr(&s->bad, 1);
    v.ray[i] = Recession(d, lp.lx[i], lp.ux[i]);
  }
  if (i < lp.m) {
    double f = source == 2 ? v.epoch_dy[i]
                           : (source == 0 ? (s->deviation ? v.y[i] : v.y[i] - v.ya[i])
                                          : v.y[i] - v.yr[i]);
    // Empty-row contradiction can be certified before the first operator evaluation.
    if (s->iterations == 0 && source != 2)
      f = v.out_ax[i] < lp.lr[i] ? 1 : (v.out_ax[i] > lp.ur[i] ? -1 : 0);
    else
      f *= lp.dr[i];
    f *= sign;
    if (!isfinite(f)) atomicOr(&s->bad, 1);
    v.f[i] = Multiplier(f, lp.lr[i], lp.ur[i]);
  }
}

__global__ void CertificateNorm(DeviceLp lp, Vectors v, Stats* partial, bool farkas) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  if (farkas) {
    if (i < lp.m) a.norm = fabs(v.f[i]);
    if (i < lp.n) {
      if (!isfinite(v.atf[i])) a.bad = 1;
      v.g[i] = Multiplier(-v.atf[i], lp.lx[i], lp.ux[i]);
      a.norm = fmax(a.norm, fabs(v.g[i]));
    }
  } else if (i < lp.n)
    a.norm = fabs(v.ray[i]);
  if (!isfinite(a.norm)) a.bad = 1;
  StoreStats(a, partial);
}

__global__ void NormalizeCertificate(DeviceLp lp, Vectors v, DeviceState* s, bool farkas) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  double norm = s->normalization;
  if (!(norm > 0) || !isfinite(norm)) return;
  if (farkas) {
    if (i < lp.m) v.f[i] /= norm;
    if (i < lp.n) v.g[i] /= norm;
  } else if (i < lp.n)
    v.ray[i] /= norm;
}

__global__ void CertificateStats(DeviceLp lp, Vectors v, Stats* partial, bool farkas) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  if (farkas) {
    if (i < lp.m) {
      double b = BoundProduct(v.f[i], lp.lr[i], lp.ur[i]);
      a.d += b;
      a.absolute += fabs(b);
    }
    if (i < lp.n) {
      double b = BoundProduct(v.g[i], lp.lx[i], lp.ux[i]);
      a.d += b;
      a.absolute += fabs(b);
      double residual = fabs(v.atf[i] + v.g[i]);
      // Judged against column j's own 1-norm (|A'f|_j <= ||A_j||_1 for ||f||_inf <= 1), so one
      // huge column cannot widen the tolerance of every other.
      a.dual = residual / fmax(1.0, lp.column_norm[i]);
      // A'f+g is left only where g was cut by an infinite bound (plus rounding). B does not
      // cover it: a feasible x needs (A'f+g)'x >= B, which large |x_j| can supply. Charge it at
      // the size of the current iterate (finite or not, it never marks the state bad).
      a.exposure = residual * fmax(1.0, fabs(v.out_x[i]));
    }
  } else {
    if (i < lp.n) {
      a.p = lp.c[i] * v.ray[i];
      a.absolute = fabs(a.p);
      a.dual = ConeViolation(v.ray[i], lp.lx[i], lp.ux[i]);
    }
    // Likewise per row: |A d|_i <= ||A_i||_1 for ||d||_inf <= 1.
    if (i < lp.m) 
      a.cone = ConeViolation(v.ad[i], lp.lr[i], lp.ur[i]) / fmax(1.0, lp.row_norm[i]);
  }
  if (!isfinite(a.d) || !isfinite(a.p) || !isfinite(a.dual) || !isfinite(a.cone) ||
      !isfinite(a.absolute))
    a.bad = 1;
  StoreStats(a, partial);
}

__global__ void SaveCertificate(DeviceLp lp, Vectors v, DeviceState* s, bool farkas) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (farkas && s->cert_source) {
    if (i < lp.m) v.saved_f[i] = v.f[i];
    if (i < lp.n) v.saved_g[i] = v.g[i];
  } else if (!farkas && s->cert_source && i < lp.n)
    v.saved_ray[i] = v.ray[i];
}

// mode: 0 KKT pass 1, 1 trial FPE, 2 reference FPE, 3 direction norm, 4 Farkas, 5 ray,
// 6 KKT pass 2 (decides optimality).
__global__ void ReduceStats(const Stats* partial, int size, DeviceState* s, int mode) {
  Stats a;
  for (int i = threadIdx.x; i < size; i += blockDim.x) a = Combine{}(a, partial[i]);
  __shared__ typename cub::BlockReduce<Stats, kBlock>::TempStorage storage;
  a = cub::BlockReduce<Stats, kBlock>(storage).Reduce(a, Combine{});
  if (threadIdx.x != 0) return;
  if (a.bad) s->bad = 1;
  if (mode == 0) {
    ++s->kkt_checks;
    s->spmv_count += 2;
    double denom = 1 + fabs(a.p) + fabs(a.d);
    s->metrics.primal = a.primal;
    s->metrics.dual = a.dual / (1 + s->cost_scale);
    s->metrics.gap = fabs(a.p - a.d) / denom;
    s->metrics.complementarity = a.comp / denom;
    s->metrics.primal_raw = a.raw;
    s->metrics.dual_raw = a.dual;
    s->metrics.complementarity_raw = a.comp;
    s->metrics.primal_objective = a.p;
    s->metrics.dual_objective = a.d;
    s->metrics.primal_index = a.primal_index;
    s->metrics.primal_raw_index = a.raw_index;
    s->metrics.dual_index = a.dual_index;
    s->metrics.complementarity_index = a.comp_index;
    if (!isfinite(denom) || !isfinite(s->metrics.gap) || !isfinite(s->metrics.dual) ||
        !isfinite(a.primal) || !isfinite(a.raw) || !isfinite(a.comp))
      s->bad = 1;
    if (!s->bad && a.primal <= s->tolerance) {
      s->feasible_valid = 1;
      s->feasible_iteration = s->iterations;
    }
    // lpKktCheck's objective error, on the objectives WITHOUT the constant term: the offset
    // cancels in p - d and would only enlarge the denominator (a huge constant would make any
    // gap pass). Sense does not matter: |p - d| and |p| + |d| are sign-free.
    double gap_model = 1 + fabs(a.p) + fabs(a.dual_pos);
    s->active_bound_norm = s->metrics.active_bound_norm = a.active_bound;
    s->active_cost_norm = s->metrics.active_cost_norm = a.active_cost;
    s->metrics.kkt_gap = fabs(a.p - a.dual_pos) / gap_model;
    s->metrics.l2_primal = sqrt(a.l2_primal_sq) / (1 + sqrt(a.b_sq));
    s->metrics.l2_dual = sqrt(a.l2_dual_sq) / (1 + sqrt(a.c_sq));
    // As cuPDLPx: the column part of the dual objective is z'x, not z times the bound chosen by
    // the sign of z, so a noise-level z on a huge bound does not dominate it.
    s->metrics.l2_gap = fabs(a.p - a.cross) / (1 + fabs(a.p) + fabs(a.cross));
    if (!isfinite(s->metrics.l2_gap)) s->bad = 1;
    if (!isfinite(s->metrics.kkt_gap) || !isfinite(a.dual_pos)) s->bad = 1;
    // The switch to the deviation form uses the strict per-row/column measures above.
    s->deviation_next |= !s->bad && a.primal <= kDeviationLevel &&
                         s->metrics.dual <= kDeviationLevel && s->metrics.gap <= kDeviationLevel &&
                         s->metrics.complementarity <= kDeviationLevel;
  } else if (mode == 6) {
    s->metrics.kkt_primal = a.rel_primal;
    s->metrics.kkt_dual = fmax(a.rel_dual, a.rel_residual);
    if (!isfinite(a.rel_primal) || !isfinite(s->metrics.kkt_dual)) s->bad = 1;
    const Metrics& k = s->metrics;
    s->optimal = !s->bad && (s->stop_strict ? k.primal <= s->tolerance &&
                                                  k.dual <= s->tolerance &&
                                                  k.gap <= s->tolerance &&
                                                  k.complementarity <= s->tolerance
                             : s->stop_l2   ? k.l2_primal <= s->tolerance &&
                                                  k.l2_dual <= s->tolerance &&
                                                  k.l2_gap <= s->tolerance
                                            : k.kkt_primal <= s->tolerance &&
                                                  k.kkt_dual <= s->tolerance &&
                                                  k.kkt_gap <= s->tolerance);
  } else if (mode == 1 || mode == 2) {
    ++s->spmv_count;
    double energy = s->omega * a.dx2 + a.dy2 / s->omega + 2 * kEta * a.cross;
    double roundoff =
        64 * 2.220446049250313e-16 * (s->omega * a.dx2 + a.dy2 / s->omega + 2 * kEta * a.absolute);
    if (!isfinite(energy) || energy < -roundoff) s->bad = 1;
    s->fpe = sqrt(fmax(energy, 0.0));
    if (mode == 2) {
      s->reference_fpe = s->fpe;
      s->previous_fpe = INFINITY;
    }
  } else if (mode == 3)
    s->normalization = a.norm;
  else {
    ++s->certificate_checks;
    if (!isfinite(a.d) || !isfinite(a.p) || !isfinite(a.absolute) || !isfinite(a.dual) ||
        !isfinite(a.cone))
      s->bad = 1;
    bool valid = !s->bad && s->normalization > 0;
    if (mode == 4) {
      s->spmv_count += 2;
      // A certificate of the original model must keep its margin after the stationarity
      // exposure at the current iterate; otherwise it does not even exclude points of the size
      // the iterates have. A reduced model's candidate only stops the reduced session.
      valid = valid && a.dual <= kCertificateTolerance &&
              a.d > 100 * kCertificateTolerance * fmax(1.0, a.absolute) &&
              (!s->exposure_check ||
               a.d - a.exposure > 100 * kCertificateTolerance * fmax(1.0, a.absolute));
      s->cert_source = valid && !s->farkas_valid;
      if (s->cert_source) {
        s->margin = a.d;
        s->stationarity = a.dual;
        s->cert_scale = s->normalization;
        s->farkas_info = CertificateInfo{s->direction_source, s->direction_sign,
                                         s->direction_iteration, s->normalization};
      }
      s->farkas_valid |= valid;
    } else {
      ++s->spmv_count;
      valid = valid && a.dual <= kCertificateTolerance &&
              a.cone <= kCertificateTolerance &&
              -a.p > 100 * kCertificateTolerance * fmax(1.0, a.absolute);
      s->cert_source = valid && !s->ray_valid;
      if (s->cert_source) {
        s->ray_margin = -a.p;
        s->ray_info = CertificateInfo{s->direction_source, s->direction_sign,
                                      s->direction_iteration, s->normalization};
        s->ray_row_violation = a.cone;
        s->ray_column_violation = a.dual;
      }
      s->ray_valid |= valid;
    }
  }
}

__global__ void WeightStats(DeviceLp lp, Vectors v, Stats* partial, const DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  Stats a;
  if (i < lp.n) {
    double d = s->deviation ? v.xp[i] : v.xp[i] - v.xa[i];
    a.dx2 = d * d;
  }
  if (i < lp.m) {
    double d = s->deviation ? v.yp[i] : v.yp[i] - v.ya[i];
    a.dy2 = d * d;
  }
  StoreStats(a, partial);
}

__global__ void UpdateWeight(const Stats* partial, int size, DeviceState* s) {
  Stats a;
  for (int i = threadIdx.x; i < size; i += blockDim.x) a = Combine{}(a, partial[i]);
  __shared__ typename cub::BlockReduce<Stats, kBlock>::TempStorage storage;
  a = cub::BlockReduce<Stats, kBlock>(storage).Reduce(a, Combine{});
  if (threadIdx.x != 0 || s->fixed_weight) return;
  s->balance_epoch_limit = 0;
  if (a.dx2 > 0 && a.dy2 > 0 && isfinite(a.dx2) && isfinite(a.dy2)) {
    double e = 0.5 * (log(a.dy2) - log(a.dx2)) - log(s->omega);
    s->integral = 0.3 * s->integral + e;
    double requested = 0.99 * e + 0.01 * s->integral;
    double change = Project(requested, -log(10.0), log(10.0));
    double next_log = Project(log(s->omega) + change, log(1e-8), log(1e8));
    // Hold the epoch length after one clipped update; shorten only when the next update is
    // clipped in the same direction. This distinguishes persistent imbalance from a transient
    // correction. The global weight bound and an unsaturated update clear this history.
    int direction =
        fabs(requested) >= log(10.0) && fabs(next_log - log(s->omega)) >= 0.99 * log(10.0)
            ? (requested > 0 ? 1 : -1)
            : 0;
    if (direction) {
      int64_t length = direction == s->balance_direction ? s->epoch / 2 : s->epoch;
      s->balance_epoch_limit = length > kMinBalanceEpoch ? length : kMinBalanceEpoch;
    }
    s->balance_direction = direction;
    s->omega = exp(next_log);
    s->best_omega = s->omega;
  } else {
    s->omega = s->best_omega;
    s->integral = 0;
    s->balance_direction = 0;
  }
}

// The new anchor is T(z); in the deviation form (from now on when requested) the deviation
// restarts at 0.
__global__ void RestartVectors(DeviceLp lp, Vectors v, const DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (s->deviation) {
    if (i < lp.n) {
      v.epoch_dx[i] = v.x[i];
      v.xa[i] += v.xp[i];
      v.x[i] = 0;
    }
    if (i < lp.m) {
      v.epoch_dy[i] = v.y[i];
      v.ya[i] += v.yp[i];
      v.y[i] = 0;
    }
    return;
  }
  if (i < lp.n) {
    v.epoch_dx[i] = v.x[i] - v.xa[i];
    v.x[i] = v.xa[i] = v.xp[i];
    if (s->deviation_next) v.x[i] = 0;
  }
  if (i < lp.m) {
    v.epoch_dy[i] = v.y[i] - v.ya[i];
    v.y[i] = v.ya[i] = v.yp[i];
    if (s->deviation_next) v.y[i] = 0;
  }
}

// The anchor terms from A'ya and A xa (`aty_anchor`, `ax_anchor`: scratch arrays).
__global__ void AnchorTerms(DeviceLp lp, Vectors v, const double* aty_anchor,
                            const double* ax_anchor, DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i == 0) s->spmv_count += 2;  // the products run at every restart; used once requested
  if (!s->deviation_next) return;
  if (i < lp.n) {
    v.ga[i] = lp.cs[i] - aty_anchor[i];
    if (!isfinite(v.ga[i])) atomicOr(&s->bad, 1);
  }
  if (i < lp.m) {
    v.lower[i] = ax_anchor[i] - lp.lrs[i];
    v.upper[i] = ax_anchor[i] - lp.urs[i];
    if (!isfinite(ax_anchor[i])) atomicOr(&s->bad, 1);
  }
}

__global__ void RestartState(DeviceState* s) {
  s->deviation |= s->deviation_next;
  ++s->restarts;
  s->epoch_iteration = s->iterations;
  s->epoch = 0;
}

__global__ void ReferenceCount(DeviceState* s) { s->spmv_count += 2; }
static __device__ void DecideImpl(DeviceState* s, Branch restart, Branch cert, bool initial,
                                  bool after_cert) {
  if (!after_cert) {
    bool do_cert =
        initial || s->iterations % 200 == 0 || s->iterations >= s->iteration_limit || Expired(s);
    SetBranch(cert, do_cert);
#ifdef HIPDLP_HOST_LOOP
    s->hold = 0;  // the check is running
#endif
    return;
  }
  // Contradictions are an error, never resolved by enumeration order.
  int conclusions = int(s->optimal) + int(s->farkas_valid) + int(s->ray_valid && s->feasible_valid);
  if (s->bad || conclusions > 1)
    s->status = Status::kNumericalError;
  else if (s->optimal)
    s->status = Status::kOptimal;
  else if (s->farkas_valid)
    s->status = Status::kInfeasible;
  else if (s->ray_valid && s->feasible_valid)
    s->status = Status::kUnbounded;
  else if (s->status == Status::kTimeLimit || Expired(s))
    s->status = Status::kTimeLimit;
  else if (s->iterations >= s->iteration_limit)
    s->status = Status::kIterationLimit;
  // A balance deadline only interrupts an epoch whose FPE has not fallen by half.
  bool trigger = !initial && s->status == Status::kRunning && s->epoch > 0 &&
                 s->iterations % kCheckInterval == 0 && s->reference_fpe > 0 &&
                 (s->iterations == kCheckInterval || s->fpe <= 0.2 * s->reference_fpe ||
                  (s->fpe <= 0.5 * s->reference_fpe && s->fpe > s->previous_fpe) ||
                  double(s->epoch) >= 0.36 * double(s->iterations) ||
                  (s->balance_epoch_limit > 0 && s->epoch >= s->balance_epoch_limit &&
                   s->fpe > 0.5 * s->reference_fpe));
  s->previous_fpe = s->fpe;
  s->restart = trigger;
  if (s->decisions && s->decision_count < kDecisionCapacity &&
      (s->iterations <= 20000 || trigger || s->iterations % 100000 == 0 ||
       s->status != Status::kRunning)) {
    s->decisions[s->decision_count++] =
        Decision{s->iterations,     s->epoch,        s->status,        trigger,
                 s->omega,          s->fpe,          s->reference_fpe, s->balance_epoch_limit,
                 s->metrics.primal, s->metrics.dual, s->metrics.gap};
  }
  SetBranch(restart, trigger);
}
__global__ void Decide(DeviceState* s, Branch restart, Branch cert, bool initial,
                       bool after_cert) {
  DecideImpl(s, restart, cert, initial, after_cert);
}
static __device__ void LoopConditionImpl(DeviceState* s, Branch loop) {
  if (s->status == Status::kRunning && s->bad) s->status = Status::kNumericalError;
  if (s->status == Status::kRunning && Expired(s)) s->status = Status::kTimeLimit;
  SetBranch(loop, s->status == Status::kRunning);
}
__global__ void LoopCondition(DeviceState* s, Branch loop) { LoopConditionImpl(s, loop); }
#ifdef HIPDLP_HOST_LOOP
// The tail of the check inside the chunk graph, skipped when the host has to launch the
// certificate or restart graph first (each of those ends with the same tail).
__global__ void DecideUnless(DeviceState* s, Branch restart, Branch cert) {
  if (!Taken(cert)) DecideImpl(s, restart, cert, false, true);
}
__global__ void LoopConditionUnless(DeviceState* s, Branch loop, Branch a, Branch b) {
  if (!Taken(a, b)) LoopConditionImpl(s, loop);
}
#endif
__global__ void Mark(DeviceState* s, int phase) {
  auto now = DeviceTime();
  if (phase == 0)
    s->scaling_end = now;
  else if (phase == 1)
    s->loop_begin = now;
  else if (phase == 2)
    s->end = now;
  else if (phase == 3)
    s->validation_ns -= now;
  else
    s->validation_ns += now;
}
// out[i] = 1-norm of line i (row of A, or row of A^T = column of A) of the unscaled matrix.
__global__ void LineNorms(int size, const int* start, const double* values, double* out,
                          DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) return;
  double norm = 0;
  for (int k = start[i]; k < start[i + 1]; ++k) norm += fabs(values[k]);
  if (!isfinite(norm)) atomicOr(&s->bad, 1);
  out[i] = norm;
}
__global__ void FillOne(double* values, int size) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) values[i] = 1;
}
__global__ void EmptyProduct(double* output, int size) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) output[i] = 0;
}

// The matrix arrays are read only (cuSPARSE takes non-const pointers). `iterative` is for
// operators executed inside the captured graph: they are preprocessed and run once before
// capture; a one-shot product needs neither.
static Spmv PrepareSpmv(Resources& r, int rows, int cols, int nnz, const int* start,
                        const int* index, const double* values, double* input, double* output,
                        double* ab, bool iterative = true) {
  Spmv op;
  op.rows = rows;
  op.empty_output = output;
  if (!nnz || !rows || !cols) return op;
  Sparse(cusparseCreateCsr(&op.matrix, rows, cols, nnz, const_cast<int*>(start),
                           const_cast<int*>(index), const_cast<double*>(values), CUSPARSE_INDEX_32I,
                           CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
  r.matrices.push_back(op.matrix);
  Sparse(cusparseCreateDnVec(&op.input, cols, input, CUDA_R_64F));
  r.dense.push_back(op.input);
  Sparse(cusparseCreateDnVec(&op.output, rows, output, CUDA_R_64F));
  r.dense.push_back(op.output);
  size_t bytes;
  Sparse(cusparseSpMV_bufferSize(r.sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, ab, op.matrix,
                                 op.input, ab + 1, op.output, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
                                 &bytes));
  // rocSPARSE rejects a buffer when it asks for none (cudaMalloc(0) gives null as well).
  op.workspace = bytes ? r.AllocateZero<char>(bytes) : nullptr;
  if (!iterative) return op;
  Sparse(cusparseSpMV_preprocess(r.sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, ab, op.matrix,
                                 op.input, ab + 1, op.output, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
                                 op.workspace));
  Sparse(cusparseSpMV(r.sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, ab, op.matrix, op.input, ab + 1,
                      op.output, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, op.workspace));
  return op;
}

static void Multiply(Resources& r, const Spmv& op, double* ab) {
  if (op.matrix)
    Sparse(cusparseSpMV(r.sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, ab, op.matrix, op.input, ab + 1,
                        op.output, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, op.workspace));
  else
    EmptyProduct<<<Blocks(op.rows), kBlock, 0, r.stream>>>(op.empty_output, op.rows);
}

#ifndef HIPDLP_HOST_LOOP
static std::vector<cudaGraphNode_t> FinishCapture(Resources& r) {
  const cudaGraphNode_t* deps;
  size_t size;
  cudaGraph_t graph;
  cudaStreamCaptureStatus status;
#if CUDART_VERSION >= 13000  // v2 was removed; the edge-data argument is unused here
  Cuda(cudaStreamGetCaptureInfo(r.stream, &status, nullptr, nullptr, &deps, nullptr, &size));
#else
  Cuda(cudaStreamGetCaptureInfo_v2(r.stream, &status, nullptr, nullptr, &deps, &size));
#endif
  std::vector<cudaGraphNode_t> result(deps, deps + size);
  Cuda(cudaGetLastError());
  Cuda(cudaStreamEndCapture(r.stream, &graph));
  return result;
}

static void StartCapture(Resources& r, cudaGraph_t graph,
                         const std::vector<cudaGraphNode_t>& deps = {}) {
  Cuda(cudaStreamBeginCaptureToGraph(r.stream, graph, deps.data(), nullptr, deps.size(),
                                     cudaStreamCaptureModeGlobal));
}

static cudaGraph_t AddConditional(Resources& r, cudaGraph_t graph,
                                  const std::vector<cudaGraphNode_t>& deps,
                                  cudaGraphConditionalHandle h, cudaGraphConditionalNodeType type,
                                  cudaGraphNode_t* node) {
  cudaGraphNodeParams p{};
  p.type = cudaGraphNodeTypeConditional;
  p.conditional.handle = h;
  p.conditional.type = type;
  p.conditional.size = 1;
#if CUDART_VERSION >= 13000
  Cuda(cudaGraphAddNode(node, graph, deps.data(), nullptr, deps.size(), &p));
#else
  Cuda(cudaGraphAddNode(node, graph, deps.data(), deps.size(), &p));
#endif
  return p.conditional.phGraph_out[0];
}
#else
// Captures what `body` enqueues on the stream into an executable graph the session owns.
template <typename Body>
static cudaGraphExec_t CaptureGraph(Resources& r, Body&& body) {
  cudaGraph_t graph;
  Cuda(cudaStreamBeginCapture(r.stream, cudaStreamCaptureModeGlobal));
  body();
  Cuda(cudaGetLastError());
  Cuda(cudaStreamEndCapture(r.stream, &graph));
  r.graphs.push_back(graph);
  cudaGraphExec_t exec;
  Cuda(cudaGraphInstantiate(&exec, graph, 0));
  r.execs.push_back(exec);
  return exec;
}
#endif

static std::string CheckOptions(const SolveOptions& options) {
  if (!std::isfinite(options.tolerance) || options.tolerance <= 0 ||
      std::isnan(options.time_limit_seconds) || options.time_limit_seconds < 0 ||
      options.iteration_limit < 0)
    return "Invalid tolerance or budget";
  return "";
}

static std::string CheckVectors(const HostLp& lp, const SolveOptions& options) {
  std::string message = CheckOptions(options);
  if (!message.empty()) return message;
  if (lp.rows < 0 || lp.columns < 0 || lp.cost.size() != size_t(lp.columns) ||
      lp.column_lower.size() != size_t(lp.columns) ||
      lp.column_upper.size() != size_t(lp.columns) || lp.row_lower.size() != size_t(lp.rows) ||
      lp.row_upper.size() != size_t(lp.rows) || (lp.sense != 1 && lp.sense != -1) ||
      !std::isfinite(lp.offset))
    return "Invalid dimensions, objective sense, or offset";
  for (double c : lp.cost)
    if (!std::isfinite(c)) return "Nonfinite cost";
  for (int64_t i = 0; i < lp.rows + int64_t(lp.columns); ++i) {
    double l = i < lp.rows ? lp.row_lower[i] : lp.column_lower[i - lp.rows];
    double u = i < lp.rows ? lp.row_upper[i] : lp.column_upper[i - lp.rows];
    if (std::isnan(l) || std::isnan(u) || l == INFINITY || u == -INFINITY || l > u)
      return "Invalid bound interval";
  }
  return "";
}

static std::string CheckMatrix(const HostLp& lp) {
  if (lp.rows < 0 || lp.columns < 0 || lp.value.size() > INT_MAX ||
      lp.row_start.size() != size_t(lp.rows) + 1 || lp.column_index.size() != lp.value.size())
    return "Invalid dimensions, objective sense, or offset";
  if (lp.row_start.front() != 0 || lp.row_start.back() != int(lp.value.size()))
    return "Invalid CSR offsets";
  for (int i = 0; i < lp.rows; ++i) {
    if (lp.row_start[i] < 0 || lp.row_start[i] > lp.row_start[i + 1] ||
        lp.row_start[i + 1] > int(lp.value.size()))
      return "Invalid CSR ordering or extent";
    int previous = -1;
    for (int k = lp.row_start[i]; k < lp.row_start[i + 1]; ++k) {
      if (lp.column_index[k] <= previous || lp.column_index[k] >= lp.columns)
        return "CSR indices must be sorted and unique";
      previous = lp.column_index[k];
    }
  }
  for (double a : lp.value)
    if (!std::isfinite(a)) return "Nonfinite matrix entry";
  return "";
}

static std::string CheckInput(const HostLp& lp, const SolveOptions& options) {
  std::string message = CheckOptions(options);
  if (message.empty()) message = CheckMatrix(lp);
  if (message.empty()) message = CheckVectors(lp, options);
  return message;
}

// --- device memory estimates ---
// Device memory one cudaMalloc takes (measured with this driver) : requests above 1 MiB are 
// rounded up to 2 MiB pages, smaller ones are carved from shared pages (64 KiB is an upper bound).
static size_t Charged(size_t bytes) {
  size_t page = bytes > (size_t(1) << 20) ? size_t(2) << 20 : size_t(64) << 10;
  return (std::max(size_t(1), bytes) + page - 1) / page * page;
}
// Bounds for the cuSPARSE 12.8 workspaces, from bufferSize queries between 300 and 90M nonzeros:
// cusparseCsr2cscEx2 (ALG1) asks for 16.13 bytes per nonzero plus < 5 KiB, cusparseSpMV
// (CSR_ALG2) for about nnz / 36.5 bytes plus 136.
static size_t TransposeWorkspace(int nnz) {
  return Charged(size_t(nnz) * 16 + size_t(nnz) / 4 + (size_t(64) << 10));
}
static size_t SpmvWorkspace(int nnz) { return Charged(size_t(nnz) / 32 + (size_t(4) << 10)); }
// A CSR with `rows` rows as UploadCsr / BuildTranspose allocate it (start, index, value).
static size_t CsrCharged(int rows, int nnz) {
  return Charged((size_t(rows) + 1) * 4) + Charged(size_t(nnz) * 4) + Charged(size_t(nnz) * 8);
}
// What the arrays above do not cover: kernel code loaded on first launch, the cuSPARSE handle and
// stream, the instantiated graph (together < 10 MiB in the sessions measured, up to Dual2_5000).
constexpr size_t kRuntimeReserve = size_t(64) << 20;
// The peak of a session given the arrays allocated before the transpose (`before`), the
// transpose's transient workspace, and everything allocated after it (`after`).
static size_t SessionPeak(size_t before, int rows, int columns, int nnz, bool transpose_borrowed,
                          size_t after) {
  if (!transpose_borrowed) before += CsrCharged(columns, nnz);
  size_t transient = transpose_borrowed || !nnz || !rows || !columns ? 0 : TransposeWorkspace(nnz);
  return before + std::max(transient, after) + kRuntimeReserve;
}

// Device memory one solve session (Run) allocates on top of the process's current use, array by
// array: the model vectors, the unscaled matrix unless borrowed, A^T and its transient conversion
// workspace unless borrowed, then the two scaled work copies, 21 column and 16 row vectors
// (scaling, iteration, anchor, check and certificate arrays, the unscaled row / column 1-norms),
// the reductions and the nine SpMV workspaces.
static size_t SolveFootprint(int rows, int columns, int nnz, bool matrix_borrowed,
                      bool transpose_borrowed) {
  size_t m = Charged(size_t(rows) * sizeof(double)), n = Charged(size_t(columns) * sizeof(double));
  size_t before = 3 * n + 2 * m + (matrix_borrowed ? 0 : CsrCharged(rows, nnz));
  size_t after = 2 * Charged(size_t(nnz) * sizeof(double)) + 21 * n + 16 * m +
                 Charged(2 * sizeof(double)) + Charged(kDecisionCapacity * sizeof(Decision)) +
                 Charged(sizeof(DeviceState)) +
                 Charged(size_t(Blocks(std::max(rows, columns))) * sizeof(Stats)) +
                 (nnz && rows && columns ? 9 * SpmvWorkspace(nnz) : 0);
  return SessionPeak(before, rows, columns, nnz, transpose_borrowed, after);
}

// What the solver needs from the device, found out at run time (not by model name).
struct DeviceInfo {
  int lanes = 32;  // threads per warp / wavefront
};
// Cooperative threads that cannot make the device-side time limit fire are useless: a stalled
// %globaltimer would let a solve run past its limit. Spin for `wait_ns` of timer time (bounded
// in reads, so a frozen timer ends the kernel too) and report how much it advanced.
__global__ void TimerProbe(unsigned long long wait_ns, unsigned long long* advanced) {
  unsigned long long t0 = DeviceTime(), t = t0;
  for (int i = 0; i < (1 << 24) && t - t0 < wait_ns; ++i) t = DeviceTime();
  *advanced = t - t0;
}
static bool TimerAdvances() {
  const unsigned long long wait_ns = 1000000;  // 1 ms
  unsigned long long* advanced;
  Cuda(cudaMalloc(&advanced, sizeof(*advanced)));
  auto begin = Clock::now();
  TimerProbe<<<1, 1>>>(wait_ns, advanced);
  unsigned long long value = 0;
  cudaError_t error = cudaMemcpy(&value, advanced, sizeof(value), cudaMemcpyDeviceToHost);
  double wall = Seconds(begin);
  cudaFree(advanced);
  Cuda(error);
  // It must advance by the wait, and not faster than the host clock says time passed.
  return value >= wait_ns && double(value) * 1e-9 <= wall;
}
// Rejects only what cannot work: no usable GPU or driver (see StatusOf), no kernel image for this
// device in the build, or a timer that does not run. Which GPUs those are is decided by the
// build (CMAKE_CUDA_ARCHITECTURES) and the machine, not by a list in the code. False fills the result.
static bool CheckDevice(SolveResult* result, DeviceInfo* info) {
  int device;
  cudaDeviceProp properties;
  Cuda(cudaGetDevice(&device));
  Cuda(cudaGetDeviceProperties(&properties, device));
  auto unsupported = [&](const std::string& message) {
    result->status = Status::kUnsupportedModel;
    result->termination_source = "DEVICE_CAPABILITY";
    result->message = message;
    return false;
  };
  std::string name = std::string(properties.name) + " (compute capability " +
                     std::to_string(properties.major) + "." + std::to_string(properties.minor) + ")";
  cudaFuncAttributes attributes;
  if (cudaFuncGetAttributes(&attributes, Begin) != cudaSuccess) {
    cudaGetLastError();
    return unsupported("This build has no code for " + name +
                       "; add its architecture to CMAKE_CUDA_ARCHITECTURES");
  }
#ifdef USE_HIP
  int wall_clock_khz = 0;
  Cuda(hipDeviceGetAttribute(&wall_clock_khz, hipDeviceAttributeWallClockRate, device));
  if (wall_clock_khz <= 0) return unsupported("No device wall clock on " + name);
  double tick_ns = 1e6 / wall_clock_khz;
  Cuda(hipMemcpyToSymbol(HIP_SYMBOL(wall_clock_ns), &tick_ns, sizeof(tick_ns)));
#endif
  if (!TimerAdvances()) return unsupported("The device timer does not advance on " + name);
  info->lanes = properties.warpSize;
  return true;
}
static void CheckSource(const SolveSource& source, const SolveOptions& options,
                        SolveResult* result) {
  const HostLp& host = *source.lp;
  result->rows = host.rows;
  result->columns = host.columns;
  if (source.matrix) {
    result->nonzeros = source.matrix->nnz;
    result->message = CheckVectors(host, options);
    if (result->message.empty() &&
        (source.matrix->rows != host.rows || source.matrix->columns != host.columns))
      result->message = "Device matrix dimensions do not match the model vectors";
  } else {
    result->nonzeros = host.value.size();
    result->message = CheckInput(host, options);
  }
  if (result->message.empty() && (source.start_primal || source.start_dual) &&
      !(source.start_primal && source.start_dual &&
        source.start_primal->size() == size_t(host.columns) &&
        source.start_dual->size() == size_t(host.rows)))
    result->message = "The start point must give a primal and a dual vector of the model size";
  if (!result->message.empty()) result->termination_source = "INPUT_VALIDATION";
}
static double SynchronizedSeconds(Resources& r, Clock::time_point since) {
  Cuda(cudaStreamSynchronize(r.stream));
  return Seconds(since);
}

static double CostScale(const HostLp& lp) {
  double scale = 0;
  for (double c : lp.cost) scale = std::max(scale, std::abs(c));
  if (lp.value.empty() || lp.columns == 0) return scale;
  std::vector<int> count(lp.columns, 0);
  for (int j : lp.column_index) ++count[j];
  auto foldable = [&](int j) {
    return count[j] == 1 && lp.cost[j] != 0 && !std::isfinite(lp.column_lower[j]) &&
           !std::isfinite(lp.column_upper[j]);
  };
  std::vector<double> folded;  // allocated at the first objective-definition row
  for (int i = 0; i < lp.rows; ++i) {
    int z = -1;
    double a = 0;
    for (int k = lp.row_start[i]; k < lp.row_start[i + 1] && z < 0; ++k)
      if (lp.value[k] != 0 && foldable(lp.column_index[k])) {
        z = lp.column_index[k];
        a = lp.value[k];
      }
    if (z < 0) continue;
    // The row binds on the side the substitution uses (FreeColumn); otherwise it is not folded.
    double y = lp.cost[z] / a, side = y > 0 ? lp.row_lower[i] : lp.row_upper[i];
    if (!std::isfinite(y) || !std::isfinite(side)) continue;
    if (folded.empty()) folded = lp.cost;
    for (int k = lp.row_start[i]; k < lp.row_start[i + 1]; ++k)
      folded[lp.column_index[k]] -= y * lp.value[k];
    folded[z] = 0;
  }
  if (folded.empty()) return scale;
  scale = 0;
  for (double c : folded)
    if (std::isfinite(c)) scale = std::max(scale, std::abs(c));
  return scale;
}

// The primal measure of the stop criterion (result.primal_valid).
static double PrimalMeasure(const Metrics& m, StopCriterion stop) {
  return stop == StopCriterion::kStrict  ? m.primal
         : stop == StopCriterion::kPdlpL2 ? m.l2_primal
                                          : m.kkt_primal;
}
// One complete GPU solve session. `source.matrix` (when given) is the unscaled constraint matrix
// of THIS problem; the scaled work copies and the transpose are made here and never alias it.
static SolveResult Run(const SolveSource& source, const SolveOptions& options, Snapshot* snapshot,
                       bool fixed_weight, MemoryMeter* shared_meter) {
  auto start = Clock::now();
  const HostLp& host = *source.lp;
  SolveResult result;
  CheckSource(source, options, &result);
  if (!result.message.empty()) return result;
  bool session_started = false;
  try {
    DeviceInfo device;
    if (!CheckDevice(&result, &device)) return result;
    const int lanes = device.lanes;
    std::optional<MemoryMeter> own_meter;
    if (!shared_meter) shared_meter = &own_meter.emplace();
    MemoryMeter& meter = *shared_meter;
    int m = host.rows, n = host.columns, nnz = int(result.nonzeros), size = std::max(m, n),
        blocks = Blocks(size);
    size_t estimate =
        SolveFootprint(m, n, nnz, source.matrix != nullptr, source.transposed != nullptr);
    size_t free_bytes, total_bytes;
    Cuda(cudaMemGetInfo(&free_bytes, &total_bytes));
    if (estimate > free_bytes) {
      result.status = Status::kOutOfMemory;
      result.termination_source = "GPU_SETUP";
      result.message = "Insufficient VRAM for FP64 model and workspaces";
      return result;
    }
    Resources r;
    auto stream = r.stream;
    Profile& profile = result.profile;
    nvtxRangePushA("UPLOAD");
    auto phase = Clock::now();
    DeviceLp lp{};
    lp.m = m;
    lp.n = n;
    lp.nnz = nnz;
    // The unscaled matrix is uploaded once (unless borrowed) and stays read-only.
    DeviceCsr uploaded;
    const DeviceCsr* unscaled = source.matrix;
    if (!unscaled) {
      uploaded = UploadCsr(host, stream, &profile);
      unscaled = &uploaded;
    }
    auto cost = host.cost;
    for (double& c : cost) c *= host.sense;
    lp.c = r.Upload(cost);
    lp.lr = r.Upload(host.row_lower);
    lp.ur = r.Upload(host.row_upper);
    lp.lx = r.Upload(host.column_lower);
    lp.ux = r.Upload(host.column_upper);
    profile.upload = SynchronizedSeconds(r, phase);
    meter.Sample();
    nvtxRangePop();
    nvtxRangePushA("TRANSPOSE");
    phase = Clock::now();
    DeviceTranspose built;
    const DeviceTranspose* transposed = source.transposed;
    if (!transposed) {
      BuildTranspose(*unscaled, r.sparse, stream, false, &built, &profile, &meter);
      transposed = &built;
    }
    profile.transpose = SynchronizedSeconds(r, phase);
    meter.Sample();
    nvtxRangePop();
    nvtxRangePushA("SCALE_COPY");
    phase = Clock::now();
    lp.row_start = const_cast<int*>(unscaled->row_start.Get());
    lp.column = const_cast<int*>(unscaled->column_index.Get());
    lp.original_a = const_cast<double*>(unscaled->value.Get());
    lp.transpose_start = const_cast<int*>(transposed->start.Get());
    lp.row = const_cast<int*>(transposed->row.Get());
    lp.original_at = const_cast<double*>(transposed->value.Get());
    // Scaling rewrites these copies in place; the copies fully overwrite the new arrays.
    lp.a = r.Allocate<double>(nnz);
    lp.at = r.Allocate<double>(nnz);
    if (nnz) {
      Cuda(cudaMemcpyAsync(lp.a, lp.original_a, size_t(nnz) * sizeof(double),
                           cudaMemcpyDeviceToDevice, stream));
      Cuda(cudaMemcpyAsync(lp.at, lp.original_at, size_t(nnz) * sizeof(double),
                           cudaMemcpyDeviceToDevice, stream));
      r.d2d_bytes += 2 * size_t(nnz) * sizeof(double);
    }
    // Initialize and the scaling passes write all of these before anything reads them.
    // LineNorms (graph prefix) writes row_norm / column_norm before any certificate reads them.
    for (double** p : {&lp.cs, &lp.lxs, &lp.uxs, &lp.dc, &lp.column_norm})
      *p = r.Allocate<double>(n);
    for (double** p : {&lp.lrs, &lp.urs, &lp.dr, &lp.row_norm}) *p = r.Allocate<double>(m);
    Vectors v{};
    // Arrays with disjoint lifetimes share storage. The per-sweep scaling factors are dead after
    // the scaling prefix; every product below is written by its SpMV right before its only
    // readers (A^T dy: FpeStats; A^T y_out: KktStats; A^T f: the Farkas kernels; A xr: DualStep;
    // A ray: the ray kernels; A^T ya, A xa: AnchorTerms), and the Farkas row multipliers f are dead
    // after SaveCertificate, before A ray is formed. Zeroed: some are SpMV warm-up inputs (see
    // below).
    double* column_scratch = r.AllocateZero<double>(n);
    double* row_scratch = r.AllocateZero<double>(m);
    lp.fc = v.atdy = v.out_aty = v.atf = column_scratch;
    lp.fr = v.axr = v.f = v.ad = row_scratch;
    // xr, y and dy are inputs of the SpMV that PrepareSpmv runs once before the graph exists
    // (a cuSPARSE warm-up), i.e. before Initialize: they are zeroed. The rest are written by
    // Initialize before anything reads them.
    for (double** p : {&v.x, &v.xp, &v.dx, &v.epoch_dx, &v.ga}) *p = r.Allocate<double>(n);
    v.xr = r.AllocateZero<double>(n);
    for (double** p : {&v.yp, &v.yr, &v.epoch_dy, &v.lower, &v.upper}) *p = r.Allocate<double>(m);
    // The anchor is the input of the anchor products' warm-up below.
    v.xa = r.AllocateZero<double>(n);
    v.ya = r.AllocateZero<double>(m);
    v.y = r.AllocateZero<double>(m);
    v.dy = r.AllocateZero<double>(m);
    for (double** p : {&v.aty, &v.out_x, &v.out_z, &v.g, &v.ray, &v.saved_g, &v.saved_ray,
                       &v.feasible})
      *p = r.AllocateZero<double>(n);
    for (double** p : {&v.out_y, &v.out_ax, &v.saved_f}) *p = r.AllocateZero<double>(m);
    std::vector<double> constants{1, 0};
    double* ab = r.Upload(constants);
    DeviceState initial{};
    initial.tolerance = options.tolerance;
    initial.cost_scale = source.cost_scale >= 0 ? source.cost_scale : CostScale(host);
    initial.iteration_limit = options.iteration_limit;
    initial.fixed_weight = fixed_weight;
    initial.exposure_check = !source.reduced;
    initial.warm_start = source.start_primal && source.start_dual;
    initial.stop_l2 = options.stop == StopCriterion::kPdlpL2;
    initial.stop_strict = options.stop == StopCriterion::kStrict;
    if (snapshot)
      initial.decisions = r.AllocateZero<Decision>(kDecisionCapacity);  // padding downloaded too
    DeviceState* s = r.Allocate<DeviceState>(1);  // overwritten by `initial` before the launch
    Stats* partial = r.AllocateZero<Stats>(blocks);
    profile.scale_copy = SynchronizedSeconds(r, phase);
    meter.Sample();
    nvtxRangePop();
    nvtxRangePushA("GPU_SETUP");
    phase = Clock::now();
    auto aty = PrepareSpmv(r, n, m, nnz, lp.transpose_start, lp.row, lp.at, v.y, v.aty, ab);
    auto axr = PrepareSpmv(r, m, n, nnz, lp.row_start, lp.column, lp.a, v.xr, v.axr, ab);
    auto atdy = PrepareSpmv(r, n, m, nnz, lp.transpose_start, lp.row, lp.at, v.dy, v.atdy, ab);
    auto out_ax =
        PrepareSpmv(r, m, n, nnz, lp.row_start, lp.column, lp.original_a, v.out_x, v.out_ax, ab);
    auto out_aty = PrepareSpmv(r, n, m, nnz, lp.transpose_start, lp.row, lp.original_at, v.out_y,
                               v.out_aty, ab);
    auto atf =
        PrepareSpmv(r, n, m, nnz, lp.transpose_start, lp.row, lp.original_at, v.f, v.atf, ab);
    auto ad = PrepareSpmv(r, m, n, nnz, lp.row_start, lp.column, lp.original_a, v.ray, v.ad, ab);
    // Anchor products, once per restart, into the scratch arrays (free between the check's last
    // reader and the reference evaluation).
    auto aty_anchor =
        PrepareSpmv(r, n, m, nnz, lp.transpose_start, lp.row, lp.at, v.ya, column_scratch, ab);
    auto ax_anchor =
        PrepareSpmv(r, m, n, nnz, lp.row_start, lp.column, lp.a, v.xa, row_scratch, ab);
    Cuda(cudaStreamSynchronize(stream));
#ifdef HIPDLP_HOST_LOOP
    // One flag per branch: written by the decision kernels, read by the host between graphs.
    enum { kLoop, kCheck, kRestart, kCert, kInitialCert, kReference, kBranches };
    unsigned* flags = r.AllocateZero<unsigned>(kBranches);
    Branch loop = flags + kLoop, check = flags + kCheck, restart = flags + kRestart,
           cert = flags + kCert, initial_cert = flags + kInitialCert,
           reference = flags + kReference;
#else
    Cuda(cudaGraphCreate(&r.graph, 0));
    cudaGraphConditionalHandle loop, check, restart, cert, initial_cert, reference;
    for (auto* h : {&loop, &check, &restart, &cert, &initial_cert, &reference})
      Cuda(cudaGraphConditionalHandleCreate(h, r.graph, 0, cudaGraphCondAssignDefault));
#endif
    auto reduce = [&](int mode) {
      ReduceStats<<<1, kBlock, 0, stream>>>(partial, blocks, s, mode);
    };
    auto kkt = [&]() {
      Mark<<<1, 1, 0, stream>>>(s, 3);
      Unscale<<<blocks, kBlock, 0, stream>>>(lp, v, s);
      Multiply(r, out_ax, ab);
      Multiply(r, out_aty, ab);
      KktStats<<<blocks, kBlock, 0, stream>>>(lp, v, partial, options.tolerance);
      reduce(0);
      KktRelativeStats<<<blocks, kBlock, 0, stream>>>(lp, v, s, partial, options.tolerance);
      reduce(6);
      SaveFeasible<<<Blocks(n), kBlock, 0, stream>>>(lp, v, s);
      Mark<<<1, 1, 0, stream>>>(s, 4);
    };
    auto fpe = [&](int mode) {
      Multiply(r, atdy, ab);
      FpeStats<<<blocks, kBlock, 0, stream>>>(lp, v, partial);
      reduce(mode);
    };
    auto anchor = [&]() {
      Multiply(r, aty_anchor, ab);
      Multiply(r, ax_anchor, ab);
      AnchorTerms<<<blocks, kBlock, 0, stream>>>(lp, v, column_scratch, row_scratch, s);
    };
    auto reference_eval = [&]() {
      ReferenceCount<<<1, 1, 0, stream>>>(s);
      Multiply(r, aty, ab);
      PrimalStep<<<Blocks(n), kBlock, 0, stream>>>(lp, v, s);
      Multiply(r, axr, ab);
      DualStep<<<blocks, kBlock, 0, stream>>>(lp, v, s, false);
      fpe(2);
    };
    auto certificates = [&]() {
      for (int source = 0; source < 3; ++source)
        for (int sign : {1, -1}) {
          CertificateDirection<<<blocks, kBlock, 0, stream>>>(lp, v, s, source, sign);
          Multiply(r, atf, ab);
          CertificateNorm<<<blocks, kBlock, 0, stream>>>(lp, v, partial, true);
          reduce(3);
          NormalizeCertificate<<<blocks, kBlock, 0, stream>>>(lp, v, s, true);
          Multiply(r, atf, ab);
          CertificateStats<<<blocks, kBlock, 0, stream>>>(lp, v, partial, true);
          reduce(4);
          SaveCertificate<<<blocks, kBlock, 0, stream>>>(lp, v, s, true);
          CertificateNorm<<<blocks, kBlock, 0, stream>>>(lp, v, partial, false);
          reduce(3);
          NormalizeCertificate<<<blocks, kBlock, 0, stream>>>(lp, v, s, false);
          Multiply(r, ad, ab);
          CertificateStats<<<blocks, kBlock, 0, stream>>>(lp, v, partial, false);
          reduce(5);
          SaveCertificate<<<blocks, kBlock, 0, stream>>>(lp, v, s, false);
        }
    };
    // Prefix includes GPU scaling and initial original-space validation.
    auto prefix = [&]() {
      Begin<<<1, 1, 0, stream>>>(s);
      FillOne<<<Blocks(m), kBlock, 0, stream>>>(lp.dr, m);
      FillOne<<<Blocks(n), kBlock, 0, stream>>>(lp.dc, n);
      for (int sweep = 0; sweep < 11; ++sweep) {
        ScalingFactors<<<WarpBlocks(m, lanes), kBlock, 0, stream>>>(
            m, lp.row_start, lp.a, lp.fr, lp.dr, sweep == 10, lanes, s);
        ScalingFactors<<<WarpBlocks(n, lanes), kBlock, 0, stream>>>(
            n, lp.transpose_start, lp.at, lp.fc, lp.dc, sweep == 10, lanes, s);
        ScaleMatrix<<<WarpBlocks(m, lanes), kBlock, 0, stream>>>(m, lp.row_start, lp.column, lp.a,
                                                                 lp.fr, lp.fc, lanes, s);
        ScaleMatrix<<<WarpBlocks(n, lanes), kBlock, 0, stream>>>(n, lp.transpose_start, lp.row, lp.at,
                                                                 lp.fc, lp.fr, lanes, s);
        ScalingGuard<<<1, 1, 0, stream>>>(s);
      }
      Initialize<<<blocks, kBlock, 0, stream>>>(lp, v, s);
      Mark<<<1, 1, 0, stream>>>(s, 0);
      LineNorms<<<Blocks(m), kBlock, 0, stream>>>(m, lp.row_start, lp.original_a, lp.row_norm, s);
      LineNorms<<<Blocks(n), kBlock, 0, stream>>>(n, lp.transpose_start, lp.original_at,
                                                  lp.column_norm, s);
      kkt();
      Decide<<<1, 1, 0, stream>>>(s, restart, initial_cert, true, false);
    };
    auto step = [&]() {
      PrimalStep<<<Blocks(n), kBlock, 0, stream>>>(lp, v, s);
      Multiply(r, axr, ab);
      DualStep<<<blocks, kBlock, 0, stream>>>(lp, v, s, true);
      Multiply(r, aty, ab);
      StepControl<<<1, 1, 0, stream>>>(s, check);
    };
    auto check_body = [&]() {
      fpe(1);
      kkt();
      Decide<<<1, 1, 0, stream>>>(s, restart, cert, false, false);
    };
    auto restart_body = [&]() {
      WeightStats<<<blocks, kBlock, 0, stream>>>(lp, v, partial, s);
      UpdateWeight<<<1, kBlock, 0, stream>>>(partial, blocks, s);
      RestartVectors<<<blocks, kBlock, 0, stream>>>(lp, v, s);
      anchor();
      RestartState<<<1, 1, 0, stream>>>(s);
      reference_eval();
    };
#ifdef HIPDLP_HOST_LOOP
    auto prefix_graph = CaptureGraph(r, prefix);
    auto initial_cert_graph = CaptureGraph(r, certificates);
    auto initial_decision = CaptureGraph(r, [&]() {
      Decide<<<1, 1, 0, stream>>>(s, restart, initial_cert, true, true);
      LoopCondition<<<1, 1, 0, stream>>>(s, reference);
    });
    auto reference_graph = CaptureGraph(r, reference_eval);
    auto loop_entry = CaptureGraph(r, [&]() {
      LoopCondition<<<1, 1, 0, stream>>>(s, loop);
      Mark<<<1, 1, 0, stream>>>(s, 1);
    });
    // Iterations are a multiple of kCheckInterval at every chunk start: a check falls due at the
    // chunk's last step, or earlier only for a terminal reason (the iterate is then Held).
    auto chunk = CaptureGraph(r, [&]() {
      for (int k = 0; k < kCheckInterval; ++k) step();
      check_body();
      DecideUnless<<<1, 1, 0, stream>>>(s, restart, cert);
      LoopConditionUnless<<<1, 1, 0, stream>>>(s, loop, cert, restart);
    });
    auto cert_graph = CaptureGraph(r, [&]() {
      certificates();
      Decide<<<1, 1, 0, stream>>>(s, restart, cert, false, true);
      LoopConditionUnless<<<1, 1, 0, stream>>>(s, loop, restart, nullptr);
    });
    auto restart_graph = CaptureGraph(r, [&]() {
      restart_body();
      LoopCondition<<<1, 1, 0, stream>>>(s, loop);
    });
    auto loop_exit = CaptureGraph(r, [&]() { Mark<<<1, 1, 0, stream>>>(s, 2); });
#else
    StartCapture(r, r.graph);
    prefix();
    auto deps = FinishCapture(r);
    cudaGraphNode_t cert_node, reference_node, loop_node, check_node, restart_node;
    auto cert_graph =
        AddConditional(r, r.graph, deps, initial_cert, cudaGraphCondTypeIf, &cert_node);
    StartCapture(r, cert_graph);
    certificates();
    FinishCapture(r);
    StartCapture(r, r.graph, {cert_node});
    Decide<<<1, 1, 0, stream>>>(s, restart, initial_cert, true, true);
    LoopCondition<<<1, 1, 0, stream>>>(s, reference);
    deps = FinishCapture(r);
    auto reference_graph =
        AddConditional(r, r.graph, deps, reference, cudaGraphCondTypeIf, &reference_node);
    StartCapture(r, reference_graph);
    reference_eval();
    FinishCapture(r);
    StartCapture(r, r.graph, {reference_node});
    LoopCondition<<<1, 1, 0, stream>>>(s, loop);
    Mark<<<1, 1, 0, stream>>>(s, 1);
    deps = FinishCapture(r);
    auto body = AddConditional(r, r.graph, deps, loop, cudaGraphCondTypeWhile, &loop_node);
    StartCapture(r, body);
    step();
    deps = FinishCapture(r);
    auto check_graph = AddConditional(r, body, deps, check, cudaGraphCondTypeIf, &check_node);
    StartCapture(r, check_graph);
    check_body();
    deps = FinishCapture(r);
    auto nested_cert = AddConditional(r, check_graph, deps, cert, cudaGraphCondTypeIf, &cert_node);
    StartCapture(r, nested_cert);
    certificates();
    FinishCapture(r);
    StartCapture(r, check_graph, {cert_node});
    Decide<<<1, 1, 0, stream>>>(s, restart, cert, false, true);
    deps = FinishCapture(r);
    auto restart_graph =
        AddConditional(r, check_graph, deps, restart, cudaGraphCondTypeIf, &restart_node);
    StartCapture(r, restart_graph);
    restart_body();
    FinishCapture(r);
    StartCapture(r, body, {check_node});
    LoopCondition<<<1, 1, 0, stream>>>(s, loop);
    FinishCapture(r);
    StartCapture(r, r.graph, {loop_node});
    Mark<<<1, 1, 0, stream>>>(s, 2);
    FinishCapture(r);
    Cuda(cudaGraphInstantiate(&r.exec, r.graph, 0));
#endif
    initial.remaining_seconds = std::max(0.0, options.time_limit_seconds - Seconds(start));
    Cuda(cudaMemcpyAsync(s, &initial, sizeof(initial), cudaMemcpyHostToDevice, stream));
    r.h2d_bytes += sizeof(initial);
    if (initial.warm_start) {
      // The start point travels through the output arrays: nothing reads them before the
      // initial check overwrites them, so a warm start allocates nothing.
      Cuda(cudaMemcpyAsync(v.out_x, source.start_primal->data(), size_t(n) * sizeof(double),
                           cudaMemcpyHostToDevice, stream));
      Cuda(cudaMemcpyAsync(v.out_y, source.start_dual->data(), size_t(m) * sizeof(double),
                           cudaMemcpyHostToDevice, stream));
      r.h2d_bytes += (size_t(n) + size_t(m)) * sizeof(double);
    }
    Cuda(cudaEventCreate(&r.event_start));
    Cuda(cudaEventCreate(&r.event_end));
    Cuda(cudaEventRecord(r.event_start, stream));
    meter.Sample();
    nvtxRangePop();
    profile.prepare = Seconds(phase);
    result.times.setup = Seconds(start);
    nvtxRangePushA("SOLVE_BEGIN/END");
    session_started = true;
#ifdef HIPDLP_HOST_LOOP
    // The host takes the branches the conditional nodes take on the device (see the #else
    // graph): one wait per check, one more after certificates or a restart.
    unsigned taken[kBranches];
    auto launch = [&](cudaGraphExec_t exec) { Cuda(cudaGraphLaunch(exec, stream)); };
    auto read = [&]() {
      Cuda(cudaMemcpyAsync(taken, flags, sizeof(taken), cudaMemcpyDeviceToHost, stream));
      Cuda(cudaStreamSynchronize(stream));
      r.d2h_bytes += sizeof(taken);
    };
    launch(prefix_graph);
    read();
    if (taken[kInitialCert]) launch(initial_cert_graph);
    launch(initial_decision);
    read();
    if (taken[kReference]) launch(reference_graph);
    launch(loop_entry);
    read();
    while (taken[kLoop]) {
      launch(chunk);
      read();
      if (taken[kCert]) {
        launch(cert_graph);
        read();
      }
      if (taken[kRestart]) {
        launch(restart_graph);
        read();
      }
    }
    launch(loop_exit);
    Cuda(cudaEventRecord(r.event_end, stream));
    Cuda(cudaStreamSynchronize(stream));
#else
    Cuda(cudaGraphLaunch(r.exec, stream));
    Cuda(cudaEventRecord(r.event_end, stream));
    // The only wait in a session. No copies, callbacks, polling or allocation above it.
    Cuda(cudaStreamSynchronize(stream));
#endif
    nvtxRangePop();
    float gpu_ms;
    Cuda(cudaEventElapsedTime(&gpu_ms, r.event_start, r.event_end));
    result.times.gpu = gpu_ms * 1e-3;
    auto download_start = Clock::now();
    nvtxRangePushA("DOWNLOAD");
    DeviceState final;
    Cuda(cudaMemcpy(&final, s, sizeof(final), cudaMemcpyDeviceToHost));
    result.status = final.status;
    result.iterations = final.iterations;
    result.restarts = final.restarts;
    result.kkt_checks = final.kkt_checks;
    result.certificate_checks = final.certificate_checks;
    result.spmv_count = final.spmv_count;
    result.metrics = final.metrics;
    result.peak_device_bytes = r.Bytes();
    result.metrics_valid = !final.bad;
    result.primal_valid =
        !final.bad && PrimalMeasure(final.metrics, options.stop) <= options.tolerance;
    result.dual_valid = final.status == Status::kOptimal;
    result.certificate_valid = final.farkas_valid;
    result.ray_valid = final.ray_valid;
    result.certificate_margin = final.margin;
    result.certificate_stationarity = final.stationarity;
    result.certificate_scale = final.cert_scale;
    result.ray_margin = final.ray_margin;
    result.ray_row_violation = final.ray_row_violation;
    result.ray_column_violation = final.ray_column_violation;
    result.farkas_info = final.farkas_info;
    result.ray_info = final.ray_info;
    result.feasible_iteration = final.feasible_iteration;
    r.d2h_bytes += sizeof(final);
    auto download = [&](double* p, int length) {
      std::vector<double> values(length);
      if (length) {
        Cuda(cudaMemcpy(values.data(), p, length * sizeof(double), cudaMemcpyDeviceToHost));
        r.d2h_bytes += length * sizeof(double);
      }
      return values;
    };
    result.primal = download(v.out_x, n);
    result.dual = download(v.out_y, m);
    result.reduced_cost = download(v.out_z, n);
    result.row_activity = download(v.out_ax, m);
    if (final.farkas_valid) {
      result.farkas_row = download(v.saved_f, m);
      result.farkas_column = download(v.saved_g, n);
    }
    if (final.ray_valid) result.ray = download(v.saved_ray, n);
    if (final.feasible_valid) result.feasible_base = download(v.feasible, n);
    if (snapshot) {
      snapshot->row_multiplier = download(lp.dr, m);
      snapshot->column_multiplier = download(lp.dc, n);
      snapshot->scaled_value = download(lp.a, nnz);
      snapshot->xa = download(v.xa, n);
      snapshot->ya = download(v.ya, m);
      // Absolute iterates, as the Snapshot documents them.
      auto absolute = [&](std::vector<double> deviation, const std::vector<double>& anchor) {
        if (final.deviation)
          for (size_t i = 0; i < deviation.size(); ++i) deviation[i] += anchor[i];
        return deviation;
      };
      snapshot->x = absolute(download(v.x, n), snapshot->xa);
      snapshot->y = absolute(download(v.y, m), snapshot->ya);
      snapshot->xp = absolute(download(v.xp, n), snapshot->xa);
      snapshot->yp = absolute(download(v.yp, m), snapshot->ya);
      snapshot->epoch_dx = download(v.epoch_dx, n);
      snapshot->epoch_dy = download(v.epoch_dy, m);
      snapshot->epoch_iteration = final.epoch_iteration;
      snapshot->omega = final.omega;
      snapshot->fpe = final.fpe;
      snapshot->decisions.resize(final.decision_count);
      if (final.decision_count)
        Cuda(cudaMemcpy(snapshot->decisions.data(), final.decisions,
                        final.decision_count * sizeof(Decision), cudaMemcpyDeviceToHost));
    }
    for (double& y : result.dual) y *= host.sense;
    for (double& z : result.reduced_cost) z *= host.sense;
    result.metrics.primal_objective = final.metrics.primal_objective * host.sense + host.offset;
    result.metrics.dual_objective = final.metrics.dual_objective * host.sense + host.offset;
    result.times.scaling = double(final.scaling_end - final.begin) * 1e-9;
    result.times.loop = double(final.end - final.loop_begin) * 1e-9;
    result.times.device = double(final.end - final.begin) * 1e-9;
    result.times.validation = double(final.validation_ns) * 1e-9;
    result.times.download = Seconds(download_start);
    profile.h2d_bytes += r.h2d_bytes;
    profile.d2h_bytes += r.d2h_bytes;
    profile.d2d_bytes += r.d2d_bytes;
    profile.peak_live_bytes = meter.Peak();
    nvtxRangePop();
    if (result.status == Status::kNumericalError)
      result.message = "Nonfinite arithmetic, invalid FPE, or contradictory conclusions";
  } catch (const DeviceFailure& e) {
    result.status = e.status;
    result.termination_source = session_started ? "CUDA_RUNTIME" : "GPU_SETUP";
    result.message = e.what();
  } catch (const std::bad_alloc&) {
    result.status = Status::kOutOfMemory;
    result.termination_source = "HOST_SETUP";
    result.message = "Host allocation failed";
  }
  result.times.total = Seconds(start);
  result.times.overshoot = std::max(0.0, result.times.total - options.time_limit_seconds);
  return result;
}

SolveResult Solve(const HostLp& lp, const SolveOptions& options) {
  return Run(SolveSource{&lp, nullptr}, options, nullptr, false, nullptr);
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
