#include <cuda_runtime.h>
#include <cusparse.h>
#include <nvtx3/nvToolsExt.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cub/block/block_reduce.cuh>
#include <optional>
#include <stdexcept>
#include <utility>

#include "device.h"
#include "solver.h"
#include "solver_internal.h"
#include "solver_test.h"

namespace hipdlp {
constexpr int kBlock = 256;
constexpr double kEta = 0.99;
constexpr double kCertificateTolerance = 1e-8;
constexpr int64_t kMinBalanceEpoch = 400; 

constexpr double kDeviationLevel = 1e-5;

using Clock = std::chrono::steady_clock;
static double Seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
static int Blocks(int size) { return std::max(1, int((int64_t(size) + kBlock - 1) / kBlock)); }

static int WarpBlocks(int rows) {
  return std::max(1, int(std::min(int64_t(4096), (int64_t(rows) + 7) / 8)));
}

__device__ unsigned long long DeviceTime() {
  unsigned long long value;
  asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(value));
  return value;
}

__device__ bool Expired(const DeviceState* s) {
  return double(DeviceTime() - s->begin) * 1e-9 >= s->remaining_seconds;
}

__device__ double Project(double v, double l, double u) { return fmin(fmax(v, l), u); }

__device__ double Multiplier(double v, double l, double u) {
  if (!isfinite(l)) v = fmin(v, 0.0);
  if (!isfinite(u)) v = fmax(v, 0.0);
  return v;
}
__device__ double Recession(double v, double l, double u) {
  if (isfinite(l)) v = fmax(v, 0.0);
  if (isfinite(u)) v = fmin(v, 0.0);
  return v;
}
__device__ double BoundProduct(double v, double l, double u) {
  return v > 0 ? v * l : (v < 0 ? v * u : 0.0);
}
__device__ double ConeViolation(double v, double l, double u) {
  double error = 0;
  if (isfinite(l)) error = fmax(error, -v);
  if (isfinite(u)) error = fmax(error, v);
  return error;
}

__device__ void BoundMetrics(Stats& s, double v, double multiplier, double l, double u, int index) {
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

__device__ void KktBounds(Stats& s, double v, double mult, double l, double u, double tol) {
  if (!isfinite(v) || !isfinite(mult)) return;
  if (isfinite(l) && (l - v) * (l - v) <= tol)
    s.active_bound = fmax(s.active_bound, fabs(l));
  else if (isfinite(u) && (v - u) * (v - u) <= tol)
    s.active_bound = fmax(s.active_bound, fabs(u));
  double bound = 1;
  if (isfinite(l) || isfinite(u)) bound = v < 0.5 * (l + u) ? l : u;
  s.dual_pos += bound * mult;
}

__device__ void StoreStats(Stats value, Stats* partial) {
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
                               double* cumulative, bool pc, DeviceState* s) {
  if (s->status != Status::kRunning) return;
  int lane = threadIdx.x % 32;
  int stride = gridDim.x * blockDim.x / 32;
  for (int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32; i < size; i += stride) {
    int active = 0;
    if (lane == 0) active = s->status == Status::kRunning && !Expired(s);
    active = __shfl_sync(0xffffffff, active, 0);
    if (!active) {
      if (lane == 0) factor[i] = 1;
      continue;
    }
    double norm = 0;
    for (int k = start[i] + lane; k < start[i + 1]; k += 32) {
      double v = fabs(a[k]);
      if (!isfinite(v)) atomicOr(&s->bad, 1);
      norm = pc ? norm + v : fmax(norm, v);
    }
    for (int offset = 16; offset > 0; offset /= 2) {
      double other = __shfl_down_sync(0xffffffff, norm, offset);
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
const double* left, const double* right, DeviceState* s) {
  if (s->status != Status::kRunning) return;
  int stride = gridDim.x * blockDim.x / 32;
  for (int i = (blockIdx.x * blockDim.x + threadIdx.x) / 32; i < size; i += stride) {
    for (int k = start[i] + threadIdx.x % 32; k < start[i + 1]; k += 32) {
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

__global__ void Initialze(DeviceLp lp, Vectors v, DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < lp.n) {
    lp.cs[i] = lp.c[i] * lp.dc[i];
    lp.lxs[i] = lp.lx[i] / lp.dc[i];
    lp.uxs[i] = lp.ux[i] / lp.dc[i];
    double start = s->warm_start ? v.out_x[i] / lp.dc[i] : 0; //out_x holds the start point
    v.x[i] = v.xp[i] = v.xr = v.xa[i] = Project(start, lp.lxs[i], lp.uxs[i]);
    v.dx[i] = v.epoch_dx[i] = 0;
    if (!isfinite(lp.cs[i]) || (isfinite(lp.lxs[i]) && !isfinite(v.x[i])) || (isfinite(lp.ux[i]) && !isfinite(lp.uxs[i])) || !isfinite(v.x[i]) ||
      !isfinite(start)) atomicOr(&s->bad, 1);
    }
    if (i < lp.m) {
      lp.lrs[i] = lp.lr[i] * lp.dr[i];
      lp.urs[i] = lp.ur[i] * lp.dr[i];
      double start = s->warm_start ? v.out_y[i] / lp.dr[i] : 0.0; // already of a legal sign
      v.y[i] = v.py[i] = v.yr[i] = v.ya[i] = start;
      v.dy[i] = v.epoch_dy[i] = 0;
      if ((isfinite(lp.lr[i]) && !isfinite(lp.lrs[i])) ||
          (isfinite(lp.ur[i]) && !isfinite(lp.urs[i])) || !isfinite(start))
        atomicOr(&s->bad, 1);
  }
}

__global__ void PrimalStep(DeviceLp lp, Vectors v, DeviceState* s) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= lp.n || s->status != Status::kRunning) return;
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
  if (s->status != Status::kRunning) return;
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

__global__ void StepControl(DeviceState* s, cudaGraphConditionalHandle check) {
  if (s->status == Status::kRunning) {
    ++s->iterations;
    ++s->epoch;
    s->spmv_count += 2;
  }
  cudaGraphSetConditional(check, s->bad || s->iterations % 40 == 0 ||
                                     s->iterations >= s->iteration_limit || Expired(s));
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
__global__ void Decide(DeviceState* s, cudaGraphConditionalHandle restart,
                       cudaGraphConditionalHandle cert, bool initial, bool after_cert) {
  if (!after_cert) {
    bool do_cert =
        initial || s->iterations % 200 == 0 || s->iterations >= s->iteration_limit || Expired(s);
    cudaGraphSetConditional(cert, do_cert);
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
                 s->iterations % 40 == 0 && s->reference_fpe > 0 &&
                 (s->iterations == 40 || s->fpe <= 0.2 * s->reference_fpe ||
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
  cudaGraphSetConditional(restart, trigger);
}
__global__ void LoopCondition(DeviceState* s, cudaGraphConditionalHandle loop) {
  if (s->status == Status::kRunning && s->bad) s->status = Status::kNumericalError;
  if (s->status == Status::kRunning && Expired(s)) s->status = Status::kTimeLimit;
  cudaGraphSetConditional(loop, s->status == Status::kRunning);
}
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
  op.workspace = r.AllocateZero<char>(bytes);
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

static std::string CheckOptions(const SolveOptions& options) {
  if (!std::isfinite(options.tolerance) || options.tolerance <= 0 ||
      std::isnan(options.time_limit_seconds) || options.time_limit_seconds < 0 ||
      options.iteration_limit < 0)
    return "Invalid tolerance or budget";
  return "";
}

std::string CheckVectors(const HostLp& lp, const SolveOptions& options) {
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

std::string CheckMatrix(const HostLp& lp) {
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

std::string CheckInput(const HostLp& lp, const SolveOptions& options) {
  std::string message = CheckOptions(options);
  if (message.empty()) message = CheckMatrix(lp);
  if (message.empty()) message = CheckVectors(lp, options);
  return message;
}


}