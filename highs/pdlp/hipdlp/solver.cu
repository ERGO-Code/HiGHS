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

__global__ void ProjectDual(DeviceLp lp, Vectors v, DeviceState* s) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < lp.m) {
    if (!isfinite(v.out_y[i])) atomicOr(&s->bad, 1);
    v.out_y[i] = Multiplier(v.out_y[i], lp.lr[i], lp.ur[i]);
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