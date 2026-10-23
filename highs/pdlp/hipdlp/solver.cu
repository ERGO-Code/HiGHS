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



}