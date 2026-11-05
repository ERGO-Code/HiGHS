#pragma once

#include "device.h"
#include "solver.h"
#include "solver_test.h"

namespace hipdlp {

// What a solve or check reads as ITS problem: the vectors and sense/offset of `lp`, and the
// unscaled constraint matrix either borrowed from the device (`matrix`, which must outlive the
// call, be complete on the device and match lp.rows/columns) or, when null, uploaded from lp.
// `matrix` is the full original model for the original-model check and certificates, and the
// reduced model for a reduced solve: the caller decides which, nothing here can tell.
struct SolveSource {
  const HostLp* lp = nullptr;
  const DeviceCsr* matrix = nullptr;
  // Optional A^T of `matrix` (same lifetime and contract); built on the device when null.
  const DeviceTranspose* transposed = nullptr;
  // A presolved model. Its INFEASIBLE never reaches a caller, it only sends the hybrid to the
  // original model, so its sessions stop on a Farkas candidate without the exposure condition
  // that an original-model certificate must also pass (see ReduceStats).
  bool reduced = false;
  // Cost norm of the dual residual (see CostScale); negative: CostScale(*lp). A presolved model
  // passes the original model's value, so both sessions use the same dual tolerance.
  double cost_scale = -1;
  // Optional start point of the session's model: its unscaled primal x and minimized-form dual
  // y, e.g. what an earlier session of the same model returned. Both or neither; null: zero.
  const std::vector<double>* start_primal = nullptr;
  const std::vector<double>* start_dual = nullptr;
};

struct DeviceLp {
  int m, n, nnz;
  int *row_start, *column, *transpose_start, *row;
  double *a, *at, *original_a, *original_at;
  double *c, *lr, *ur, *lx, *ux, *cs, *lrs, *urs, *lxs, *uxs;
  double *dr, *dc, *fr, *fc;
  // 1-norms of the rows / columns of the unscaled matrix: the scale of each certificate residual.
  double *row_norm, *column_norm;
};

// Below kDeviationLevel the iterate is kept as the anchor (xa, ya) plus a deviation: x, y, xp, yp,
// xr, yr are then deviations from it (DeviceState::deviation). Near a solution the steps are far
// below the precision of the anchor's entries; formed from absolute values they lose those digits
// to cancellation at every step, and under the reflection that error accumulates over an epoch
// (2026-10-14 tail experiments). In the deviation form a step only involves small numbers and the
// anchor terms ga = c - A'ya and the row bound residuals lower = A xa - l, upper = A xa - u, which
// are rounded once per epoch.
struct Vectors {
  double *x, *y, *xp, *yp, *xr, *yr, *xa, *ya, *dx, *dy;
  double *ga, *lower, *upper;
  double *epoch_dx, *epoch_dy;
  double *aty, *axr, *atdy;
  double *out_x, *out_y, *out_z, *out_ax, *out_aty;
  double *f, *g, *atf, *ray, *ad;
  double *saved_f, *saved_g, *saved_ray, *feasible;
};

struct DeviceState {
  Status status;
  int bad, feasible_valid, optimal, farkas_valid, ray_valid, restart;
  int fixed_weight, cert_source, exposure_check;
  int deviation, deviation_next;  // deviation form active / requested from the next restart
  int hold;  // HIPDLP_HOST_LOOP: a check is due, the iterate waits for the end of the chunk
  int64_t iterations, epoch, restarts, kkt_checks, certificate_checks, spmv_count, iteration_limit;
  double tolerance, remaining_seconds, omega, best_omega, integral;
  double cost_scale;  // dual residual denominator 1 + cost_scale (CostScale)
  int warm_start;  // Initialize starts from out_x / out_y (scaled by 1/dc, 1/dr), not from zero
  int stop_l2;     // OPTIMAL is decided by the L2 measures (StopCriterion::kPdlpL2)
  int stop_strict; // ... by the strict four measures (StopCriterion::kStrict)
  double active_bound_norm, active_cost_norm;  // of the last check (lpKktCheck definitions)
  double fpe, reference_fpe, previous_fpe, normalization;
  double margin, stationarity, cert_scale;
  double ray_margin, ray_row_violation, ray_column_violation;
  Metrics metrics;
  unsigned long long begin, scaling_end, loop_begin, end, validation_ns;
  Decision* decisions;
  int decision_count;
  int direction_source, direction_sign;
  int64_t direction_iteration, epoch_iteration;
  int64_t balance_epoch_limit;  // 0: no extra deadline; otherwise local to the next epoch.
  int balance_direction;
  CertificateInfo farkas_info, ray_info;
  int64_t feasible_iteration;
};

// lpKktCheck pass 1 and the PDLP / cuPDLPx L2 measures (Metrics::l2_*): squared violation and
// residual sums and the squared norms of b and c they are divided by. `cross` carries cuPDLPx's
// dual objective z'x + b-part (l2_gap).
struct KktSums {
  double p = 0, d = 0, dual_pos = 0, cross = 0;
  double l2_primal_sq = 0, l2_dual_sq = 0, b_sq = 0, c_sq = 0;
  double active_bound = 0, active_cost = 0;
  double primal = 0, dual = 0, comp = 0, raw = 0;
  int primal_index = -1, raw_index = -1, dual_index = -1, comp_index = -1;
  int bad = 0;
  __device__ void Merge(const KktSums& b) {
    p += b.p;
    d += b.d;
    dual_pos += b.dual_pos;
    cross += b.cross;
    l2_primal_sq += b.l2_primal_sq;
    l2_dual_sq += b.l2_dual_sq;
    b_sq += b.b_sq;
    c_sq += b.c_sq;
    active_bound = fmax(active_bound, b.active_bound);
    active_cost = fmax(active_cost, b.active_cost);
    if (b.primal > primal) {
      primal = b.primal;
      primal_index = b.primal_index;
    }
    if (b.dual > dual) {
      dual = b.dual;
      dual_index = b.dual_index;
    }
    if (b.comp > comp) {
      comp = b.comp;
      comp_index = b.comp_index;
    }
    if (b.raw > raw) {
      raw = b.raw;
      raw_index = b.raw_index;
    }
    bad |= b.bad;
  }
};

// lpKktCheck pass 2 (the relative failures).
struct KktRelSums {
  double rel_primal = 0, rel_dual = 0, rel_residual = 0;
  int bad = 0;
  __device__ void Merge(const KktRelSums& b) {
    rel_primal = fmax(rel_primal, b.rel_primal);
    rel_dual = fmax(rel_dual, b.rel_dual);
    rel_residual = fmax(rel_residual, b.rel_residual);
    bad |= b.bad;
  }
};

// Fixed-point error of a step (trial and reference).
struct FpeSums {
  double dx2 = 0, dy2 = 0, cross = 0, absolute = 0;
  int bad = 0;
  __device__ void Merge(const FpeSums& b) {
    dx2 += b.dx2;
    dy2 += b.dy2;
    cross += b.cross;
    absolute += b.absolute;
    bad |= b.bad;
  }
};

// Squared movement of the primal and dual iterates, for the primal weight.
struct WeightSums {
  double dx2 = 0, dy2 = 0;
  __device__ void Merge(const WeightSums& b) {
    dx2 += b.dx2;
    dy2 += b.dy2;
  }
};

// Infinity norm of a certificate direction.
struct NormMax {
  double norm = 0;
  int bad = 0;
  __device__ void Merge(const NormMax& b) {
    norm = fmax(norm, b.norm);
    bad |= b.bad;
  }
};

struct FarkasSums {
  double d = 0, absolute = 0, dual = 0;
  double exposure = 0;  // sum |A'f+g|_j * max(1, |x_j|) over the current iterate x
  int bad = 0;
  __device__ void Merge(const FarkasSums& b) {
    d += b.d;
    absolute += b.absolute;
    exposure += b.exposure;
    dual = fmax(dual, b.dual);
    bad |= b.bad;
  }
};

struct RaySums {
  double p = 0, absolute = 0, dual = 0, cone = 0;
  int bad = 0;
  __device__ void Merge(const RaySums& b) {
    p += b.p;
    absolute += b.absolute;
    dual = fmax(dual, b.dual);
    cone = fmax(cone, b.cone);
    bad |= b.bad;
  }
};

struct Combine {
  template <class T>
  __device__ T operator()(T a, const T& b) const {
    a.Merge(b);
    return a;
  }
};

// The partial results of one block-per-slot reduction share one buffer; it is sized for the widest reduction (the KKT sums)
constexpr size_t kPartialBytes =
    std::max({sizeof(KktSums), sizeof(KktRelSums), sizeof(FpeSums), sizeof(WeightSums),
              sizeof(NormMax), sizeof(FarkasSums), sizeof(RaySums)});

struct Resources {
 public:
  Resources() {
    try {
      allocations.reserve(80);
      dense.reserve(18);
      matrices.reserve(9);
      Cuda(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
      Sparse(cusparseCreate(&sparse));
      Sparse(cusparseSetStream(sparse, stream));
      Sparse(cusparseSetPointerMode(sparse, CUSPARSE_POINTER_MODE_DEVICE));
    } catch (...) {
      if (sparse) cusparseDestroy(sparse);
      if (stream) cudaStreamDestroy(stream);
      throw;
    }
  }
  Resources(const Resources&) = delete;
  Resources& operator=(const Resources&) = delete;
  ~Resources() {
    cudaStreamCaptureStatus capture;
    if (stream && cudaStreamIsCapturing(stream, &capture) == cudaSuccess &&
        capture != cudaStreamCaptureStatusNone) {
      cudaGraph_t abandoned;
      cudaStreamEndCapture(stream, &abandoned);
    }
    if (event_start) cudaEventDestroy(event_start);
    if (event_end) cudaEventDestroy(event_end);
    if (exec) cudaGraphExecDestroy(exec);
    if (graph) cudaGraphDestroy(graph);
    for (auto e : execs) cudaGraphExecDestroy(e);
    for (auto g : graphs) cudaGraphDestroy(g);
    for (auto p : dense) cusparseDestroyDnVec(p);
    for (auto p : matrices) cusparseDestroySpMat(p);
    if (sparse) cusparseDestroy(sparse);
    for (void* p : allocations) cudaFree(p);
    if (stream) cudaStreamDestroy(stream);
  }
  // Uninitialized: only for arrays that a copy or kernel overwrites completely before anything
  // reads them. Everything else uses AllocateZero.
  template <typename T>
  T* Allocate(size_t count) {
    T* p;
    size_t bytes = std::max(size_t(1), count) * sizeof(T);
    Cuda(cudaMalloc(&p, bytes));
    try {
      allocations.push_back(p);
    } catch (...) {
      cudaFree(p);
      throw;
    }
    allocated_bytes += bytes;
    return p;
  }
  template <typename T>
  T* AllocateZero(size_t count) {
    T* p = Allocate<T>(count);
    Cuda(cudaMemsetAsync(p, 0, std::max(size_t(1), count) * sizeof(T), stream));
    return p;
  }
  template <typename T>
  T* Upload(const std::vector<T>& values) {
    T* p = Allocate<T>(values.size());
    if (!values.empty()) {
      Cuda(cudaMemcpyAsync(p, values.data(), values.size() * sizeof(T), cudaMemcpyHostToDevice,
                           stream));
      h2d_bytes += values.size() * sizeof(T);
    }
    return p;
  }
  size_t Bytes() const { return allocated_bytes; }
  uint64_t h2d_bytes = 0, d2h_bytes = 0, d2d_bytes = 0;
  cudaStream_t stream = nullptr;
  cusparseHandle_t sparse = nullptr;
  cudaGraph_t graph = nullptr;
  cudaGraphExec_t exec = nullptr;
  std::vector<cudaGraph_t> graphs;  // HIPDLP_HOST_LOOP: one graph per branch
  std::vector<cudaGraphExec_t> execs;
  cudaEvent_t event_start = nullptr, event_end = nullptr;
  std::vector<cusparseDnVecDescr_t> dense;
  std::vector<cusparseSpMatDescr_t> matrices;

 private:
  std::vector<void*> allocations;
  size_t allocated_bytes = 0;
};
struct Spmv {
  cusparseSpMatDescr_t matrix = nullptr;
  cusparseDnVecDescr_t input = nullptr, output = nullptr;
  void* workspace = nullptr;
  int rows;
  double* empty_output;
};

}  // namespace hipdlp
