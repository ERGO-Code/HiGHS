#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hipdlp {

enum class Status {
  kRunning,
  kOptimal,
  kInfeasible,
  kUnbounded,
  kTimeLimit,
  kIterationLimit,
  kInvalidModel,
  kUnsupportedModel,
  kOutOfMemory,
  kNumericalError,
  kDeviceError
};
const char* StatusName(Status status);

// Immutable input to Solve: CSR with explicit IEEE infinities for absent bounds.
// Duplicate entries must be combined by the caller. sense is +1 (min), -1 (max).
struct HostLp {
  int rows = 0, columns = 0;
  std::vector<int> row_start, column_index;
  std::vector<double> value, cost, row_lower, row_upper, column_lower, column_upper;
  int sense = 1;
  double offset = 0;
};
// What decides OPTIMAL (SolveOptions::stop). kPdlpL2 (default): the standard PDLP / cuPDLPx
// relative-L2 measures, Metrics::l2_*, so a run stops where cuPDLPx stops. kLpKktCheck: the HiGHS
// lpKktCheck measures, Metrics::kkt_*, the per-component checks HiGHS applies to a solution
// itself; they are stricter than L2 when a few components violate their bounds by a lot.
// kStrict: the strict element-wise four measures (primal, dual, gap, complementarity) that decided
// OPTIMAL before 2026-10-14, for runs comparable with that history.
// All are evaluated on the same point and always reported.
enum class StopCriterion { kPdlpL2, kLpKktCheck, kStrict };
struct SolveOptions {
  double tolerance = 1e-4;
  double time_limit_seconds = std::numeric_limits<double>::infinity();
  int64_t iteration_limit = std::numeric_limits<int64_t>::max();
  StopCriterion stop = StopCriterion::kPdlpL2;
};
struct Metrics {
  // Strict element-wise report: each bound violation over 1 + |violated bound|, the dual residual
  // over 1 + max|c|, gap and complementarity over 1 + |p| + |d|. It decides OPTIMAL only
  // with StopCriterion::kStrict (the certificates' feasible base still uses `primal`).
  double primal = 0, dual = 0, gap = 0, complementarity = 0;
  double primal_raw = 0, dual_raw = 0, complementarity_raw = 0;
  double primal_objective = 0, dual_objective = 0;
  int primal_index = -1, primal_raw_index = -1, dual_index = -1, complementarity_index = -1;
  // primal/complementarity: rows first, then columns. dual: column index.
  // The optimality gate of StopCriterion::kLpKktCheck: the measures of HiGHS lpKktCheck for a
  // solution without a basis, at the session tolerance (relative primal and dual infeasibility,
  // relative dual residual, relative primal-dual objective error) with the same active bound /
  // active cost norms. With that criterion OPTIMAL iff kkt_primal, kkt_dual and kkt_gap are all
  // <= tolerance (lpKktCheck itself keeps "Optimal" up to 100x that).
  double kkt_primal = 0, kkt_dual = 0, kkt_gap = 0;
  double active_bound_norm = 0, active_cost_norm = 0;
  // The standard PDLP / cuPDLPx termination measures of the same point, the gate of
  // StopCriterion::kPdlpL2 and reported for the others, all in unscaled coordinates without the
  // objective offset: l2_primal is the L2 norm of
  // the row and column bound violations over 1 + ||b||_2 (b: each row's lower bound when finite
  // and different from the upper one, its upper bound when finite), l2_dual the L2 norm of
  // c - A^T y - z over 1 + ||c||_2 and l2_gap is |p - d| / (1 + |p| + |d|) with cuPDLPx's dual
  // objective d: each row's bound by the sign of its multiplier plus z'x for the columns. As in
  // cuPDLPx (whose z is the primal projection's slack) z is position-based here: nonzero only on a
  // bound the value is within sqrt(tolerance) of, with the sign that bound allows; the rest of
  // c - A^T y counts as dual residual. cuPDLPx stops when all three are <= its tolerance.
  double l2_primal = 0, l2_dual = 0, l2_gap = 0;
};
struct Timings {
  double read = 0, presolve = 0, setup = 0, scaling = 0, loop = 0;
  double validation = 0, download = 0, postsolve = 0, total = 0, overshoot = 0;
  double gpu = 0, device = 0;
  // Input work after `read` that is charged to the solve budget and included in `total`.
  double pack = 0;
  // Wall clock of the dedicated original-model KKT check of a hybrid solve (setup + one bounded
  // GPU execution + download); 0 when no check ran. Included in `total`, not added to another
  // field.
  double final_check = 0;
};
// Where the GPU preparation time went and what crossed the PCIe bus or was copied on the
// device. Times are wall clocks of stream-synchronized phases (they do not add up to `total`:
// presolve scans, solve sessions and validation are reported elsewhere). Bytes count every
// copy made by the solver, including presolve batches and the final downloads.
struct Profile {
  double upload = 0;      // H2D of matrices and vectors, including device allocation
  double transpose = 0;   // CSR -> CSC on the device (cusparseCsr2cscEx2)
  double scale_copy = 0;  // D2D work copies of the unscaled matrices and zero fills
  double prepare = 0;     // cuSPARSE buffers/preprocess and CUDA graph construction
  double pack = 0;        // reduced CSR construction on the device
  uint64_t h2d_bytes = 0, d2h_bytes = 0, d2d_bytes = 0;
  // Reduced matrices built on the CPU and uploaded because the device had no staging memory for
  // the device construction (never for device faults, which are reported as errors).
  int pack_cpu_fallbacks = 0;
  size_t transpose_workspace_bytes = 0;  // largest cusparseCsr2cscEx2 workspace queried
  // Device memory in use by this process (cudaMemGetInfo) above the level at entry, sampled
  // at the end of every setup step; allocations of all stages that were alive together.
  size_t peak_live_bytes = 0;
};
struct CertificateInfo {
  // 0: current-anchor; 1: current-reflected; 2: saved epoch displacement;
  // 3: initial row violation.
  int source = -1, sign = 0;
  int64_t iteration = 0;  // Local to the GPU session that generated it.
  double scale = 0;       // Original-coordinate infinity norm before normalization.
};
// Where the input time went; unused fields stay zero. See docs for the stage boundaries.
struct InputReport {
  std::string reader;  // native
  bool pipelined = false;
  uint64_t file_bytes = 0, text_bytes = 0;
  double scan = 0, drain = 0, finalize = 0, inflate = 0, wait = 0;
  uint64_t model_hash = 0;  // LpFingerprint of the original model when requested
};
struct SolveResult {
  Status status = Status::kInvalidModel;
  std::string message, termination_source = "GPU", fallback_reason;
  // Hybrid: reduced sessions rerun at a tighter tolerance after the original model rejected the
  // recovered solution (the original is solved instead only when that does not help).
  int reduced_retries = 0;
  std::string metrics_model = "ORIGINAL";
  bool metrics_valid = false;
  bool reduced_metrics_valid = false;
  int rows = 0, columns = 0, nonzeros = 0;
  int reduced_rows = 0, reduced_columns = 0, reduced_nonzeros = 0;
  int64_t iterations = 0, restarts = 0, kkt_checks = 0, certificate_checks = 0, spmv_count = 0;
  bool primal_valid = false, dual_valid = false, certificate_valid = false, ray_valid = false;
  std::vector<double> primal, dual, reduced_cost, row_activity;
  std::vector<double> farkas_row, farkas_column, ray, feasible_base;
  CertificateInfo farkas_info, ray_info;
  int64_t feasible_iteration = 0;
  // Certificate residuals are relative: stationarity is max_j |A'f+g|_j / max(1, ||A_j||_1) and
  // ray_row_violation is max_i cone violation of (A d)_i / max(1, ||A_i||_1), for the normalized
  // certificate; a certificate is accepted when each is at most 1e-8.
  double certificate_margin = 0, certificate_stationarity = 0, certificate_scale = 0;
  double ray_margin = 0, ray_row_violation = 0, ray_column_violation = 0;
  Metrics metrics, reduced_metrics;
  Timings times;
  InputReport input;
  Profile profile;
  size_t peak_device_bytes = 0;  // Bytes allocated by the largest single solve/validation
};
SolveResult Solve(const HostLp& lp, const SolveOptions& options = {});

}  // namespace hipdlp
