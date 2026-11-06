/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.cpp
 * @brief
 * @author Julian Hall, Yanyu Zhou
 */
#include "pdlp/HiPdlpWrapper.h"

#include <algorithm>
#include <limits>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "pdlp/hipdlp/solver.h"

HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object) {
  return solveLpHiPdlp(solver_object.options_, solver_object.timer_,
                       solver_object.lp_, solver_object.basis_,
                       solver_object.solution_, solver_object.model_status_,
                       solver_object.highs_info_, solver_object.callback_);
}

namespace {

// HighsLp (column-wise) -> hipdlp::HostLp (CSR, sorted unique columns)
hipdlp::HostLp toHostLp(const HighsLp& lp) {
  // Work on a copy: lp is const. Going through column-wise first means the
  // row-wise copy always has ascending column indices within each row,
  // whatever format lp.a_matrix_ is in.
  HighsSparseMatrix matrix = lp.a_matrix_;
  matrix.ensureColwise();
  matrix.ensureRowwise();
  const HighsInt num_nz = matrix.start_[lp.num_row_];

  hipdlp::HostLp host;
  host.rows = lp.num_row_;
  host.columns = lp.num_col_;
  host.row_start.assign(matrix.start_.begin(),
                        matrix.start_.begin() + lp.num_row_ + 1);
  host.column_index.assign(matrix.index_.begin(),
                           matrix.index_.begin() + num_nz);
  host.value.assign(matrix.value_.begin(), matrix.value_.begin() + num_nz);
  host.cost = lp.col_cost_;
  host.row_lower = lp.row_lower_;
  host.row_upper = lp.row_upper_;
  host.column_lower = lp.col_lower_;
  host.column_upper = lp.col_upper_;
  // ObjSense::kMinimize = 1, kMaximize = -1: the same convention as HostLp
  host.sense = static_cast<int>(lp.sense_);
  host.offset = lp.offset_;
  return host;
}

hipdlp::SolveOptions toSolveOptions(const HighsOptions& options,
                                    HighsTimer& timer) {
  hipdlp::SolveOptions solve_options;
  // The solver has a single tolerance. A user-set kkt_tolerance wins.
  solve_options.tolerance = options.kkt_tolerance != kDefaultKktTolerance
                                ? options.kkt_tolerance
                                : options.pdlp_optimality_tolerance;
  // What decides "optimal" (see hipdlp::StopCriterion in hipdlp/solver.h):
  //   kPdlpL2      the PDLP / cuPDLPx relative-L2 measures (the library
  //   default) kLpKktCheck  the per-component measures of HiGHS lpKktCheck
  //   kStrict      the strict element-wise measures
  // Change this line to use another one.
  solve_options.stop = hipdlp::StopCriterion::kPdlpL2;
  // Remaining share of the global time limit
  solve_options.time_limit_seconds =
      options.time_limit >= kHighsInf
          ? std::numeric_limits<double>::infinity()
          : std::max(0.0, options.time_limit - timer.read());
  solve_options.iteration_limit = options.pdlp_iteration_limit >= kHighsIInf
                                      ? std::numeric_limits<int64_t>::max()
                                      : options.pdlp_iteration_limit;
  return solve_options;
}

HighsStatus solveError(HighsModelStatus& model_status) {
  model_status = HighsModelStatus::kSolveError;
  return HighsStatus::kError;
}

hipdlp::SolveResult solveHost(const HighsOptions& options,
                              const hipdlp::HostLp& host,
                              const hipdlp::SolveOptions& solve_options) {
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (options.hipdlp_backend != kHiPdlpBackendCpuString)
    return hipdlp::Solve(host, solve_options);
#else
  if (options.hipdlp_backend == kHiPdlpBackendCudaString ||
      options.hipdlp_backend == kHiPdlpBackendRocmString)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "HiPDLP: hipdlp_backend=%s needs a GPU build; using the CPU "
                 "solver\n",
                 options.hipdlp_backend.c_str());
#endif
  return hipdlp::SolveCpu(host, solve_options);
}

}  // namespace

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  resetModelStatusAndHighsInfo(model_status, highs_info);
  highs_basis.valid = false;
  highs_solution.clear();

  const hipdlp::HostLp host = toHostLp(lp);
  const hipdlp::SolveResult result =
      solveHost(options, host, toSolveOptions(options, timer));

  highs_info.pdlp_iteration_count =
      static_cast<HighsInt>(std::min<int64_t>(result.iterations, kHighsIInf));
  // Name the backend, and report GPU time only when there was a GPU
  const long long iterations = static_cast<long long>(result.iterations);
  if (result.backend == "CPU")
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "HiPDLP (CPU): %s after %lld iterations; total %.2fs\n",
                 hipdlp::StatusName(result.status), iterations,
                 result.times.total);
  else
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "HiPDLP (%s): %s after %lld iterations; GPU time %.2fs, "
                 "total %.2fs\n",
                 result.backend.c_str(), hipdlp::StatusName(result.status),
                 iterations, result.times.gpu, result.times.total);

  switch (result.status) {
    case hipdlp::Status::kOptimal:
      model_status = HighsModelStatus::kOptimal;
      break;
    case hipdlp::Status::kInfeasible:
      model_status = HighsModelStatus::kInfeasible;
      break;
    case hipdlp::Status::kUnbounded:
      model_status = HighsModelStatus::kUnbounded;
      break;
    case hipdlp::Status::kTimeLimit:
      model_status = HighsModelStatus::kTimeLimit;
      break;
    case hipdlp::Status::kIterationLimit:
      model_status = HighsModelStatus::kIterationLimit;
      break;
    default:
      // kInvalidModel, kUnsupportedModel (no usable GPU, or none the build
      // has code for), kOutOfMemory, kNumericalError, kDeviceError
      highsLogUser(options.log_options, HighsLogType::kError,
                   "HiPDLP failed with status %s: %s\n",
                   hipdlp::StatusName(result.status), result.message.c_str());
      return solveError(model_status);
  }

  if (result.primal.size() != static_cast<size_t>(lp.num_col_) ||
      result.row_activity.size() != static_cast<size_t>(lp.num_row_) ||
      result.reduced_cost.size() != static_cast<size_t>(lp.num_col_) ||
      result.dual.size() != static_cast<size_t>(lp.num_row_)) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiPDLP returned a solution of the wrong size\n");
    return solveError(model_status);
  }
  // Duals are already in the HiGHS convention (original objective sense,
  // reduced cost = cost - A^T y); the point is the last checked iterate.
  highs_solution.col_value = result.primal;
  highs_solution.row_value = result.row_activity;
  highs_solution.col_dual = result.reduced_cost;
  highs_solution.row_dual = result.dual;
  highs_solution.value_valid = true;
  highs_solution.dual_valid = result.status == hipdlp::Status::kOptimal;
  return HighsStatus::kOk;
}
