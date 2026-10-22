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

/*
HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  return hipdlp::api::solve(options, timer, lp, highs_basis, highs_solution,
                            model_status, highs_info, callback,
                            hipdlp::api::defaultMigrationOptions());
}
*/

#ifdef HIGHS_HIPDLP_CUDA_DETECTED

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

hipdlp::SolveOptions toSolveOptions(const HighsOptions& options, HighsTimer& timer) {
  hipdlp::SolveOptions solve_options;
  // The solver has a single tolerance, applied to the HiGHS lpKktCheck
  // measures. A user-set kkt_tolerance wins.
solve_options.tolerance = options.kkt_tolerance != kDefaultKktTolerance
                                ? options.kkt_tolerance
                                : options.pdlp_optimality_tolerance;
  // Remaining share of the global time limit
  solve_options.time_limit_seconds =
      options.time_limit >= kHighsInf
          ? std::numeric_limits<double>::infinity()
          : std::max(0.0, options.time_limit - timer.read());
  solve_options.iteration_limit = options.pdlp_iteration_limit >= kHighsIInf
                                      ? std::numeric_limits<int64_t>::max()
                                      : options.pdlp_iteration_limit;
  // HiGHS has already presolved: Solve() does not presolve anyway
  solve_options.presolve = false;
  return solve_options;
}

HighsStatus solveError(HighsModelStatus& model_status) {
  model_status = HighsModelStatus::kSolveError;
  return HighsStatus::kError;
}

} // namespace

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  resetModelStatusAndHighsInfo(model_status, highs_info);
  highs_basis.valid = false;
  highs_solution.clear();

  const hipdlp::HostLp host = toHostLp(lp);
  const hipdlp::SolveResult result = hipdlp::Solve(host, toSolveOptions(options, timer));

  highs_info.pdlp_iteration_count = static_cast<HighsInt>(
      std::min<int64_t>(result.iterations, kHighsIInf));
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "HiPDLP: %s after %lld iterations; GPU time %.2fs, total "
               "%.2fs\n",
               hipdlp::StatusName(result.status),
               static_cast<long long>(result.iterations), result.times.gpu,
               result.times.total);

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
      // kInvalidModel, kUnsupportedModel (not an sm_89 GPU), kOutOfMemory,
      // kNumericalError, kDeviceError
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

#else  // !HIGHS_HIPDLP_CUDA_DETECTED

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  resetModelStatusAndHighsInfo(model_status, highs_info);
  highs_basis.valid = false;
  highs_solution.clear();
  highsLogUser(options.log_options, HighsLogType::kError,
               "HiPDLP requires a CUDA build: configure with "
               "-DHIGHS_HIPDLP_CUDA=ON\n");
  model_status = HighsModelStatus::kSolveError;
  return HighsStatus::kError;
}

#endif

