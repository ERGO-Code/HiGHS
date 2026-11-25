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
                       solver_object.highs_info_, solver_object.callback_,
                       nullptr, solver_object.hipdlp_omega_);
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
                                    HighsTimer& timer,
                                    const HiPdlpRestart* restart) {
  hipdlp::SolveOptions solve_options;
  // Stop on the tolerances that HiGHS then checks the solution with
  // (lpKktCheck): primal and dual feasibility, dual residual and optimality.
  // A user-set kkt_tolerance replaces all of them, as it does there. The L2 gap
  // keeps pdlp_optimality_tolerance, as cuPDLP-C does. (The primal residual
  // needs no test: the row activities returned are A x.)
  const bool uniform = options.kkt_tolerance != kDefaultKktTolerance;
  // A restart scales each of them (see HiPdlpRestart)
  const HiPdlpRestart none;
  const HiPdlpRestart& scales = restart ? *restart : none;
  auto tolerance = [&](double option, double scale) {
    return scale * (uniform ? options.kkt_tolerance : option);
  };
  solve_options.primal_tolerance =
      tolerance(options.primal_feasibility_tolerance, scales.primal_scale);
  solve_options.dual_tolerance =
      tolerance(options.dual_feasibility_tolerance, scales.dual_scale);
  solve_options.dual_residual_tolerance =
      tolerance(options.dual_residual_tolerance, scales.residual_scale);
  solve_options.l2_gap_tolerance =
      tolerance(options.pdlp_optimality_tolerance, scales.objective_scale);
  solve_options.kkt_gap_tolerance =
      tolerance(options.optimality_tolerance, scales.objective_scale);
  // The single tolerance of the strict criterion, the certificates and the
  // CPU solver, as before (a restart: at the tightest scale)
  solve_options.tolerance =
      tolerance(options.pdlp_optimality_tolerance,
                std::min({scales.primal_scale, scales.dual_scale,
                          scales.residual_scale, scales.objective_scale}));
  // What decides "optimal" (see hipdlp::StopCriterion in hipdlp/solver.h):
  // kPdlpL2, the PDLP / cuPDLPx relative-L2 measures (the library default);
  // kLpKktCheck, the per-component measures of HiGHS lpKktCheck; kStrict, the
  // strict element-wise measures. Change this line to use another one.
  solve_options.stop = hipdlp::StopCriterion::kPdlpL2;
  // The L2 stop also needs the lpKktCheck measures within 10 times their
  // tolerances, so that it cannot end far outside what lpKktCheck accepts
  solve_options.kkt_guard = 10;
  // A restart aims at what lpKktCheck measures
  if (restart) solve_options.stop = hipdlp::StopCriterion::kLpKktCheck;
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
                              const hipdlp::SolveOptions& solve_options,
                              const hipdlp::StartPoint* start) {
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (options.hipdlp_backend != kHiPdlpBackendCpuString)
    return start ? hipdlp::Solve(host, solve_options, *start)
                 : hipdlp::Solve(host, solve_options);
#else
  if (options.hipdlp_backend == kHiPdlpBackendCudaString ||
      options.hipdlp_backend == kHiPdlpBackendRocmString)
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "HiPDLP: hipdlp_backend=%s needs a GPU build; using the CPU "
                 "solver\n",
                 options.hipdlp_backend.c_str());
#endif

  return start ? hipdlp::SolveCpu(host, solve_options, *start)
               : hipdlp::SolveCpu(host, solve_options);
}

}  // namespace

HiPdlpKktViolations hiPdlpKktViolations(const HighsOptions& options,
                                        const HighsInfo& info) {
  // The largest KKT error relative to its tolerance, from the HighsInfo of
  // the last lpKktCheck of a solution without a basis: lpKktCheck changes
  // "Optimal" to "Unknown" when this exceeds 100
  const bool uniform = options.kkt_tolerance != kDefaultKktTolerance;
  auto relative = [&](const double error, const double tolerance) {
    return error / (uniform ? options.kkt_tolerance : tolerance);
  };
  HiPdlpKktViolations violations;
  // HiPDLP has one tolerance for primal infeasibilities and residuals
  violations.primal = std::max(relative(info.max_relative_primal_infeasibility,
                                        options.primal_feasibility_tolerance),
                               relative(info.max_relative_primal_residual_error,
                                        options.primal_residual_tolerance));
  violations.dual = relative(info.max_relative_dual_infeasibility,
                             options.dual_feasibility_tolerance);
  violations.dual_residual = relative(info.max_relative_dual_residual_error,
                                      options.dual_residual_tolerance);
  violations.objective =
      relative(info.primal_dual_objective_error, options.optimality_tolerance);
  return violations;
}

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback, const HiPdlpRestart* restart,
                          double* final_omega) {
  resetModelStatusAndHighsInfo(model_status, highs_info);
  highs_basis.valid = false;
  if (final_omega) *final_omega = 0;
  // Copy the start before highs_solution is cleared: it may be the same object.
  // HiPDLP's dual is the minimized-form y: the HiGHS row dual times the sense.
  std::vector<double> start_primal, start_dual;
  hipdlp::StartPoint start;
  if (restart && restart->start && restart->start->value_valid &&
      restart->start->col_value.size() == static_cast<size_t>(lp.num_col_)) {
    const HighsSolution& from = *restart->start;
    start_primal = from.col_value;
    start_dual.assign(lp.num_row_, 0.0);
    if (from.dual_valid &&
        from.row_dual.size() == static_cast<size_t>(lp.num_row_)) {
      const double sense = static_cast<double>(lp.sense_);
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
        start_dual[iRow] = sense * from.row_dual[iRow];
    }
    start.primal = &start_primal;
    start.dual = &start_dual;
    start.omega = restart->omega;
  }
  highs_solution.clear();

  const hipdlp::HostLp host = toHostLp(lp);
  const hipdlp::SolveOptions solve_options =
      toSolveOptions(options, timer, restart);
  if (restart)
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "HiPDLP: restart %s with tolerances %g (primal), %g (dual), "
                 "%g (dual residual), %g (objective) and the KKT-check stop\n",
                 start.primal ? "from the previous solution" : "from zero",
                 solve_options.primal_tolerance, solve_options.dual_tolerance,
                 solve_options.dual_residual_tolerance,
                 solve_options.kkt_gap_tolerance);
  const hipdlp::SolveResult result =
      solveHost(options, host, solve_options, start.primal ? &start : nullptr);
  if (final_omega) *final_omega = result.omega;

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
