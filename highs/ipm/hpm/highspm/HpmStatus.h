#ifndef HIGHSPM_IPM_STATUS_H
#define HIGHSPM_IPM_STATUS_H

#include <map>
#include <string>

#include "ipm/ipx/ipx_status.h"

namespace highspm {

enum IpmStatus {
  // With these statuses, the solver can proceed with refining, if requested
  kIpmStatusNotRun,
  kIpmStatusMaxIter,
  kIpmStatusNoProgress,
  kIpmStatusImprecise,

  // With these statuses, the solver should stop and not attempt refining
  kIpmStatusStop,
  kIpmStatusError,
  kIpmStatusTimeLimit,
  kIpmStatusUserInterrupt,
  kIpmStatusPrimalInfeasible,
  kIpmStatusDualInfeasible,
  kIpmStatusUnknown,

  // Solver is optimal
  kIpmStatusOptimal,
  kIpmStatusPDFeas,
  kIpmStatusBasic
};

enum LinearSolverStatus {
  kLinearSolverStatusOk = 0,
  kLinearSolverStatusErrorOom,
  kLinearSolverStatusErrorAnalyse,
  kLinearSolverStatusErrorFactorise,
  kLinearSolverStatusErrorSolve
};

inline IpmStatus IpxToHpmStatus(Int ipx_status) {
  static const std::map<Int, IpmStatus> status_map{
      {IPX_STATUS_not_run, kIpmStatusNotRun},
      {IPX_STATUS_iter_limit, kIpmStatusMaxIter},
      {IPX_STATUS_no_progress, kIpmStatusNoProgress},
      {IPX_STATUS_failed, kIpmStatusError},
      {IPX_STATUS_time_limit, kIpmStatusTimeLimit},
      {IPX_STATUS_user_interrupt, kIpmStatusUserInterrupt},
      {IPX_STATUS_primal_infeas, kIpmStatusPrimalInfeasible},
      {IPX_STATUS_dual_infeas, kIpmStatusDualInfeasible},
      {IPX_STATUS_optimal, kIpmStatusPDFeas},
      {IPX_STATUS_imprecise, kIpmStatusImprecise}};

  auto found = status_map.find(ipx_status);
  if (found != status_map.end()) return found->second;
  return kIpmStatusUnknown;
}

inline std::string statusString(IpmStatus status) {
  static const std::map<IpmStatus, std::string> status_map{
      {kIpmStatusNotRun, "not run"},
      {kIpmStatusMaxIter, "max iterations"},
      {kIpmStatusNoProgress, "no progress"},
      {kIpmStatusError, "internal error"},
      {kIpmStatusTimeLimit, "time limit"},
      {kIpmStatusUserInterrupt, "user interrupt"},
      {kIpmStatusPrimalInfeasible, "primal infeasible"},
      {kIpmStatusDualInfeasible, "dual infeasible"},
      {kIpmStatusPDFeas, "primal-dual feasible"},
      {kIpmStatusBasic, "crossover optimal"}};

  auto found = status_map.find(status);
  if (found != status_map.end()) return found->second;
  return "unknown";
}

}  // namespace highspm

#endif