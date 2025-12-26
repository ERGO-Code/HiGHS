#ifndef HIPO_STATUS_H
#define HIPO_STATUS_H

#include <map>
#include <string>

#include "ipm/ipx/ipx_status.h"

namespace hipo {

// Status is used both as return value for intermediate functions and as final
// status of the solver.

enum Status {
  // used only to return positive status
  kStatusOk,

  // Stopped status: solver did not converge and does not have errors
  kStatusNotRun,
  kStatusNoProgress,
  kStatusMaxIter,
  kStatusTimeLimit,
  kStatusUserInterrupt,

  // Failed status: solver has some error or interrupt
  kStatusFailed,
  kStatusError,
  kStatusOoM,
  kStatusErrorAnalyse,
  kStatusErrorFactorise,
  kStatusErrorSolve,
  kStatusBadModel,
  kStatusUnknown,

  // Solved status: solver found a solution
  kStatusSolved,
  kStatusPrimalInfeasible,
  kStatusDualInfeasible,
  kStatusImprecise,
  kStatusPDFeas,
  kStatusBasic
};

inline Status IpxToHipoStatus(Int ipx_status) {
  static const std::map<Int, Status> status_map{
      {IPX_STATUS_not_run, kStatusNotRun},
      {IPX_STATUS_iter_limit, kStatusMaxIter},
      {IPX_STATUS_no_progress, kStatusNoProgress},
      {IPX_STATUS_failed, kStatusError},
      {IPX_STATUS_time_limit, kStatusTimeLimit},
      {IPX_STATUS_user_interrupt, kStatusUserInterrupt},
      {IPX_STATUS_primal_infeas, kStatusPrimalInfeasible},
      {IPX_STATUS_dual_infeas, kStatusDualInfeasible},
      {IPX_STATUS_optimal, kStatusPDFeas},
      {IPX_STATUS_imprecise, kStatusImprecise}};

  auto found = status_map.find(ipx_status);
  if (found != status_map.end()) return found->second;
  return kStatusUnknown;
}

inline std::string statusString(Status status) {
  static const std::map<Status, std::string> status_map{
      {kStatusNotRun, "not run"},
      {kStatusMaxIter, "max iterations"},
      {kStatusNoProgress, "no progress"},
      {kStatusImprecise, "imprecise"},
      {kStatusError, "internal error"},
      {kStatusOoM, "out of memory"},
      {kStatusErrorAnalyse, "error in analyse phase"},
      {kStatusErrorFactorise, "error in factorise phase"},
      {kStatusErrorSolve, "error in solve phase"},
      {kStatusBadModel, "invalid model"},
      {kStatusTimeLimit, "time limit"},
      {kStatusUserInterrupt, "user interrupt"},
      {kStatusPrimalInfeasible, "primal infeasible"},
      {kStatusDualInfeasible, "dual infeasible"},
      {kStatusPDFeas, "primal-dual feasible"},
      {kStatusBasic, "crossover optimal"}};

  auto found = status_map.find(status);
  if (found != status_map.end()) return found->second;
  return "unknown";
}

}  // namespace hipo

#endif