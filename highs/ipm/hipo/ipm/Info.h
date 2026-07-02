#ifndef HIPO_INFO_H
#define HIPO_INFO_H

#include "Options.h"
#include "Status.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/ipx/info.h"

namespace hipo {

struct Info {
  // Size of problem, as seen by the solver
  Int m_solver, n_solver;

  // Size of original problem
  Int m_original, n_original;

  // Status of solver, see Status.h
  Status status = kStatusNotSet;
  Int error = kOk;

  // Status of ipm iterations up to pd feas solution found (potentially using
  // ipx to refine, in which case the ipx status_ipm is converted to hipo
  // status).
  Status status_phase1 = kStatusNotSet;

  // Status of ipm iterations after pd feas solution found, or status of
  // crossover if it is run, in which case ipx status_crossover is converted to
  // hipo status.
  Status status_phase2 = kStatusNotSet;

  // residuals and objectives of final solution
  double p_res_rel, p_res_abs, d_res_rel, d_res_abs, p_obj, d_obj, pd_gap;

  // Number of ipm iterations performed
  Int ipm_iter = 0;

  // True if ipx was used to refine solution or for crossover.
  bool ipx_used_refine = false;
  bool ipx_used_crossover = false;

  // Info from ipx
  ipx::Info ipx_info;

  // Number of correctors used
  Int correctors;

  // Total times to form matrix, factorise and solve linear systems
  double analyse_NE_time{};
  double analyse_AS_time{};
  double matrix_time{};
  double AS_structure_time{};
  double NE_structure_time{};
  double factor_time{};
  double solve_time{};
  double residual_time{};
  double omega_time{};

  // Counters
  Int factor_number{};
  Int solve_number{};
};

}  // namespace hipo

#endif