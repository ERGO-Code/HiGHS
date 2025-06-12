#ifndef HIPO_INFO_H
#define HIPO_INFO_H

#include "Options.h"
#include "Status.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "ipm/ipx/info.h"

namespace hipo {

struct Info {
  // Size of problem, as seen by the solver
  Int m_solver, n_solver;

  // Size of original problem
  Int m_original, n_original;

  // Status of solver, see IpmStatus.h
  Status status = kStatusNotRun;

  // residuals and objectives of final solution
  double p_res_rel, p_res_abs, d_res_rel, d_res_abs, p_obj, d_obj, pd_gap;

  // Number of ipm iterations performed
  Int ipm_iter = 0;

  // True if ipx was invoked, whether to refine solution or for crossover
  bool ipx_used = false;

  // Info from ipx
  ipx::Info ipx_info;

  // Number of correctors used
  Int correctors;

  // Nla option used
  OptionNla option_nla;

  // Parallel option used
  OptionParallel option_par;

  // Total times to form matrix, factorise and solve linear systems
  double analyse_NE_time{};
  double analyse_AS_time{};
  double matrix_time{};
  double factor_time{};
  double solve_time{};

  // Counters
  Int factor_number{};
  Int solve_number{};

  // Information on dense columns
  Int num_dense_cols{};
  double max_col_density{};
};

}  // namespace hipo

#endif