/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.h
 * @brief
 */
#ifndef IPM_IPX_WRAPPER_H_
#define IPM_IPX_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "ipm/IpxSolution.h"
#include "ipm/ipx/ipx_status.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsSolution.h"

#ifdef HPM
#include "ipm/hpm/highspm/HpmSolver.h"
#endif

HighsStatus solveLpIpx(HighsLpSolverObject& solver_object);

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, HighsBasis& highs_basis,
                       HighsSolution& highs_solution,
                       HighsModelStatus& model_status, HighsInfo& highs_info,
                       HighsCallback& callback);

#ifdef HPM
HighsStatus solveLpHpm(HighsLpSolverObject& solver_object);

HighsStatus solveLpHpm(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, HighsBasis& highs_basis,
                       HighsSolution& highs_solution,
                       HighsModelStatus& model_status, HighsInfo& highs_info,
                       HighsCallback& callback);

#endif

void fillInIpxData(const HighsLp& lp, ipx::Int& num_col, ipx::Int& num_row,
                   double& offset, std::vector<double>& obj,
                   std::vector<double>& col_lb, std::vector<double>& col_ub,
                   std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                   std::vector<double>& Ax, std::vector<double>& rhs,
                   std::vector<char>& constraint_type);

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag);

HighsStatus reportIpxIpmCrossoverStatus(const HighsOptions& options,
                                        const ipx::Int status,
                                        const bool ipm_status);

HighsStatus reportHpmStatus(const HighsOptions& options,
                            const highspm::Int status,
                            const highspm::HpmSolver& hpm);

HighsStatus reportHpmCrossoverStatus(const HighsOptions& options,
                                     const ipx::Int status);

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    std::string solver, std::string message,
                    const int value = -1);

bool illegalIpxSolvedStatus(const ipx::Info& ipx_info,
                            const HighsOptions& options);

bool illegalIpxStoppedIpmStatus(const ipx::Info& ipx_info,
                                const HighsOptions& options);

bool illegalIpxStoppedCrossoverStatus(const ipx::Info& ipx_info,
                                      const HighsOptions& options);

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info);

void reportHpmNoProgress(const HighsOptions& options,
                         const highspm::HpmInfo& hpm_info);

void getHighsNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
                               const ipx::Int num_col, const ipx::Int num_row,
                               const std::vector<double>& rhs,
                               const std::vector<char>& constraint_type,
                               const ipx::LpSolver& lps,
                               const HighsModelStatus model_status,
                               HighsSolution& highs_solution);

void getHpmNonVertexSolution(
    const HighsOptions& options, const HighsLp& lp, const highspm::Int num_col,
    const highspm::Int num_row, const std::vector<double>& rhs,
    const std::vector<char>& constraint_type, const highspm::HpmSolver& hpm,
    const HighsModelStatus model_status, HighsSolution& highs_solution);

void reportSolveData(const HighsLogOptions& log_options,
                     const ipx::Info& ipx_info);
#endif
