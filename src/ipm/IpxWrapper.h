/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IPM_IPX_WRAPPER_H_
#define IPM_IPX_WRAPPER_H_

#include <algorithm>

#include "ipm/IpxSolution.h"
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsOptions.h"

IpxStatus fillInIpxData(const HighsLp& lp, ipx::Int& num_col,
                        std::vector<double>& obj, std::vector<double>& col_lb,
                        std::vector<double>& col_ub, ipx::Int& num_row,
                        std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                        std::vector<double>& Ax, std::vector<double>& rhs,
                        std::vector<char>& constraint_type);

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag);

HighsStatus reportIpxIpmCrossoverStatus(const HighsOptions& options,
                                        const ipx::Int status,
                                        const bool ipm_status);

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    std::string message, const int value = -1);

bool illegalIpxSolvedStatus(ipx::Info& ipx_info, const HighsOptions& options);

bool illegalIpxStoppedIpmStatus(ipx::Info& ipx_info,
                                const HighsOptions& options);

bool illegalIpxStoppedCrossoverStatus(ipx::Info& ipx_info,
                                      const HighsOptions& options);

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info);

HighsStatus analyseIpmNoProgress(const ipx::Info& ipx_info,
                                 const ipx::Parameters& parameters,
                                 HighsModelStatus& unscaled_model_status);

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, bool& imprecise_solution,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsIterationCounts& iteration_counts,
                       HighsModelStatus& unscaled_model_status,
                       HighsSolutionParams& unscaled_solution_params);

#endif
