/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapperEmpty.h
 * @brief
 */
#ifndef IPM_IPX_WRAPPER_EMPTY_H_
#define IPM_IPX_WRAPPER_EMPTY_H_

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, bool& imprecise_solution,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsIterationCounts& iteration_counts,
                       HighsModelStatus& unscaled_model_status,
                       HighsSolutionParams& unscaled_solution_params) {
  unscaled_model_status = HighsModelStatus::kNotset;
  return HighsStatus::kError;
}

#endif
