/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/CupdlpWrapper.h
 * @brief
 */
#ifndef PDLP_CUPDLP_WRAPPER_H_
#define PDLP_CUPDLP_WRAPPER_H_

#include <algorithm>
#include <cassert>

//#include "ipm/IpxSolution.h"
//#include "ipm/ipx/ipx_status.h"
//#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsSolution.h"

HighsStatus solveLpCupdlp(HighsLpSolverObject& solver_object);

HighsStatus solveLpCupdlp(const HighsOptions& options,
			  HighsTimer& timer,
			  const HighsLp& lp, 
			  HighsBasis& highs_basis,
			  HighsSolution& highs_solution,
			  HighsModelStatus& model_status,
			  HighsInfo& highs_info,
			  HighsCallback& callback);
#endif
