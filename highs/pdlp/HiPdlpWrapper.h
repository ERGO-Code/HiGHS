/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.h
 * @brief
 */
#ifndef PDLP_HIPDLP_WRAPPER_H_
#define PDLP_HIPDLP_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "lp_data/HighsSolution.h"
//#include "pdlp/cupdlp/cupdlp.h"

//typedef enum CONSTRAINT_TYPE { EQ = 0, LEQ, GEQ, BOUND } constraint_type;

HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object);

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback);

#endif
