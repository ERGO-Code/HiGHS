/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSOLUTIONDEBUG_H_
#define SIMPLEX_HIGHSSOLUTIONDEBUG_H_

#include "lp_data/HighsLp.h"
bool equalSolutionParams(const HighsSolutionParams& solution_params0,
                         const HighsSolutionParams& solution_params1);
bool equalSolutionObjectiveParams(const HighsSolutionParams& solution_params0,
                                  const HighsSolutionParams& solution_params1);
bool equalSolutionStatusParams(const HighsSolutionParams& solution_params0,
                               const HighsSolutionParams& solution_params1);
bool equalSolutionInfeasibilityParams(
    const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);



#endif  // SIMPLEX_HIGHSSOLUTIONDEBUG_H_
