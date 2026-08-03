/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.cpp
 * @brief
 */
#include "util/HighsSolverSelect.h"

HighsSolverSelect selectSolver(const HighsLp& lp) {

    return HighsSolverSelect::kDualSimplex;
}