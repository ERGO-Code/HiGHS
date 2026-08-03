/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.h
 * @brief
 */
#ifndef UTIL_HIGHS_SOLVER_SELECT_H_
#define UTIL_HIGHS_SOLVER_SELECT_H_

class HighsLp;

enum class HighsSolverSelect {
    kDualSimplex,
    kPrimalSimplex,
    kIpx,
    kHipo,
    kCupdlp,
    kHipdlp
};

HighsSolverSelect selectSolver(const HighsLp& lp);

#endif /* UTIL_HIGHS_SOLVER_SELECT_H_ */