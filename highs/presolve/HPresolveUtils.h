/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_HIGHS_PRESOLVE_UTILS_H_
#define PRESOLVE_HIGHS_PRESOLVE_UTILS_H_

#include "lp_data/HConst.h"

namespace presolve {

enum class SingletonRowResult {
  kRedundant,
  kPrimalInfeasible,
  kBoundsTightened,
};

SingletonRowResult computeSingletonRowBounds(
    double val, double rowLower, double rowUpper, double colLower,
    double colUpper, double primalFeastol, double maxAbsColVal, bool isIntegral,
    double& lb, double& ub, bool& lowerTightened, bool& upperTightened);

}  // namespace presolve
#endif
