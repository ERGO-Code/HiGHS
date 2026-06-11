/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_RKO_H_
#define HIGHS_RKO_H_

#include "lp_data/HighsLp.h"

bool rkoHeuristic(const HighsLp* lp, std::vector<double>& solution);

#endif
