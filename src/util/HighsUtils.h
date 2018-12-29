/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSUTILS_H_
#define UTIL_HIGHSUTILS_H_

#include "HConfig.h"
#include <vector>

/**
 * @brief Logical check of double being +Infinity
 */
bool highs_isInfinity(
		      double val //!< Value being tested against +Infinity
		      );
#ifdef HiGHSDEV
void util_anVecV(const char* message, int vecDim, std::vector<double>& vec, bool anVLs);
#endif
#endif // UTIL_HIGHSUTILS_H_
