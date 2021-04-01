/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsInt.h
 * @brief The definition for the integer type to use
 * @author Leona Gottwald
 */

#ifndef UTIL_HIGHS_INT_H_
#define UTIL_HIGHS_INT_H_

#include <stdint.h>

#include "HConfig.h"

#ifdef HIGHSINT64
typedef int64_t HighsInt;
typedef uint64_t HighsUInt;
#else
typedef int HighsInt;
typedef unsigned int HighsUInt;
#endif

#endif