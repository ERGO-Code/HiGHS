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
/**@file HighsMemoryAllocation.h
 * @brief Utilities for memory allocation that return true if successful
 */

#ifndef UTIL_HIGHS_MEMORY_ALLOCATION_H_
#define UTIL_HIGHS_MEMORY_ALLOCATION_H_

#include <set>
#include <vector>

#include "util/HighsInt.h"

bool okUint8Resize(std::vector<uint8_t>& use_vector, const HighsInt dimension,
                   const bool value = false);

bool okHighsIntResize(std::vector<HighsInt>& use_vector,
                      const HighsInt dimension, const HighsInt value = 0);

bool okHighsIntReserve(std::vector<HighsInt>& use_vector,
                       const HighsInt dimension);

bool okHighsIntAssign(std::vector<HighsInt>& use_vector,
                      const HighsInt dimension, const HighsInt value = 0);

bool okHighsIntSetResize(std::vector<std::set<HighsInt>>& use_vector,
                         const HighsInt dimension);

bool okDoubleResize(std::vector<double>& use_vector, const HighsInt dimension,
                    const double value = 0);

#endif
