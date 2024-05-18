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

#include "util/HighsMemoryAllocation.h"

#include <cstdio>

bool okUint8Resize(std::vector<uint8_t>& use_vector, const HighsInt dimension,
                   const bool value) {
  try {
    use_vector.resize(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::HighsUint8Resize fails with %s\n", e.what());
    return false;
  }
  return true;
}

bool okHighsIntResize(std::vector<HighsInt>& use_vector,
                      const HighsInt dimension, const HighsInt value) {
  printf("HighsMemoryAllocation::HighsIntAssign %d values\n", int(dimension));
  try {
    use_vector.resize(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::HighsIntResize fails with %s\n", e.what());
    return false;
  }
  if (dimension >= 46) return false;
  return true;
}

bool okHighsIntReserve(std::vector<HighsInt>& use_vector,
                       const HighsInt dimension) {
  try {
    use_vector.reserve(dimension);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::HighsIntReserve fails with %s\n", e.what());
    return false;
  }
  return true;
}

bool okHighsIntAssign(std::vector<HighsInt>& use_vector,
                      const HighsInt dimension, const HighsInt value) {
  try {
    use_vector.assign(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::HighsIntAssign fails with %s\n", e.what());
    return false;
  }
  return true;
}

bool okHighsIntSetResize(std::vector<std::set<HighsInt>>& use_vector,
                         const HighsInt dimension) {
  try {
    use_vector.resize(dimension, std::set<HighsInt>());
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::HighsIntSetResize fails with %s\n",
           e.what());
    return false;
  }
  return true;
}

bool okDoubleResize(std::vector<double>& use_vector, const HighsInt dimension,
                    const double value) {
  try {
    use_vector.resize(dimension, value);
  } catch (const std::bad_alloc& e) {
    printf("HighsMemoryAllocation::DoubleResize fails with %s\n", e.what());
    return false;
  }
  return true;
}
