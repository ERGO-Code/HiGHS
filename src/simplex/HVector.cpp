/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HVector.cpp
 * @brief
 */
#include "simplex/HVector.h"

#include <cassert>
#include <cmath>

#include "lp_data/HConst.h"
#include "stdio.h"  //Just for temporary printf

void HVector::setup(HighsInt size_) {
  /*
   * Initialise an HVector instance
   */
  size = size_;
  count = 0;
  index.resize(size);
  array.assign(size, 0);
  cwork.assign(size + 6400, 0);  // MAX invert
  iwork.assign(size * 4, 0);

  packCount = 0;
  packIndex.resize(size);
  packValue.resize(size);

  // Initialise three values that are initialised in clear(), but
  // weren't originally initialised in setup(). Probably doesn't
  // matter, since clear() is usually called before a vector is
  // (re-)used.
  packFlag = false;
  syntheticTick = 0;
  next = 0;
}

void HVector::clear() {
  /*
   * Clear an HVector instance
   */
  // Standard HVector to clear
  HighsInt clearVector_inDense = count < 0 || count > size * 0.3;
  if (clearVector_inDense) {
    // Treat the array as full if there are no indices or too many indices
    array.assign(size, 0);
  } else {
    // Zero according to the indices of (possible) nonzeros
    for (HighsInt i = 0; i < count; i++) {
      array[index[i]] = 0;
    }
  }
  // Reset the flag to indicate when to pack
  packFlag = false;
  // Zero the number of stored indices
  count = 0;
  // Zero the synthetic clock for operations with this vector
  syntheticTick = 0;
  // Initialise the next value
  next = 0;
}

void HVector::tight() {
  /*
   * Packing: Zero values in Vector.array which exceed kHighsTiny in
   * magnitude
   */
  HighsInt totalCount = 0;
  for (HighsInt i = 0; i < count; i++) {
    const HighsInt my_index = index[i];
    const double value = array[my_index];
    if (fabs(value) > kHighsTiny) {
      index[totalCount++] = my_index;
    } else {
      array[my_index] = 0;
    }
  }
  count = totalCount;
}

void HVector::pack() {
  /*
   * Packing (if packFlag set): Pack values/indices in Vector.array
   * into packValue/Index
   */
  if (packFlag) {
    packFlag = false;
    packCount = 0;
    for (HighsInt i = 0; i < count; i++) {
      const HighsInt ipack = index[i];
      packIndex[packCount] = ipack;
      packValue[packCount] = array[ipack];
      packCount++;
    }
  }
}

void HVector::copy(const HVector* from) {
  /*
   * Copy from another HVector structure to this instance
   */
  clear();
  syntheticTick = from->syntheticTick;
  const HighsInt fromCount = count = from->count;
  const HighsInt* fromIndex = &from->index[0];
  const double* fromArray = &from->array[0];
  for (HighsInt i = 0; i < fromCount; i++) {
    const HighsInt iFrom = fromIndex[i];
    const double xFrom = fromArray[iFrom];
    index[i] = iFrom;
    array[iFrom] = xFrom;
  }
}

double HVector::norm2() {
  /*
   * Compute the squared 2-norm of the vector
   */
  const HighsInt workCount = count;
  const HighsInt* workIndex = &index[0];
  const double* workArray = &array[0];

  double result = 0;
  for (HighsInt i = 0; i < workCount; i++) {
    double value = workArray[workIndex[i]];
    result += value * value;
  }
  return result;
}

void HVector::saxpy(const double pivotX, const HVector* pivot) {
  /*
   * Add a multiple pivotX of *pivot into this vector, maintaining
   * indices of nonzeros but not tracking cancellation
   */
  HighsInt workCount = count;
  HighsInt* workIndex = &index[0];
  double* workArray = &array[0];

  const HighsInt pivotCount = pivot->count;
  const HighsInt* pivotIndex = &pivot->index[0];
  const double* pivotArray = &pivot->array[0];

  for (HighsInt k = 0; k < pivotCount; k++) {
    const HighsInt iRow = pivotIndex[k];
    const double x0 = workArray[iRow];
    const double x1 = x0 + pivotX * pivotArray[iRow];
    if (x0 == 0) workIndex[workCount++] = iRow;
    workArray[iRow] = (fabs(x1) < kHighsTiny) ? kHighsZero : x1;
  }
  count = workCount;
}
