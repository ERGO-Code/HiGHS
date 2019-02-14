/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HVector.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HVector.h"
#include "lp_data/HConst.h"

#include <cassert>
#include <cmath>
#include "stdio.h"  //Just for temporary printf

void HVector::setup(int size_) {
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

  pWd = dfSparseDaStr;
  packMap.clear();
  valueP1.assign(size, ilP1);
  valueP2.assign(size, ilP2);
}

void HVector::clear() {
  /*
   * Clear an HVector instance
   */
  if (pWd == dfSparseDaStr) {
    // Standard HVector to clear
    int clearVector_inDense = count < 0 || count > size * 0.3;
    if (clearVector_inDense) {
      // Treat the array as full if there are no indices or too many indices
      array.assign(size, 0);
    } else {
      // Zero according to the indices of (possible) nonzeros
      for (int i = 0; i < count; i++) {
        array[index[i]] = 0;
      }
    }
  } else if (pWd == p0SparseDaStr) {
    // Clear a value-index map data structure
    packMap.clear();
  } else if (pWd == p1SparseDaStr) {
    // Clear the 1-byte pointers to packed values
    for (int i = 0; i < count; i++) {
      valueP1[index[i]] = ilP1;
    }
  } else if (pWd == p2SparseDaStr) {
    // Clear the 2-byte pointers to packed values
    for (int i = 0; i < count; i++) {
      valueP2[index[i]] = ilP2;
    }
  }
  // Possibly check that the vector is cleared
  bool ckClear = false;
  if (ckClear) {
    for (int i = 0; i < size; i++) {
      if (array[i] != 0) {
        printf("Error: cleared array[%d]=%g\n", i, array[i]);
        fflush(stdout);
      }
      if (valueP1[i] != ilP1) {
        printf("Error: cleared valueP1[%d]=%d\n", i, valueP1[i]);
        fflush(stdout);
      }
      if (valueP2[i] != ilP2) {
        printf("Error: cleared valueP2[%d]=%d\n", i, valueP2[i]);
        fflush(stdout);
      }
    }
  }
  // Reset the flag to indicate when to pack
  packFlag = false;
  // Zero the number of stored indices
  count = 0;

  // fakeTick is soon to be deleted
  //    fakeTick = 0;

  // Zero the synthetic clock for operations with this vector
  syntheticTick = 0;

  // Initialise the next value
  next = 0;

  // Set the data structure type to the default value
  pWd = dfSparseDaStr;
}

void HVector::tight() {
  /*
   * Packing: Zero values in Vector.array which exceed HIGHS_CONST_TINY in
   * magnitude
   */
  if (pWd != dfSparseDaStr) {
    printf("ERROR: HVector::tight() not implemented for pWd=%d\n", pWd);
  }
  int totalCount = 0;
  for (int i = 0; i < count; i++) {
    const int my_index = index[i];
    const double value = array[my_index];
    if (fabs(value) > HIGHS_CONST_TINY) {
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
   * into packValue/Index when using default data structure, just
   * indices when using the 1-byte pointer ultra-sparse representation
   */
  if (packFlag) {
    packFlag = false;
    packCount = 0;
    if (pWd == dfSparseDaStr) {
      for (int i = 0; i < count; i++) {
        const int ipack = index[i];
        packIndex[packCount] = ipack;
        packValue[packCount] = array[ipack];
        packCount++;
      }
    } else if (pWd >= p1SparseDaStr) {
      for (int i = 0; i < count; i++) {
        packIndex[packCount] = index[i];
        packCount++;
      }
    }
  }
}

void HVector::copy(const HVector *from) {
  /*
   * Copy from another HVector structure to this instance
   */
  clear();
  //    fakeTick = from->fakeTick;
  syntheticTick = from->syntheticTick;
  const int fromCount = count = from->count;
  const int *fromIndex = &from->index[0];
  const double *fromArray = &from->array[0];
  const int frompWd = from->pWd;
  if (frompWd == dfSparseDaStr) {
    for (int i = 0; i < fromCount; i++) {
      const int iFrom = fromIndex[i];
      const double xFrom = fromArray[iFrom];
      index[i] = iFrom;
      array[iFrom] = xFrom;
    }
  } else if (frompWd == p1SparseDaStr) {
    const unsigned char *fromValueP1 = &from->valueP1[0];
    const double *fromPackValue = &from->packValue[0];
    for (int i = 0; i < fromCount; i++) {
      const int iFrom = fromIndex[i];
      const int valueP = fromValueP1[iFrom];
      const double xFrom = fromPackValue[valueP];
      index[i] = iFrom;
      array[iFrom] = xFrom;
    }
  } else if (frompWd == p2SparseDaStr) {
    const unsigned short *fromValueP2 = &from->valueP2[0];
    const double *fromPackValue = &from->packValue[0];
    for (int i = 0; i < fromCount; i++) {
      const int iFrom = fromIndex[i];
      const int valueP = fromValueP2[iFrom];
      const double xFrom = fromPackValue[valueP];
      index[i] = iFrom;
      array[iFrom] = xFrom;
    }
  }
}

double HVector::norm2() {
  /*
   * Compute the squared 2-norm of the vector
   */
  const int workCount = count;
  const int *workIndex = &index[0];
  const double *workArray = &array[0];

  if (pWd != dfSparseDaStr) {
    printf("ERROR: HVector::norm2() not implemented for pWd=%d\n", pWd);
  }
  double result = 0;
  for (int i = 0; i < workCount; i++) {
    double value = workArray[workIndex[i]];
    result += value * value;
  }
  return result;
}

void HVector::saxpy(const double pivotX, const HVector *pivot) {
  /*
   * Add a multiple pivotX of *pivot into this vector, maintaining
   * indices of nonzeros but not tracking cancellation
   */
  int workCount = count;
  int *workIndex = &index[0];
  double *workArray = &array[0];

  const int pivotCount = pivot->count;
  const int *pivotIndex = &pivot->index[0];
  const double *pivotArray = &pivot->array[0];

  if (pWd != dfSparseDaStr) {
    printf(
        "ERROR: HVector::saxpy(const double pivotX, const HVector *pivot) not "
        "implemented for pWd=%d\n",
        pWd);
  }
  for (int k = 0; k < pivotCount; k++) {
    const int iRow = pivotIndex[k];
    const double x0 = workArray[iRow];
    const double x1 = x0 + pivotX * pivotArray[iRow];
    if (x0 == 0) workIndex[workCount++] = iRow;
    workArray[iRow] = (fabs(x1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : x1;
  }
  count = workCount;
}
