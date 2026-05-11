/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HVector.h
 * @brief Vector structure for HiGHS
 */
#ifndef UTIL_HVECTOR_BASE_H_
#define UTIL_HVECTOR_BASE_H_

#include <cmath>
#include <vector>

#include "util/HPreFetch.h"
#include "util/HighsInt.h"

using std::vector;

template <typename Real>
class HVectorBase {
 public:
  void setup(HighsInt size_) {
    size = size_;
    count = 0;
    index.resize(size);
    array.assign(size, Real{0});
    cwork.assign(size + 6400, 0);
    iwork.assign(size * 4, 0);
    packCount = 0;
    packIndex.resize(size);
    packValue.resize(size);
    packFlag = false;
    synthetic_tick = 0;
    next = 0;
  }

  void clear() {
    HighsInt dense_clear = count < 0 || count > size * 0.3;
    if (dense_clear) {
      array.assign(size, Real{0});
    } else {
      for (HighsInt i = 0; i < count; i++) {
        array[index[i]] = 0;
      }
    }
    this->clearScalars();
  }

  void clearScalars() {
    this->packFlag = false;
    this->count = 0;
    this->synthetic_tick = 0;
    this->next = 0;
  }

  HighsInt size;
  HighsInt count;
  vector<HighsInt> index;
  vector<Real> array;

  double synthetic_tick;

  vector<char> cwork;
  vector<HighsInt> iwork;
  HVectorBase<Real>* next;

  void tight() {
    HighsInt totalCount = 0;
    using std::abs;
    if (count < 0) {
      for (auto& val : array)
        if (abs(val) < 1e-14) val = 0;
    } else {
      for (HighsInt i = 0; i < count; i++) {
        const HighsInt my_index = index[i];
        const Real& value = array[my_index];
        if (abs(value) >= 1e-14) {
          index[totalCount++] = my_index;
        } else {
          array[my_index] = Real{0};
        }
      }
      count = totalCount;
    }
  }

  void pack() {
    if (!packFlag) return;
    packFlag = false;
    packCount = 0;
    for (HighsInt i = 0; i < count; i++) {
      const HighsInt ipack = index[i];
      packIndex[packCount] = ipack;
      packValue[packCount] = array[ipack];
      packCount++;
    }
  }

  void reIndex() {
    if (count >= 0 && count <= size * 0.1) return;
    count = 0;
    for (HighsInt i = 0; i < size; i++)
      if ((double)array[i]) index[count++] = i;
  }

  bool packFlag;
  HighsInt packCount;
  vector<HighsInt> packIndex;
  vector<Real> packValue;

  template <typename FromReal>
  void copy(const HVectorBase<FromReal>* from) {
    clear();
    synthetic_tick = from->synthetic_tick;
    const HighsInt fromCount = count = from->count;
    const HighsInt* fromIndex = &from->index[0];
    const FromReal* fromArray = &from->array[0];
    for (HighsInt i = 0; i < fromCount; i++) {
      const HighsInt iFrom = fromIndex[i];
      const FromReal xFrom = fromArray[iFrom];
      index[i] = iFrom;
      array[iFrom] = Real(xFrom);
    }
  }

  Real norm2() const {
    const HighsInt workCount = count;
    const HighsInt* workIndex = &index[0];
    const Real* workArray = &array[0];
    Real result = Real{0};
    for (HighsInt i = 0; i < workCount; i++) {
      Real value = workArray[workIndex[i]];
      result += value * value;
    }
    return result;
  }

  template <typename RealPivX, typename RealPiv>
  void saxpy(const RealPivX pivotX, const HVectorBase<RealPiv>* pivot) {
    HighsInt workCount = count;
    HighsInt* workIndex = &index[0];
    Real* workArray = &array[0];

    const HighsInt pivotCount = pivot->count;
    const HighsInt* pivotIndex = &pivot->index[0];
    const RealPiv* pivotArray = &pivot->array[0];

    using std::abs;
    HighsInt k = 0;
    for (; k + HPC_PREFETCH_DIST < pivotCount; k++) {
      HighsInt preRow = pivotIndex[k + HPC_PREFETCH_DIST];
      HPC_PREFETCH_RD(&pivotArray[preRow]);
      HPC_PREFETCH_WR(&workArray[preRow]);
      const HighsInt iRow = pivotIndex[k];
      const Real x0 = workArray[iRow];
      const Real x1 = Real(x0 + pivotX * pivotArray[iRow]);
      if (HIGHS_UNLIKELY(x0 == Real{0})) workIndex[workCount++] = iRow;
      workArray[iRow] = (HIGHS_UNLIKELY(abs(x1) < 1e-14)) ? 1e-50 : x1;
    }
    for (; k < pivotCount; k++) {
      const HighsInt iRow = pivotIndex[k];
      const Real x0 = workArray[iRow];
      const Real x1 = Real(x0 + pivotX * pivotArray[iRow]);
      if (HIGHS_UNLIKELY(x0 == Real{0})) workIndex[workCount++] = iRow;
      workArray[iRow] = (HIGHS_UNLIKELY(abs(x1) < 1e-14)) ? 1e-50 : x1;
    }
    count = workCount;
  }

  bool isEqual(const HVectorBase<Real>& v0) {
    if (this->size != v0.size) return false;
    if (this->count != v0.count) return false;
    if (this->index != v0.index) return false;
    if (this->array != v0.array) return false;
    if (this->synthetic_tick != v0.synthetic_tick) return false;
    return true;
  }
};

#endif /* UTIL_HVECTOR_BASE_H_ */
