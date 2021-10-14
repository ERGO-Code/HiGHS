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
/**@file simplex/HVector.h
 * @brief Vector structure for HiGHS
 */
#ifndef SIMPLEX_HVECTOR_BASE_H_
#define SIMPLEX_HVECTOR_BASE_H_

#include <cmath>
#include <map>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsInt.h"

// using std::map;
using std::vector;

/**
 * @brief Class for the vector structure for HiGHS
 */
template <typename Real>
class HVectorBase {
 public:
  /**
   * @brief Initialise a vector
   */
  void setup(HighsInt size_  //!< Dimension of the vector to be initialised
  ) {
    /*
     * Initialise an HVector instance
     */
    size = size_;
    count = 0;
    index.resize(size);
    array.assign(size, Real{0});
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
    synthetic_tick = 0;
    next = 0;
  }

  /**
   * @brief Clear the vector
   *
   */
  void clear() {
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
    synthetic_tick = 0;
    // Initialise the next value
    next = 0;
  }

  HighsInt size;           //!< Dimension of the vector
  HighsInt count;          //!< Number of nonzeros
  vector<HighsInt> index;  //!< Packed indices of nonzeros
  vector<Real> array;      //!< Full-length array of values

  double synthetic_tick;  //!< Synthetic clock for operations with this vector

  // For update
  vector<char> cwork;       //!< char working buffer for UPDATE
  vector<HighsInt> iwork;   //!< integer working buffer for UPDATE
  HVectorBase<Real>* next;  //!< Allows vectors to be linked for PAMI

  /*
   * Zero values in Vector.array that exceed kHighsTiny in magnitude
   */
  void tight() {
    /*
     * Zero values in Vector.array that do not exceed kHighsTiny in magnitude
     */
    HighsInt totalCount = 0;
    using std::abs;
    for (HighsInt i = 0; i < count; i++) {
      const HighsInt my_index = index[i];
      const double value = array[my_index];
      if (abs(value) >= kHighsTiny) {
        index[totalCount++] = my_index;
      } else {
        array[my_index] = 0;
      }
    }
    count = totalCount;
  }

  /**
   * @brief Packing (if packFlag set): Pack values/indices in Vector.array into
   * packValue/Index
   *
   */
  void pack() {
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

  bool packFlag;               //!< Flag to indicate whether to pack or not
  HighsInt packCount;          //!< Number of nonzeros packed
  vector<HighsInt> packIndex;  //!< Packed indices
  vector<Real> packValue;      //!< Packed values

  /**
   * @brief Copy from another HVector structure to this instanc
   */
  template <typename FromReal>
  void copy(const HVectorBase<FromReal>*
                from  //!< Source of HVector structure to be copied
  ) {
    /*
     * Copy from another HVector structure to this instance
     * The real type of "from" does not need to be the same, but must be
     * convertible to this HVector's real type.
     */
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

  /**
   * @brief Compute the squared 2-norm of the vector
   */
  Real norm2() {
    /*
     * Compute the squared 2-norm of the vector
     */
    const HighsInt workCount = count;
    const HighsInt* workIndex = &index[0];
    const double* workArray = &array[0];

    Real result = Real{0};
    for (HighsInt i = 0; i < workCount; i++) {
      Real value = workArray[workIndex[i]];
      result += value * value;
    }
    return result;
  }

  /**
   * @brief Add a multiple pivotX of *pivot into this vector,
   * maintaining indices of nonzeros but not tracking cancellation
   */
  template <typename RealPivX, typename RealPiv>
  void saxpy(const RealPivX pivotX,  //!< The multiple of *pivot to be added
             const HVectorBase<RealPiv>*
                 pivot  //!< The vector whose multiple is to be added
  ) {
    /*
     * Add a multiple pivotX of *pivot into this vector, maintaining
     * indices of nonzeros but not tracking cancellation.
     * The real types may all be different but must mix in operations and be
     * convertible to this HVector's real type.
     */
    HighsInt workCount = count;
    HighsInt* workIndex = &index[0];
    Real* workArray = &array[0];

    const HighsInt pivotCount = pivot->count;
    const HighsInt* pivotIndex = &pivot->index[0];
    const RealPiv* pivotArray = &pivot->array[0];

    using std::abs;
    for (HighsInt k = 0; k < pivotCount; k++) {
      const HighsInt iRow = pivotIndex[k];
      const Real x0 = workArray[iRow];
      const Real x1 = Real(x0 + pivotX * pivotArray[iRow]);
      if (x0 == Real{0}) workIndex[workCount++] = iRow;
      workArray[iRow] = (abs(x1) < kHighsTiny) ? kHighsZero : x1;
    }
    count = workCount;
  }

  bool isEqual(HVectorBase<Real>& v0) {
    if (this->size != v0.size) return false;
    if (this->count != v0.count) return false;
    if (this->index != v0.index) return false;
    if (this->array != v0.array) return false;
    //  if (this->index.size() != v0.index.size()) return false;
    //  for (HighsInt el = 0; el < (HighsInt)this->index.size(); el++)
    //    if (this->index[el] != v0.index[el]) return false;
    if (this->synthetic_tick != v0.synthetic_tick) return false;
    return true;
  }
};

#endif /* SIMPLEX_HVECTOR_H_ */
