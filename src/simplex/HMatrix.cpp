/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HMatrix.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HMatrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "lp_data/HConst.h"
#include "simplex/HVector.h"

using std::fabs;
using std::max;
using std::swap;

void HMatrix::setup(HighsInt numCol_, HighsInt numRow_, const HighsInt* Astart_,
                    const HighsInt* Aindex_, const double* Avalue_,
                    const int8_t* nonbasicFlag_) {
  // Copy the A matrix and setup row-wise matrix with the nonbasic
  // columns before the basic columns for a general set of nonbasic
  // variables
  //
  // Copy A
  numCol = numCol_;
  numRow = numRow_;
  Astart.assign(Astart_, Astart_ + numCol_ + 1);

  HighsInt AcountX = Astart_[numCol_];
  Aindex.assign(Aindex_, Aindex_ + AcountX);
  Avalue.assign(Avalue_, Avalue_ + AcountX);

  // Build row copy - pointers
  std::vector<HighsInt> AR_Bend;
  ARstart.resize(numRow + 1);
  AR_Nend.assign(numRow, 0);
  AR_Bend.assign(numRow, 0);
  // Count the nonzeros of nonbasic and basic columns in each row
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    if (nonbasicFlag_[iCol]) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        AR_Nend[iRow]++;
      }
    } else {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        AR_Bend[iRow]++;
      }
    }
  }
  ARstart[0] = 0;
  for (HighsInt i = 0; i < numRow; i++)
    ARstart[i + 1] = ARstart[i] + AR_Nend[i] + AR_Bend[i];
  for (HighsInt i = 0; i < numRow; i++) {
    AR_Bend[i] = ARstart[i] + AR_Nend[i];
    AR_Nend[i] = ARstart[i];
  }
  // Build row copy - elements
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    if (nonbasicFlag_[iCol]) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        HighsInt iPut = AR_Nend[iRow]++;
        ARindex[iPut] = iCol;
        ARvalue[iPut] = Avalue[k];
      }
    } else {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        HighsInt iPut = AR_Bend[iRow]++;
        ARindex[iPut] = iCol;
        ARvalue[iPut] = Avalue[k];
      }
    }
  }
  // Initialise the density of the PRICE result
  //  row_apDensity = 0;
}

void HMatrix::setup_lgBs(HighsInt numCol_, HighsInt numRow_,
                         const HighsInt* Astart_, const HighsInt* Aindex_,
                         const double* Avalue_) {
  // Copy the A matrix and setup row-wise matrix with the nonbasic
  // columns before the basic columns for a logical basis
  //
  // Copy A
  numCol = numCol_;
  numRow = numRow_;
  Astart.assign(Astart_, Astart_ + numCol_ + 1);

  HighsInt AcountX = Astart_[numCol_];
  Aindex.assign(Aindex_, Aindex_ + AcountX);
  Avalue.assign(Avalue_, Avalue_ + AcountX);

  // Build row copy - pointers
  ARstart.resize(numRow + 1);
  AR_Nend.assign(numRow, 0);
  for (HighsInt k = 0; k < AcountX; k++) AR_Nend[Aindex[k]]++;
  ARstart[0] = 0;
  for (HighsInt i = 1; i <= numRow; i++)
    ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
  for (HighsInt i = 0; i < numRow; i++) AR_Nend[i] = ARstart[i];

  // Build row copy - elements
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      HighsInt iRow = Aindex[k];
      HighsInt iPut = AR_Nend[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }
  // Initialise the density of the PRICE result
  //  row_apDensity = 0;
}

void HMatrix::update(HighsInt variable_in, HighsInt variable_out) {
  if (variable_in < numCol) {
    for (HighsInt k = Astart[variable_in]; k < Astart[variable_in + 1]; k++) {
      HighsInt iRow = Aindex[k];
      HighsInt iFind = ARstart[iRow];
      HighsInt iSwap = --AR_Nend[iRow];
      while (ARindex[iFind] != variable_in) iFind++;
      // todo @ Julian : this assert can fail
      assert(iFind >= 0 && iFind < int(ARindex.size()));
      assert(iSwap >= 0 && iSwap < int(ARindex.size()));
      swap(ARindex[iFind], ARindex[iSwap]);
      swap(ARvalue[iFind], ARvalue[iSwap]);
    }
  }

  if (variable_out < numCol) {
    for (HighsInt k = Astart[variable_out]; k < Astart[variable_out + 1]; k++) {
      HighsInt iRow = Aindex[k];
      HighsInt iFind = AR_Nend[iRow];
      HighsInt iSwap = AR_Nend[iRow]++;
      while (ARindex[iFind] != variable_out) iFind++;
      swap(ARindex[iFind], ARindex[iSwap]);
      swap(ARvalue[iFind], ARvalue[iSwap]);
    }
  }
}

double HMatrix::compute_dot(HVector& vector, HighsInt iCol) const {
  double result = 0;
  if (iCol < numCol) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      result += vector.array[Aindex[k]] * Avalue[k];
  } else {
    result = vector.array[iCol - numCol];
  }
  return result;
}

void HMatrix::collect_aj(HVector& vector, HighsInt iCol,
                         double multiplier) const {
  if (iCol < numCol) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      HighsInt index = Aindex[k];
      double value0 = vector.array[index];
      double value1 = value0 + multiplier * Avalue[k];
      if (value0 == 0) vector.index[vector.count++] = index;
      vector.array[index] =
          (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
    }
  } else {
    HighsInt index = iCol - numCol;
    double value0 = vector.array[index];
    double value1 = value0 + multiplier;
    if (value0 == 0) vector.index[vector.count++] = index;
    vector.array[index] =
        (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
  }
}

void HMatrix::priceByColumn(HVector& row_ap, const HVector& row_ep) const {
  // Alias
  HighsInt ap_count = 0;
  HighsInt* ap_index = &row_ap.index[0];
  double* ap_array = &row_ap.array[0];
  const double* ep_array = &row_ep.array[0];
  // Computation
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double value = 0;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      value += ep_array[Aindex[k]] * Avalue[k];
    }
    if (fabs(value) > HIGHS_CONST_TINY) {
      ap_array[iCol] = value;
      ap_index[ap_count++] = iCol;
    }
  }
  row_ap.count = ap_count;
}

void HMatrix::priceByRowSparseResult(HVector& row_ap,
                                     const HVector& row_ep) const {
  // Vanilla hyper-sparse row-wise PRICE
  // Set up parameters so that priceByRowSparseResultWithSwitch runs as vanilla
  // hyper-sparse PRICE
  const double historical_density =
      -0.1;           // Historical density always forces hyper-sparse PRICE
  HighsInt fm_i = 0;  // Always start from first index of row_ep
  const double switch_density = 1.1;  // Never switch to standard row-wise PRICE
  priceByRowSparseResultWithSwitch(row_ap, row_ep, historical_density, fm_i,
                                   switch_density);
}

void HMatrix::priceByRowSparseResultWithSwitch(HVector& row_ap,
                                               const HVector& row_ep,
                                               double historical_density,
                                               HighsInt from_i,
                                               double switch_density) const {
  // (Continue) hyper-sparse row-wise PRICE with possible switches to
  // standard row-wise PRICE either immediately based on historical
  // density or during hyper-sparse PRICE if there is too much fill-in
  // Alias
  HighsInt ap_count = row_ap.count;
  HighsInt* ap_index = &row_ap.index[0];
  double* ap_array = &row_ap.array[0];
  const HighsInt ep_count = row_ep.count;
  const HighsInt* ep_index = &row_ep.index[0];
  const double* ep_array = &row_ep.array[0];
  // Computation

  HighsInt nx_i = from_i;
  // Possibly don't perform hyper-sparse PRICE based on historical density
  if (historical_density <= hyperPRICE) {
    for (HighsInt i = nx_i; i < ep_count; i++) {
      HighsInt iRow = ep_index[i];
      // Possibly switch to standard row-wise price
      HighsInt iRowNNz = AR_Nend[iRow] - ARstart[iRow];
      double lc_dsty = (1.0 * ap_count) / numCol;
      bool price_by_row_sw =
          ap_count + iRowNNz >= numCol || lc_dsty > switch_density;
      if (price_by_row_sw) break;
      double multiplier = ep_array[iRow];
      for (HighsInt k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
        HighsInt index = ARindex[k];
        double value0 = ap_array[index];
        double value1 = value0 + multiplier * ARvalue[k];
        if (value0 == 0) ap_index[ap_count++] = index;
        ap_array[index] =
            (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
      }
      nx_i = i + 1;
    }
    row_ap.count = ap_count;
  }
  from_i = nx_i;
  if (from_i < ep_count) {
    // PRICE is not complete: finish without maintaining nonzeros of result
    priceByRowDenseResult(row_ap, row_ep, from_i);
  } else {
    // PRICE is complete maintaining nonzeros of result
    // Try to remove cancellation
    priceByRowSparseResultRemoveCancellation(row_ap);
  }
}

void HMatrix::priceByRowDenseResult(HVector& row_ap, const HVector& row_ep,
                                    HighsInt from_i) const {
  // (Continue) standard row-wise PRICE
  // Alias
  HighsInt* ap_index = &row_ap.index[0];
  double* ap_array = &row_ap.array[0];
  const HighsInt ep_count = row_ep.count;
  const HighsInt* ep_index = &row_ep.index[0];
  const double* ep_array = &row_ep.array[0];
  // Computation
  for (HighsInt i = from_i; i < ep_count; i++) {
    HighsInt iRow = ep_index[i];
    double multiplier = ep_array[iRow];
    for (HighsInt k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
      HighsInt index = ARindex[k];
      double value0 = ap_array[index];
      double value1 = value0 + multiplier * ARvalue[k];
      ap_array[index] =
          (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
    }
  }
  // Determine indices of nonzeros in PRICE result
  HighsInt ap_count = 0;
  for (HighsInt index = 0; index < numCol; index++) {
    double value1 = ap_array[index];
    if (fabs(value1) < HIGHS_CONST_TINY) {
      ap_array[index] = 0;
    } else {
      ap_index[ap_count++] = index;
    }
  }
  row_ap.count = ap_count;
}

void HMatrix::priceByRowSparseResultRemoveCancellation(HVector& row_ap) const {
  // Alias
  HighsInt* ap_index = &row_ap.index[0];
  double* ap_array = &row_ap.array[0];
  // Try to remove cancellation
  HighsInt ap_count = 0;
  ap_count = row_ap.count;
  const HighsInt apcount1 = ap_count;
  ap_count = 0;
  for (HighsInt i = 0; i < apcount1; i++) {
    const HighsInt index = ap_index[i];
    const double value = ap_array[index];
    if (fabs(value) > HIGHS_CONST_TINY) {
      ap_index[ap_count++] = index;
    } else {
      ap_array[index] = 0;
    }
  }
  row_ap.count = ap_count;
}
