/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPreData.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/HPreData.h"

#include "util/HighsUtils.h"

using std::cout;
using std::endl;
using std::setw;

namespace presolve {

HPreData::HPreData() {}

double HPreData::getRowValue(HighsInt i) {
  double sum = 0;
  for (HighsInt k = ARstart[i]; k < ARstart[i + 1]; k++)
    if (flagRow[ARindex[k]]) sum += ARvalue[k] * valuePrimal[ARindex[k]];
  return sum;
}

double HPreData::getaij(HighsInt i, HighsInt j) {
  HighsInt k = ARstart[i];
  while (j != ARindex[k] && k <= ARstart[i + 1]) k++;
  if (k == ARstart[i + 1]) {
    // cout<<"Error: No such element in A: row "<<i<<", column "<<j<<endl;
    // exit(0);
  }
  return ARvalue[k];
}

bool HPreData::isZeroA(HighsInt i, HighsInt j) {
  HighsInt k = ARstart[i];
  while (k < ARstart[i + 1] && j != ARindex[k]) k++;
  if (k == ARstart[i + 1]) {
    return true;
  }
  return false;
}

void HPreData::makeARCopy() {
  // Make a AR copy
  highsSparseTranspose(numRow, numCol, Astart, Aindex, Avalue, ARstart, ARindex,
                       ARvalue);
}

void HPreData::makeACopy() {
  // Make a A copy

  vector<HighsInt> iwork(numColOriginal, 0);
  Astart.assign(numColOriginal + 1, 0);
  const HighsInt AcountX = ARindex.size();
  Aindex.resize(AcountX);
  Avalue.resize(AcountX);
  for (HighsInt k = 0; k < AcountX; k++)
    if (ARindex[k] < numColOriginal) iwork[ARindex[k]]++;
  for (HighsInt i = 1; i <= numColOriginal; i++)
    Astart[i] = Astart[i - 1] + iwork[i - 1];
  for (HighsInt i = 0; i < numColOriginal; i++) iwork[i] = Astart[i];
  for (HighsInt iRow = 0; iRow < numRowOriginal; iRow++) {
    for (HighsInt k = ARstart[iRow]; k < ARstart[iRow + 1]; k++) {
      const HighsInt iColumn = ARindex[k];
      if (iColumn != numColOriginal) {
        HighsInt iPut = iwork[iColumn]++;
        Aindex[iPut] = iRow;
        Avalue[iPut] = ARvalue[k];
      }
    }
  }

  Aend.resize(numColOriginal + 1, 0);
  for (HighsInt i = 0; i < numColOriginal; i++) Aend[i] = Astart[i + 1];
}

void initPresolve(PresolveStats& stats) {
  std::cout << "Init Presolve form HiGHS" << std::endl;
}

}  // namespace presolve
