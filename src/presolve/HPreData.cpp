/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
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

double HPreData::getRowValue(int i) {
  double sum = 0;
  for (int k = ARstart[i]; k < ARstart[i + 1]; k++)
    if (flagRow[ARindex[k]]) sum += ARvalue[k] * valuePrimal[ARindex[k]];
  return sum;
}

double HPreData::getaij(int i, int j) {
  int k = ARstart[i];
  while (j != ARindex[k] && k <= ARstart[i + 1]) k++;
  if (k == ARstart[i + 1]) {
    // cout<<"Error: No such element in A: row "<<i<<", column "<<j<<endl;
    // exit(0);
  }
  return ARvalue[k];
}

bool HPreData::isZeroA(int i, int j) {
  int k = ARstart[i];
  while (k < ARstart[i + 1] && j != ARindex[k]) k++;
  if (k == ARstart[i + 1]) {
    return true;
  }
  return false;
}

void HPreData::makeARCopy() {
  // Make a AR copy
  highsSparseTranspose(numRow, numCol, Astart, Aindex, Avalue, ARstart, ARindex, ARvalue);
}

void HPreData::makeACopy() {
  // Make a A copy
  highsSparseTranspose(numColOriginal, numRowOriginal, ARstart, ARindex, ARvalue, Astart, Aindex, Avalue);

  Aend.resize(numColOriginal + 1, 0);
  for (int i = 0; i < numColOriginal; i++) Aend[i] = Astart[i + 1];
}

void initPresolve(PresolveStats& stats) {
  std::cout << "Init Presolve form HiGHS" << std::endl;
}

}  // namespace presolve