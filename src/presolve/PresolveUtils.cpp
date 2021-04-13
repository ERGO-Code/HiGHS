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
#include "presolve/PresolveUtils.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lp_data/HConst.h"

namespace presolve {

using std::setw;

void printRowOneLine(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue) {
  assert(row >= 0 && row < numRow);

  // go over row and sum
  // col
  // if flagCol[]
  // a_ij * value_j

  double sum = 0.0;
  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++) {
    const HighsInt col = ARindex[k];
    assert(col >= 0 && col <= numCol);
    sum += ARvalue[k] * values[col];
  }

  std::cout << "row " << row << ": " << flagRow[row] << "   " << rowLower[row]
            << " <= " << sum << " <= " << rowUpper[row] << std::endl;
}

void printRow(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue) {
  assert(row >= 0 && row < numRow);

  std::cout << "row " << row << ": " << flagRow[row] << "   " << rowLower[row]
            << " <= ... <= " << rowUpper[row] << std::endl
            << "..." << std::endl;
  // go over row and print
  // col
  // flagCol[] ..... next col in row
  // a_ij
  // x_j
  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++) {
    const HighsInt col = ARindex[k];
    assert(col >= 0 && col <= numCol);
    (void)col;
  }

  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << ARindex[k] << " ";

  std::cout << std::endl;

  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << flagCol[ARindex[k]] << " ";

  std::cout << std::endl;
  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << ARvalue[k] << " ";

  std::cout << std::endl;
  for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << values[ARindex[k]] << " ";

  std::cout << std::endl;
}

void printCol(
    const HighsInt col, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& row_values, const std::vector<HighsInt>& Astart,
    const std::vector<HighsInt>& Aend, const std::vector<HighsInt>& Aindex,
    const std::vector<double>& Avalue) {
  assert(col >= 0 && col < numCol);

  std::cout << "col" << col << ": " << flagCol[col] << "   " << colLower[col]
            << " <= ... <= " << colUpper[col] << std::endl
            << "..." << std::endl;

  // go over col and print
  // row flagRow[] a_ij x_j
  // ...
  // next row in column

  for (HighsInt k = Astart[col]; k < Aend[col]; k++) {
    const HighsInt row = Aindex[k];
    assert(row >= 0 && row <= numRow);
    std::cout << setw(3) << row << " ";
    std::cout << setw(3) << flagRow[row] << " ";
    std::cout << setw(3) << Avalue[k] << " ";
    std::cout << setw(3) << row_values[row] << " ";
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

void printRowWise(
    const HighsInt numRow, const HighsInt numCol,
    const std::vector<double>& colCost, const std::vector<double>& colLower,
    const std::vector<double>& colUpper, const std::vector<double>& rowLower,
    const std::vector<double>& rowUpper, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue) {
  const HighsInt rows = numRow;
  const HighsInt cols = numCol;

  std::cout << "\n-----cost-----\n";

  for (HighsUInt i = 0; i < colCost.size(); i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------AR-|-L-U-----\n";
  for (HighsInt i = 0; i < rows; i++) {
    for (HighsInt j = 0; j < cols; j++) {
      HighsInt ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1])
        std::cout << ARvalue[ind];
      else
        std::cout << "   ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }
  std::cout << "------l------\n";
  for (HighsInt i = 0; i < cols; i++) {
    if (colLower[i] > -kHighsInf)
      std::cout << colLower[i] << " ";
    else
      std::cout << "-inf";
  }
  std::cout << std::endl;
  std::cout << "------u------\n";
  for (HighsInt i = 0; i < cols; i++) {
    if (colUpper[i] < kHighsInf)
      std::cout << colUpper[i] << " ";
    else
      std::cout << "inf ";
  }
  std::cout << std::endl;
}

void printA(const HighsInt numRow, const HighsInt numCol,
            const std::vector<double>& colCost,
            const std::vector<double>& rowLower,
            const std::vector<double>& rowUpper,
            const std::vector<double>& colLower,
            const std::vector<double>& colUpper,
            const std::vector<HighsInt>& Astart,
            const std::vector<HighsInt>& Aindex, std::vector<double>& Avalue) {
  char buff[7];
  std::cout << "\n-----cost-----\n";

  for (HighsInt i = 0; i < numCol; i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------A-|-b-----\n";
  for (HighsInt i = 0; i < numRow; i++) {
    for (HighsInt j = 0; j < numCol; j++) {
      HighsInt ind = Astart[j];
      while (Aindex[ind] != i && ind < Astart[j + 1]) ind++;
      // if a_ij is nonzero print
      if (Aindex[ind] == i && ind < Astart[j + 1]) {
        std::cout << Avalue[ind] << " ";
      } else
        std::cout << " ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }
  std::cout << "------l------\n";
  for (HighsInt i = 0; i < numCol; i++) {
    if (colLower[i] > -kHighsInf)
      std::cout << colLower[i] << " ";
    else
      std::cout << "-inf ";
    std::cout << setw(9) << buff;
  }
  std::cout << std::endl;
  std::cout << "------u------\n";
  for (HighsInt i = 0; i < numCol; i++) {
    if (colUpper[i] < kHighsInf)
      std::cout << colUpper[i] << " ";
    else
      std::cout << "inf ";
  }
  std::cout << std::endl;
}

void printAR(const HighsInt numRow, const HighsInt numCol,
             const std::vector<double>& colCost,
             const std::vector<double>& rowLower,
             const std::vector<double>& rowUpper,
             const std::vector<HighsInt>& ARstart,
             const std::vector<HighsInt>& ARindex,
             std::vector<double>& ARvalue) {
  std::cout << "\n-----cost-----\n";

  for (HighsInt i = 0; i < numCol; i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------AR-|-b-----\n";
  for (HighsInt i = 0; i < numRow; i++) {
    for (HighsInt j = 0; j < numCol; j++) {
      HighsInt ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1]) {
        std::cout << ARvalue[ind] << " ";
      } else
        std::cout << " ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }

  std::cout << std::endl;
}

}  // namespace presolve
