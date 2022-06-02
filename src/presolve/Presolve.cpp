/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/Presolve.cpp
 * @brief
 */
#include "presolve/Presolve.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "presolve/HighsLpPropagator.h"
#include "presolve/PresolveUtils.h"

namespace presolve {

using std::cout;
using std::endl;
using std::flush;
using std::get;
using std::ios;
using std::list;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::setprecision;
using std::setw;
using std::stringstream;

void Presolve::removeFixedCol(HighsInt j) {
  assert(std::isfinite(colUpper[j]));
  setPrimalValue(j, colUpper.at(j));
  addChange(kFixedCol, 0, j);
  if (iPrint > 0)
    cout << "PR: Fixed variable " << j << " = " << colUpper.at(j)
         << ". Column eliminated." << endl;

  countRemovedCols(kFixedCol);

  for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k) {
    if (flagRow.at(Aindex.at(k))) {
      HighsInt i = Aindex.at(k);

      if (nzRow.at(i) == 0) {
        removeEmptyRow(i);
        if (status == Stat::kInfeasible) return;
        countRemovedRows(kFixedCol);
      }
    }
  }
}

void Presolve::removeEmptyRow(HighsInt i) {
  // Analyse dependency on numerical tolerance
  double value = min(rowLower.at(i), -rowUpper.at(i));
  timer.updateNumericsRecord(kNumericsEmptyRowBound, value);
  if (rowLower.at(i) <= empty_row_bound_tolerance &&
      rowUpper.at(i) >= -empty_row_bound_tolerance) {
    if (iPrint > 0) cout << "PR: Empty row " << i << " removed. " << endl;
    flagRow.at(i) = 0;
    valueRowDual.at(i) = 0;
    addChange(kEmptyRow, i, 0);
  } else {
    if (iPrint > 0) cout << "PR: Problem infeasible." << endl;
    status = kInfeasible;
    return;
  }
}

void Presolve::addChange(PresolveRule type, HighsInt row, HighsInt col) {
  change ch;
  ch.type = type;
  ch.row = row;
  ch.col = col;
  chng.push(ch);

  if (type < kPresolveRulesCount) timer.addChange(type);
}

// when setting a value to a primal variable and eliminating row update b,
// singleton Rows linked list, number of nonzeros in rows
void Presolve::setPrimalValue(const HighsInt j, const double value) {
  flagCol.at(j) = 0;
  if (!hasChange) hasChange = true;
  valuePrimal.at(j) = value;

  // update nonzeros
  for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k) {
    HighsInt row = Aindex.at(k);
    if (flagRow.at(row)) {
      nzRow.at(row)--;

      // update singleton row list
      if (nzRow.at(row) == 1) singRow.push_back(row);
    }
  }

  // update values if necessary
  if (fabs(value) > 0) {
    // RHS
    vector<pair<HighsInt, double>> bndsL, bndsU;

    for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k)
      if (flagRow.at(Aindex.at(k))) {
        const HighsInt row = Aindex[k];
        // std::cout << row << " " << rowLower[row] << " " << rowUpper[row] <<
        // std::endl;

        if (iKKTcheck == 1) {
          bndsL.push_back(make_pair(row, rowLower.at(row)));
          bndsU.push_back(make_pair(row, rowUpper.at(row)));
        }
        if (rowLower.at(row) > -kHighsInf)
          rowLower.at(row) -= Avalue.at(k) * value;
        if (rowUpper.at(row) < kHighsInf)
          rowUpper.at(row) -= Avalue.at(k) * value;

        if (implRowValueLower.at(row) > -kHighsInf)
          implRowValueLower.at(row) -= Avalue.at(k) * value;
        if (implRowValueUpper.at(row) < kHighsInf)
          implRowValueUpper.at(row) -= Avalue.at(k) * value;

        if (nzRow.at(row) == 0) {
          if (rowLower[row] - rowUpper[row] > tol) {
            status = kInfeasible;
            return;
          }
          if (rowLower[row] > tol || rowUpper[row] < -tol) {
            status = kInfeasible;
            return;
          }

          flagRow[row] = 0;
          addChange(PresolveRule::kEmptyRow, row, j);
        }
      }

    if (iKKTcheck == 1) {
      chk2.rLowers.push(bndsL);
      chk2.rUppers.push(bndsU);
    }

    // shift objective
    if (colCost.at(j) != 0) objShift += colCost.at(j) * value;
  }
}

void Presolve::resizeImpliedBounds() {
  // implied bounds for crashes
  // row duals
  vector<double> temp = implRowDualLower;
  vector<double> teup = implRowDualUpper;
  implRowDualLower.resize(numRow);
  implRowDualUpper.resize(numRow);

  HighsInt k = 0;
  for (HighsInt i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      implRowDualLower.at(k) = temp.at(i);
      implRowDualUpper.at(k) = teup.at(i);
      k++;
    }

  // row value
  temp = implRowValueLower;
  teup = implRowValueUpper;
  implRowValueLower.resize(numRow);
  implRowValueUpper.resize(numRow);
  k = 0;
  for (HighsInt i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      if (temp.at(i) < rowLower.at(i)) temp.at(i) = rowLower.at(i);
      implRowValueLower.at(k) = temp.at(i);
      if (teup.at(i) > rowUpper.at(i)) teup.at(i) = rowUpper.at(i);
      implRowValueUpper.at(k) = teup.at(i);
      k++;
    }

  // column dual
  temp = implColDualLower;
  teup = implColDualUpper;
  implColDualLower.resize(numCol);
  implColDualUpper.resize(numCol);

  k = 0;
  for (HighsInt i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      implColDualLower.at(k) = temp.at(i);
      implColDualUpper.at(k) = teup.at(i);
      k++;
    }

  // column value
  temp = implColLower;
  teup = implColUpper;
  implColLower.resize(numCol);
  implColUpper.resize(numCol);

  k = 0;
  for (HighsInt i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      if (temp.at(i) < colLower.at(i)) temp.at(i) = colLower.at(i);
      implColLower.at(k) = temp.at(i);
      if (teup.at(i) > colUpper.at(i)) teup.at(i) = colUpper.at(i);
      implColUpper.at(k) = teup.at(i);
      k++;
    }
}

void Presolve::checkKkt(bool final) {
  // final = true or intermediate = true
  if (!iKKTcheck) return;

  // update row value done in initState below.

  std::cout << "~~~~~~~~ " << std::endl;
  bool intermediate = !final;
  dev_kkt_check::State state = initState(intermediate);

  dev_kkt_check::KktInfo info = dev_kkt_check::initInfo();

  bool pass = dev_kkt_check::checkKkt(state, info);
  if (final) {
    if (pass)
      std::cout << "KKT PASS" << std::endl;
    else
      std::cout << "KKT FAIL" << std::endl;
  }
  std::cout << "~~~~~~~~ " << std::endl;
}

/***
 * A'y + c = z
 *
 * returns y_row = -(A'y      +   c   - z )/a_rowcol
 *               (except row)  (at el)
 */
void Presolve::countRemovedRows(PresolveRule rule) {
  timer.increaseCount(true, rule);
}

void Presolve::countRemovedCols(PresolveRule rule) {
  timer.increaseCount(false, rule);
  if (timer.time_limit > 0 &&
      timer.timer_.readRunHighsClock() > timer.time_limit)
    status = Stat::kTimeout;
}

dev_kkt_check::State Presolve::initState(const bool intermediate) {
  // update row value
  rowValue.assign(numRowOriginal, 0);
  for (HighsInt i = 0; i < numRowOriginal; ++i) {
    if (flagRow[i])
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        const HighsInt col = ARindex[k];
        if (flagCol[col]) rowValue.at(i) += valuePrimal.at(col) * ARvalue.at(k);
      }
  }

  if (!intermediate)
    return dev_kkt_check::State(
        numCol, numRow, Astart, Aend, Aindex, Avalue, ARstart, ARindex, ARvalue,
        colCost, colLower, colUpper, rowLower, rowUpper, flagCol, flagRow,
        colValue, colDual, rowValue, rowDual, col_status, row_status);

  // if intermediate step use checker's row and col bounds and cost
  return chk2.initState(numColOriginal, numRowOriginal, Astart, Aend, Aindex,
                        Avalue, ARstart, ARindex, ARvalue, flagCol, flagRow,
                        valuePrimal, valueColDual, rowValue, valueRowDual,
                        col_status, row_status);
}

}  // namespace presolve
