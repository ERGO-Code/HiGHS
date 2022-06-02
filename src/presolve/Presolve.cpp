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

void Presolve::load(const HighsLp& lp, bool mip) {
  timer.recordStart(kMatrixCopy);
  numCol = lp.num_col_;
  numRow = lp.num_row_;
  numTot = numTot;
  Astart = lp.a_matrix_.start_;
  Aindex = lp.a_matrix_.index_;
  Avalue = lp.a_matrix_.value_;
  this->mip = mip;

  colCost = lp.col_cost_;
  objShift = lp.offset_;
  if (lp.sense_ == ObjSense::kMaximize) {
    objShift = -objShift;
    for (HighsUInt col = 0; col < lp.col_cost_.size(); col++)
      colCost[col] = -colCost[col];
  }

  integrality = lp.integrality_;
  colLower = lp.col_lower_;
  colUpper = lp.col_upper_;
  rowLower = lp.row_lower_;
  rowUpper = lp.row_upper_;

  modelName = lp.model_name_;
  timer.recordFinish(kMatrixCopy);
}

// printing with cout goes here.
void reportDev(const string& message) {
  std::cout << message << std::flush;
  return;
}

void printMainLoop(const MainLoop& l) {
  std::cout << "    loop : " << l.rows << "," << l.cols << "," << l.nnz << "   "
            << std::endl;
}

void printDevStats(const DevStats& stats) {
  assert(stats.n_loops == (HighsInt)stats.loops.size());

  std::cout << "dev-presolve-stats::" << std::endl;
  std::cout << "  n_loops = " << stats.n_loops << std::endl;
  std::cout << "    loop : rows, cols, nnz " << std::endl;
  for (const MainLoop l : stats.loops) printMainLoop(l);
  return;
}

void getRowsColsNnz(const std::vector<HighsInt>& flagRow,
                    const std::vector<HighsInt>& flagCol,
                    const std::vector<HighsInt>& nzRow,
                    const std::vector<HighsInt>& nzCol, HighsInt& _rows,
                    HighsInt& _cols, HighsInt& _nnz) {
  HighsInt numCol = flagCol.size();
  HighsInt numRow = flagRow.size();
  HighsInt rows = 0;
  HighsInt cols = 0;

  std::vector<HighsInt> nnz_rows(numRow, 0);
  std::vector<HighsInt> nnz_cols(numCol, 0);

  HighsInt total_rows = 0;
  HighsInt total_cols = 0;

  for (HighsInt i = 0; i < numRow; i++)
    if (flagRow.at(i)) {
      rows++;
      nnz_rows[i] += nzRow[i];
      total_rows += nzRow[i];
    }

  for (HighsInt j = 0; j < numCol; j++)
    if (flagCol.at(j)) {
      cols++;
      nnz_cols[j] += nzCol[j];
      total_cols += nzCol[j];
    }

  // Nonzeros.
  assert(total_cols == total_rows);

  _rows = rows;
  _cols = cols;
  _nnz = total_cols;
}

void Presolve::trimA() {
  //  HighsInt cntEl = 0;
  //  for (HighsInt j = 0; j < numCol; ++j)
  //    if (flagCol.at(j)) cntEl += nzCol.at(j);

  vector<pair<HighsInt, size_t>> vp;
  vp.reserve(numCol);

  for (HighsInt i = 0; i != numCol; ++i) {
    vp.push_back(make_pair(Astart.at(i), i));
  }

  // Sorting will put lower values ahead of larger ones,
  // resolving ties using the original index
  sort(vp.begin(), vp.end());

  vector<HighsInt> Aendtmp;
  Aendtmp = Aend;

  HighsInt iPut = 0;
  for (size_t i = 0; i != vp.size(); ++i) {
    HighsInt col = vp.at(i).second;
    if (flagCol.at(col)) {
      HighsInt k = vp.at(i).first;
      Astart.at(col) = iPut;
      while (k < Aendtmp.at(col)) {
        if (flagRow.at(Aindex.at(k))) {
          Avalue[iPut] = Avalue.at(k);
          Aindex[iPut] = Aindex.at(k);
          iPut++;
        }
        k++;
      }
      Aend.at(col) = iPut;
    }
  }
  Avalue.resize(iPut);
  Aindex.resize(iPut);
}

void Presolve::resizeProblem() {
  HighsInt nz = 0;
  HighsInt nR = 0;
  HighsInt nC = 0;

  // arrays to keep track of indices
  rIndex.assign(numRow, -1);
  cIndex.assign(numCol, -1);

  for (HighsInt i = 0; i < numRow; ++i)
    if (flagRow.at(i)) {
      nz += nzRow.at(i);
      rIndex.at(i) = nR;
      nR++;
    }

  for (HighsInt i = 0; i < numCol; ++i)
    if (flagCol.at(i)) {
      cIndex.at(i) = nC;
      nC++;
    }

  // counts
  numRowOriginal = numRow;
  numColOriginal = numCol;
  numRow = nR;
  numCol = nC;
  numTot = nR + nC;

  if (iPrint < 0) {
    stringstream ss;
    ss << ",  Reduced : " << numRow << ",  " << numCol << ",  ";
    reportDev(ss.str());
  }

  chk2.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);

  // This is where status = kEmpty was set if nR + nC == 0
  assert(nR + nC > 0);
  if (nR + nC == 0) return;

  // matrix
  vector<HighsInt> iwork(numCol, 0);
  Astart.assign(numCol + 1, 0);
  Aend.assign(numCol + 1, 0);
  Aindex.resize(nz);
  Avalue.resize(nz);

  for (HighsInt i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i))
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        const HighsInt j = ARindex.at(k);
        if (flagCol.at(j)) iwork.at(cIndex.at(j))++;
      }

  for (HighsInt i = 1; i <= numCol; ++i)
    Astart.at(i) = Astart.at(i - 1) + iwork.at(i - 1);
  for (HighsInt i = 0; i < numCol; ++i) iwork.at(i) = Aend.at(i) = Astart.at(i);
  for (HighsInt i = 0; i < numRowOriginal; ++i) {
    if (flagRow.at(i)) {
      HighsInt iRow = rIndex.at(i);
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        const HighsInt j = ARindex.at(k);
        if (flagCol.at(j)) {
          HighsInt iCol = cIndex.at(j);
          HighsInt iPut = iwork.at(iCol)++;
          Aindex.at(iPut) = iRow;
          Avalue.at(iPut) = ARvalue.at(k);
        }
      }
    }
  }

  if (iPrint < 0) {
    stringstream ss;
    ss << Avalue.size() << ", ";
    reportDev(ss.str());
  }

  // also call before trimming
  resizeImpliedBounds();

  // cost, bounds
  colCostAtEl = colCost;
  vector<double> tempCost = colCost;
  vector<double> temp = colLower;
  vector<double> teup = colUpper;
  vector<HighsVarType> tempintegrality = integrality;

  colCost.resize(numCol);
  colLower.resize(numCol);
  colUpper.resize(numCol);
  integrality.resize(numCol);

  HighsInt k = 0;
  for (HighsInt i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      colCost.at(k) = tempCost.at(i);
      colLower.at(k) = temp.at(i);
      colUpper.at(k) = teup.at(i);
      integrality.at(k) = tempintegrality.at(i);
      k++;
    }

  // RHS and bounds
  temp = rowLower;
  teup = rowUpper;
  rowLower.resize(numRow);
  rowUpper.resize(numRow);
  k = 0;
  for (HighsInt i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      rowLower.at(k) = temp.at(i);
      rowUpper.at(k) = teup.at(i);
      k++;
    }
}

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

void Presolve::detectImpliedIntegers() {
  std::vector<HighsInt> numcont(numRow);
  std::vector<HighsInt> equations;
  equations.reserve(numRow);

  for (HighsInt i = 0; i != numRow; ++i) {
    if (!flagRow[i]) continue;
    if (rowLower[i] != rowUpper[i]) continue;

    const HighsInt start = ARstart[i];
    const HighsInt end = ARstart[i + 1];

    for (HighsInt j = start; j != end; ++j) {
      if (!flagCol[ARindex[j]]) continue;
      if (integrality[ARindex[j]] == HighsVarType::kContinuous) ++numcont[i];
    }

    if (numcont[i] == 1) equations.push_back(i);
  }

  HighsInt numimplint = 0;
  HighsInt primalimplint;

  for (size_t k = 0; k != equations.size(); ++k) {
    HighsInt i = equations[k];
    if (numcont[i] != 1) continue;

    const HighsInt start = ARstart[i];
    const HighsInt end = ARstart[i + 1];

    HighsInt cont = -1;
    for (HighsInt j = start; j != end; ++j) {
      if (!flagCol[ARindex[j]]) continue;
      if (integrality[ARindex[j]] == HighsVarType::kContinuous) {
        cont = j;
        break;
      }
    }

    assert(cont != -1);
    double b = rowUpper[i] / ARvalue[cont];
    if (std::abs(b - std::floor(b + 0.5)) > 1e-9) continue;

    bool impliedint = true;
    for (HighsInt j = start; j != end; ++j) {
      if (j == cont) continue;
      if (!flagCol[ARindex[j]]) continue;

      double val = ARvalue[j] / ARvalue[cont];
      if (std::abs(val - std::floor(val + 0.5)) > 1e-9) {
        impliedint = false;
        break;
      }
    }

    if (!impliedint) continue;

    HighsInt col = ARindex[cont];
    const HighsInt colstart = Astart[col];
    const HighsInt colend = Aend[col];
    integrality[col] = HighsVarType::kImplicitInteger;
    roundIntegerBounds(col);
    ++numimplint;

    for (HighsInt j = colstart; j != colend; ++j) {
      if (--numcont[Aindex[j]] == 1) {
        assert(rowLower[Aindex[j]] == rowUpper[Aindex[j]]);
        equations.push_back(Aindex[j]);
      }
    }
  }

  highsLogDev(log_options, HighsLogType::kVerbose,
              "found %" HIGHSINT_FORMAT
              " implied integers with primal detection method\n",
              numimplint);

  primalimplint = numimplint;

  for (HighsInt i = 0; i != numCol; ++i) {
    if (!flagCol[i]) continue;
    if (integrality[i] != HighsVarType::kContinuous) continue;

    const HighsInt colstart = Astart[i];
    const HighsInt colend = Aend[i];
    bool haseq = false;
    for (HighsInt j = colstart; j != colend; ++j) {
      HighsInt row = Aindex[j];
      if (!flagRow[row]) continue;
      if (rowLower[row] == rowUpper[row]) {
        haseq = true;
        break;
      }
    }

    if (haseq) continue;

    bool impliedinteger = true;
    for (HighsInt j = colstart; j != colend; ++j) {
      HighsInt row = Aindex[j];
      if (!flagRow[row]) continue;

      if (rowUpper[row] != kHighsInf) {
        double val = rowUpper[row] / Avalue[j];

        if (std::abs(val - std::floor(val + 0.5)) > 1e-9) {
          impliedinteger = false;
          break;
        }
      }

      if (rowLower[row] != -kHighsInf) {
        double val = rowLower[row] / Avalue[j];

        if (std::abs(val - std::floor(val + 0.5)) > 1e-9) {
          impliedinteger = false;
          break;
        }
      }

      const HighsInt start = ARstart[row];
      const HighsInt end = ARstart[row + 1];
      for (HighsInt k = start; k != end; ++k) {
        if (ARindex[k] == i) continue;
        if (!flagCol[ARindex[k]]) continue;

        if (integrality[ARindex[k]] == HighsVarType::kContinuous) {
          impliedinteger = false;
          break;
        }

        double val = ARvalue[k] / Avalue[j];
        if (std::abs(val - std::floor(val + 0.5)) > 1e-9) {
          impliedinteger = false;
          break;
        }
      }

      if (!impliedinteger) break;
    }

    if (!impliedinteger) continue;
    integrality[i] = HighsVarType::kImplicitInteger;
    roundIntegerBounds(i);
    ++numimplint;
  }

  highsLogDev(log_options, HighsLogType::kVerbose,
              "found %" HIGHSINT_FORMAT
              " implied integers with dual detection method\n",
              numimplint - primalimplint);

  highsLogDev(log_options, HighsLogType::kVerbose,
              "implint detection found %" HIGHSINT_FORMAT " implied integers\n",
              numimplint);
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

void Presolve::roundIntegerBounds(const HighsInt col) {
  // for mip we check if the bounds can be rounded
  if (mip && integrality[col] != HighsVarType::kContinuous) {
    if (colLower[col] != -kHighsInf)
      colLower[col] =
          ceil(colLower[col] - default_primal_feasiblility_tolerance);

    if (colUpper[col] != kHighsInf)
      colUpper[col] =
          floor(colUpper[col] + default_primal_feasiblility_tolerance);
  }
}

void Presolve::fillStackRowBounds(HighsInt row) {
  postValue.push(rowUpper.at(row));
  postValue.push(rowLower.at(row));
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

HighsPostsolveStatus Presolve::postsolve(const HighsSolution& reduced_solution,
                                         const HighsBasis& reduced_basis,
                                         HighsSolution& recovered_solution,
                                         HighsBasis& recovered_basis) {
  colValue = reduced_solution.col_value;
  colDual = reduced_solution.col_dual;
  rowDual = reduced_solution.row_dual;

  col_status = reduced_basis.col_status;
  row_status = reduced_basis.row_status;

  makeACopy();  // so we can efficiently calculate primal and dual values

  //	iKKTcheck = false;
  // set corresponding parts of solution vectors:
  HighsInt j_index = 0;
  vector<HighsInt> eqIndexOfReduced(numCol, -1);
  vector<HighsInt> eqIndexOfReduROW(numRow, -1);
  for (HighsInt i = 0; i < numColOriginal; ++i)
    if (cIndex.at(i) > -1) {
      eqIndexOfReduced.at(j_index) = i;
      ++j_index;
    }
  j_index = 0;
  for (HighsInt i = 0; i < numRowOriginal; ++i)
    if (rIndex.at(i) > -1) {
      eqIndexOfReduROW.at(j_index) = i;
      ++j_index;
    }

  vector<HighsBasisStatus> temp_col_status = col_status;
  vector<HighsBasisStatus> temp_row_status = row_status;

  nonbasicFlag.assign(numColOriginal + numRowOriginal, 1);
  col_status.assign(numColOriginal, HighsBasisStatus::kNonbasic);  // Was LOWER
  row_status.assign(numRowOriginal, HighsBasisStatus::kNonbasic);  // Was LOWER

  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt iCol = eqIndexOfReduced.at(i);
    assert(iCol < (HighsInt)valuePrimal.size());
    assert(iCol < (HighsInt)valueColDual.size());
    assert(iCol >= 0);
    valuePrimal[iCol] = colValue.at(i);
    valueColDual[iCol] = colDual.at(i);
    col_status.at(iCol) = temp_col_status.at(i);
  }

  for (HighsInt i = 0; i < numRow; ++i) {
    HighsInt iRow = eqIndexOfReduROW.at(i);
    valueRowDual[iRow] = rowDual.at(i);
    row_status.at(iRow) = temp_row_status.at(i);
  }

  // cmpNBF(-1, -1);
  // testBasisMatrixSingularity();

  if (iKKTcheck) {
    cout << std::endl << "~~~~~ KKT check on HiGHS solution ~~~~~\n";
    checkKkt();
  }

  vector<HighsInt> fRjs;
  while (!chng.empty()) {
    change c = chng.top();
    chng.pop();
    // cout<<"chng.pop:       "<<c.col<<"       "<<c.row << endl;

    setBasisElement(c);
    switch (c.type) {
      case kAggregator: {
        // restore solution, basis, flags, and colCostAtEl
        aggregatorStack.back().postsolveStack.undo(
            flagCol, flagRow, valuePrimal, valueColDual, valueRowDual,
            col_status, row_status);

        // restore AR to state before the aggregator got called
        ARstart = std::move(aggregatorStack.back().ARstartAtCall);
        ARindex = std::move(aggregatorStack.back().ARindexAtCall);
        ARvalue = std::move(aggregatorStack.back().ARvalueAtCall);
        colCostAtEl = std::move(aggregatorStack.back().colCostAtCall);
        aggregatorStack.pop_back();

        // restore A from AR
        makeACopy();
        break;
      }
      case kTwoColSingTrivial: {
        // WIP
        HighsInt y = (HighsInt)postValue.top();
        postValue.pop();
        HighsInt x = (HighsInt)postValue.top();
        postValue.pop();
        assert(x == c.col);
        flagRow[c.row] = true;
        flagCol[x] = true;
        flagCol[y] = true;
        row_status.at(c.row) = HighsBasisStatus::kBasic;
        break;
      }
      case kDoubletonEquation: {  // Doubleton equation row
        getDualsDoubletonEquation(c.row, c.col);

        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout
                << "----KKT check after doubleton equation re-introduced. Row: "
                << c.row << ", column " << c.col << " -----\n";
          chk2.addChange(17, c.row, c.col, valuePrimal[c.col],
                         valueColDual[c.col], valueRowDual[c.row]);
          checkKkt();
        }
        // exit(2);
        break;
      }
      case kDoubletonEquationRowBoundsUpdate: {
        // new bounds from doubleton equation, retrieve old ones
        // just for KKT check, not called otherwise
        chk2.addChange(171, c.row, c.col, 0, 0, 0);
        break;
      }
      case kDoubletonEquationNewXNonzero: {
        // matrix transformation from doubleton equation, case x still there
        // case new x is not 0
        // just change value of entry in row for x
        HighsInt indi;
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == c.col) break;
        ARvalue.at(indi) = postValue.top();
        for (indi = Astart[c.col]; indi < Aend[c.col]; ++indi)
          if (Aindex.at(indi) == c.row) break;
        Avalue.at(indi) = postValue.top();

        if (iKKTcheck == 1)
          chk2.addChange(172, c.row, c.col, postValue.top(), 0, 0);
        postValue.pop();

        break;
      }
      case kDoubletonEquationXZeroInitially: {
        // matrix transformation from doubleton equation, retrieve old value
        // case when row does not have x initially: entries for row i swap x and
        // y cols

        const HighsInt yindex = (HighsInt)postValue.top();
        postValue.pop();

        // reverse AR for case when x is zero and y entry has moved
        HighsInt indi;
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == c.col) break;
        ARvalue.at(indi) = postValue.top();
        ARindex.at(indi) = yindex;

        // reverse A for case when x is zero and y entry has moved
        for (indi = Astart[c.col]; indi < Aend[c.col]; ++indi)
          if (Aindex.at(indi) == c.row) break;

        // recover x: column decreases by 1
        // if indi is not Aend-1 swap elements indi and Aend-1
        if (indi != Aend[c.col] - 1) {
          double tmp = Avalue[Aend[c.col] - 1];
          HighsInt tmpi = Aindex[Aend[c.col] - 1];
          Avalue[Aend[c.col] - 1] = Avalue.at(indi);
          Aindex[Aend[c.col] - 1] = Aindex.at(indi);
          Avalue.at(indi) = tmp;
          Aindex.at(indi) = tmpi;
        }
        Aend[c.col]--;

        // recover y: column increases by 1
        // update A: append X column to end of array
        HighsInt st = Avalue.size();
        for (HighsInt ind = Astart[yindex]; ind < Aend[yindex]; ++ind) {
          Avalue.push_back(Avalue.at(ind));
          Aindex.push_back(Aindex.at(ind));
        }
        Avalue.push_back(postValue.top());
        Aindex.push_back(c.row);
        Astart[yindex] = st;
        Aend[yindex] = Avalue.size();

        double topp = postValue.top();
        postValue.pop();
        if (iKKTcheck == 1) {
          chk2.addChange(173, c.row, c.col, topp, (double)yindex, 0);
        }

        break;
      }
      case kDoubletonEquationNewXZeroArUpdate: {
        // sp case x disappears row representation change
        HighsInt indi;
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == numColOriginal) break;
        ARindex.at(indi) = c.col;
        ARvalue.at(indi) = postValue.top();

        postValue.pop();

        break;
      }
      case kDoubletonEquationNewXZeroAUpdate: {
        // sp case x disappears column representation change
        // here A is copied from AR array at end of presolve so need to expand x
        // column  Aend[c.col]++; wouldn't do because old value is overriden
        double oldXvalue = postValue.top();
        postValue.pop();
        HighsInt x = c.col;

        // update A: append X column to end of array
        HighsInt st = Avalue.size();
        for (HighsInt ind = Astart.at(x); ind < Aend.at(x); ++ind) {
          Avalue.push_back(Avalue.at(ind));
          Aindex.push_back(Aindex.at(ind));
        }
        Avalue.push_back(oldXvalue);
        Aindex.push_back(c.row);
        Astart.at(x) = st;
        Aend.at(x) = Avalue.size();

        break;
      }
      case kEmptyRow: {
        valueRowDual[c.row] = 0;
        flagRow[c.row] = 1;
        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after empty row " << c.row
                 << " re-introduced-----\n";
          chk2.addChange(0, c.row, 0, 0, 0, 0);
          checkKkt();
        }
        break;
      }
      case kSingRow: {
        // valuePrimal is already set for this one, colDual also, we need
        // rowDual. AR copy keeps full matrix.  col dual maybe infeasible, we
        // need to check.  recover old bounds and see
        getDualsSingletonRow(c.row, c.col);

        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after singleton row " << c.row
                 << " re-introduced. Variable: " << c.col << " -----\n";
          chk2.addChange(1, c.row, c.col, valuePrimal[c.col],
                         valueColDual[c.col], valueRowDual[c.row]);
          checkKkt();
        }
        break;
      }
      case kForcingRowVariable:
        fRjs.push_back(c.col);
        flagCol[c.col] = 1;
        if (iKKTcheck == 1 && valuePrimal[c.col] != 0)
          chk2.addChange(22, c.row, c.col, 0, 0, 0);
        break;
      case kForcingRow: {
        string str = getDualsForcingRow(c.row, fRjs);

        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after forcing row " << c.row
                 << " re-introduced. Variable(s): " << str << " -----\n";
          chk2.addChange(3, c.row, 0, 0, 0, valueRowDual[c.row]);
          checkKkt();
        }
        fRjs.clear();
        break;
      }
      case kRedundantRow: {
        // this is not zero if the row bounds got relaxed and transferred to a
        // column which then had a nonzero dual.
        valueRowDual[c.row] = 0;

        flagRow[c.row] = 1;

        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after redundant row " << c.row
                 << " re-introduced.----------------\n";
          checkKkt();
        }
        break;
      }
      case kFreeSingCol:
      case kImpliedFreeSingCol: {
        // colDual rowDual already set.
        // calculate row value without xj
        double aij = getaij(c.row, c.col);
        double sum = 0;
        for (HighsInt k = ARstart[c.row]; k < ARstart[c.row + 1]; ++k)
          if (flagCol.at(ARindex.at(k)))
            sum += valuePrimal.at(ARindex.at(k)) * ARvalue.at(k);

        double rowlb = postValue.top();
        postValue.pop();
        double rowub = postValue.top();
        postValue.pop();

        // calculate xj
        if (valueRowDual[c.row] < 0) {
          // row is at lower bound
          valuePrimal[c.col] = (rowlb - sum) / aij;
        } else if (valueRowDual[c.row] > 0) {
          // row is at upper bound
          valuePrimal[c.col] = (rowub - sum) / aij;
        } else if (rowlb == rowub)
          valuePrimal[c.col] = (rowlb - sum) / aij;
        else if (colCostAtEl[c.col] > 0) {
          // we are interested in the lowest possible value of x:
          // max { l_j, bound implied by row i }
          double bndL;
          if (aij > 0)
            bndL = (rowlb - sum) / aij;
          else
            bndL = (rowub - sum) / aij;
          valuePrimal[c.col] = max(colLowerOriginal[c.col], bndL);
        } else if (colCostAtEl[c.col] < 0) {
          // we are interested in the highest possible value of x:
          // min { u_j, bound implied by row i }
          double bndU;
          if (aij < 0)
            bndU = (rowlb - sum) / aij;
          else
            bndU = (rowub - sum) / aij;
          valuePrimal[c.col] = min(colUpperOriginal[c.col], bndU);
        } else {  // cost is zero
          double bndL, bndU;
          if (aij > 0) {
            bndL = (rowlb - sum) / aij;
            bndU = (rowub - sum) / aij;
          } else {
            bndL = (rowub - sum) / aij;
            bndU = (rowlb - sum) / aij;
          }
          double valuePrimalUB = min(colUpperOriginal[c.col], bndU);
          double valuePrimalLB = max(colLowerOriginal[c.col], bndL);
          if (valuePrimalUB < valuePrimalLB - tol) {
            cout << "Postsolve error: inconsistent bounds for implied free "
                    "column singleton "
                 << c.col << endl;
          }

          if (fabs(valuePrimalLB) < fabs(valuePrimalUB))
            valuePrimal[c.col] = valuePrimalLB;
          else
            valuePrimal[c.col] = valuePrimalUB;
        }
        sum = sum + valuePrimal[c.col] * aij;

        double costAtTimeOfElimination = postValue.top();
        postValue.pop();
        objShift += (costAtTimeOfElimination * sum) / aij;

        flagRow[c.row] = 1;
        flagCol[c.col] = 1;
        // valueRowDual[c.row] = 0;

        if (iKKTcheck == 1) {
          chk2.addCost(c.col, costAtTimeOfElimination);
          if (c.type == kFreeSingCol && chk2.print == 1)
            cout << "----KKT check after free col singleton " << c.col
                 << " re-introduced. Row: " << c.row << " -----\n";
          else if (c.type == kImpliedFreeSingCol && chk2.print == 1)
            cout << "----KKT check after implied free col singleton " << c.col
                 << " re-introduced. Row: " << c.row << " -----\n";
          chk2.addChange(4, c.row, c.col, valuePrimal[c.col],
                         valueColDual[c.col], valueRowDual[c.row]);
          checkKkt();
        }
        break;
      }
      case kSingColDoubletonIneq: {
        // column singleton in a doubleton equation.
        // colDual already set. need valuePrimal from stack. maybe change
        // rowDual depending on bounds. old bounds kept in oldBounds. variables
        // j,k : we eliminated j and are left with changed bounds on k and no
        // row. c.col is column COL (K) - eliminated, j is with new bounds
        pair<HighsInt, vector<double>> p = oldBounds.top();
        oldBounds.pop();
        const HighsInt j = p.first;
        vector<double> v = p.second;
        // double lbNew = v[0];
        // double ubNew = v[1];
        double cjNew = v[2];

        p = oldBounds.top();
        oldBounds.pop();
        v = p.second;
        double ubOld = v[1];
        double lbOld = v[0];
        double cjOld = v[2];

        p = oldBounds.top();
        oldBounds.pop();
        v = p.second;
        double ubCOL = v[1];
        double lbCOL = v[0];
        double ck = v[2];

        double rowlb = postValue.top();
        postValue.pop();
        double rowub = postValue.top();
        postValue.pop();
        double aik = postValue.top();
        postValue.pop();
        double aij = postValue.top();
        postValue.pop();
        double xj = valuePrimal.at(j);

        // calculate xk, depending on signs of coeff and cost
        double upp = kHighsInf;
        double low = -kHighsInf;

        if ((aij > 0 && aik > 0) || (aij < 0 && aik < 0)) {
          if (rowub < kHighsInf) upp = (rowub - aij * xj) / aik;
          if (rowlb > -kHighsInf) low = (rowlb - aij * xj) / aik;
        } else {
          if (rowub < kHighsInf) upp = (rowub - aij * xj) / aik;
          if (rowlb > -kHighsInf) low = (rowlb - aij * xj) / aik;
        }

        double xkValue = 0;
        if (ck == 0) {
          if (low < 0 && upp > 0)
            xkValue = 0;
          else if (fabs(low) < fabs(upp))
            xkValue = low;
          else
            xkValue = upp;
        }

        else if ((ck > 0 && aik > 0) || (ck < 0 && aik < 0)) {
          assert(low > -kHighsInf);
          xkValue = low;
        } else if ((ck > 0 && aik < 0) || (ck < 0 && aik > 0)) {
          assert(low < kHighsInf);
          xkValue = upp;
        }

        // primal value and objective shift
        valuePrimal[c.col] = xkValue;
        objShift += -cjNew * xj + cjOld * xj + ck * xkValue;

        // fix duals
        double rowVal = aij * xj + aik * xkValue;

        // If row is strictly between bounds:
        // Row is basic and column is non basic.
        if ((rowub == kHighsInf || (rowub - rowVal > tol)) &&
            (rowlb == -kHighsInf || (rowVal - rowlb > tol))) {
          row_status.at(c.row) = HighsBasisStatus::kBasic;
          col_status.at(c.col) = HighsBasisStatus::kNonbasic;
          valueRowDual[c.row] = 0;
          flagRow[c.row] = 1;
          valueColDual[c.col] = getColumnDualPost(c.col);
        } else {
          // row is at a bound
          // case fabs(rowlb - rowub) < tol
          double lo = -kHighsInf;
          double up = kHighsInf;

          if (fabs(rowub - rowVal) <= tol) {
            lo = 0;
            up = kHighsInf;
          } else if (fabs(rowlb - rowVal) <= tol) {
            lo = -kHighsInf;
            up = 0;
          }

          colCostAtEl.at(j) = cjOld;  // revert cost before calculating duals
          getBoundOnLByZj(c.row, j, &lo, &up, lbOld, ubOld);
          getBoundOnLByZj(c.row, c.col, &lo, &up, lbCOL, ubCOL);

          // calculate yi
          if (lo - up > tol)
            cout << "PR: Error in postsolving doubleton inequality " << c.row
                 << " : inconsistent bounds for its dual value." << std::endl;

          // WARNING: bound_row_dual not used. commented out to surpress warning
          // but maybe this causes trouble. Look into when you do dual postsolve
          // again (todo)
          //
          //
          // double bound_row_dual = 0;
          // if (lo > 0) {
          //   bound_row_dual = lo;
          // } else if (up < 0) {
          //   bound_row_dual = up;
          // }

          // kxx
          // if (lo > 0 || up < 0)
          if (lo > 0 || up < 0 || ck != 0) {
            // row is nonbasic
            // since either dual value zero for it is infeasible
            // or the column cost has changed for col j hence the row dual has
            // to be nonzero to balance out the Stationarity of Lagrangian.
            row_status.at(c.row) = HighsBasisStatus::kNonbasic;
            col_status.at(c.col) = HighsBasisStatus::kBasic;
            valueColDual[c.col] = 0;
            flagRow[c.row] = 1;
            valueRowDual[c.row] = getRowDualPost(c.row, c.col);
            valueColDual[j] = getColumnDualPost(j);
          } else {
            // zero row dual is feasible, set row to basic and column to
            // nonbasic.
            row_status.at(c.row) = HighsBasisStatus::kBasic;
            col_status.at(c.col) = HighsBasisStatus::kNonbasic;
            valueRowDual[c.row] = 0;
            flagRow[c.row] = 1;
            valueColDual[c.col] = getColumnDualPost(c.col);
          }
        }

        flagCol[c.col] = 1;

        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after col singleton " << c.col
                 << " in doubleton ineq re-introduced. Row: " << c.row
                 << " -----\n";

          chk2.addChange(5, c.row, c.col, valuePrimal[c.col],
                         valueColDual[c.col], valueRowDual[c.row]);
          checkKkt();
        }
        // exit(2);
        break;
      }
      case kEmptyCol:
      case kDominatedCols:
      case kWeaklyDominatedCols: {
        // got valuePrimal, need colDual
        if (c.type != kEmptyCol) {
          double z = colCostAtEl[c.col];
          for (HighsInt k = Astart[c.col]; k < Astart[c.col + 1]; ++k)
            if (flagRow.at(Aindex.at(k)))
              z = z + valueRowDual.at(Aindex.at(k)) * Avalue.at(k);
          valueColDual[c.col] = z;
        }

        flagCol[c.col] = 1;
        if (iKKTcheck == 1) {
          if (c.type == kEmptyCol && chk2.print == 1)
            cout << "----KKT check after empty column " << c.col
                 << " re-introduced.-----------\n";
          else if (c.type == kDominatedCols && chk2.print == 1)
            cout << "----KKT check after dominated column " << c.col
                 << " re-introduced.-----------\n";
          else if (c.type == kWeaklyDominatedCols && chk2.print == 1)
            cout << "----KKT check after weakly dominated column " << c.col
                 << " re-introduced.-----------\n";

          chk2.addChange(6, 0, c.col, valuePrimal[c.col], valueColDual[c.col],
                         0);
          checkKkt();
        }
        break;
      }

      case kFixedCol: {
        // got valuePrimal, need colDual
        valueColDual[c.col] = getColumnDualPost(c.col);

        flagCol[c.col] = 1;
        if (iKKTcheck == 1) {
          if (chk2.print == 1)
            cout << "----KKT check after fixed variable " << c.col
                 << " re-introduced.-----------\n";
          chk2.addChange(7, 0, c.col, valuePrimal[c.col], valueColDual[c.col],
                         0);
          checkKkt();
        }
        break;
      }
    }
    // cmpNBF(c.row, c.col);
  }

  // cmpNBF();

  // Check number of basic variables
  HighsInt num_basic_var = 0;
  for (HighsInt iCol = 0; iCol < numColOriginal; iCol++) {
    if (col_status[iCol] == HighsBasisStatus::kBasic) {
      assert(num_basic_var < numRowOriginal);
      if (num_basic_var == numRowOriginal) {
        printf("Error in postsolve: more basic variables than rows\n");
        break;
      }
      num_basic_var++;
    }
  }
  for (HighsInt iRow = 0; iRow < numRowOriginal; iRow++) {
    // HighsInt iVar = numColOriginal + iRow;
    if (row_status[iRow] == HighsBasisStatus::kBasic) {
      assert(num_basic_var < numRowOriginal);
      if (num_basic_var == numRowOriginal) {
        printf("Error from postsolve: more basic variables than rows\n");
        break;
      }
      num_basic_var++;
    }
  }
  // Return error if the number of basic variables does not equal the
  // number of rows in the original LP
  assert(num_basic_var == numRowOriginal);
  if (num_basic_var != numRowOriginal) {
    printf("Error from postsolve: number of basic variables = %" HIGHSINT_FORMAT
           " != %" HIGHSINT_FORMAT
           " = number "
           "of rows\n",
           num_basic_var, numRowOriginal);
    return HighsPostsolveStatus::kBasisError;
  }

  // now recover original model data to pass back to HiGHS
  // A is already recovered!
  // however, A is expressed in terms of Astart, Aend and columns are in
  // different order so
  makeACopy();

  numRow = numRowOriginal;
  numCol = numColOriginal;
  numTot = numRow + numCol;

  rowUpper = rowUpperOriginal;
  rowLower = rowLowerOriginal;

  colUpper = colUpperOriginal;
  colLower = colLowerOriginal;

  colCost = colCostOriginal;

  colValue = valuePrimal;
  colDual = valueColDual;
  rowDual = valueRowDual;

  rowValue.assign(numRow, 0);
  for (HighsInt i = 0; i < numRowOriginal; ++i) {
    for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
      rowValue.at(i) += valuePrimal.at(ARindex.at(k)) * ARvalue.at(k);
  }

  // cout<<"Singularity check at end of postsolve: ";
  // testBasisMatrixSingularity();

  if (iKKTcheck != 0) {
    cout << "~~~~~ KKT check of postsolved solution with DevKkt checker ~~~~~"
         << std::endl;

    checkKkt(true);
  }

  // Save solution to PresolveComponentData.
  recovered_solution.col_value = colValue;
  recovered_solution.col_dual = colDual;
  recovered_solution.row_value = rowValue;
  recovered_solution.row_dual = rowDual;

  recovered_basis.col_status = col_status;
  recovered_basis.row_status = row_status;

  return HighsPostsolveStatus::kSolutionRecovered;
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

void Presolve::setBasisElement(change c) {
  // col_status starts off as [numCol] and has already been increased to
  // [numColOriginal] and row_status starts off as [numRow] and has already been
  // increased to [numRowOriginal] so fill fill in gaps in both

  switch (c.type) {
    case kEmptyRow: {
      if (report_postsolve) {
        printf("2.1 : Recover row %3" HIGHSINT_FORMAT " as %3" HIGHSINT_FORMAT
               " (basic): empty row\n",
               c.row, numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::kBasic;
      break;
    }
    case kRedundantRow: {
      if (report_postsolve) {
        printf("2.3 : Recover row %3" HIGHSINT_FORMAT " as %3" HIGHSINT_FORMAT
               " (basic): redundant\n",
               c.row, numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::kBasic;
      break;
    }
    case kFreeSingCol:
    case kImpliedFreeSingCol: {
      if (report_postsolve) {
        printf("2.4a: Recover col %3" HIGHSINT_FORMAT " as %3" HIGHSINT_FORMAT
               " (basic): implied free singleton "
               "column\n",
               c.col, numColOriginal + c.row);
      }
      col_status.at(c.col) = HighsBasisStatus::kBasic;

      if (report_postsolve) {
        printf("2.5b: Recover row %3" HIGHSINT_FORMAT " as %3" HIGHSINT_FORMAT
               " (nonbasic): implied free singleton "
               "column\n",
               c.row, numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::kNonbasic;  // Was LOWER
      break;
    }
    case kEmptyCol:
    case kDominatedCols:
    case kWeaklyDominatedCols: {
      if (report_postsolve) {
        printf("2.7 : Recover column %3" HIGHSINT_FORMAT
               " (nonbasic): weakly dominated column\n",
               c.col);
      }
      col_status.at(c.col) = HighsBasisStatus::kNonbasic;  // Was LOWER
      break;
    }
    case kFixedCol: {  // fixed variable:
      // check if it was NOT after singRow
      if (chng.size() > 0)
        if (chng.top().type != kSingRow) {
          if (report_postsolve) {
            printf("2.8 : Recover column %3" HIGHSINT_FORMAT
                   " (nonbasic): weakly dominated "
                   "column\n",
                   c.col);
          }
          col_status.at(c.col) = HighsBasisStatus::kNonbasic;  // Was LOWER
        }
      break;
    }
    default:
      break;
  }
}

/***
 * lo and up refer to the place storing the current bounds on y_row
 *
 */
void Presolve::getBoundOnLByZj(HighsInt row, HighsInt j, double* lo, double* up,
                               double colLow, double colUpp) {
  double cost = colCostAtEl.at(j);  // valueColDual.at(j);
  double x = -cost;

  double sum = 0;
  for (HighsInt kk = Astart.at(j); kk < Aend.at(j); ++kk)
    if (flagRow.at(Aindex.at(kk))) {
      sum = sum + Avalue.at(kk) * valueRowDual.at(Aindex.at(kk));
    }
  x = x - sum;

  double aij = getaij(row, j);
  x = x / aij;

  if (fabs(colLow - colUpp) < tol)
    return;  // here there is no restriction on zj so no bound on y

  if ((valuePrimal.at(j) - colLow) > tol &&
      (colUpp - valuePrimal.at(j)) > tol) {
    // set both bounds
    if (x < *up) *up = x;
    if (x > *lo) *lo = x;
  }

  else if ((valuePrimal.at(j) == colLow && aij < 0) ||
           (valuePrimal.at(j) == colUpp && aij > 0)) {
    if (x < *up) *up = x;
  } else if ((valuePrimal.at(j) == colLow && aij > 0) ||
             (valuePrimal.at(j) == colUpp && aij < 0)) {
    if (x > *lo) *lo = x;
  }
}

/**
 * returns z_col
 * z = A'y + c
 */
double Presolve::getColumnDualPost(HighsInt col) {
  HighsInt row;
  double z;
  double sum = 0;
  for (HighsInt cnt = Astart.at(col); cnt < Aend.at(col); cnt++)
    if (flagRow.at(Aindex.at(cnt))) {
      row = Aindex.at(cnt);
      sum = sum + valueRowDual.at(row) * Avalue.at(cnt);
    }
  z = sum + colCostAtEl.at(col);
  return z;
}

/***
 * A'y + c = z
 *
 * returns y_row = -(A'y      +   c   - z )/a_rowcol
 *               (except row)  (at el)
 */
double Presolve::getRowDualPost(HighsInt row, HighsInt col) {
  double x = 0;

  for (HighsInt kk = Astart.at(col); kk < Aend.at(col); ++kk)
    if (flagRow.at(Aindex.at(kk)) && Aindex.at(kk) != row)
      x = x + Avalue.at(kk) * valueRowDual.at(Aindex.at(kk));

  x = x + colCostAtEl.at(col) - valueColDual.at(col);

  double y = getaij(row, col);
  return -x / y;
}

string Presolve::getDualsForcingRow(HighsInt row, vector<HighsInt>& fRjs) {
  double z;
  stringstream ss;
  HighsInt j;

  double lo = -kHighsInf;
  double up = kHighsInf;
  HighsInt lo_col = -1;
  HighsInt up_col = -1;

  double cost, sum;

  for (size_t jj = 0; jj < fRjs.size(); ++jj) {
    j = fRjs[jj];

    pair<HighsInt, vector<double>> p = oldBounds.top();
    vector<double> v = get<1>(p);
    oldBounds.pop();
    double colLow = v[0];
    double colUpp = v[1];

    // calculate bound x imposed by zj
    double save_lo = lo;
    double save_up = up;
    getBoundOnLByZj(row, j, &lo, &up, colLow, colUpp);
    if (lo > save_lo) lo_col = j;
    if (up < save_up) up_col = j;
  }

  // calculate yi
  if (lo > up)
    cout << "PR: Error in postsolving forcing row " << row
         << " : inconsistent bounds for its dual value.\n";

  if (lo <= 0 && up >= 0) {
    valueRowDual.at(row) = 0;
    row_status[row] = HighsBasisStatus::kBasic;
  } else if (lo > 0) {
    // row is set to basic and column to non-basic but that should change
    row_status[row] = HighsBasisStatus::kNonbasic;
    col_status.at(lo_col) = HighsBasisStatus::kBasic;
    valueRowDual.at(row) = lo;
    valueColDual.at(lo_col) = 0;
    // valueColDual[lo_col] should be zero since it imposed the lower bound.
  } else if (up < 0) {
    // row is set to basic and column to non-basic but that should change
    row_status[row] = HighsBasisStatus::kNonbasic;
    col_status.at(up_col) = HighsBasisStatus::kBasic;
    valueRowDual.at(row) = up;
    valueColDual.at(up_col) = 0;
  }

  flagRow.at(row) = 1;

  for (size_t jj = 0; jj < fRjs.size(); ++jj) {
    j = fRjs[jj];
    if (lo > 0 && j == lo_col) continue;
    if (up < 0 && j == up_col) continue;

    col_status[j] = HighsBasisStatus::kNonbasic;

    cost = valueColDual.at(j);
    sum = 0;
    for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k)
      if (flagRow.at(Aindex.at(k))) {
        sum = sum + valueRowDual.at(Aindex.at(k)) * Avalue.at(k);
        // cout<<" row "<<Aindex.at(k)<<" dual
        // "<<valueRowDual.at(Aindex.at(k))<<" a_"<<Aindex.at(k)<<"_"<<j<<"\n";
      }
    z = cost + sum;

    valueColDual.at(j) = z;

    if (iKKTcheck == 1) {
      ss << j;
      ss << " ";
      chk2.addChange(2, 0, j, valuePrimal.at(j), valueColDual.at(j), cost);
    }
  }

  return ss.str();
}

void Presolve::getDualsSingletonRow(const HighsInt row, const HighsInt col) {
  pair<HighsInt, vector<double>> bnd = oldBounds.top();
  oldBounds.pop();

  valueRowDual.at(row) = 0;

  const double cost = postValue.top();
  postValue.pop();
  colCostAtEl[col] = cost;

  const double aij = getaij(row, col);
  const double l = (get<1>(bnd))[0];
  const double u = (get<1>(bnd))[1];
  const double lrow = (get<1>(bnd))[2];
  const double urow = (get<1>(bnd))[3];

  flagRow.at(row) = 1;

  HighsBasisStatus local_status = col_status.at(col);
  if (local_status != HighsBasisStatus::kBasic) {
    // x was not basic but is now
    // if x is strictly between original bounds or a_ij*x_j is at a bound.
    if (fabs(valuePrimal.at(col) - l) > tol &&
        fabs(valuePrimal.at(col) - u) > tol) {
      if (report_postsolve) {
        printf("3.1 : Make column %3" HIGHSINT_FORMAT
               " basic and row %3" HIGHSINT_FORMAT " nonbasic\n",
               col, row);
      }
      col_status.at(col) = HighsBasisStatus::kBasic;
      row_status.at(row) = HighsBasisStatus::kNonbasic;  // Was LOWER
      valueColDual[col] = 0;
      valueRowDual[row] = getRowDualPost(row, col);
    } else {
      // column is at bound
      const bool isRowAtLB = fabs(aij * valuePrimal[col] - lrow) < tol;
      const bool isRowAtUB = fabs(aij * valuePrimal[col] - urow) < tol;

      const double save_dual = valueColDual[col];
      valueColDual[col] = 0;
      const double row_dual = getRowDualPost(row, col);

      if ((isRowAtLB && !isRowAtUB && row_dual > 0) ||
          (!isRowAtLB && isRowAtUB && row_dual < 0) ||
          (!isRowAtLB && !isRowAtUB)) {
        // make row basic
        row_status.at(row) = HighsBasisStatus::kBasic;
        valueRowDual[row] = 0;
        valueColDual[col] = save_dual;
      } else {
        // column is basic
        col_status.at(col) = HighsBasisStatus::kBasic;
        row_status.at(row) = HighsBasisStatus::kNonbasic;
        valueColDual[col] = 0;
        valueRowDual[row] = getRowDualPost(row, col);
      }
    }
  } else {
    // x is basic
    if (report_postsolve) {
      printf("3.3 : Make row %3" HIGHSINT_FORMAT " basic\n", row);
    }
    row_status.at(row) = HighsBasisStatus::kBasic;
    valueRowDual[row] = 0;
    // if the row dual is zero it does not contribute to the column dual.
  }
}

void Presolve::getDualsDoubletonEquation(const HighsInt row,
                                         const HighsInt col) {
  // colDual already set. need valuePrimal from stack. maybe change rowDual
  // depending on bounds. old bounds kept in oldBounds. variables j,k : we
  // eliminated col(k)(c.col) and are left with changed bounds on j and no row.
  //                               y x

  constexpr bool report = false;

  pair<HighsInt, vector<double>> p = oldBounds.top();
  oldBounds.pop();
  vector<double> v = get<1>(p);
  const HighsInt x = get<0>(p);
  assert(x >= 0 && x <= numColOriginal);
  const double ubxNew = v[1];
  const double lbxNew = v[0];
  const double cxNew = v[2];

  p = oldBounds.top();
  oldBounds.pop();
  v = get<1>(p);
  const double ubxOld = v[1];
  const double lbxOld = v[0];
  const double cxOld = v[2];

  p = oldBounds.top();
  oldBounds.pop();
  v = get<1>(p);
  const double uby = v[1];
  const double lby = v[0];
  const double cy = v[2];

  const HighsInt y = col;
  assert(y >= 0 && y <= numColOriginal);

  const double b = postValue.top();
  postValue.pop();
  const double aky = postValue.top();
  postValue.pop();
  const double akx = postValue.top();
  postValue.pop();
  const double valueX = valuePrimal.at(x);

  // primal value and objective shift
  valuePrimal.at(y) = (b - akx * valueX) / aky;
  objShift += -cxNew * valueX + cxOld * valueX + cy * valuePrimal.at(y);

  // column cost of x
  colCostAtEl.at(x) = cxOld;

  flagRow.at(row) = 1;
  flagCol.at(y) = 1;

  if (mip) return;

  const HighsBasisStatus x_status_reduced = col_status.at(x);
  bool x_make_basic = false;
  bool y_make_basic = false;
  bool row_basic = false;
  // x stayed, y was removed
  if (valuePrimal.at(y) - lby > tol && uby - valuePrimal.at(y) > tol) {
    // If column y has value between bounds set it to basic.
    col_status.at(y) = HighsBasisStatus::kBasic;
    row_status.at(row) = HighsBasisStatus::kNonbasic;

    // makeYBasic();
    valueColDual.at(y) = 0;
    valueRowDual.at(row) = getRowDualPost(row, y);
    if (report) printf("4.2 : Make column %3" HIGHSINT_FORMAT " basic\n", y);
    return;
  }

  if (((x_status_reduced == HighsBasisStatus::kNonbasic ||
        x_status_reduced == HighsBasisStatus::kUpper) &&
       fabs(valueX - ubxNew) < tol && ubxNew < ubxOld) ||
      ((x_status_reduced == HighsBasisStatus::kNonbasic ||
        x_status_reduced == HighsBasisStatus::kLower) &&
       fabs(valueX - lbxNew) < tol && lbxNew > lbxOld) ||
      (fabs(valueX - lbxNew) < tol && fabs(lbxOld - lbxNew) < tol &&
       (x_status_reduced == HighsBasisStatus::kUpper ||
        x_status_reduced == HighsBasisStatus::kLower))) {
    if (ubxNew > lbxNew) {
      // Column x is nonbasic at reduced solution at a reduced bound but needs
      // to be changed to basic since this bound is expanding.
      assert(col_status.at(y) == HighsBasisStatus::kNonbasic);

      x_make_basic = true;
      // makeXBasic()
    }

    if (ubxNew == lbxNew) {
      if (report) {
        printf(
            "4.5 : Maybe dual restriction on column x : no longer on feas "
            "side.\n");
        if (ubxNew == lbxOld && ubxNew < ubxOld)
          std::cout << " u.  " << valueColDual[x] << std::endl;
        if (lbxNew == ubxOld && lbxNew > lbxOld)
          std::cout << " l.  " << valueColDual[x] << std::endl;

        if ((ubxNew == lbxOld && ubxNew < ubxOld && valueColDual[x] < tol) ||
            (lbxNew == ubxOld && lbxNew > lbxOld && valueColDual[x] > tol)) {
          if (report) printf("4.6 : Change dual of x to zero.\n");
        }

        // make x basic.
        valueColDual.at(x) = 0;
        valueRowDual.at(row) = getRowDualPost(row, x);
        valueColDual.at(y) = getColumnDualPost(y);
        col_status.at(x) = HighsBasisStatus::kBasic;
        row_status.at(row) = HighsBasisStatus::kNonbasic;
        if (report)
          printf("4.77 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
        return;
      }
    }

    if (x_make_basic) {
      // transfer dual of x to dual of row
      valueColDual.at(x) = 0;
      valueRowDual.at(row) = getRowDualPost(row, x);
      valueColDual.at(y) = getColumnDualPost(y);

      if (lby == -kHighsInf || uby == kHighsInf || fabs(lby - uby) > tol) {
        // Make sure y is at a bound
        assert(fabs(valuePrimal[y] - lby) < tol ||
               fabs(valuePrimal[y] - uby) < tol);

        // Check that y is dual feasible
        bool feasible = true;
        if (fabs(valuePrimal[y] - lby) < tol && valueColDual[y] < -tol)
          feasible = false;
        if (fabs(valuePrimal[y] - uby) < tol && valueColDual[y] > tol)
          feasible = false;

        if (feasible) {
          col_status.at(x) = HighsBasisStatus::kBasic;
          row_status.at(row) = HighsBasisStatus::kNonbasic;
          if (report)
            printf("4.1 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
          return;
        }
        // Y not dual feasible
        y_make_basic = true;
      }
      // If not feasble X will remail nonbasic and we will make the
    }
  }

  if (!y_make_basic) {
    if (lby != uby) {
      assert(fabs(lby - valuePrimal[y]) < tol ||
             fabs(uby - valuePrimal[y]) < tol);
      // If postsolved column y is at a bound. If lby != uby we have a
      // restriction on the dual sign of y.
      // col_status.at(y) = HighsBasisStatus::kBasic;
      // row_status.at(row) = HighsBasisStatus::kNonbasic;

      // valueColDual.at(y) = 0;
      // valueRowDual.at(row) = getRowDualPost(row, y);

      // if (report) printf("4.2 : Make column %3" HIGHSINT_FORMAT " basic\n",
      // y);
    } else {
      // Column y is at a bound.
      assert(fabs(uby - valuePrimal[y]) < tol ||
             fabs(lby - valuePrimal[y]) < tol);

      // Check if tight.
      if (fabs(lby - uby) < tol) {
        // assert(fabs(lby - valuePrimal[y]) < tol);
        // no restriction so can make row basic but we don't have to always
        // since it can make the dual of X infeasible (check below)
        row_basic = true;
      }  // Else Will need to check dual feasibility of y dual.

      if (x_status_reduced != HighsBasisStatus::kBasic) {
        // make x basic.
        valueColDual.at(x) = 0;
        valueRowDual.at(row) = getRowDualPost(row, x);
        valueColDual.at(y) = getColumnDualPost(y);
        col_status.at(x) = HighsBasisStatus::kBasic;
        row_status.at(row) = HighsBasisStatus::kNonbasic;
        if (report)
          printf("4.778 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
        return;
      }
    }
  }

  // Print & check some info.
  // if (x_status_reduced == HighsBasisStatus::kBasic)
  //   std::cout << "BASIC" << std::endl;
  // else
  //   std::cout << "NOT BASIC" << std::endl;

  // std::cout << "valueXDual " << valueColDual[x] << std::endl;

  // see if X at a bound
  if (fabs(valueX - ubxNew) < tol || fabs(valueX - lbxNew) < tol) {
    if ((fabs(valueX - ubxNew) < tol && ubxNew < ubxOld) ||
        (fabs(valueX - lbxNew) < tol && lbxNew > lbxOld)) {
      // std::cout << "     4.122" << std::endl;
      // std::cout << lbxOld << " lbxOld " << std::endl;
      // std::cout << ubxOld << " ubxOld " << std::endl;
      // std::cout << lbxNew << " lbxNew " << std::endl;
      // std::cout << ubxNew << " ubxNew " << std::endl;
      // std::cout << valueX << " val  " << std::endl;
      // if X was non basic make it basic
      if (x_status_reduced != HighsBasisStatus::kBasic) {
        valueColDual.at(x) = 0;
        valueRowDual.at(row) = getRowDualPost(row, x);
        valueColDual.at(y) = getColumnDualPost(y);

        // Check dual feasibility of y.
        bool feasible = true;
        if (lby > -kHighsInf && lby < uby && fabs(lby - valuePrimal[y]) < tol)
          if (valueColDual[y] < 0) feasible = false;
        if (uby < kHighsInf && lby < uby && fabs(uby - valuePrimal[y]) < tol)
          if (valueColDual[y] > 0) feasible = false;

        if (feasible) {
          col_status.at(x) = HighsBasisStatus::kBasic;
          row_status.at(row) = HighsBasisStatus::kNonbasic;
          if (report)
            printf("4.122778 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
          return;
        } else {
          if (report)
            printf("4.1227785 : Make column %3" HIGHSINT_FORMAT " basic\n", y);
          // dual of y needs to change. make y basic by proceeding below.
        }
      }

      // if X was basic make y basic.
      valueColDual.at(y) = 0;
      valueRowDual.at(row) = getRowDualPost(row, y);
      valueColDual.at(x) = getColumnDualPost(x);
      col_status.at(y) = HighsBasisStatus::kBasic;
      row_status.at(row) = HighsBasisStatus::kNonbasic;
      if (report)
        printf("4.122779 : Make column %3" HIGHSINT_FORMAT " basic\n", y);
      return;
    } else {
      // std::cout << "     4.002" << std::endl;
      row_basic = false;
    }
  } else {
    // X strictly between bounds
    assert(x_status_reduced == HighsBasisStatus::kBasic);
    assert(valuePrimal[x] - lbxNew > tol && ubxNew - valuePrimal[x] > tol);
  }

  if (row_basic) {
    assert(col_status.at(y) == HighsBasisStatus::kNonbasic);
    row_status.at(row) = HighsBasisStatus::kBasic;

    valueRowDual.at(row) = 0;
    valueColDual.at(y) = getColumnDualPost(y);

    if (report) printf("4.1 : Make row    %3" HIGHSINT_FORMAT " basic\n", row);
  } else {
    // Try Y Basic.

    col_status.at(y) = HighsBasisStatus::kBasic;
    row_status.at(row) = HighsBasisStatus::kNonbasic;

    valueColDual.at(y) = 0;
    valueRowDual.at(row) = getRowDualPost(row, y);

    if (report) printf("4.4 : Make column %3" HIGHSINT_FORMAT " basic\n", y);

    // Check complementary slackness on x.
    if ((valueColDual[x] < -tol && fabs(valuePrimal[x] - lbxOld) > tol) ||
        (valueColDual[x] > tol && fabs(ubxOld - valuePrimal[x]) > tol)) {
      if (x_status_reduced != HighsBasisStatus::kBasic) {
        // make X basic.
        valueColDual.at(x) = 0;
        valueRowDual.at(row) = getRowDualPost(row, x);
        valueColDual.at(y) = getColumnDualPost(y);
        col_status.at(x) = HighsBasisStatus::kBasic;
        col_status.at(y) = HighsBasisStatus::kNonbasic;
        row_status.at(row) = HighsBasisStatus::kNonbasic;
        if (report)
          printf("4.779 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
        return;
      }
      // If X already basic and y can not be feasibly made basic then the row
      // remains as the only option. X not working out

      // row_status.at(row) = HighsBasisStatus::kBasic;
      // col_status.at(y) = HighsBasisStatus::kNonbasic;

      // valueRowDual.at(row) = 0;
      // valueColDual.at(y) = getColumnDualPost(y);

      // if (report) printf("4.7791 : Make row    %3" HIGHSINT_FORMAT "
      // basic\n", row);
      if (report)
        printf("??? 4.7791 : Make row    %3" HIGHSINT_FORMAT " basic\n", row);
    }

    // Check dual feasibility of y.
    bool feasible = true;
    if (lby > -kHighsInf && lby < uby && fabs(lby - valuePrimal[y]) < tol)
      if (valueColDual[y] < 0) feasible = false;
    if (uby < kHighsInf && lby < uby && fabs(uby - valuePrimal[y]) < tol)
      if (valueColDual[y] > 0) feasible = false;

    if (!feasible) {
      // make X basic.
      valueColDual.at(x) = 0;
      valueRowDual.at(row) = getRowDualPost(row, x);
      valueColDual.at(y) = getColumnDualPost(y);
      col_status.at(x) = HighsBasisStatus::kBasic;
      col_status.at(y) = HighsBasisStatus::kNonbasic;
      row_status.at(row) = HighsBasisStatus::kNonbasic;
      if (report)
        printf("4.879 : Make column %3" HIGHSINT_FORMAT " basic\n", x);
    }

    // y is at a bound with infeasible dual
    // : attempt to make basic
  }
}

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
