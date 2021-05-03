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
/**@file presolve/Presolve.cpp
 * @brief
 */
#include "presolve/Presolve.h"

//#include "simplex/HFactor.h"
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
  numCol = lp.numCol_;
  numRow = lp.numRow_;
  numTot = numTot;
  Astart = lp.Astart_;
  Aindex = lp.Aindex_;
  Avalue = lp.Avalue_;
  this->mip = mip;

  colCost = lp.colCost_;
  objShift = lp.offset_;
  if (lp.sense_ == ObjSense::kMaximize) {
    objShift = -objShift;
    for (HighsUInt col = 0; col < lp.colCost_.size(); col++)
      colCost[col] = -colCost[col];
  }

  integrality = lp.integrality_;
  colLower = lp.colLower_;
  colUpper = lp.colUpper_;
  rowLower = lp.rowLower_;
  rowUpper = lp.rowUpper_;

  modelName = lp.model_name_;
  timer.recordFinish(kMatrixCopy);
}

void Presolve::setNumericalTolerances() {
  const bool use_original_tol = false;
  const double zero_tolerance = 1e-16;
  if (use_original_tol) {
    inconsistent_bounds_tolerance = tol;
    fixed_column_tolerance =
        zero_tolerance;  // Since exact equality is currently used
    doubleton_equation_bound_tolerance = tol;
    doubleton_inequality_bound_tolerance = tol;
    presolve_small_matrix_value = tol;
    empty_row_bound_tolerance = tol;
    dominated_column_tolerance = tol;
    weakly_dominated_column_tolerance = tol;
  } else {
    // Tolerance on bounds being inconsistent: should be twice
    // primal_feasibility_tolerance since bounds inconsistent by this
    // value can be satisfied to within the primal feasibility tolerance
    // by a primal vlaue at their midpoint. The following is twice the
    // default primal_feasibility_tolerance.
    inconsistent_bounds_tolerance = 2 * default_primal_feasiblility_tolerance;
    // Tolerance on column bounds differences being considered to be
    // zero, allowing a column to be fixed
    fixed_column_tolerance =
        zero_tolerance;  // Since exact equality is currently used
    //        2 * default_primal_feasiblility_tolerance;
    // Tolerance on bound differences being considered to be zero,
    // allowing a doubleton to be treated as an equation. What value
    // this should have is unclear. It could depend on the coefficients
    // of the two variables and the values of the bounds, as there's an
    // implicit infeasibility created when the optimal value for one
    // variable is substituted to deduce the optimal value of the other.
    doubleton_equation_bound_tolerance =
        2 * default_primal_feasiblility_tolerance;
    doubleton_inequality_bound_tolerance = doubleton_equation_bound_tolerance;
    // Need to decide when a matrix coefficient changed by substitution
    // is zeroed: should be the small_matrix_value, for which the
    // following is the default value
    presolve_small_matrix_value = default_small_matrix_value;
    // Tolerance on the lower and upper bound being sufficiently close
    // to zero to allowing an empty row to be removed, rather than have
    // the LP deduced as infeasible. This should be t =
    // primal_feasibility_tolerance since the row activity of zero
    // satisfies a lower bound of at most t, and an upper bound of at
    // least -t. The following is the default
    // primal_feasibility_tolerance.
    empty_row_bound_tolerance = default_primal_feasiblility_tolerance;
    dominated_column_tolerance = default_dual_feasiblility_tolerance;
    weakly_dominated_column_tolerance = default_dual_feasiblility_tolerance;
  }
  timer.model_name = modelName;
  // Initialise the numerics records. JAJH thinks that this has to be
  // done here, as the tolerances are only known in Presolve.h/cpp so
  // have to be passed in
  timer.presolve_numerics.resize(kPresolveNumericsCount);
  timer.initialiseNumericsRecord(kNumericsInconsistentBounds,
                                 "Inconsistent bounds",
                                 inconsistent_bounds_tolerance);
  timer.initialiseNumericsRecord(kNumericsFixedColumn, "Fixed column",
                                 fixed_column_tolerance);
  timer.initialiseNumericsRecord(kNumericsDoubletonEquationBound,
                                 "Doubleton equation bound",
                                 doubleton_equation_bound_tolerance);
  timer.initialiseNumericsRecord(kNumericsDoubletonInequalityBound,
                                 "Doubleton inequality bound",
                                 doubleton_inequality_bound_tolerance);
  timer.initialiseNumericsRecord(kNumericsSmallMatrixValue,
                                 "Small matrix value",
                                 presolve_small_matrix_value);
  timer.initialiseNumericsRecord(kNumericsEmptyRowBound, "Empty row bounds",
                                 empty_row_bound_tolerance);
  timer.initialiseNumericsRecord(kNumericsDominatedColumn, "Dominated column",
                                 dominated_column_tolerance);
  timer.initialiseNumericsRecord(kNumericsWeaklyDominatedColumn,
                                 "Weakly dominated column",
                                 weakly_dominated_column_tolerance);
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

void Presolve::reportDevMidMainLoop() {
  if (iPrint == 0) return;

  HighsInt rows = 0;
  HighsInt cols = 0;
  HighsInt nnz = 0;
  getRowsColsNnz(flagRow, flagCol, nzRow, nzCol, rows, cols, nnz);

  std::cout << "                                             counts " << rows
            << ",  " << cols << ", " << nnz << std::endl;
}

void Presolve::reportDevMainLoop() {
  if (iPrint == 0) {
    if (timer.getTime() > 10)
      highsLogDev(log_options, HighsLogType::kVerbose,
                  "Presolve finished main loop %" HIGHSINT_FORMAT "... ",
                  stats.dev.n_loops + 1);
  } else {
    HighsInt rows = 0;
    HighsInt cols = 0;
    HighsInt nnz = 0;

    getRowsColsNnz(flagRow, flagCol, nzRow, nzCol, rows, cols, nnz);

    stats.dev.n_loops++;
    stats.dev.loops.push_back(MainLoop{rows, cols, nnz});

    std::cout << "Starting loop " << stats.dev.n_loops;

    printMainLoop(stats.dev.loops[stats.dev.n_loops - 1]);
  }
  return;
}

void Presolve::removeEmpty() {
  // cols
  for (HighsInt col = 0; col < numCol; col++) {
    if (flagCol[col])
      if (nzCol[col] == 0) {
        removeEmptyColumn(col);
      }
  }

  // rows
  for (HighsInt row = 0; row < numRow; row++) {
    if (flagRow[row])
      if (nzRow[row] == 0) {
        removeEmptyRow(row);
      }
  }
}

HighsInt Presolve::runPresolvers(const std::vector<Presolver>& order) {
  //***************** main loop ******************

  checkBoundsAreConsistent();
  if (status) return status;

  if (iPrint) std::cout << "----> fixed cols" << std::endl;

  for (Presolver main_loop_presolver : order) {
    double time_start = timer.timer_.readRunHighsClock();
    if (iPrint) std::cout << "----> ";
    auto it = kPresolverNames.find(main_loop_presolver);
    assert(it != kPresolverNames.end());
    if (iPrint) std::cout << (*it).second << std::endl;

    switch (main_loop_presolver) {
      case Presolver::kMainEmpty:
        removeEmpty();
        removeFixed();
        break;
      case Presolver::kMainRowSingletons:
        timer.recordStart(kRemoveRowSingletons);
        removeRowSingletons();
        timer.recordFinish(kRemoveRowSingletons);
        break;
      case Presolver::kMainForcing:
        timer.recordStart(kRemoveForcingConstraints);
        removeForcingConstraints();
        timer.recordFinish(kRemoveForcingConstraints);
        break;
      case Presolver::kMainColSingletons:
        timer.recordStart(kRemoveColumnSingletons);
        removeColumnSingletons();
        timer.recordFinish(kRemoveColumnSingletons);
        break;
      case Presolver::kMainDoubletonEq:
        timer.recordStart(kRemoveDoubletonEquations);
        removeDoubletonEquations();
        timer.recordFinish(kRemoveDoubletonEquations);
        break;
      case Presolver::kMainDominatedCols:
        timer.recordStart(kRemoveDominatedColumns);
        removeDominatedColumns();
        timer.recordFinish(kRemoveDominatedColumns);
        break;
      case Presolver::kMainSingletonsOnly:
        // To implement
        // timer.recordStart(SING_ONLY);
        removeSingletonsOnly();
        // timer.recordFinish(SING_ONLY);
        break;
      case Presolver::kMainMipDualFixing:
        timer.recordStart(kMipDualFixing);
        applyMipDualFixing();
        timer.recordFinish(kMipDualFixing);
        break;
    }

    double time_end = timer.timer_.readRunHighsClock();
    if (iPrint)
      std::cout << (*it).second << " time: " << time_end - time_start
                << std::endl;
    reportDevMidMainLoop();
    if (status) return status;
  }

  //***************** main loop ******************
  return status;
}

// void Presolve::removeSingletonsOnly() {
//   for (HighsInt row = 0; row < numRow; row++) {
//     if (!flagRow[row]) continue;
//     bool valid = true;
//     HighsInt nz_col = 0;
//     for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++) {
//       const HighsInt col = ARindex[k];
//       if (!flagCol[col]) continue;
//       if (nzCol[col] != 1) {
//         valid = false;
//         break;
//       }
//       nz_col++;
//     }
//     if (!valid) continue;
//     if (nz_col == 0) {
//       flagRow[row] = false;
//       continue;
//     }

//     std::cout << "Singletons only row found! nzcol = " << nz_col << " L = "
//     << rowLower[row] << " U = " << rowUpper[row] << std::endl;
//   }
// }

void Presolve::removeFixed() {
  timer.recordStart(kFixedCol);
  for (HighsInt j = 0; j < numCol; ++j)
    if (flagCol.at(j)) {
      // Analyse dependency on numerical tolerance
      timer.updateNumericsRecord(kNumericsFixedColumn,
                                 fabs(colUpper.at(j) - colLower.at(j)));
      roundIntegerBounds(j);
      if (fabs(colUpper.at(j) - colLower.at(j)) > fixed_column_tolerance)
        continue;
      removeFixedCol(j);
      if (status) {
        timer.recordFinish(kFixedCol);
        return;
      }
    }
  timer.recordFinish(kFixedCol);
}

HighsInt Presolve::presolve(HighsInt print) {
  timer.start_time = timer.getTime();
  bool aggregatorCalled = false;

  if (iPrint > 0) {
    cout << "Presolve started ..." << endl;
    cout << "Original problem ... N=" << numCol << "  M=" << numRow << endl;
  }

  if (iPrint < 0) {
    stringstream ss;
    ss << "dev-presolve: model:      rows, colx, nnz , " << modelName << ":  "
       << numRow << ",  " << numCol << ",  " << (HighsInt)Avalue.size()
       << std::endl;
    reportDev(ss.str());
  }

  initializeVectors();
  if (status) return status;

  // removeFixed();
  // if (status) return status;

  HighsInt iter = 1;
  if (order.size() == 0) {
    // pre_release_order:
    order.push_back(Presolver::kMainEmpty);
    order.push_back(Presolver::kMainRowSingletons);
    order.push_back(Presolver::kMainForcing);
    order.push_back(Presolver::kMainRowSingletons);
    order.push_back(Presolver::kMainDoubletonEq);
    order.push_back(Presolver::kMainRowSingletons);
    order.push_back(Presolver::kMainColSingletons);
    order.push_back(Presolver::kMainDominatedCols);
    if (mip) order.push_back(Presolver::kMainMipDualFixing);
    // wip
    // order.push_back(Presolver::kMainSingletonsOnly);
  }

  const double reduction_pct_for_further_presolve_iteration = 0.05;
  HighsInt model_cols_rows = numCol + numRow;
  HighsInt prev_cols_rows;
  HighsInt current_cols_rows = model_cols_rows;
  //  double prev_diff = 0;
  // Else: The order has been modified for experiments
  while (hasChange == 1) {
    if (max_iterations > 0 && iter > max_iterations) break;
    hasChange = false;
    //    printf("presolve iteration %" HIGHSINT_FORMAT " (max=%"
    //    HIGHSINT_FORMAT ")\n", iter, max_iterations);
    reportDevMainLoop();
    timer.recordStart(kRunPresolvers);
    HighsInt run_status = runPresolvers(order);
    timer.recordFinish(kRunPresolvers);
    assert(run_status == status);
    if (run_status) return status;

    // Exit check
    prev_cols_rows = current_cols_rows;
    current_cols_rows = 0;
    HighsInt current_num_col = 0;
    HighsInt current_num_row = 0;
    for (HighsInt i = 0; i < numRow; i++)
      if (flagRow[i]) current_cols_rows++;
    current_num_row = current_cols_rows;
    for (HighsInt i = 0; i < numCol; i++)
      if (flagCol[i]) current_cols_rows++;
    current_num_col = current_cols_rows - current_num_row;
    HighsInt diff = prev_cols_rows - current_cols_rows;
    double iteration_reduction_pct =
        100 * (1.0 * diff) / (1.0 * model_cols_rows);
    highsLogDev(log_options, HighsLogType::kVerbose,
                "Iteration %2" HIGHSINT_FORMAT
                " (Presolve)   Current number rows = %9" HIGHSINT_FORMAT
                "; cols = %9" HIGHSINT_FORMAT
                ": "
                "Reduction this iteration (%9" HIGHSINT_FORMAT ") is %5.2f%%\n",
                iter, current_num_row, current_num_col, diff,
                iteration_reduction_pct);
    if (current_cols_rows == 0) break;
    iter++;

    if (hasChange) {
      if (iteration_reduction_pct <
          reduction_pct_for_further_presolve_iteration)
        hasChange = false;
    }

    if (!hasChange && !aggregatorCalled) {
      aggregatorCalled = true;
      runAggregator();
      if (mip) detectImpliedIntegers();
      runPropagator();
      prev_cols_rows = current_cols_rows;
      current_cols_rows = 0;
      HighsInt current_num_col = 0;
      HighsInt current_num_row = 0;
      for (HighsInt i = 0; i < numRow; i++)
        if (flagRow[i]) current_num_row++;
      for (HighsInt i = 0; i < numCol; i++)
        if (flagCol[i]) current_num_col++;
      current_cols_rows = current_num_col + current_num_row;
      HighsInt diff = prev_cols_rows - current_cols_rows;
      double iteration_reduction_pct =
          100 * (1.0 * diff) / (1.0 * model_cols_rows);
      highsLogDev(log_options, HighsLogType::kVerbose,
                  "Iteration %2" HIGHSINT_FORMAT
                  " (Aggregator) Current number rows = %9" HIGHSINT_FORMAT
                  "; cols = %9" HIGHSINT_FORMAT
                  ": "
                  "Reduction this iteration (%9" HIGHSINT_FORMAT
                  ") is %5.2f%%\n",
                  iter, current_num_row, current_num_col, diff,
                  iteration_reduction_pct);
      iter++;
    }
  }

  // std::cout << "   MAIN LOOP ITER = " << iter << std::endl;

  reportDevMainLoop();

  timer.recordStart(kResizeMatrix);
  checkForChanges(iter);
  timer.recordFinish(kResizeMatrix);

  timer.updateInfo();

  if (iPrint != 0) printDevStats(stats.dev);

  return status;
}

HighsPresolveStatus Presolve::presolve() {
  timer.recordStart(kTotalPresolveTime);
  HighsPresolveStatus presolve_status = HighsPresolveStatus::kNotReduced;
  HighsInt result = presolve(0);
  switch (result) {
    case Stat::kUnboundedOrInfeasible:
      presolve_status = HighsPresolveStatus::kUnboundedOrInfeasible;
      break;
    case Stat::kInfeasible:
      presolve_status = HighsPresolveStatus::kInfeasible;
      break;
    case Stat::kReduced:
      if (numCol > 0 || numRow > 0)
        presolve_status = HighsPresolveStatus::kReduced;
      else
        presolve_status = HighsPresolveStatus::kReducedToEmpty;
      break;
    case Stat::kEmpty:
      presolve_status = HighsPresolveStatus::kEmpty;
      break;
    case Stat::kOptimal:
      // reduced problem solution indicated as optimal by
      // the solver.
      break;
    case Stat::kTimeout:
      presolve_status = HighsPresolveStatus::kTimeout;
  }
  timer.recordFinish(kTotalPresolveTime);
  if (iPrint > 0) {
    timer.reportClocks();
    timer.reportNumericsRecords();
  }
  return presolve_status;
}

void Presolve::checkBoundsAreConsistent() {
  for (HighsInt col = 0; col < numCol; col++) {
    if (flagCol[col]) {
      // Analyse dependency on numerical tolerance
      timer.updateNumericsRecord(kNumericsInconsistentBounds,
                                 colLower[col] - colUpper[col]);
      roundIntegerBounds(col);
      if (colLower[col] - colUpper[col] > inconsistent_bounds_tolerance) {
        status = kInfeasible;
        return;
      }
    }
  }

  for (HighsInt row = 0; row < numRow; row++) {
    if (flagRow[row]) {
      // Analyse dependency on numerical tolerance
      timer.updateNumericsRecord(kNumericsInconsistentBounds,
                                 rowLower[row] - rowUpper[row]);
      if (rowLower[row] - rowUpper[row] > inconsistent_bounds_tolerance) {
        status = kInfeasible;
        return;
      }
    }
  }
}

/**
 * returns <x, y>
 * 		   <x, -1> if we need to skip row
 *
 * 		   row is of form akx_x + aky_y = b,
 */
pair<HighsInt, HighsInt> Presolve::getXYDoubletonEquations(const HighsInt row) {
  // todo, for mip presolve also check integrality of right hand side value to
  // detect integer feasibility
  pair<HighsInt, HighsInt> colIndex;
  // row is of form akx_x + aky_y = b, where k=row and y is present in fewer
  // constraints

  HighsInt col1 = -1;
  HighsInt col2 = -1;
  double val1 = 0.0;
  double val2 = 0.0;
  HighsInt kk = ARstart.at(row);
  while (kk < ARstart.at(row + 1)) {
    if (flagCol.at(ARindex.at(kk))) {
      if (col1 == -1) {
        col1 = ARindex.at(kk);
        val1 = std::abs(ARvalue[kk]);
      } else if (col2 == -1) {
        col2 = ARindex.at(kk);
        val2 = std::abs(ARvalue[kk]);
      } else {
        cout << "ERROR: doubleton eq row" << row
             << " has more than two variables. \n";
        col2 = -2;
        break;
      }
      ++kk;
    } else
      ++kk;
  }
  if (col2 == -1)
    cout << "ERROR: doubleton eq row" << row
         << " has less than two variables. \n";
  if (col2 < 0) {
    colIndex.second = -1;
    return colIndex;
  }

  HighsInt x, y;
  if (mip && (integrality[col1] == HighsVarType::kInteger ||
              integrality[col2] == HighsVarType::kInteger)) {
    if (integrality[col1] != integrality[col2]) {
      // only one of the columns is integral, select the non-integral column to
      // be substituted out
      if (integrality[col1] == HighsVarType::kInteger) {
        // printf("column %" HIGHSINT_FORMAT " integral, column %"
        // HIGHSINT_FORMAT " is not\n", col1, col2);
        assert(integrality[col2] != HighsVarType::kInteger);
        y = col2;
        x = col1;
      } else {
        // printf("column %" HIGHSINT_FORMAT " integral, column %"
        // HIGHSINT_FORMAT " is not\n", col2, col1);
        assert(integrality[col2] == HighsVarType::kInteger);
        assert(integrality[col1] != HighsVarType::kInteger);
        y = col1;
        x = col2;
      }
    } else {
      // printf(
      //     "column %" HIGHSINT_FORMAT " with coefficient %f is integral and
      //     column %" HIGHSINT_FORMAT " with " "coefficient %f too\n", col1,
      //     val1, col2, val2);
      // both columns are integral, need to choose column with smaller
      // coefficient value, if coefficients are the same choose the column with
      // fewer nonzeros
      if (val1 < val2) {
        y = col1;
        x = col2;
      } else if (val2 < val1) {
        y = col2;
        x = col1;
      } else if (nzCol.at(col1) <= nzCol.at(col2)) {
        y = col1;
        x = col2;
      } else {
        x = col1;
        y = col2;
      }
    }
  } else if (nzCol.at(col1) <= nzCol.at(col2)) {
    y = col1;
    x = col2;
  } else {
    x = col1;
    y = col2;
  }

  colIndex.first = x;
  colIndex.second = y;
  return colIndex;
}

void Presolve::processRowDoubletonEquation(const HighsInt row, const HighsInt x,
                                           const HighsInt y, const double akx,
                                           const double aky, const double b) {
  // std::cout << "col 2... c = " << colCost.at(2)<< std::endl;
  // presolve::printCol(2, numRow, numCol, flagRow, flagCol, colLower,
  //                    colUpper, valueRowDual, Astart, Aend, Aindex, Avalue);

  postValue.push(akx);
  postValue.push(aky);
  postValue.push(b);

  // modify bounds on variable x (j), variable y (col,k) is substituted out
  // double aik = Avalue.at(k);
  // double aij = Avalue.at(kk);
  pair<double, double> p = getNewBoundsDoubletonConstraint(row, y, x, aky, akx);
  double low = p.first;
  double upp = p.second;

  // add old bounds of x to checker and for postsolve
  if (iKKTcheck == 1) {
    vector<pair<HighsInt, double>> bndsL, bndsU, costS;
    bndsL.push_back(make_pair(x, colLower.at(x)));
    bndsU.push_back(make_pair(x, colUpper.at(x)));
    costS.push_back(make_pair(x, colCost.at(x)));

    chk2.cLowers.push(bndsL);
    chk2.cUppers.push(bndsU);
    chk2.costs.push(costS);
  }

  vector<double> bnds({colLower.at(y), colUpper.at(y), colCost.at(y)});
  vector<double> bnds2({colLower.at(x), colUpper.at(x), colCost.at(x)});
  oldBounds.push(make_pair(y, bnds));
  oldBounds.push(make_pair(x, bnds2));

  if (low > colLower.at(x)) colLower.at(x) = low;
  if (upp < colUpper.at(x)) colUpper.at(x) = upp;

  roundIntegerBounds(x);

  // modify cost of xj
  colCost.at(x) = colCost.at(x) - colCost.at(y) * akx / aky;
  objShift += colCost.at(y) * b / aky;

  // for postsolve: need the new bounds too
  assert(x >= 0 && x < numCol);
  vector<double> bnds3({colLower.at(x), colUpper.at(x), colCost.at(x)});
  oldBounds.push(make_pair(x, bnds3));

  addChange(kDoubletonEquation, row, y);

  // remove y (col) and the row
  if (iPrint > 0)
    cout << "PR: Doubleton equation removed. Row " << row << ", column " << y
         << ", column left is " << x << "    nzy=" << nzCol.at(y) << endl;

  flagRow.at(row) = 0;
  nzCol.at(x)--;

  countRemovedRows(kDoubletonEquation);
  countRemovedCols(kDoubletonEquation);

  //----------------------------
  flagCol.at(y) = 0;
  if (!hasChange) hasChange = true;
}

void Presolve::caseTwoSingletonsDoubletonInequality(const HighsInt row,
                                                    const HighsInt x,
                                                    const HighsInt y) {
  // std::cout << "Call caseTwoSing..." << std::endl;

  // std::cout << "Two column singletons: row " << row << ", x = " << x << ", y
  // = " << y << std::endl; std::cout << "                     cx = " <<
  // colCost[x] << "  cy = " << colCost[y] << std::endl; std::cout << " ax = "
  // << getaij(row, x) << "  ay = " << getaij(row, y) << std::endl; std::cout <<
  // "   L = " << rowLower[row] << "  U = " << rowUpper[row] << std::endl;
  // std::cout << "   lx = " << colLower[x] << "  ux = " << colUpper[x] <<
  // std::endl; std::cout << "   ly = " << colLower[y] << "  uy = " <<
  // colUpper[y] << std::endl;

  assert(nzRow[row] = 2);
  assert(nzCol[x] = 1);
  assert(nzCol[y] = 1);

  assert(flagCol[x]);
  assert(flagCol[y]);
  assert(flagRow[row]);

  // trivial case
  // if(rowLower[row] == 0 && rowUpper[row] == 0) {
  //   if ((colLower[x] <= 0 && colUpper[x] >= 0) &&
  //       (colLower[y] <= 0 && colUpper[y] >= 0)) {

  //     // primal and dual values set to 0 already. just flagRow
  //     flagRow[row] = false;
  //     flagCol[x] = false;
  //     flagCol[y] = false;
  //     postValue.push((double)y);
  //     addChange(PresolveRule::kTwoColSingTrivial, row, x);
  //     std::cout << "Trivial case row " << row << std::endl;
  //   }
  // }
}

void Presolve::removeDoubletonEquations() {
  if (timer.reachLimit()) {
    status = Stat::kTimeout;
    return;
  }
  timer.recordStart(kDoubletonEquation);
  // flagCol should have one more element at end which is zero
  // needed for AR matrix manipulation
  if ((HighsInt)flagCol.size() == numCol) flagCol.push_back(0);

  HighsInt iter = 0;

  for (HighsInt row = 0; row < numRow; row++) {
    if (flagRow.at(row)) {
      // Analyse dependency on numerical tolerance
      if (nzRow.at(row) == 2 && rowLower[row] > -kHighsInf &&
          rowUpper[row] < kHighsInf) {
        // Possible doubleton equation
        timer.updateNumericsRecord(kNumericsDoubletonEquationBound,
                                   fabs(rowLower[row] - rowUpper[row]));
      }
      if (nzRow.at(row) == 2 && rowLower[row] > -kHighsInf &&
          rowUpper[row] < kHighsInf &&
          fabs(rowLower[row] - rowUpper[row]) <=
              doubleton_equation_bound_tolerance) {
        // row is of form akx_x + aky_y = b, where k=row and y is present in
        // fewer constraints
        const double b = rowLower.at(row);
        pair<HighsInt, int> colIndex = getXYDoubletonEquations(row);
        const HighsInt x = colIndex.first;
        const HighsInt y = colIndex.second;

        if (x >= 0 && y == -1) {
          // no second variable
          nzRow[row]--;
          continue;
        }

        // two singletons case handled elsewhere
        if (y < 0 || ((nzCol.at(y) == 1 && nzCol.at(x) == 1))) {
          caseTwoSingletonsDoubletonInequality(row, x, y);
          continue;
        }

        // singleton rows only in y column which is present in fewer constraints
        // and eliminated. bool rs_only = true; for (HighsInt k = Astart.at(y);
        // k < Aend.at(y); ++k)
        //   if (flagRow.at(Aindex.at(k)) && Aindex.at(k) != row) {
        //     if (nzRow[row]  > 1) {
        //       rs_only = false;
        //       break;
        //     }
        //   }
        // if (rs_only) continue;

        const double akx = getaij(row, x);
        const double aky = getaij(row, y);
        processRowDoubletonEquation(row, x, y, akx, aky, b);
        if (status) {
          timer.recordFinish(kDoubletonEquation);
          return;
        }

        // printRow(row, numRow, numCol, flagRow, flagCol, rowLower, rowUpper,
        //          valuePrimal, ARstart, ARindex, ARvalue);

        for (HighsInt k = Astart.at(y); k < Aend.at(y); ++k)
          if (flagRow.at(Aindex.at(k)) && Aindex.at(k) != row) {
            const HighsInt i = Aindex.at(k);
            const double aiy = Avalue.at(k);

            // update row bounds
            if (iKKTcheck == 1) {
              vector<pair<HighsInt, double>> bndsL, bndsU;
              bndsL.push_back(make_pair(i, rowLower.at(i)));
              bndsU.push_back(make_pair(i, rowUpper.at(i)));
              chk2.rLowers.push(bndsL);
              chk2.rUppers.push(bndsU);
              addChange(kDoubletonEquationRowBoundsUpdate, i, y);
            }

            if (rowLower.at(i) > -kHighsInf) rowLower.at(i) -= b * aiy / aky;
            if (rowUpper.at(i) < kHighsInf) rowUpper.at(i) -= b * aiy / aky;

            if (implRowValueLower.at(i) > -kHighsInf)
              implRowValueLower.at(i) -= b * aiy / aky;
            if (implRowValueUpper.at(i) < kHighsInf)
              implRowValueUpper.at(i) -= b * aiy / aky;

            // update matrix coefficients
            if (isZeroA(i, x)) {
              UpdateMatrixCoeffDoubletonEquationXzero(i, x, y, aiy, akx, aky);
              // std::cout << "   . row " << i << " zero " << std::endl;
            } else {
              UpdateMatrixCoeffDoubletonEquationXnonZero(i, x, y, aiy, akx,
                                                         aky);
              // std::cout << "   . row " << i << " zero " << std::endl;
            }
          }
        if (Avalue.size() > 40000000) {
          trimA();
        }

        iter++;
      }
    }
  }
  timer.recordFinish(kDoubletonEquation);
}

void Presolve::UpdateMatrixCoeffDoubletonEquationXzero(
    const HighsInt i, const HighsInt x, const HighsInt y, const double aiy,
    const double akx, const double aky) {
  // case x is zero initially
  // row nonzero count doesn't change here
  // cout<<"case: x not present "<<i<<" "<<endl;

  // update AR
  HighsInt ind;
  for (ind = ARstart.at(i); ind < ARstart.at(i + 1); ++ind)
    if (ARindex.at(ind) == y) {
      break;
    }

  assert(ARvalue.at(ind) == aiy);

  postValue.push(aiy);
  postValue.push(y);
  addChange(kDoubletonEquationXZeroInitially, i, x);

  ARindex.at(ind) = x;
  ARvalue.at(ind) = -aiy * akx / aky;

  // update A: append X column to end of array
  const HighsInt st = Avalue.size();
  for (HighsInt ind = Astart.at(x); ind < Aend.at(x); ++ind) {
    Avalue.push_back(Avalue.at(ind));
    Aindex.push_back(Aindex.at(ind));
  }
  Avalue.push_back(-aiy * akx / aky);
  Aindex.push_back(i);
  Astart.at(x) = st;
  Aend.at(x) = Avalue.size();

  nzCol.at(x)++;
  // nzRow does not change here.
}

void Presolve::UpdateMatrixCoeffDoubletonEquationXnonZero(
    const HighsInt i, const HighsInt x, const HighsInt y, const double aiy,
    const double akx, const double aky) {
  HighsInt ind;

  // update nonzeros: for removal of
  nzRow.at(i)--;
  if (nzRow.at(i) == 1) singRow.push_back(i);

  if (nzRow.at(i) == 0) {
    // singRow.remove(i);
    removeEmptyRow(i);
    countRemovedRows(kDoubletonEquation);
  }

  double xNew;
  for (ind = ARstart.at(i); ind < ARstart.at(i + 1); ++ind)
    if (ARindex.at(ind) == x) break;

  xNew = ARvalue.at(ind) - (aiy * akx) / aky;
  // Analyse dependency on numerical tolerance
  timer.updateNumericsRecord(kNumericsSmallMatrixValue, fabs(xNew));
  if (fabs(xNew) > presolve_small_matrix_value) {
    // case new x != 0
    // cout<<"case: x still there row "<<i<<" "<<endl;

    postValue.push(ARvalue.at(ind));
    addChange(kDoubletonEquationNewXNonzero, i, x);
    ARvalue.at(ind) = xNew;

    // update A:
    for (ind = Astart.at(x); ind < Aend.at(x); ++ind)
      if (Aindex.at(ind) == i) {
        break;
      }
    Avalue.at(ind) = xNew;
  } else {
    // case new x == 0
    // cout<<"case: x also disappears from row "<<i<<" "<<endl;
    // update nz row
    nzRow.at(i)--;
    // update singleton row list
    if (nzRow.at(i) == 1) singRow.push_back(i);

    if (nzRow.at(i) == 0) {
      removeEmptyRow(i);
      countRemovedRows(kDoubletonEquation);
    }

    if (nzRow.at(i) > 0) {
      // AR update
      // set ARindex of element for x to numCol
      // flagCol[numCol] = false
      // mind when resizing: should be OK
      postValue.push(ARvalue.at(ind));

      ARindex.at(ind) = numCol;

      addChange(kDoubletonEquationNewXZeroArUpdate, i, x);
    }

    if (nzCol.at(x) > 0) {
      // A update for case when x is zero: move x entry to end and set
      // Aend to be Aend - 1;
      HighsInt indi;
      for (indi = Astart.at(x); indi < Aend.at(x); ++indi)
        if (Aindex.at(indi) == i) break;

      postValue.push(Avalue.at(indi));

      // if indi is not Aend-1 swap elements indi and Aend-1
      if (indi != Aend.at(x) - 1) {
        double tmp = Avalue.at(Aend.at(x) - 1);
        HighsInt tmpi = Aindex.at(Aend.at(x) - 1);
        Avalue.at(Aend.at(x) - 1) = Avalue.at(indi);
        Aindex.at(Aend.at(x) - 1) = Aindex.at(indi);
        Avalue.at(indi) = tmp;
        Aindex.at(indi) = tmpi;
      }
      Aend.at(x)--;
      addChange(kDoubletonEquationNewXZeroAUpdate, i, x);
    }

    // update nz col
    nzCol.at(x)--;
    // update singleton col list
    if (nzCol.at(x) == 1) singCol.push_back(x);
    if (nzCol.at(x) == 0) {
      removeEmptyColumn(x);
    }
  }
}

void Presolve::trimA() {
  HighsInt cntEl = 0;
  for (HighsInt j = 0; j < numCol; ++j)
    if (flagCol.at(j)) cntEl += nzCol.at(j);

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

  if (nR + nC == 0) {
    status = kEmpty;
    return;
  }

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

void Presolve::initializeVectors() {
  // copy original bounds
  colCostOriginal = colCost;
  rowUpperOriginal = rowUpper;
  rowLowerOriginal = rowLower;
  colUpperOriginal = colUpper;
  colLowerOriginal = colLower;

  makeARCopy();

  valueRowDual.resize(numRow);
  valuePrimal.resize(numCol);
  valueColDual.resize(numCol);

  flagCol.assign(numCol, 1);
  flagRow.assign(numRow, 1);

  if (iKKTcheck) setKKTcheckerData();

  nzCol.assign(numCol, 0);
  nzRow.assign(numRow, 0);

  for (HighsInt i = 0; i < numRow; ++i) {
    nzRow.at(i) = ARstart.at(i + 1) - ARstart.at(i);
    if (nzRow.at(i) == 1) singRow.push_back(i);
    if (nzRow.at(i) == 0) {
      timer.recordStart(kEmptyRow);
      removeEmptyRow(i);
      countRemovedRows(kEmptyRow);
      timer.recordFinish(kEmptyRow);
    }
  }

  Aend.resize(numCol + 1);
  for (HighsInt i = 0; i < numCol; ++i) {
    Aend.at(i) = Astart.at(i + 1);
    nzCol.at(i) = Aend.at(i) - Astart.at(i);
    if (nzCol.at(i) == 1) singCol.push_back(i);
  }

  implColUpper = colUpper;  // working copies of primal variable bounds
  implColLower = colLower;
  implColLowerRowIndex.assign(numCol, -1);
  implColUpperRowIndex.assign(numCol, -1);

  implRowDualLowerSingColRowIndex.assign(numRow, -1);
  implRowDualUpperSingColRowIndex.assign(numRow, -1);
  implRowDualLower.assign(numRow, -kHighsInf);
  implRowDualUpper.assign(numRow, kHighsInf);

  implColDualLower.assign(numCol, -kHighsInf);
  implColDualUpper.assign(numCol, kHighsInf);
  implRowValueLower = rowLower;
  implRowValueUpper = rowUpper;

  for (HighsInt i = 0; i < numRow; ++i) {
    if (rowLower.at(i) <= -kHighsInf) implRowDualUpper.at(i) = 0;
    if (rowUpper.at(i) >= kHighsInf) implRowDualLower.at(i) = 0;
  }

  for (HighsInt i = 0; i < numCol; ++i) {
    if (colLower.at(i) <= -kHighsInf) implColDualUpper.at(i) = 0;
    if (colUpper.at(i) >= kHighsInf) implColDualLower.at(i) = 0;
  }

  colCostAtEl = colCost;

  // initialize integrality information
  if (!mip) integrality.assign(numCol, HighsVarType::kContinuous);
}

void Presolve::runAggregator() {
  // run the aggregator and store back the modified matrix
  timer.recordStart(kAggregator);
  aggregatorStack.emplace_back();

  aggregatorStack.back().colCostAtCall = colCost;

  {
    HAggregator aggregator(rowLower, rowUpper, colCost, objShift, integrality,
                           colLower, colUpper);

    aggregator.fromDynamicCSC(Avalue, Aindex, Astart, Aend, flagRow, flagCol);
    aggregatorStack.back().postsolveStack = aggregator.run();

    if (aggregatorStack.back().postsolveStack.empty()) {
      if (mip) {
        // todo, we need to store the state of the matrix before
        // doing mip presolve steps that alter the coefficients,
        // currently only coefficient tightening.
        // this is a bit hacky but works for now as those steps
        // run after the first call of the aggregator
        aggregatorStack.back().ARvalueAtCall = std::move(ARvalue);
        aggregatorStack.back().ARindexAtCall = std::move(ARindex);
        aggregatorStack.back().ARstartAtCall = std::move(ARstart);
        aggregator.toCSC(Avalue, Aindex, Astart);
        makeARCopy();
        Aend.resize(numCol + 1);
        std::copy(Astart.begin() + 1, Astart.end(), Aend.begin());
        chng.emplace(change{kAggregator, 0, 0});
      } else
        aggregatorStack.pop_back();
      timer.recordFinish(kAggregator);
      return;
    }
    aggregator.toCSC(Avalue, Aindex, Astart);
  }

  hasChange = true;
  chng.emplace(change{kAggregator, 0, 0});

  aggregatorStack.back().postsolveStack.unsetFlags(flagRow, flagCol);

  // store old AR and create updated AR matrix
  aggregatorStack.back().ARvalueAtCall = std::move(ARvalue);
  aggregatorStack.back().ARindexAtCall = std::move(ARindex);
  aggregatorStack.back().ARstartAtCall = std::move(ARstart);
  makeARCopy();

  // set nzRow and nzCol, remove empty rows
  nzCol.assign(numCol, 0);
  nzRow.assign(numRow, 0);

  for (HighsInt i = 0; i < numRow; ++i) {
    if (!flagRow[i]) continue;
    nzRow.at(i) = ARstart.at(i + 1) - ARstart.at(i);
    if (nzRow.at(i) == 1) singRow.push_back(i);
    if (nzRow.at(i) == 0) {
      timer.recordStart(kEmptyRow);
      removeEmptyRow(i);
      countRemovedRows(kEmptyRow);
      timer.recordFinish(kEmptyRow);
    }
  }

  Aend.resize(numCol + 1);
  for (HighsInt i = 0; i < numCol; ++i) {
    if (!flagCol[i]) continue;
    Aend.at(i) = Astart.at(i + 1);
    nzCol.at(i) = Aend.at(i) - Astart.at(i);
    if (nzCol.at(i) == 1) singCol.push_back(i);
  }

  implColUpper = colUpper;  // working copies of primal variable bounds
  implColLower = colLower;
  implColLowerRowIndex.assign(numCol, -1);
  implColUpperRowIndex.assign(numCol, -1);

  implRowDualLowerSingColRowIndex.assign(numRow, -1);
  implRowDualUpperSingColRowIndex.assign(numRow, -1);
  implRowDualLower.assign(numRow, -kHighsInf);
  implRowDualUpper.assign(numRow, kHighsInf);

  implColDualLower.assign(numCol, -kHighsInf);
  implColDualUpper.assign(numCol, kHighsInf);
  implRowValueLower = rowLower;
  implRowValueUpper = rowUpper;

  for (HighsInt i = 0; i < numRow; ++i) {
    if (rowLower.at(i) <= -kHighsInf) implRowDualUpper.at(i) = 0;
    if (rowUpper.at(i) >= kHighsInf) implRowDualLower.at(i) = 0;
  }

  for (HighsInt i = 0; i < numCol; ++i) {
    if (colLower.at(i) <= -kHighsInf) implColDualUpper.at(i) = 0;
    if (colUpper.at(i) >= kHighsInf) implColDualLower.at(i) = 0;
  }

  timer.recordFinish(kAggregator);
}

void Presolve::runPropagator() {
  // run the propagator to compute tightened bounds
  // which can help with other reductions like redundant rows or forcing
  // constraints and also decrease the number of total iterations required due
  // to their exploitation with the long-step rule for the dual simplex ratio
  // test
  HighsLpPropagator propagator(colLower, colUpper, integrality, Avalue, Aindex,
                               Astart, Aend, ARvalue, ARindex, ARstart, flagRow,
                               flagCol, rowLower, rowUpper);
  propagator.computeRowActivities();
  HighsInt nboundchgs = propagator.propagate();
  highsLogDev(log_options, HighsLogType::kVerbose,
              "Propagation found %" HIGHSINT_FORMAT " bound changes\n",
              nboundchgs);
  // propagation found nothing, so we can stop here. Only for mip we also try
  // coefficient tightening
  if (!mip && nboundchgs == 0) return;

  if (mip) {
    HighsInt ntotalcoeffchgs = 0;
    while (true) {
      HighsInt ncoeffchgs = propagator.tightenCoefficients();
      highsLogDev(log_options, HighsLogType::kVerbose,
                  "tightened %" HIGHSINT_FORMAT " coefficients\n", ncoeffchgs);
      // if no coefficients where tightened we can stop
      if (ncoeffchgs == 0) break;
      ntotalcoeffchgs += ncoeffchgs;

      hasChange = true;
      nboundchgs = propagator.propagate();

      if (propagator.infeasible()) {
        status = kInfeasible;
        return;
      }
      // if no further bounds where changed we can stop
      if (nboundchgs == 0) break;
    }

    if (ntotalcoeffchgs != 0) {
      // row sides might have changed during coefficient tightening
      implRowValueLower = rowLower;
      implRowValueUpper = rowUpper;
    }

    // no bounds where tightened, so we can return here
    if (propagator.getNumChangedBounds() == 0) return;
  }

  HighsInt ntightened = 0;
  // we cannot use the tightest bounds that we obtained by propagation
  // as then the dual postsolve step might not work anymore. Instead
  // we relax the bounds by a wide enough margin so that they cannot
  // be used in a basic feasible solution
  for (HighsInt i = 0; i != numCol; ++i) {
    if (!flagCol[i]) continue;

    if (colLower[i] >= propagator.colLower_[i] &&
        colUpper[i] <= propagator.colUpper_[i])
      continue;

    if (mip) {
      if (colLower[i] < propagator.colLower_[i]) {
        colLower[i] = propagator.colLower_[i];
        ++ntightened;
      }

      if (colUpper[i] > propagator.colUpper_[i]) {
        colUpper[i] = propagator.colUpper_[i];
        ++ntightened;
      }

      roundIntegerBounds(i);
      if (std::abs(colUpper[i] - colLower[i]) <= fixed_column_tolerance)
        removeFixedCol(i);

      continue;
    }

    HighsInt start = Astart[i];
    HighsInt end = Aend[i];
    double minabs = 1.0;

    for (HighsInt j = start; j != end; ++j) {
      HighsInt row = Aindex[j];
      if (!flagRow[row]) continue;

      minabs = std::min(minabs, std::abs(Avalue[j]));
    }

    // use a big enough factor of the feasibility tolerance and further increase
    // it if the column has coefficient values below 1.0 so that it is ensured
    // that no primal feasible solution can have this variable sitting at this
    // artifical bound
    double margin = 128.0 * default_primal_feasiblility_tolerance / minabs;

    // now widen the bounds and check if they are tighter than the previous
    // bounds
    double abslb = std::abs(propagator.colLower_[i]);

    if (abslb <= 1e8) {
      // use a relative minimal margin to ensure large bounds are relaxed enough
      double minlowermargin = default_primal_feasiblility_tolerance * abslb;
      propagator.colLower_[i] -= std::max(margin, minlowermargin);
      if (propagator.colLower_[i] > colLower[i]) {
        colLower[i] = propagator.colLower_[i];
        ++ntightened;
      }
    }

    double absub = std::abs(propagator.colUpper_[i]);

    if (absub <= 1e8) {
      // use a relative minimal margin to ensure large bounds are relaxed enough
      double minuppermargin = default_primal_feasiblility_tolerance * absub;
      propagator.colUpper_[i] += std::max(margin, minuppermargin);
      if (propagator.colUpper_[i] < colUpper[i]) {
        colUpper[i] = propagator.colUpper_[i];
        ++ntightened;
      }
    }
  }

  implColLower = colLower;
  implColUpper = colUpper;

  highsLogDev(log_options, HighsLogType::kVerbose,
              "Tightened %" HIGHSINT_FORMAT " bounds\n", ntightened);
  if (ntightened != 0) hasChange = true;
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

void Presolve::applyMipDualFixing() {
  for (HighsInt i = 0; i != numCol; ++i) {
    if (!flagCol[i] || integrality[i] != HighsVarType::kInteger) continue;

    HighsInt start = Astart[i];
    HighsInt end = Aend[i];
    HighsInt nuplocks = 0;
    HighsInt ndownlocks = 0;

    if (colCost[i] > 0 || colUpper[i] == kHighsInf) ++nuplocks;

    if (colCost[i] < 0 || colLower[i] == -kHighsInf) ++ndownlocks;

    if (ndownlocks != 0 && nuplocks != 0) continue;

    for (HighsInt j = start; j != end; ++j) {
      HighsInt row = Aindex[j];
      if (!flagRow[row]) continue;
      double lower;
      double upper;

      if (Avalue[j] < 0) {
        lower = -rowUpper[row];
        upper = -rowLower[row];
      } else {
        lower = rowLower[row];
        upper = rowUpper[row];
      }

      if (lower != -kHighsInf) ++ndownlocks;
      if (upper != kHighsInf) ++nuplocks;

      if (ndownlocks != 0 && nuplocks != 0) break;
    }

    if (ndownlocks == 0) {
      colUpper[i] = colLower[i];
      removeFixedCol(i);
      timer.increaseCount(false, kMipDualFixing);
    } else if (nuplocks == 0) {
      colLower[i] = colUpper[i];
      removeFixedCol(i);
      timer.increaseCount(false, kMipDualFixing);
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

void Presolve::removeEmptyColumn(HighsInt j) {
  flagCol.at(j) = 0;
  // singCol.remove(j);
  double value;
  if ((colCost.at(j) < 0 && colUpper.at(j) >= kHighsInf) ||
      (colCost.at(j) > 0 && colLower.at(j) <= -kHighsInf)) {
    if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
    status = kUnboundedOrInfeasible;
    return;
  }

  if (colCost.at(j) > 0)
    value = colLower.at(j);
  else if (colCost.at(j) < 0)
    value = colUpper.at(j);
  else if (colUpper.at(j) >= 0 && colLower.at(j) <= 0)
    value = 0;
  else if (colUpper.at(j) < 0)
    value = colUpper.at(j);
  else
    value = colLower.at(j);

  setPrimalValue(j, value);
  valueColDual.at(j) = colCost.at(j);

  addChange(kEmptyCol, 0, j);

  if (iPrint > 0)
    cout << "PR: Column: " << j
         << " eliminated: all nonzero rows have been removed. Cost = "
         << colCost.at(j) << ", value = " << value << endl;

  countRemovedCols(kEmptyCol);
}

void Presolve::rowDualBoundsDominatedColumns() {
  HighsInt col, i, k;

  // todo, do not use integer variables to derive bounds on the row duals
  // using continous and implied integer variables is valid as the complementary
  // slackness condition holds when all integer columns are fixed to any integer
  // value

  // for each row calc yihat and yibar and store in implRowDualLower and
  // implRowDualUpper
  for (list<HighsInt>::iterator it = singCol.begin(); it != singCol.end(); ++it)
    if (flagCol.at(*it)) {
      col = *it;
      if (mip && integrality[col] == HighsVarType::kInteger) continue;
      k = getSingColElementIndexInA(col);
      if (k < 0) continue;
      assert(k < (HighsInt)Aindex.size());
      i = Aindex.at(k);

      if (!flagRow.at(i)) {
        cout << "ERROR: column singleton " << col << " is in row " << i
             << " which is already mapped off\n";
        exit(-1);
      }

      if (colLower.at(col) <= -kHighsInf || colUpper.at(col) >= kHighsInf) {
        if (colLower.at(col) > -kHighsInf && colUpper.at(col) >= kHighsInf) {
          if (Avalue.at(k) > 0)
            if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
              implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
          if (Avalue.at(k) < 0)
            if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
              implRowDualLower.at(i) = colCost.at(col) / Avalue.at(k);
        } else if (colLower.at(col) <= -kHighsInf &&
                   colUpper.at(col) < kHighsInf) {
          if (Avalue.at(k) > 0)
            if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
              implRowDualUpper.at(i) = -colCost.at(col) / Avalue.at(k);
          if (Avalue.at(k) < 0)
            if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
              implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
        } else if (colLower.at(col) <= -kHighsInf &&
                   colUpper.at(col) >= kHighsInf) {
          // all should be removed earlier but use them
          if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
            implRowDualLower.at(i) = colCost.at(col) / Avalue.at(k);
          if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
            implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
        }

        if (implRowDualLower.at(i) > implRowDualUpper.at(i)) {
          cout << "Error: inconstistent bounds for Lagrange multiplier for row "
               << i << " detected after column singleton " << col
               << ". In presolve::dominatedColumns" << endl;
          exit(0);
        }
      }
    }
}

pair<double, double> Presolve::getImpliedColumnBounds(HighsInt j) {
  pair<double, double> out;
  double e = 0;
  double d = 0;

  HighsInt i;
  for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k) {
    i = Aindex.at(k);
    if (flagRow.at(i)) {
      if (Avalue.at(k) < 0) {
        if (implRowDualUpper.at(i) < kHighsInf)
          e += Avalue.at(k) * implRowDualUpper.at(i);
        else {
          e = -kHighsInf;
          break;
        }
      } else {
        if (implRowDualLower.at(i) > -kHighsInf)
          e += Avalue.at(k) * implRowDualLower.at(i);
        else {
          e = -kHighsInf;
          break;
        }
      }
    }
  }

  for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k) {
    i = Aindex.at(k);
    if (flagRow.at(i)) {
      if (Avalue.at(k) < 0) {
        if (implRowDualLower.at(i) > -kHighsInf)
          d += Avalue.at(k) * implRowDualLower.at(i);
        else {
          d = kHighsInf;
          break;
        }
      } else {
        if (implRowDualUpper.at(i) < kHighsInf)
          d += Avalue.at(k) * implRowDualUpper.at(i);
        else {
          d = kHighsInf;
          break;
        }
      }
    }
  }

  if (e > d) {
    cout << "Error: inconstistent bounds for Lagrange multipliers for column "
         << j << ": e>d. In presolve::dominatedColumns" << endl;
    exit(-1);
  }
  out.first = d;
  out.second = e;
  return out;
}

void Presolve::removeDominatedColumns() {
  // todo, use only continuous and implied integer columns for calculating row
  // dual upper and lower bounds then the domination (and also weak domination
  // check) should be correct for integers

  // for each column j calculate e and d and check:
  double e, d;
  pair<double, double> p;
  if (timer.reachLimit()) {
    status = Stat::kTimeout;
    return;
  }
  for (HighsInt j = 0; j < numCol; ++j)
    if (flagCol.at(j)) {
      p = getImpliedColumnBounds(j);
      d = p.first;
      e = p.second;

      // Analyse dependency on numerical tolerance
      bool dominated = colCost.at(j) - d > tol;
      timer.updateNumericsRecord(kNumericsDominatedColumn, colCost.at(j) - d);
      if (!dominated) {
        timer.updateNumericsRecord(kNumericsDominatedColumn, e - colCost.at(j));
      }

      // check if it is dominated
      if (colCost.at(j) - d > tol) {
        if (colLower.at(j) <= -kHighsInf) {
          if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
          status = kUnboundedOrInfeasible;
          return;
        }
        setPrimalValue(j, colLower.at(j));
        addChange(kDominatedCols, 0, j);
        if (iPrint > 0)
          cout << "PR: Dominated column " << j
               << " removed. Value := " << valuePrimal.at(j) << endl;
        countRemovedCols(kDominatedCols);
      } else if (colCost.at(j) - e < -tol) {
        if (colUpper.at(j) >= kHighsInf) {
          if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
          status = kUnboundedOrInfeasible;
          return;
        }
        setPrimalValue(j, colUpper.at(j));
        addChange(kDominatedCols, 0, j);
        if (iPrint > 0)
          cout << "PR: Dominated column " << j
               << " removed. Value := " << valuePrimal.at(j) << endl;
        countRemovedCols(kDominatedCols);
      } else {
        // update implied bounds
        if (implColDualLower.at(j) < (colCost.at(j) - d))
          implColDualLower.at(j) = colCost.at(j) - d;
        if (implColDualUpper.at(j) > (colCost.at(j) - e))
          implColDualUpper.at(j) = colCost.at(j) - e;
        if (implColDualLower.at(j) > implColDualUpper.at(j))
          cout << "INCONSISTENT\n";

        removeIfWeaklyDominated(j, d, e);
        continue;
      }
      if (status) return;
    }
}

void Presolve::removeIfWeaklyDominated(const HighsInt j, const double d,
                                       const double e) {
  HighsInt i;
  // check if it is weakly dominated: Excluding singletons!
  if (nzCol.at(j) > 1) {
    // Analyse dependency on numerical tolerance
    bool possible = d < kHighsInf && colLower.at(j) > -kHighsInf;
    timer.updateNumericsRecord(kNumericsWeaklyDominatedColumn,
                               fabs(colCost.at(j) - d));
    if (possible &&
        fabs(colCost.at(j) - d) < weakly_dominated_column_tolerance) {
      if (e > -kHighsInf && colUpper.at(j) < kHighsInf)
        timer.updateNumericsRecord(kNumericsWeaklyDominatedColumn,
                                   fabs(colCost.at(j) - e));
    }

    if (d < kHighsInf &&
        fabs(colCost.at(j) - d) < weakly_dominated_column_tolerance &&
        colLower.at(j) > -kHighsInf) {
      setPrimalValue(j, colLower.at(j));
      addChange(kWeaklyDominatedCols, 0, j);
      if (iPrint > 0)
        cout << "PR: Weakly Dominated column " << j
             << " removed. Value := " << valuePrimal.at(j) << endl;

      countRemovedCols(kWeaklyDominatedCols);
    } else if (e > -kHighsInf &&
               fabs(colCost.at(j) - e) < weakly_dominated_column_tolerance &&
               colUpper.at(j) < kHighsInf) {
      setPrimalValue(j, colUpper.at(j));
      addChange(kWeaklyDominatedCols, 0, j);
      if (iPrint > 0)
        cout << "PR: Weakly Dominated column " << j
             << " removed. Value := " << valuePrimal.at(j) << endl;

      countRemovedCols(kWeaklyDominatedCols);
    } else {
      double bnd;

      // calculate new bounds
      if (!mip || integrality[j] != HighsVarType::kInteger) {
        if (colLower.at(j) > -kHighsInf || colUpper.at(j) >= kHighsInf)
          for (HighsInt kk = Astart.at(j); kk < Aend.at(j); ++kk)
            if (flagRow.at(Aindex.at(kk)) && d < kHighsInf) {
              i = Aindex.at(kk);
              if (Avalue.at(kk) > 0 && implRowDualLower.at(i) > -kHighsInf) {
                bnd = -(colCost.at(j) + d) / Avalue.at(kk) +
                      implRowDualLower.at(i);
                if (bnd < implRowDualUpper.at(i) &&
                    !(bnd < implRowDualLower.at(i)))
                  implRowDualUpper.at(i) = bnd;
              } else if (Avalue.at(kk) < 0 &&
                         implRowDualUpper.at(i) < kHighsInf) {
                bnd = -(colCost.at(j) + d) / Avalue.at(kk) +
                      implRowDualUpper.at(i);
                if (bnd > implRowDualLower.at(i) &&
                    !(bnd > implRowDualUpper.at(i)))
                  implRowDualLower.at(i) = bnd;
              }
            }

        if (colLower.at(j) <= -kHighsInf || colUpper.at(j) < kHighsInf)
          for (HighsInt kk = Astart.at(j); kk < Aend.at(j); ++kk)
            if (flagRow.at(Aindex.at(kk)) && e > -kHighsInf) {
              i = Aindex.at(kk);
              if (Avalue.at(kk) > 0 && implRowDualUpper.at(i) < kHighsInf) {
                bnd = -(colCost.at(j) + e) / Avalue.at(kk) +
                      implRowDualUpper.at(i);
                if (bnd > implRowDualLower.at(i) &&
                    !(bnd > implRowDualUpper.at(i)))
                  implRowDualLower.at(i) = bnd;
              } else if (Avalue.at(kk) < 0 &&
                         implRowDualLower.at(i) > -kHighsInf) {
                bnd = -(colCost.at(j) + e) / Avalue.at(kk) +
                      implRowDualLower.at(i);
                if (bnd < implRowDualUpper.at(i) &&
                    !(bnd < implRowDualLower.at(i)))
                  implRowDualUpper.at(i) = bnd;
              }
            }
      }
    }
  }
}

void Presolve::setProblemStatus(const HighsInt s) {
  if (s == kInfeasible)
    cout << "NOT-OPT status = 1, returned from solver after presolve: Problem "
            "infeasible.\n";
  else if (s == kUnboundedOrInfeasible)
    cout << "NOT-OPT status = 2, returned from solver after presolve: Problem "
            "unbounded.\n";
  else if (s == 0) {
    status = kOptimal;
    return;
  } else
    cout << "unknown problem status returned from solver after presolve: " << s
         << endl;
  status = s;
}

void Presolve::setKKTcheckerData() {
  // after initializing equations.
  chk2.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);
}

pair<double, double> Presolve::getNewBoundsDoubletonConstraint(
    const HighsInt row, const HighsInt col, const HighsInt j, const double aik,
    const double aij) {
  HighsInt i = row;

  double upp = kHighsInf;
  double low = -kHighsInf;

  roundIntegerBounds(col);

  if (aij > 0 && aik > 0) {
    if (colLower.at(col) > -kHighsInf && rowUpper.at(i) < kHighsInf)
      upp = (rowUpper.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < kHighsInf && rowLower.at(i) > -kHighsInf)
      low = (rowLower.at(i) - aik * colUpper.at(col)) / aij;
  } else if (aij > 0 && aik < 0) {
    if (colLower.at(col) > -kHighsInf && rowLower.at(i) > -kHighsInf)
      low = (rowLower.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < kHighsInf && rowUpper.at(i) < kHighsInf)
      upp = (rowUpper.at(i) - aik * colUpper.at(col)) / aij;
  } else if (aij < 0 && aik > 0) {
    if (colLower.at(col) > -kHighsInf && rowUpper.at(i) < kHighsInf)
      low = (rowUpper.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < kHighsInf && rowLower.at(i) > -kHighsInf)
      upp = (rowLower.at(i) - aik * colUpper.at(col)) / aij;
  } else {
    if (colLower.at(col) > -kHighsInf && rowLower.at(i) > -kHighsInf)
      upp = (rowLower.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < kHighsInf && rowUpper.at(i) < kHighsInf)
      low = (rowUpper.at(i) - aik * colUpper.at(col)) / aij;
  }

  if (upp - low < -inconsistent_bounds_tolerance) {
    if (iPrint > 0)
      std::cout << "Presolve warning: inconsistent bounds in doubleton "
                   "constraint row "
                << row << std::endl;
  }

  return make_pair(low, upp);
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

void Presolve::removeFreeColumnSingleton(const HighsInt col, const HighsInt row,
                                         const HighsInt k) {
  if (iPrint > 0)
    cout << "PR: Free column singleton " << col << " removed. Row " << row
         << " removed." << endl;

  // modify costs
  vector<pair<HighsInt, double>> newCosts;
  HighsInt j;
  for (HighsInt kk = ARstart.at(row); kk < ARstart.at(row + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      newCosts.push_back(make_pair(j, colCost.at(j)));
      colCost.at(j) =
          colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);
    }
  }
  if (iKKTcheck == 1) chk2.costs.push(newCosts);

  flagCol.at(col) = 0;
  postValue.push(colCost.at(col));
  fillStackRowBounds(row);

  valueColDual.at(col) = 0;
  valueRowDual.at(row) = -colCost.at(col) / Avalue.at(k);

  double b = valueRowDual[row] < 0 ? rowLower[row] : rowUpper[row];
  objShift += colCost.at(col) * b / Avalue.at(k);

  addChange(kFreeSingCol, row, col);
  removeRow(row);

  countRemovedCols(kFreeSingCol);
  countRemovedRows(kFreeSingCol);
}

bool Presolve::removeColumnSingletonInDoubletonInequality(const HighsInt col,
                                                          const HighsInt i,
                                                          const HighsInt k) {
  // second column index j
  // second column row array index kk
  HighsInt j = -1;

  // count
  HighsInt kk = ARstart.at(i);
  while (kk < ARstart.at(i + 1)) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col)
      break;
    else
      ++kk;
  }
  if (kk == ARstart.at(i + 1))
    cout << "ERROR: nzRow[" << i << "]=2, but no second variable in row. \n";

  // only inequality case and case two singletons here,
  // others handled in doubleton equation
  // Analyse dependency on numerical tolerance
  if (nzCol.at(j) > 1)
    timer.updateNumericsRecord(kNumericsDoubletonInequalityBound,
                               fabs(rowLower.at(i) - rowUpper.at(i)));
  if ((fabs(rowLower.at(i) - rowUpper.at(i)) <
       doubleton_inequality_bound_tolerance) &&
      (nzCol.at(j) > 1))
    return false;

  // additional check if it is indeed implied free
  // needed since we handle inequalities and it may not be true
  // low and upp to be tighter than original bounds for variable col
  // so it is indeed implied free and we can remove it
  pair<double, double> p =
      getNewBoundsDoubletonConstraint(i, j, col, ARvalue.at(kk), Avalue.at(k));
  if (!(colLower.at(col) <= p.first && colUpper.at(col) >= p.second)) {
    return false;
  }

  postValue.push(ARvalue.at(kk));
  postValue.push(Avalue.at(k));

  // modify bounds on variable j, variable col (k) is substituted out
  // double aik = Avalue.at(k);
  // double aij = Avalue.at(kk);
  p = getNewBoundsDoubletonConstraint(i, col, j, Avalue.at(k), ARvalue.at(kk));
  double low = p.first;
  double upp = p.second;

  // add old bounds of xj to checker and for postsolve
  if (iKKTcheck == 1) {
    vector<pair<HighsInt, double>> bndsL, bndsU, costS;
    bndsL.push_back(make_pair(j, colLower.at(j)));
    bndsU.push_back(make_pair(j, colUpper.at(j)));
    costS.push_back(make_pair(j, colCost.at(j)));
    chk2.cLowers.push(bndsL);
    chk2.cUppers.push(bndsU);
    chk2.costs.push(costS);
  }

  vector<double> bndsCol({colLower.at(col), colUpper.at(col), colCost.at(col)});
  vector<double> bndsJ({colLower.at(j), colUpper.at(j), colCost.at(j)});
  oldBounds.push(make_pair(col, bndsCol));
  oldBounds.push(make_pair(j, bndsJ));

  // modify bounds of xj
  if (low > colLower.at(j)) colLower.at(j) = low;
  if (upp < colUpper.at(j)) colUpper.at(j) = upp;

  // modify cost of xj
  colCost.at(j) =
      colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);

  // for postsolve: need the new bounds too
  // oldBounds.push_back(colLower.at(j)); oldBounds.push_back(colUpper.at(j));
  bndsJ.at(0) = (colLower.at(j));
  bndsJ.at(1) = (colUpper.at(j));
  bndsJ.at(2) = (colCost.at(j));
  oldBounds.push(make_pair(j, bndsJ));

  // remove col as free column singleton
  if (iPrint > 0)
    cout << "PR: Column singleton " << col
         << " in a doubleton inequality constraint removed. Row " << i
         << " removed. variable left is " << j << endl;

  flagCol.at(col) = 0;
  fillStackRowBounds(i);
  countRemovedCols(kSingColDoubletonIneq);
  countRemovedRows(kSingColDoubletonIneq);

  valueColDual.at(col) = 0;
  valueRowDual.at(i) =
      -colCost.at(col) /
      Avalue.at(k);  // may be changed later, depending on bounds.
  addChange(kSingColDoubletonIneq, i, col);

  // if not special case two column singletons
  if (nzCol.at(j) > 1)
    removeRow(i);
  else if (nzCol.at(j) == 1)
    removeSecondColumnSingletonInDoubletonRow(j, i);

  return true;
}

void Presolve::removeSecondColumnSingletonInDoubletonRow(const HighsInt j,
                                                         const HighsInt i) {
  // case two singleton columns
  // when we get here bounds on xj are updated so we can choose low/upper one
  // depending on the cost of xj
  flagRow.at(i) = 0;
  double value;
  if (colCost.at(j) > 0) {
    if (colLower.at(j) <= -kHighsInf) {
      if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
      status = kUnboundedOrInfeasible;
      return;
    }
    value = colLower.at(j);
  } else if (colCost.at(j) < 0) {
    if (colUpper.at(j) >= kHighsInf) {
      if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
      status = kUnboundedOrInfeasible;
      return;
    }
    value = colUpper.at(j);
  } else {  //(colCost.at(j) == 0)
    if (colUpper.at(j) >= 0 && colLower.at(j) <= 0)
      value = 0;
    else if (fabs(colUpper.at(j)) < fabs(colLower.at(j)))
      value = colUpper.at(j);
    else
      value = colLower.at(j);
  }
  setPrimalValue(j, value);
  addChange(kSingColDoubletonIneqSecondSingCol, 0, j);
  if (iPrint > 0)
    cout << "PR: Second singleton column " << j << " in doubleton row " << i
         << " removed.\n";
  countRemovedCols(kSingColDoubletonIneq);
  // singCol.remove(j);
}

void Presolve::removeZeroCostColumnSingleton(const HighsInt col,
                                             const HighsInt row,
                                             const HighsInt k) {
  assert(Aindex[k] == row);
  assert(fabs(colCost[col]) < tol);
  std::cout << "Zero cost column singleton: col = " << col << ", row " << row
            << ", coeff = " << Avalue[k] << ", cost = " << colCost[col]
            << std::endl;
  std::cout << "   L = " << rowLower[row] << "  U = " << rowUpper[row]
            << std::endl;
  std::cout << "   l = " << colLower[col] << "  u = " << colUpper[col]
            << std::endl;
}

void Presolve::removeColumnSingletons() {
  list<HighsInt>::iterator it = singCol.begin();

  if (timer.reachLimit()) {
    status = Stat::kTimeout;
    return;
  }

  while (it != singCol.end()) {
    if (flagCol[*it]) {
      const HighsInt col = *it;
      assert(0 <= col && col <= numCol);

      const HighsInt k = getSingColElementIndexInA(col);
      if (k < 0) {
        it = singCol.erase(it);
        if (k == -2) flagCol[col] = 0;
        continue;
      }
      assert(k < (HighsInt)Aindex.size());
      const HighsInt i = Aindex.at(k);
      const double ai = Avalue[k];

      // todo if the variable type of column is HighsVarType::INTEGRAL check the
      // integrality of all coefficients divided by the coefficient of the
      // integral singleton variable coefficients before doing a substitution,
      // for now we skip the reductions
      if (mip && integrality[col] == HighsVarType::kInteger) {
        // for integral columns only handle equality rows
        if (rowLower[i] != rowUpper[i]) {
          ++it;
          continue;
        }

        bool suitable = false;
        for (HighsInt kk = ARstart[i]; kk < ARstart[i + 1]; ++kk) {
          HighsInt j = ARindex[kk];
          if (flagCol[j] && j != col) {
            if (integrality[col] != HighsVarType::kInteger) {
              suitable = false;
              break;
            }

            double substval = ARvalue[kk] / ai;
            double intval = std::floor(substval + 0.5);
            if (std::abs(substval - intval) > tol) {
              suitable = false;
              break;
            }
          }
        }

        if (!suitable) {
          ++it;
          continue;
        }

        printf("integer singleton can be subsituted\n");
      }

      // zero cost
      bool on_zero_cost = false && !mip;
      if (on_zero_cost && fabs(colCost.at(col)) < tol) {
        removeZeroCostColumnSingleton(col, i, k);
        it = singCol.erase(it);
        continue;
      }

      if (colLower.at(col) <= -kHighsInf && colUpper.at(col) >= kHighsInf) {
        removeFreeColumnSingleton(col, i, k);
        it = singCol.erase(it);
        continue;
      }

      // implied free
      const bool result = removeIfImpliedFree(col, i, k);
      if (result) {
        it = singCol.erase(it);
        continue;
      }

      // todo, I think this case might not work for integral variables
      // singleton column in a doubleton inequality
      // case two column singletons
      if (!mip || integrality[col] != HighsVarType::kInteger)
        if (nzRow.at(i) == 2) {
          const bool result_di =
              removeColumnSingletonInDoubletonInequality(col, i, k);
          if (result_di) {
            if (status) return;
            it = singCol.erase(it);
            continue;
          }
        }
      it++;

      if (status) return;
    } else
      it = singCol.erase(it);
  }
}

void Presolve::removeSingletonsOnly() {
  for (HighsInt row = 0; row < numRow; row++) {
    if (!flagRow[row]) continue;
    bool valid = true;
    HighsInt nz_col = 0;
    for (HighsInt k = ARstart[row]; k < ARstart[row + 1]; k++) {
      const HighsInt col = ARindex[k];
      if (!flagCol[col]) continue;
      if (nzCol[col] != 1) {
        valid = false;
        break;
      }
      nz_col++;
    }
    if (!valid) continue;
    if (nz_col == 0) {
      flagRow[row] = false;
      continue;
    }

    std::cout << "Singletons only row found! nzcol = " << nz_col
              << " L = " << rowLower[row] << " U = " << rowUpper[row]
              << std::endl;
  }
  // timer.recordStart(KNAPSACK);

  list<HighsInt>::iterator it = singCol.begin();
  while (it != singCol.end()) {
    const HighsInt col = *it;
    if (!flagCol[col]) {
      it = singCol.erase(it);
      continue;
    }

    bool remove = isKnapsack(col);
    if (remove) {
      removeKnapsack(col);
      it = singCol.erase(it);
      continue;
    }
    it++;
  }

  // timer.recordFinish(KNAPSACK);
}

void Presolve::removeKnapsack(const HighsInt col) {
  for (HighsInt k = Astart[col]; k < Aend[col]; k++) {
    assert(Aindex[k] >= 0 && Aindex[k] <= numRow);
    // todo:
  }

  return;
}

bool Presolve::isKnapsack(const HighsInt col) const {
  for (HighsInt k = Astart[col]; k < Aend[col]; k++) {
    assert(Aindex[k] >= 0 && Aindex[k] <= numRow);
    if (flagRow[Aindex[k]]) {
      if (nzCol[Aindex[k]] != 1) return false;
    }
  }
  return true;
}

pair<double, double> Presolve::getBoundsImpliedFree(double lowInit,
                                                    double uppInit,
                                                    const HighsInt col,
                                                    const HighsInt i,
                                                    const HighsInt k) {
  double low = lowInit;
  double upp = uppInit;

  // use implied bounds with original bounds
  HighsInt j;
  double l, u;
  // if at any stage low becomes  or upp becomes inf break loop
  // can't use bounds for variables generated by the same row.
  // low
  for (HighsInt kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      // check if new bounds are precisely implied bounds from same row
      if (i != implColLowerRowIndex.at(j))
        l = max(colLower.at(j), implColLower.at(j));
      else
        l = colLower.at(j);
      if (i != implColUpperRowIndex.at(j))
        u = min(colUpper.at(j), implColUpper.at(j));
      else
        u = colUpper.at(j);

      if ((Avalue.at(k) < 0 && ARvalue.at(kk) > 0) ||
          (Avalue.at(k) > 0 && ARvalue.at(kk) < 0))
        if (l <= -kHighsInf) {
          low = -kHighsInf;
          break;
        } else
          low -= ARvalue.at(kk) * l;
      else if (u >= kHighsInf) {
        low = -kHighsInf;
        break;
      } else
        low -= ARvalue.at(kk) * u;
    }
  }
  // upp
  for (HighsInt kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      // check if new bounds are precisely implied bounds from same row
      if (i != implColLowerRowIndex.at(j))
        l = max(colLower.at(j), implColLower.at(j));
      else
        l = colLower.at(j);
      if (i != implColUpperRowIndex.at(j))
        u = min(colUpper.at(j), implColUpper.at(j));
      else
        u = colUpper.at(j);
      // if at any stage low becomes  or upp becomes inf it's not implied free
      // low::
      if ((Avalue.at(k) < 0 && ARvalue.at(kk) > 0) ||
          (Avalue.at(k) > 0 && ARvalue.at(kk) < 0))
        if (u >= kHighsInf) {
          upp = kHighsInf;
          break;
        } else
          upp -= ARvalue.at(kk) * u;
      else if (l <= -kHighsInf) {
        upp = kHighsInf;
        break;
      } else
        upp -= ARvalue.at(kk) * l;
    }
  }
  return make_pair(low, upp);
}

void Presolve::removeImpliedFreeColumn(const HighsInt col, const HighsInt i,
                                       const HighsInt k) {
  if (iPrint > 0)
    cout << "PR: Implied free column singleton " << col << " removed.  Row "
         << i << " removed." << endl;

  countRemovedCols(kImpliedFreeSingCol);
  countRemovedRows(kImpliedFreeSingCol);

  // modify costs
  HighsInt j;
  vector<pair<HighsInt, double>> newCosts;
  for (HighsInt kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      newCosts.push_back(make_pair(j, colCost.at(j)));
      colCost.at(j) =
          colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);
    }
  }
  if (iKKTcheck == 1) chk2.costs.push(newCosts);

  flagCol.at(col) = 0;
  postValue.push(colCost.at(col));
  fillStackRowBounds(i);

  valueColDual.at(col) = 0;
  valueRowDual.at(i) = -colCost.at(col) / Avalue.at(k);
  double b = valueRowDual[i] < 0 || rowUpper[i] == kHighsInf ? rowLower[i]
                                                             : rowUpper[i];
  assert(std::isfinite(b));
  objShift += colCost.at(col) * b / Avalue.at(k);
  addChange(kImpliedFreeSingCol, i, col);
  removeRow(i);
}

bool Presolve::removeIfImpliedFree(HighsInt col, HighsInt i, HighsInt k) {
  // first find which bound is active for row i
  // A'y + c = z so yi = -ci/aij
  double aij = getaij(i, col);
  if (aij != Avalue.at(k)) cout << "ERROR during implied free";
  double yi = -colCost.at(col) / aij;
  double low, upp;

  if (yi > 0) {
    if (rowUpper.at(i) >= kHighsInf) return false;
    low = rowUpper.at(i);
    upp = rowUpper.at(i);
  } else if (yi < 0) {
    if (rowLower.at(i) <= -kHighsInf) return false;
    low = rowLower.at(i);
    upp = rowLower.at(i);
  } else {
    low = rowLower.at(i);
    upp = rowUpper.at(i);
  }

  pair<double, double> p = getBoundsImpliedFree(low, upp, col, i, k);
  low = p.first;
  upp = p.second;

  if (low > -kHighsInf) low = low / Avalue.at(k);
  if (upp < kHighsInf) upp = upp / Avalue.at(k);

  // if implied free
  if (colLower.at(col) <= low && low <= upp && upp <= colUpper.at(col)) {
    removeImpliedFreeColumn(col, i, k);
    return true;
  }
  // else calculate implied bounds
  else if (colLower.at(col) <= low && low <= upp) {
    if (implColLower.at(col) < low) {
      implColLower.at(col) = low;
      implColUpperRowIndex.at(col) = i;
    }
    implColDualUpper.at(col) = 0;
  } else if (low <= upp && upp <= colUpper.at(col)) {
    if (implColUpper.at(col) > upp) {
      implColUpper.at(col) = upp;
      implColUpperRowIndex.at(col) = i;
    }
    implColDualLower.at(col) = 0;
  }

  return false;
}

// used to remove column too, now possible to just modify bounds
void Presolve::removeRow(HighsInt i) {
  hasChange = true;
  flagRow.at(i) = 0;
  for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
    HighsInt j = ARindex.at(k);
    if (flagCol.at(j)) {
      nzCol.at(j)--;
      // if now singleton add to list
      if (nzCol.at(j) == 1) {
        HighsInt index = getSingColElementIndexInA(j);
        if (index >= 0)
          singCol.push_back(j);
        else
          cout << "Warning: Column " << j
               << " with 1 nz but not in singCol or? Row removing of " << i
               << ". Ignored.\n";
      }
      // if it was a singleton column remove from list and problem
      if (nzCol.at(j) == 0) removeEmptyColumn(j);
    }
  }
}

void Presolve::fillStackRowBounds(HighsInt row) {
  postValue.push(rowUpper.at(row));
  postValue.push(rowLower.at(row));
}

pair<double, double> Presolve::getImpliedRowBounds(HighsInt row) {
  double g = 0;
  double h = 0;

  HighsInt col;
  for (HighsInt k = ARstart.at(row); k < ARstart.at(row + 1); ++k) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      if (ARvalue.at(k) < 0) {
        if (colUpper.at(col) < kHighsInf)
          g += ARvalue.at(k) * colUpper.at(col);
        else {
          g = -kHighsInf;
          break;
        }
      } else {
        if (colLower.at(col) > -kHighsInf)
          g += ARvalue.at(k) * colLower.at(col);
        else {
          g = -kHighsInf;
          break;
        }
      }
    }
  }

  for (HighsInt k = ARstart.at(row); k < ARstart.at(row + 1); ++k) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      if (ARvalue.at(k) < 0) {
        if (colLower.at(col) > -kHighsInf)
          h += ARvalue.at(k) * colLower.at(col);
        else {
          h = kHighsInf;
          break;
        }
      } else {
        if (colUpper.at(col) < kHighsInf)
          h += ARvalue.at(k) * colUpper.at(col);
        else {
          h = kHighsInf;
          break;
        }
      }
    }
  }
  return make_pair(g, h);
}

void Presolve::setVariablesToBoundForForcingRow(const HighsInt row,
                                                const bool isLower) {
  HighsInt k, col;
  if (iPrint > 0)
    cout << "PR: Forcing row " << row
         << " removed. Following variables too:   nzRow=" << nzRow.at(row)
         << endl;

  flagRow.at(row) = 0;
  addChange(kForcingRow, row, 0);
  k = ARstart.at(row);
  while (k < ARstart.at(row + 1)) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      double value;
      if ((ARvalue.at(k) < 0 && isLower) || (ARvalue.at(k) > 0 && !isLower))
        value = colUpper.at(col);
      else
        value = colLower.at(col);

      setPrimalValue(col, value);
      valueColDual.at(col) = colCost.at(col);
      vector<double> bnds({colLower.at(col), colUpper.at(col)});
      oldBounds.push(make_pair(col, bnds));
      addChange(kForcingRowVariable, 0, col);

      if (iPrint > 0)
        cout << "PR:      Variable  " << col << " := " << value << endl;
      countRemovedCols(kForcingRow);
    }
    ++k;
  }

  countRemovedRows(kForcingRow);
}

void Presolve::dominatedConstraintProcedure(const HighsInt i, const double g,
                                            const double h) {
  HighsInt j;
  double val;
  if (h < kHighsInf) {
    // fill in implied bounds arrays
    if (h < implRowValueUpper.at(i)) {
      implRowValueUpper.at(i) = h;
    }
    if (h <= rowUpper.at(i)) implRowDualLower.at(i) = 0;

    // calculate implied bounds for discovering free column singletons
    for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
      j = ARindex.at(k);
      if (flagCol.at(j)) {
        if (ARvalue.at(k) < 0 && colLower.at(j) > -kHighsInf) {
          val = (rowLower.at(i) - h) / ARvalue.at(k) + colLower.at(j);
          if (val < implColUpper.at(j)) {
            implColUpper.at(j) = val;
            implColUpperRowIndex.at(j) = i;
          }
        } else if (ARvalue.at(k) > 0 && colUpper.at(j) < kHighsInf) {
          val = (rowLower.at(i) - h) / ARvalue.at(k) + colUpper.at(j);
          if (val > implColLower.at(j)) {
            implColLower.at(j) = val;
            implColLowerRowIndex.at(j) = i;
          }
        }
      }
    }
  }
  if (g > -kHighsInf) {
    // fill in implied bounds arrays
    if (g > implRowValueLower.at(i)) {
      implRowValueLower.at(i) = g;
    }
    if (g >= rowLower.at(i)) implRowDualUpper.at(i) = 0;

    // calculate implied bounds for discovering free column singletons
    for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
      HighsInt j = ARindex.at(k);
      if (flagCol.at(j)) {
        if (ARvalue.at(k) < 0 && colUpper.at(j) < kHighsInf) {
          val = (rowUpper.at(i) - g) / ARvalue.at(k) + colUpper.at(j);
          if (val > implColLower.at(j)) {
            implColLower.at(j) = val;
            implColLowerRowIndex.at(j) = i;
          }
        } else if (ARvalue.at(k) > 0 && colLower.at(j) > -kHighsInf) {
          val = (rowUpper.at(i) - g) / ARvalue.at(k) + colLower.at(j);
          if (val < implColUpper.at(j)) {
            implColUpper.at(j) = val;
            implColUpperRowIndex.at(j) = i;
          }
        }
      }
    }
  }
}

void Presolve::removeForcingConstraints() {
  double g, h;
  pair<double, double> implBounds;

  if (timer.reachLimit()) {
    status = Stat::kTimeout;
    return;
  }
  for (HighsInt i = 0; i < numRow; ++i)
    if (flagRow.at(i)) {
      if (status) return;
      if (nzRow.at(i) == 0) {
        removeEmptyRow(i);
        countRemovedRows(kEmptyRow);
        continue;
      }

      // removeRowSingletons will handle just after removeForcingConstraints
      if (nzRow.at(i) == 1) continue;

      implBounds = getImpliedRowBounds(i);

      g = implBounds.first;
      h = implBounds.second;

      // Infeasible row
      if (g > rowUpper.at(i) || h < rowLower.at(i)) {
        if (iPrint > 0) cout << "PR: Problem infeasible." << endl;
        status = kInfeasible;
        return;
      }
      // Forcing row
      else if (g == rowUpper.at(i)) {
        setVariablesToBoundForForcingRow(i, true);
      } else if (h == rowLower.at(i)) {
        setVariablesToBoundForForcingRow(i, false);
      }
      // Redundant row
      else if (g >= rowLower.at(i) && h <= rowUpper.at(i)) {
        removeRow(i);
        addChange(kRedundantRow, i, 0);
        if (iPrint > 0)
          cout << "PR: Redundant row " << i << " removed." << endl;
        countRemovedRows(kRedundantRow);
      }
      // Dominated constraints
      else {
        dominatedConstraintProcedure(i, g, h);
        continue;
      }
    }
}

void Presolve::removeRowSingletons() {
  if (timer.reachLimit()) {
    status = Stat::kTimeout;
    return;
  }
  timer.recordStart(kSingRow);

  list<HighsInt>::iterator it = singRow.begin();
  while (it != singRow.end()) {
    if (flagRow[*it]) {
      const HighsInt i = *it;
      assert(i >= 0 && i < numRow);

      const HighsInt k = getSingRowElementIndexInAR(i);
      if (k < 0) {
        it = singRow.erase(it);
        // kxx
        continue;
      }

      const HighsInt j = ARindex.at(k);

      // add old bounds OF X to checker and for postsolve
      if (iKKTcheck == 1) {
        vector<pair<HighsInt, double>> bndsL, bndsU, costS;
        bndsL.push_back(make_pair(j, colLower.at(j)));
        bndsU.push_back(make_pair(j, colUpper.at(j)));
        costS.push_back(make_pair(j, colCost.at(j)));

        chk2.cLowers.push(bndsL);
        chk2.cUppers.push(bndsU);
        chk2.costs.push(costS);
      }

      vector<double> bnds(
          {colLower.at(j), colUpper.at(j), rowLower.at(i), rowUpper.at(i)});
      oldBounds.push(make_pair(j, bnds));

      double aij = ARvalue.at(k);
      /*		//before update bounds of x take it out of rows with
      implied row bounds for (HighsInt r = Astart.at(j); r<Aend.at(j); r++) { if
      (flagRow[Aindex[r]]) { HighsInt rr = Aindex[r]; if (implRowValueLower[rr]
      > -kHighsInf) { if (aij > 0) implRowValueLower[rr] =
      implRowValueLower[rr] - aij*colLower.at(j); else implRowValueLower[rr] =
      implRowValueLower[rr] - aij*colUpper.at(j);
                      }
                      if (implRowValueUpper[rr] < kHighsInf) {
                              if (aij > 0)
                                      implRowValueUpper[rr] =
      implRowValueUpper[rr] - aij*colUpper.at(j); else implRowValueUpper[rr] =
      implRowValueUpper[rr] - aij*colLower.at(j);
                      }
              }
      }*/

      // todo, round bounds for integer and implied integer variables, be
      // careful with tolerances, i.e. use floor(ub + eps) for rounding down an
      // upper bound update bounds of X
      if (aij > 0) {
        if (rowLower.at(i) != -kHighsInf)
          colLower.at(j) =
              max(max(rowLower.at(i) / aij, -kHighsInf), colLower.at(j));
        if (rowUpper.at(i) != kHighsInf)
          colUpper.at(j) =
              min(min(rowUpper.at(i) / aij, kHighsInf), colUpper.at(j));
      } else if (aij < 0) {
        if (rowLower.at(i) != -kHighsInf)
          colUpper.at(j) =
              min(min(rowLower.at(i) / aij, kHighsInf), colUpper.at(j));
        if (rowUpper.at(i) != kHighsInf)
          colLower.at(j) =
              max(max(rowUpper.at(i) / aij, -kHighsInf), colLower.at(j));
      }

      /*		//after update bounds of x add to rows with implied row
      bounds for (HighsInt r = Astart.at(j); r<Aend.at(j); r++) { if
      (flagRow[r]) { HighsInt rr = Aindex[r]; if (implRowValueLower[rr] >
      -kHighsInf) { if (aij > 0) implRowValueLower[rr] =
      implRowValueLower[rr] + aij*colLower.at(j); else implRowValueLower[rr] =
      implRowValueLower[rr] + aij*colUpper.at(j);
                      }
                      if (implRowValueUpper[rr] < kHighsInf) {
                              if (aij > 0)
                                      implRowValueUpper[rr] =
      implRowValueUpper[rr] + aij*colUpper.at(j); else implRowValueUpper[rr] =
      implRowValueUpper[rr] + aij*colLower.at(j);
                      }
              }
      }*/

      roundIntegerBounds(j);

      // check for feasibility
      // Analyse dependency on numerical tolerance
      timer.updateNumericsRecord(kNumericsInconsistentBounds,
                                 colLower.at(j) - colUpper.at(j));
      if (colLower.at(j) - colUpper.at(j) > inconsistent_bounds_tolerance) {
        status = kInfeasible;
        timer.recordFinish(kSingRow);
        return;
      }

      if (iPrint > 0)
        cout << "PR: Singleton row " << i << " removed. Bounds of variable  "
             << j << " modified: l= " << colLower.at(j)
             << " u=" << colUpper.at(j) << ", aij = " << aij << endl;

      addChange(kSingRow, i, j);
      postValue.push(colCost.at(j));
      removeRow(i);

      if (flagCol.at(j)) {
        // Analyse dependency on numerical tolerance
        timer.updateNumericsRecord(kNumericsFixedColumn,
                                   fabs(colUpper.at(j) - colLower.at(j)));
        if (fabs(colUpper.at(j) - colLower.at(j)) <= fixed_column_tolerance)
          removeFixedCol(j);
      }
      countRemovedRows(kSingRow);

      if (status) {
        timer.recordFinish(kSingRow);
        return;
      }
      it = singRow.erase(it);
    } else {
      it++;
    }
  }
  timer.recordFinish(kSingRow);
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

void Presolve::checkForChanges(HighsInt iteration) {
  if (iteration <= 2) {
    // flagCol has one more element at end which is zero
    // from removeDoubletonEquatoins, needed for AR matrix manipulation
    if (none_of(flagCol.begin(), flagCol.begin() + numCol,
                [](HighsInt i) { return i == 0; }) &&
        none_of(flagRow.begin(), flagRow.begin() + numRow,
                [](HighsInt i) { return i == 0; })) {
      if (iPrint > 0)
        cout << "PR: No variables were eliminated at presolve." << endl;
      noPostSolve = true;
      return;
    }
  }
  resizeProblem();
  status = Stat::kReduced;
}

// void Presolve::reportTimes() {
//   HighsInt reportList[] = {kEmptyRow,
//                       kFixedCol,
//                       kSingRow,
//                       kDoubletonEquation,
//                       kForcingRow,
//                       kRedundantRow,
//                       kFreeSingCol,
//                       kSingColDoubletonIneq,
//                       kImpliedFreeSingCol,
//                       kDominatedCols,
//                       kWeaklyDominatedCols};
//   HighsInt reportCount = sizeof(reportList) / sizeof(HighsInt);

//   printf("Presolve rules ");
//   for (HighsInt i = 0; i < reportCount; ++i) {
//     printf(" %s", timer.itemNames[reportList[i]].c_str());
//     cout << flush;
//   }

//   printf("\n");
//   cout << "Time spent     " << flush;
//   for (HighsInt i = 0; i < reportCount; ++i) {
//     float f = (float)timer.itemTicks[reportList[i]];
//     if (f < 0.01)
//       cout << setw(4) << " <.01 ";
//     else
//       printf(" %3.2f ", f);
//   }
//   printf("\n");
// }

// void Presolve::recordCounts(const string fileName) {
//   ofstream myfile;
//   myfile.open(fileName.c_str(), ios::app);
//   HighsInt reportList[] = {kEmptyRow,
//                       kFixedCol,
//                       kSingRow,
//                       kDoubletonEquation,
//                       kForcingRow,
//                       kRedundantRow,
//                       kFreeSingCol,
//                       kSingColDoubletonIneq,
//                       kImpliedFreeSingCol,
//                       kDominatedCols,
//                       kWeaklyDominatedCols,
//                       kEmptyCol};
//   HighsInt reportCount = sizeof(reportList) / sizeof(HighsInt);

//   myfile << "Problem " << modelName << ":\n";
//   myfile << "Rule   , removed rows , removed cols , time  \n";

//   HighsInt cRows = 0, cCols = 0;
//   for (HighsInt i = 0; i < reportCount; ++i) {
//     float f = (float)timer.itemTicks[reportList[i]];

//     myfile << setw(7) << timer.itemNames[reportList[i]].c_str() << ", "
//            << setw(7) << countRemovedRows[reportList[i]] << ", " << setw(7)
//            << countRemovedCols[reportList[i]] << ", ";
//     if (f < 0.001)
//       myfile << setw(7) << " <.001 ";
//     else
//       myfile << setw(7) << setprecision(3) << f;
//     myfile << endl;

//     cRows += countRemovedRows[reportList[i]];
//     cCols += countRemovedCols[reportList[i]];
//   }

//   if (!noPostSolve) {
//     if (cRows != numRowOriginal - numRow) cout << "Wrong row reduction
//     count\n"; if (cCols != numColOriginal - numCol) cout << "Wrong col
//     reduction count\n";

//     myfile << setw(7) << "Total "
//            << ", " << setw(7) << numRowOriginal - numRow << ", " << setw(7)
//            << numColOriginal - numCol;
//   } else {
//     myfile << setw(7) << "Total "
//            << ", " << setw(7) << 0 << ", " << setw(7) << 0;
//   }
//   myfile << endl << " \\\\ " << endl;
//   myfile.close();
// }

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

HighsInt Presolve::getSingRowElementIndexInAR(HighsInt i) {
  assert(i >= 0 && i < numRow);
  HighsInt k = ARstart.at(i);
  while (k < ARstart[i + 1] && !flagCol.at(ARindex.at(k))) k++;
  if (k >= ARstart.at(i + 1)) {
    return -1;
  }
  HighsInt rest = k + 1;
  while (rest < ARstart.at(i + 1) && !flagCol.at(ARindex.at(rest))) ++rest;
  if (rest < ARstart.at(i + 1)) {
    return -1;
  }
  return k;
}

HighsInt Presolve::getSingColElementIndexInA(HighsInt j) {
  HighsInt k = Astart.at(j);
  assert(k >= 0 && k < (HighsInt)Aindex.size());
  assert(Aindex[k] >= 0 && Aindex[k] < numRow);
  assert(flagRow.size() == (unsigned int)numRow);

  while (!flagRow.at(Aindex.at(k))) ++k;
  if (k >= Aend.at(j)) {
    assert(nzCol[j] == 0);
    return -2;
  }
  HighsInt rest = k + 1;
  while (rest < Aend.at(j) && !flagRow.at(Aindex.at(rest))) ++rest;
  if (rest < Aend.at(j)) {
    // Occurs if a singleton column is no longer singleton.
    return -1;
  }
  return k;
}

void Presolve::testAnAR(HighsInt post) {
  HighsInt rows = numRow;
  HighsInt cols = numCol;

  double valueA = 0;
  double valueAR = 0;
  bool hasValueA, hasValueAR;

  if (post) {
    rows = numRowOriginal;
    cols = numColOriginal;
  }

  // check that A = AR
  for (HighsInt i = 0; i < rows; ++i) {
    for (HighsInt j = 0; j < cols; ++j) {
      if (post == 0)
        if (!flagRow.at(i) || !flagCol.at(j)) continue;
      hasValueA = false;
      for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k)
        if (Aindex.at(k) == i) {
          hasValueA = true;
          valueA = Avalue.at(k);
        }

      hasValueAR = false;
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
        if (ARindex.at(k) == j) {
          hasValueAR = true;
          valueAR = ARvalue.at(k);
        }

      if (hasValueA != hasValueAR)
        cout << "    MATRIX is0   DIFF row=" << i << " col=" << j
             << "           ------------A: " << hasValueA
             << "  AR: " << hasValueAR << endl;
      else if (hasValueA && valueA != valueAR)
        cout << "    MATRIX VAL  DIFF row=" << i << " col=" << j
             << "           ------------A: " << valueA << "  AR: " << valueAR
             << endl;
    }
  }

  if (post == 0) {
    // check nz
    HighsInt nz = 0;
    for (HighsInt i = 0; i < rows; ++i) {
      if (!flagRow.at(i)) continue;
      nz = 0;
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
        if (flagCol.at(ARindex.at(k))) nz++;
      if (nz != nzRow.at(i))
        cout << "    NZ ROW      DIFF row=" << i << " nzRow=" << nzRow.at(i)
             << " actually " << nz << "------------" << endl;
    }

    for (HighsInt j = 0; j < cols; ++j) {
      if (!flagCol.at(j)) continue;
      nz = 0;
      for (HighsInt k = Astart.at(j); k < Aend.at(j); ++k)
        if (flagRow.at(Aindex.at(k))) nz++;
      if (nz != nzCol.at(j))
        cout << "    NZ COL      DIFF col=" << j << " nzCol=" << nzCol.at(j)
             << " actually " << nz << "------------" << endl;
    }
  }
}

HighsPostsolveStatus Presolve::primalPostsolve(
    const std::vector<double>& reduced_solution,
    HighsSolution& recovered_solution) {
  colValue = reduced_solution;

  makeACopy();  // so we can efficiently calculate primal and dual
  // values

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

  for (HighsInt i = 0; i < numCol; ++i) {
    HighsInt iCol = eqIndexOfReduced.at(i);
    assert(iCol < (HighsInt)valuePrimal.size());
    assert(iCol >= 0);
    valuePrimal[iCol] = colValue.at(i);
    // valueColDual[iCol] = colDual.at(i);
    // col_status.at(iCol) = temp_col_status.at(i);
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

    switch (c.type) {
      case kAggregator: {
        // restore solution, basis, flags, and colCostAtEl
        aggregatorStack.back().postsolveStack.undo(flagCol, flagRow,
                                                   valuePrimal);

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
        break;
      }
      case kDoubletonEquation: {  // Doubleton equation row
        getDualsDoubletonEquation(c.row, c.col);

        // exit(2);
        break;
      }
      case kDoubletonEquationRowBoundsUpdate: {
        // new bounds from doubleton equation, retrieve old ones
        // just for KKT check, not called otherwise
        // chk2.addChange(171, c.row, c.col, 0, 0, 0);
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

        postValue.pop();

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
        flagRow[c.row] = 1;
        break;
      }
      case kSingRow: {
        // valuePrimal is already set for this one, colDual also, we need
        // rowDual. AR copy keeps full matrix.  col dual maybe infeasible, we
        // need to check.  recover old bounds and see
        // getDualsSingletonRow(c.row, c.col);
        oldBounds.pop();
        postValue.pop();

        break;
      }
      case kForcingRowVariable:
        oldBounds.pop();
        flagCol[c.col] = 1;
        break;
      case kForcingRow: {
        flagRow[c.row] = 1;
        break;
      }
      case kRedundantRow: {
        flagRow[c.row] = 1;
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
        if (rowlb == rowub)
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
        break;
      }
      case kSingColDoubletonIneq: {
        assert(false);
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

        // exit(2);
        break;
      }
      case kEmptyCol:
      case kDominatedCols:
      case kWeaklyDominatedCols: {
        // got valuePrimal, need colDual
        flagCol[c.col] = 1;
        break;
      }

      case kFixedCol: {
        // got valuePrimal, need colDual
        flagCol[c.col] = 1;
        break;
      }
    }
    // cmpNBF(c.row, c.col);
  }

  // cmpNBF();

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

  rowValue.assign(numRow, 0);
  for (HighsInt i = 0; i < numRowOriginal; ++i) {
    for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
      rowValue.at(i) += valuePrimal.at(ARindex.at(k)) * ARvalue.at(k);
  }

  // cout<<"Singularity check at end of postsolve: ";
  // testBasisMatrixSingularity();

  // Save solution to PresolveComponentData.
  recovered_solution.col_value = colValue;
  recovered_solution.row_value = rowValue;

  return HighsPostsolveStatus::kSolutionRecovered;
}
// todo: error reporting.
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

/* testing and dev
HighsInt Presolve::testBasisMatrixSingularity() {

        HFactor factor;

        //resize matrix in M so we can pass to factor
        HighsInt i, j, k;
        HighsInt nz = 0;
        HighsInt nR = 0;
        HighsInt nC = 0;

        numRowOriginal = rowLowerOriginal.size();
        numColOriginal = colLowerOriginal.size();
        //arrays to keep track of indices
        vector<HighsInt> rIndex_(numRowOriginal, -1);
        vector<HighsInt> cIndex_(numColOriginal, -1);

        for (i=0;i<numRowOriginal;++i)
                if (flagRow.at(i)) {
                        for (j = ARstart.at(i); j<ARstart.at(i+1); ++j)
                                if (flagCol[ARindex.at(j)])
                                        nz ++;
                        rIndex_.at(i) = nR;
                        nR++;
                        }

        for (i=0;i<numColOriginal;++i)
                if (flagCol.at(i)) {
                        cIndex_.at(i) = nC;
                        nC++;
                }


        //matrix
        vector<HighsInt>    Mstart(nC + 1, 0);
        vector<HighsInt>    Mindex(nz);
        vector<double> Mvalue(nz);

    vector<HighsInt> iwork(nC, 0);

    for (i = 0;i<numRowOriginal; ++i)
        if (flagRow.at(i))
            for (HighsInt k = ARstart.at(i); k < ARstart.at(i+1);++k ) {
                j = ARindex.at(k);
                if (flagCol.at(j))
                                iwork[cIndex_.at(j)]++;
                        }
    for (i = 1; i <= nC; ++i)
        Mstart.at(i) = Mstart[i - 1] + iwork[i - 1];
   for (i = 0; i < numColOriginal; ++i)
        iwork.at(i) = Mstart.at(i);

   for (i = 0; i < numRowOriginal; ++i) {
        if (flagRow.at(i)) {
                        HighsInt iRow = rIndex_.at(i);
                    for (k = ARstart.at(i); k < ARstart[i + 1];++k ) {
                        j = ARindex.at(k);
                        if (flagCol.at(j)) {
                                HighsInt iCol = cIndex_.at(j);
                                    HighsInt iPut = iwork[iCol]++;
                                    Mindex[iPut] = iRow;
                                    Mvalue[iPut] = ARvalue.at(k);
                                }
                    }
                }
    }

    vector<HighsInt>  bindex(nR);
    HighsInt countBasic=0;

    printf("To recover this test need to use col/row_status\n");
     for (HighsInt i=0; i< nonbasicFlag.size();++i) {
         if (nonbasicFlag.at(i) == 0)
                         countBasic++;
     }

     if (countBasic != nR)
         cout<<" Wrong count of basic variables: != numRow"<<endl;

     HighsInt c=0;
     for (HighsInt i=0; i< nonbasicFlag.size();++i) {
         if (nonbasicFlag.at(i) == 0) {
                        if (i < numColOriginal)
                                bindex[c] = cIndex_.at(i);
                        else
                                bindex[c] = nC + rIndex_[i - numColOriginal];
                        c++;
         }
    }

        factor.setup(nC, nR, &Mstart[0], &Mindex[0], &Mvalue[0],  &bindex[0]);
/ *	if (1) // for this check both A and M are the full matrix again
        {
                if (nC - numColOriginal != 0)
                        cout<<"columns\n";
                if (nR - numRowOriginal != 0)
                        cout<<"rows\n";
                for (HighsInt i=0; i< Mstart.size();++i)
                        if (Mstart.at(i) - Astart.at(i) != 0)
                                cout<<"Mstart "<<i<<"\n";
                for (HighsInt i=0; i< Mindex.size();++i)
                        if (Mindex.at(i) - Aindex.at(i) != 0)
                                cout<<"Mindex "<<i<<"\n";
                for (HighsInt i=0; i< Mvalue.size();++i)
                        if (Mvalue.at(i) - Avalue.at(i) != 0)
                                cout<<"Mvalue "<<i<<"\n";
                for (HighsInt i=0; i< bindex.size();++i)
                        if (nonbasicFlag.at(i) - nbffull.at(i) != 0)
                                cout<<"nbf "<<i<<"\n";
        } * /

        try {
        factor.build();
    } catch (runtime_error& error) {
        cout << error.what() << endl;
        cout << "Postsolve: could not factorize basis matrix." << endl;
        return 0;
    }
    cout << "Postsolve: basis matrix successfully factorized." << endl;

    return 1;
}*/

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
