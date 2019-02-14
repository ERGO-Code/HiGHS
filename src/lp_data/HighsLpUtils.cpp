/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "util/HighsUtils.h"
#include "lp_data/HighsStatus.h"

void getLpCosts(
		const HighsLp& lp,
		int firstcol,
		int lastcol,
		double* XcolCost
		) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) XcolCost[col - firstcol] = lp.colCost_[col];
}

void getLpColBounds(const HighsLp& lp,
		    int firstcol,
		    int lastcol,
		    double* XcolLower,
		    double* XcolUpper) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) {
    if (XcolLower != NULL) XcolLower[col - firstcol] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - firstcol] = lp.colUpper_[col];
  }
}

void getLpRowBounds(const HighsLp& lp,
		    int firstrow,
		    int lastrow,
		    double* XrowLower,
		    double* XrowUpper) {
  assert(0 <= firstrow);
  assert(firstrow <= lastrow);
  assert(lastrow < lp.numRow_);
  for (int row = firstrow; row <= lastrow; ++row) {
    if (XrowLower != NULL) XrowLower[row - firstrow] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - firstrow] = lp.rowUpper_[row];
  }
}

// Get a single coefficient from the matrix
void getLpMatrixCoefficient(const HighsLp& lp, int row, int col, double *val) {
  assert(row >= 0 && row < lp.numRow_);
  assert(col >= 0 && col < lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);
#endif
  //  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);

  int get_el = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    //  printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //  lp.Aindex_[el], row);
    if (lp.Aindex_[el] == row) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    //  printf("model.util_getCoeff: Cannot find row %d in column %d\n", row, col);
    *val = 0;
  } else {
    //  printf("model.util_getCoeff: Found row %d in column %d as element %d:
    //  value %g\n", row, col, get_el, lp.Avalue_[get_el]);
    *val = lp.Avalue_[get_el];
  }
}

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void reportLp(const HighsLp &lp) {
  reportLpBrief(lp);
  reportLpColVec(lp);
  reportLpRowVec(lp);
  reportLpColMtx(lp);
}

// Report the LP briefly
void reportLpBrief(const HighsLp &lp) {
  reportLpDimensions(lp);
  reportLpObjSense(lp);
}

// Report the LP dimensions
void reportLpDimensions(const HighsLp &lp) {
  HighsPrintMessage(ML_MINIMAL,
                    "LP has %d columns, %d rows and %d nonzeros\n",
                    lp.numCol_, lp.numRow_, lp.Astart_[lp.numCol_]);
}

// Report the LP objective sense
void reportLpObjSense(const HighsLp &lp) {
  if (lp.sense_ == OBJSENSE_MINIMIZE)
    HighsPrintMessage(ML_MINIMAL, "Objective sense is minimize\n");
  else if (lp.sense_ == OBJSENSE_MAXIMIZE)
    HighsPrintMessage(ML_MINIMAL, "Objective sense is maximize\n");
  else
    HighsPrintMessage(ML_MINIMAL,
                      "Objective sense is ill-defined as %d\n", lp.sense_);
}

// Report the vectors of LP column data
void reportLpColVec(const HighsLp &lp) {
  if (lp.numCol_ <= 0) return;
  HighsPrintMessage(ML_VERBOSE,
                    "  Column        Lower        Upper         Cost\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(ML_VERBOSE, "%8d %12g %12g %12g\n", iCol,
                      lp.colLower_[iCol], lp.colUpper_[iCol], lp.colCost_[iCol]);
  }
}

// Report the vectors of LP row data
void reportLpRowVec(const HighsLp &lp) {
  if (lp.numRow_ <= 0) return;
  HighsPrintMessage(ML_VERBOSE,
                    "     Row        Lower        Upper\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsPrintMessage(ML_VERBOSE, "%8d %12g %12g\n", iRow,
                      lp.rowLower_[iRow], lp.rowUpper_[iRow]);
  }
}

// Report the LP column-wise matrix
void reportLpColMtx(const HighsLp &lp) {
  if (lp.numCol_ <= 0) return;
  HighsPrintMessage(ML_VERBOSE,
                    "Column Index              Value\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(ML_VERBOSE, "    %8d Start   %10d\n", iCol,
                      lp.Astart_[iCol]);
    for (int el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
      HighsPrintMessage(ML_VERBOSE, "          %8d %12g\n",
                        lp.Aindex_[el], lp.Avalue_[el]);
    }
  }
  HighsPrintMessage(ML_VERBOSE, "             Start   %10d\n",
                    lp.Astart_[lp.numCol_]);
}

/*
void reportLpSolution(HighsModelObject &highs_model) {
  HighsLp lp = highs_model.solver_lp_;
  reportLpBrief(lp);
  //  simplex_interface.report_simplex_solution_status();
  assert(lp.numCol_ > 0);
  assert(lp.numRow_ > 0);
  vector<double> colPrimal(lp.numCol_);
  vector<double> colDual(lp.numCol_);
  vector<int> colStatus(lp.numCol_);
  vector<double> rowPrimal(lp.numRow_);
  vector<double> rowDual(lp.numRow_);
  vector<int> rowStatus(lp.numRow_);
  //  util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  //  if (util_convertWorkingToBaseStat(&colStatus[0], &rowStatus[0])) return;
  //  util_reportColVecSol(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_, colPrimal, colDual, colStatus);
  //  util_reportRowVecSol(lp.numRow_, lp.rowLower_, lp.rowUpper_, rowPrimal, rowDual, rowStatus);
}
*/



HighsStatus checkLp(const HighsLp& lp) {
  // Check dimensions.
  if (lp.numCol_ <= 0 || lp.numRow_ <= 0)
    return HighsStatus::LpError;

  // Check vectors.
  if ((int)lp.colCost_.size() != lp.numCol_)
    return HighsStatus::LpError;

  if ((int)lp.colLower_.size() != lp.numCol_ ||
      (int)lp.colUpper_.size() != lp.numCol_)
    return HighsStatus::LpError;
  if ((int)lp.rowLower_.size() != lp.numRow_ ||
      (int)lp.rowUpper_.size() != lp.numRow_)
    return HighsStatus::LpError;

  for (int i = 0; i < lp.numRow_; i++)
    if (lp.rowLower_[i] < -HIGHS_CONST_INF || lp.rowUpper_[i] > HIGHS_CONST_INF)
      return HighsStatus::LpError;

  for (int j = 0; j < lp.numCol_; j++) {
    if (lp.colCost_[j] < -HIGHS_CONST_INF || lp.colCost_[j] > HIGHS_CONST_INF)
      return HighsStatus::LpError;

    if (lp.colLower_[j] < -HIGHS_CONST_INF || lp.colUpper_[j] > HIGHS_CONST_INF)
      return HighsStatus::LpError;
    if (lp.colLower_[j] > lp.colUpper_[j] + kBoundTolerance)
      return HighsStatus::LpError;
  }

  // Check matrix.
  if ((size_t)lp.nnz_ != lp.Avalue_.size())
    return HighsStatus::LpError;
  if (lp.nnz_ <= 0) return HighsStatus::LpError;
  if ((int)lp.Aindex_.size() != lp.nnz_)
    return HighsStatus::LpError;

  if ((int)lp.Astart_.size() != lp.numCol_ + 1)
    return HighsStatus::LpError;
  // Was lp.Astart_[i] >= lp.nnz_ (below), but this is wrong when the
  // last column is empty. Need to check as follows, and also check
  // the entry lp.Astart_[lp.numCol_] > lp.nnz_
  for (int i = 0; i < lp.numCol_; i++) {
    if (lp.Astart_[i] > lp.Astart_[i + 1] || lp.Astart_[i] > lp.nnz_ ||
        lp.Astart_[i] < 0)
      return HighsStatus::LpError;
  }
  if (lp.Astart_[lp.numCol_] > lp.nnz_ ||
      lp.Astart_[lp.numCol_] < 0)
      return HighsStatus::LpError;

  for (int k = 0; k < lp.nnz_; k++) {
    if (lp.Aindex_[k] < 0 || lp.Aindex_[k] >= lp.numRow_)
      return HighsStatus::LpError;
    if (lp.Avalue_[k] < -HIGHS_CONST_INF || lp.Avalue_[k] > HIGHS_CONST_INF)
      return HighsStatus::LpError;
  }

  return HighsStatus::OK;
}


#ifdef HiGHSDEV
void util_analyseLp(const HighsLp &lp, const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_analyseVectorValues("Column costs", lp.numCol_, lp.colCost_, false);
  util_analyseVectorValues("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_analyseVectorValues("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_analyseVectorValues("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_analyseVectorValues("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_analyseVectorValues("Matrix entries", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  util_analyseModelBounds("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_analyseModelBounds("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
}
#endif
