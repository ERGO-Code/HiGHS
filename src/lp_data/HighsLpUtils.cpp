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
  HighsLp lp = highs_model.simplex_lp_;
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

int validate_col_bounds(int XnumCol, const double* XcolLower, const double* XcolUpper) {
  assert(XnumCol >= 0);
  if (XnumCol == 0) return 0;
  for (int col = 0; col < XnumCol; col++) {
    if (highs_isInfinity(XcolLower[col])) return col + 1;
    if (XcolLower[col] > XcolUpper[col]) return XnumCol + col + 1;
    if (highs_isInfinity(-XcolUpper[col])) return -(col + 1);
  }
}

void filter_col_bounds(HighsLp& lp, int XfromCol, int XtoCol, const double infinite_bound) {
  assert(XfromCol >= 0);
  assert(XtoCol < lp.numCol_);
  assert(XfromCol <= XtoCol);
  int numChangedLowerBounds = 0;
  int numChangedUpperBounds = 0;
  for (int col = XfromCol; col <= XtoCol; col++) {
    if (lp.colLower_[col] <= -infinite_bound) {
      lp.colLower_[col] = -HIGHS_CONST_INF;
      numChangedLowerBounds++;
    }
    if (lp.colUpper_[col] >= infinite_bound) {
      lp.colUpper_[col] = HIGHS_CONST_INF;
      numChangedUpperBounds++;
    }
  }
  if (numChangedLowerBounds)
    HighsLogMessage(HighsMessageType::WARNING, "%12d col lower bounds below %12g interpreted as -Infinity",
		    numChangedLowerBounds, -infinite_bound);
  if (numChangedUpperBounds)
    HighsLogMessage(HighsMessageType::WARNING, "%12d col upper bounds above %12g interpreted as +Infinity",
		    numChangedUpperBounds, infinite_bound);
}

int validate_row_bounds(int XnumRow, const double* XrowLower, const double* XrowUpper) {
  assert(XnumRow >= 0);
  if (XnumRow == 0) return 0;
  for (int row = 0; row < XnumRow; row++) {
    if (highs_isInfinity(XrowLower[row])) return row + 1;
    if (XrowLower[row] > XrowUpper[row]) return XnumRow + row + 1;
    if (highs_isInfinity(-XrowUpper[row])) return -(row + 1);
  }
}

void filter_row_bounds(HighsLp& lp, int XfromRow, int XtoRow, double infinite_bound) {
  assert(XfromRow >= 0);
  assert(XtoRow < lp.numRow_);
  assert(XfromRow <= XtoRow);
  int numChangedLowerBounds = 0;
  int numChangedUpperBounds = 0;
  for (int row = XfromRow; row <= XtoRow; row++) {
    if (lp.rowLower_[row] <= -infinite_bound) {
      lp.rowLower_[row] = -HIGHS_CONST_INF;
      numChangedLowerBounds++;
    }
    if (lp.rowUpper_[row] >= infinite_bound) {
      lp.rowUpper_[row] = HIGHS_CONST_INF;
      numChangedUpperBounds++;
    }
  }
  if (numChangedLowerBounds)
    HighsLogMessage(HighsMessageType::WARNING, "%12d row lower bounds below %12g interpreted as -Infinity",
		    numChangedLowerBounds, -infinite_bound);
  if (numChangedUpperBounds)
    HighsLogMessage(HighsMessageType::WARNING, "%12d row upper bounds above %12g interpreted as +Infinity",
		    numChangedUpperBounds, infinite_bound);
}

int validate_matrix_indices(int XnumRow, int XnumCol, int XnumNZ, const int* XAstart, const int* XAindex) {
  assert(XnumRow >= 0);
  assert(XnumCol >= 0);
  assert(XnumNZ >= 0);
  if (XnumRow == 0) return 0;
  if (XnumCol == 0) return 0;
  if (XnumNZ == 0) return 0;
  vector <int> colVec;
  colVec.assign(XnumRow, 0);
  for (int col = 0; col < XnumCol; col++) {
    int fromEl = XAstart[col];
    int toEl;
    if (col < XnumCol-1) {
      toEl = XAstart[col+1]-1;
    } else {
      toEl = XnumNZ-1;
    }
    if (fromEl > toEl+1 || fromEl > XnumNZ || fromEl < 0) return -(XnumCol + col + 1);
    for (int el = fromEl; el <= toEl; el++) {
      int row = XAindex[el];
      if (row < 0) return col + 1;
      if (row >= XnumRow) return -(col + 1);
      if (colVec[row]) return XnumCol + col + 1;
    }
    for (int el = XAstart[col]; el <= toEl; el++) colVec[XAindex[el]] = 0;
  }
  int col = XnumCol;
  int fromEl = XAstart[col];
  if (fromEl > XnumNZ || fromEl < 0) return -(XnumCol + col + 1);
}

void filter_matrix_values(HighsLp& lp, int XfromCol, int XtoCol, double small_matrix_value) {
  assert(XfromCol >= 0);
  assert(XtoCol < lp.numCol_);
  assert(XfromCol <= XtoCol);
  assert(small_matrix_value >=0);
  int numRemovedValues = 0;
  int numTrueNZ = lp.Astart_[XfromCol];
  for (int col = XfromCol; col <= XtoCol; col++) {
    int fromEl = lp.Astart_[col];
    lp.Astart_[col] = numTrueNZ;
    for (int el = fromEl; el < lp.Astart_[col+1]; el++) {
      if (fabs(lp.Avalue_[el]) <= small_matrix_value) {
	numRemovedValues++;
      } else {
	lp.Aindex_[numTrueNZ] = lp.Aindex_[el];
	lp.Avalue_[numTrueNZ] = lp.Avalue_[el];
	numTrueNZ++;
      }
    }
  }
  if (numRemovedValues) {
    for (int col = XtoCol+1; col < lp.numCol_; col++) {
      int fromEl = lp.Astart_[col];
      lp.Astart_[col] = numTrueNZ;
	for (int el = fromEl; el < lp.Astart_[col+1]; el++) {
	  lp.Aindex_[numTrueNZ] = lp.Aindex_[el];
	  lp.Avalue_[numTrueNZ] = lp.Avalue_[el];
	  numTrueNZ++;
	}
    }
    lp.Astart_[lp.numCol_] = numTrueNZ;
    lp.nnz_ = numTrueNZ;
    HighsLogMessage(HighsMessageType::WARNING, "%12d matrix values less than %12g removed",
		    numRemovedValues, small_matrix_value);
  }
}

void filter_row_matrix_values(int XnumRow, int XnumNZ, int* XARstart, int* XARindex, double* XARvalue, double small_matrix_value) {
  assert(XnumRow >= 0);
  assert(XnumNZ >= 0);
  assert(small_matrix_value >=0);
  if (XnumRow == 0) return;
  if (XnumNZ == 0) return;
  int numRemovedValues = 0;
  int numTrueNZ = 0;
  for (int row = 0; row < XnumRow; row++) {
    int fromEl = XARstart[row];
    XARstart[row] = numTrueNZ;
    for (int el = fromEl; el < XARstart[row+1]; el++) {
      if (fabs(XARvalue[el]) <= small_matrix_value) {
	numRemovedValues++;
      } else {
	XARindex[numTrueNZ] = XARindex[el];
	XARvalue[numTrueNZ] = XARvalue[el];
	numTrueNZ++;
      }
    }
  }
  XARstart[XnumRow] = numTrueNZ;
  if (numRemovedValues) 
    HighsLogMessage(HighsMessageType::WARNING, "%12d matrix values less than %12g removed",
		    numRemovedValues, small_matrix_value);
}


#ifdef HiGHSDEV
void util_analyseLp(const HighsLp &lp, const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_analyseVectorValues("Column costs", lp.numCol_, lp.colCost_, false);
  util_analyseVectorValues("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_analyseVectorValues("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_analyseVectorValues("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_analyseVectorValues("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_analyseVectorValues("Matrix sparsity", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  util_analyseMatrixSparsity("Constraint matrix", lp.numCol_, lp.numRow_, lp.Astart_, lp.Aindex_);
  util_analyseModelBounds("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_analyseModelBounds("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
}
#endif
