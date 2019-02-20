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

int assess_lp(const HighsLp& lp, const HighsOptions& options) {
  // Assess the bounds and matrix indices, returning on error unless addition is forced
  int return_code;
  bool error_found = false;
  return_code = assess_col_bounds(lp.numCol_, &lp.colLower_[0], &lp.colUpper_[0], options.infinite_bound);
  error_found = error_found || return_code;
  return_code = assess_row_bounds(lp.numRow_, &lp.rowLower_[0], &lp.rowUpper_[0], options.infinite_bound);
  error_found = error_found || return_code;
  return_code = assess_matrix_indices(lp.numRow_, lp.numCol_, lp.nnz_, &lp.Astart_[0], &lp.Aindex_[0]);
  error_found = error_found || return_code;
}

int add_lp_cols(HighsLp& lp,
		int XnumNewCol, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
		int XnumNewNZ, const int *XAstart, const int *XAindex, const double *XAvalue,
		const HighsOptions& options, const bool force) {
  int return_code = append_lp_cols(lp, XnumNewCol, XcolCost, XcolLower, XcolUpper,
				  XnumNewNZ, XAstart, XAindex, XAvalue,
				  options, force);
  lp.numCol_ += XnumNewCol;
  return return_code;
}

int append_lp_cols(HighsLp& lp,
		   int XnumNewCol, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
		   int XnumNewNZ, const int *XAstart, const int *XAindex, const double *XAvalue,
		   const HighsOptions& options, const bool force) {
  int newNumCol = lp.numCol_ + XnumNewCol;

  // Assess the bounds and matrix indices, returning on error unless addition is forced
  int return_code = assess_col_bounds(XnumNewCol, XcolLower, XcolUpper, options.infinite_bound);
  if (return_code && !force) return return_code;
  return_code = assess_matrix_indices(lp.numRow_, XnumNewCol, XnumNewNZ, XAstart, XAindex);
  if (return_code && !force) return return_code;

  // Append the columns to the LP vectors and matrix
  append_cols_to_lp_vectors(lp, XnumNewCol, XcolCost, XcolLower, XcolUpper);
  append_cols_to_lp_matrix(lp, XnumNewCol, XnumNewNZ, XAstart, XAindex, XAvalue);

  // Normalise the new LP column bounds and matrix columns
  int numChangedBounds = normalise_col_bounds(lp, lp.numCol_, newNumCol, options.infinite_bound);
  int numRemovedEntries = normalise_matrix_entries(lp, lp.numCol_, newNumCol, options.small_matrix_value);
  return return_code;
}

int assess_col_bounds(int XnumCol, const double* XcolLower, const double* XcolUpper, double infinite_bound) {
  assert(XnumCol >= 0);
  if (XnumCol == 0) return 0;

  int return_code = 0;
  for (int col = 0; col < XnumCol; col++) {
    bool legalLowerBound = XcolLower[col] < infinite_bound;
    assert(legalLowerBound);
    if (!legalLowerBound) {
      HighsLogMessage(HighsMessageType::ERROR, "Col %12d has lower bound of %12g >= %12g = Infinity",
		      col, XcolLower[col], infinite_bound);
      if (!return_code) return_code = col + 1;
    }
    bool legalLowerUpperBound = XcolLower[col] <= XcolUpper[col];
    if (!legalLowerUpperBound) {
      HighsLogMessage(HighsMessageType::WARNING, "Col %12d has inconsistent bounds [%12g, %12g]", col, XcolLower[col], XcolUpper[col]);
      if (!return_code) return_code = XnumCol + col + 1;
    }
    bool legalUpperBound = XcolUpper[col] > -infinite_bound;
    assert(legalUpperBound);
    if (!legalUpperBound) {
      HighsLogMessage(HighsMessageType::ERROR, "Col %12d has upper bound of %12g <= %12g = -Infinity",
		      col, XcolUpper[col], -infinite_bound);
      if (!return_code) return_code = -(col + 1);
    }
  }
  return return_code;
}

void append_cols_to_lp_vectors(HighsLp &lp, int XnumNewCol,
			       const double*XcolCost, const double *XcolLower, const double *XcolUpper) {
  assert(XnumNewCol >= 0);
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  lp.colCost_.resize(newNumCol);
  lp.colLower_.resize(newNumCol);
  lp.colUpper_.resize(newNumCol);

  for (int col = 0; col < XnumNewCol; col++) {
    lp.colCost_[lp.numCol_ + col] = XcolCost[col];
    lp.colLower_[lp.numCol_ + col] = XcolLower[col];
    lp.colUpper_[lp.numCol_ + col] = XcolUpper[col];
  }
}

int normalise_col_bounds(HighsLp& lp, int XfromCol, int XtoCol, const double infinite_bound) {
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
  int numChangedBounds = numChangedLowerBounds + numChangedUpperBounds;
  return numChangedBounds;
}

int add_lp_rows(HighsLp& lp,
		int XnumNewRow, const double *XrowLower,  const double *XrowUpper,
		int XnumNewNZ, const int *XARstart, const int *XARindex, const double *XARvalue,
		const HighsOptions& options, const bool force) {
  int return_code = append_lp_rows(lp, XnumNewRow, XrowLower, XrowUpper,
				  XnumNewNZ, XARstart, XARindex, XARvalue,
				  options, force);
  lp.numRow_ += XnumNewRow;
  return return_code;
}

int append_lp_rows(HighsLp& lp,
		   int XnumNewRow, const double *XrowLower,  const double *XrowUpper,
		   int XnumNewNZ, const int *XARstart, const int *XARindex, const double *XARvalue,
		   const HighsOptions& options, const bool force) {
  int newNumRow = lp.numRow_ + XnumNewRow;
  int return_code = assess_row_bounds(XnumNewRow, XrowLower, XrowUpper, options.infinite_bound);
  if (return_code && !force) return return_code;
  append_rows_to_lp_vectors(lp, XnumNewRow, XrowLower, XrowUpper);
  int numChangedBounds = normalise_row_bounds(lp, lp.numRow_, newNumRow, options.infinite_bound);
  return return_code;
}

int assess_row_bounds(int XnumRow, const double* XrowLower, const double* XrowUpper, double infinite_bound) {
  assert(XnumRow >= 0);
  if (XnumRow == 0) return 0;

  int return_code = 0;
  for (int row = 0; row < XnumRow; row++) {
    bool legalLowerBound = XrowLower[row] < infinite_bound;
    assert(legalLowerBound);
    if (!legalLowerBound) {
      HighsLogMessage(HighsMessageType::ERROR, "Row %12d has lower bound of %12g >= %12g = Infinity",
		      row, XrowLower[row], infinite_bound);
      if (!return_code) return_code = row + 1;
    }
    bool legalLowerUpperBound = XrowLower[row] <= XrowUpper[row];
    if (!legalLowerUpperBound) {
      HighsLogMessage(HighsMessageType::WARNING, "Row %12d has inconsistent bounds [%12g, %12g]", row, XrowLower[row], XrowUpper[row]);
      if (!return_code) return_code = XnumRow + row + 1;
    }
    bool legalUpperBound = XrowUpper[row] > -infinite_bound;
    assert(legalUpperBound);
    if (!legalUpperBound) {
      HighsLogMessage(HighsMessageType::ERROR, "Row %12d has upper bound of %12g <= %12g = -Infinity",
		      row, XrowUpper[row], -infinite_bound);
      if (!return_code) return_code = -(row + 1);
    }
  }
  return return_code;
}

void append_rows_to_lp_vectors(HighsLp &lp, int XnumNewRow,
			       const double *XrowLower, const double *XrowUpper) {
  assert(XnumNewRow >= 0);
  if (XnumNewRow == 0) return;
  int newNumRow = lp.numRow_ + XnumNewRow;
  lp.rowLower_.resize(newNumRow);
  lp.rowUpper_.resize(newNumRow);

  for (int row = 0; row < XnumNewRow; row++) {
    lp.rowLower_[lp.numRow_ + row] = XrowLower[row];
    lp.rowUpper_[lp.numRow_ + row] = XrowUpper[row];
  }
}

int normalise_row_bounds(HighsLp& lp, int XfromRow, int XtoRow, double infinite_bound) {
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
  int numChangedBounds = numChangedLowerBounds + numChangedUpperBounds;
  return numChangedBounds;
}

int assess_matrix_indices(int XnumRow, int XnumCol, int XnumNZ, const int* XAstart, const int* XAindex) {
  assert(XnumRow >= 0);
  assert(XnumCol >= 0);
  assert(XnumNZ >= 0);
  if (XnumRow == 0) return 0;
  if (XnumCol == 0) return 0;
  if (XnumNZ == 0) return 0;

  int return_code = 0;
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
    bool legalFromEl = fromEl >=0 && fromEl <= toEl+1 && fromEl <= XnumNZ;
    assert(legalFromEl);
    if (!legalFromEl) {
      HighsLogMessage(HighsMessageType::ERROR, "Col %12d has illegal start %12d wrt 0, next start %12d or number of nonzeros %12d",
		      col, fromEl, toEl, XnumNZ);
      if (!return_code) return_code = -(XnumCol + col + 1);
    }
    for (int el = fromEl; el <= toEl; el++) {
      int row = XAindex[el];
      bool legalRow = row >= 0;
      assert(legalRow);
      if (!legalRow) {
	HighsLogMessage(HighsMessageType::ERROR, "Col %12d has illegal index %12d", col, row);
	if (!return_code) return_code = col + 1;
      }
      legalRow = row < XnumRow;
      assert(legalRow);
      if (!legalRow) {
	HighsLogMessage(HighsMessageType::ERROR, "Col %12d has illegal index %12d >= %12g", col, row, XnumRow);
	if (!return_code) return_code = -(col + 1);
      }
      legalRow = colVec[row] == 0;
      assert(legalRow);
      if (!legalRow) {
	  HighsLogMessage(HighsMessageType::ERROR, "Col %12d has duplicate index %12d", col, row);	  
	  if (!return_code) return_code = XnumCol + col + 1;
	}
      }
    for (int el = XAstart[col]; el <= toEl; el++) colVec[XAindex[el]] = 0;
  }
  int col = XnumCol;
  int fromEl = XAstart[col];
  bool legalFromEl = fromEl >=0 && fromEl <= XnumNZ;
  assert(legalFromEl);
  if (!legalFromEl) {
    HighsLogMessage(HighsMessageType::ERROR, "Col %12d has illegal start %12d wrt 0 or number of nonzeros %12d",
		    col, fromEl, XnumNZ);
    if (!return_code) return_code = -(XnumCol + col + 1);
  }
  return return_code;

}


void append_cols_to_lp_matrix(HighsLp &lp, int XnumNewCol,
			      int XnumNewNZ, const int *XAstart, const int *XAindex, const double *XAvalue) {
  assert(XnumNewCol >= 0);
  assert(XnumNewNZ >= 0);
  if (XnumNewCol == 0) return;
  int newNumCol = lp.numCol_ + XnumNewCol;
  lp.Astart_.resize(newNumCol + 1);
  for (int col = 0; col < XnumNewCol; col++) {
    lp.Astart_[lp.numCol_ + col + 1] = lp.Astart_[lp.numCol_];
  }
  
  // Determine the current number of nonzeros
  int currentNumNZ = lp.Astart_[lp.numCol_];
  
  // Determine the new number of nonzeros and resize the column-wise matrix
  // arrays
  int newNumNZ = currentNumNZ + XnumNewNZ;
  lp.Aindex_.resize(newNumNZ);
  lp.Avalue_.resize(newNumNZ);
  
  // Append the new columns
  for (int col = 0; col < XnumNewCol; col++) {
    lp.Astart_[lp.numCol_ + col] = XAstart[col] + currentNumNZ;
  }
  lp.Astart_[lp.numCol_ + XnumNewCol] = newNumNZ;
  
  for (int el = 0; el < XnumNewNZ; el++) {
    int row = XAindex[el];
    assert(row >= 0);
    assert(row < lp.numRow_);
    lp.Aindex_[currentNumNZ + el] = row;
    lp.Avalue_[currentNumNZ + el] = XAvalue[el];
  }
}

void append_rows_to_lp_matrix(HighsLp &lp, int XnumNewRow,
			      int XnumNewNZ, const int *XARstart, const int *XARindex, const double *XARvalue) {
  assert(XnumNewRow >= 0);
  assert(XnumNewNZ >= 0);
  // Check that nonzeros aren't being appended to a matrix with no columns
  assert(XnumNewNZ == 0 || lp.numCol_ > 0);
  if (XnumNewRow == 0) return;
  int newNumRow = lp.numRow_ + XnumNewRow;

  // NB SCIP doesn't have XARstart[XnumNewRow] defined, so have to use XnumNewNZ for last
  // entry
  if (XnumNewNZ == 0) return;
  int currentNumNZ = lp.Astart_[lp.numCol_];
  vector<int> Alength;
  Alength.assign(lp.numCol_, 0);
  for (int el = 0; el < XnumNewNZ; el++) {
    int col = XARindex[el];
    //      printf("El %2d: appending entry in column %2d\n", el, col); 
    assert(col >= 0);
    assert(col < lp.numCol_);
    Alength[col]++;
  }
  // Determine the new number of nonzeros and resize the column-wise matrix arrays
  int newNumNZ = currentNumNZ + XnumNewNZ;
  lp.Aindex_.resize(newNumNZ);
  lp.Avalue_.resize(newNumNZ);

  // Append the new rows
  // Shift the existing columns to make space for the new entries
  int nwEl = newNumNZ;
  for (int col = lp.numCol_ - 1; col >= 0; col--) {
    // printf("Column %2d has additional length %2d\n", col, Alength[col]);
    int Astart_Colp1 = nwEl;
    nwEl -= Alength[col];
    // printf("Shift: nwEl = %2d\n", nwEl);
    for (int el = lp.Astart_[col + 1] - 1; el >= lp.Astart_[col]; el--) {
      nwEl--;
      // printf("Shift: Over-writing lp.Aindex_[%2d] with lp.Aindex_[%2d]=%2d\n",
      // nwEl, el, lp.Aindex_[el]);
      lp.Aindex_[nwEl] = lp.Aindex_[el];
      lp.Avalue_[nwEl] = lp.Avalue_[el];
    }
    lp.Astart_[col + 1] = Astart_Colp1;
  }
  // printf("After shift: nwEl = %2d\n", nwEl);
  assert(nwEl == 0);
  // util_reportColMtx(lp.numCol_, lp.Astart_, lp.Aindex_, lp.Avalue_);

  // Insert the new entries
  for (int row = 0; row < XnumNewRow; row++) {
    int fEl = XARstart[row];
    int lEl = (row < XnumNewRow - 1 ? XARstart[row + 1] : XnumNewNZ) - 1;
    for (int el = fEl; el <= lEl; el++) {
      int col = XARindex[el];
      nwEl = lp.Astart_[col + 1] - Alength[col];
      Alength[col]--;
      // printf("Insert: row = %2d; col = %2d; lp.Astart_[col+1]-Alength[col] =
      // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
      // lp.Astart_[col+1]-Alength[col], Alength[col], nwEl);
      assert(nwEl >= 0);
      assert(el >= 0);
      // printf("Insert: Over-writing lp.Aindex_[%2d] with lp.Aindex_[%2d]=%2d\n",
      // nwEl, el, lp.Aindex_[el]);
      lp.Aindex_[nwEl] = lp.numRow_ + row;
      lp.Avalue_[nwEl] = XARvalue[el];
    }
  }
}

int normalise_matrix_entries(HighsLp& lp, int XfromCol, int XtoCol, double small_matrix_value) {
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
  return numRemovedValues;
}

int normalise_row_matrix_entries(int XnumCol, int XnumRow, int XnumNZ, int* XARstart, int* XARindex, double* XARvalue, double small_matrix_value) {
  assert(XnumCol >= 0);
  assert(XnumRow >= 0);
  assert(XnumNZ >= 0);
  assert(small_matrix_value >=0);
  if (XnumRow == 0) return 0;
  if (XnumNZ == 0) return 0;
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
  return numRemovedValues;
}


void del_cols_from_lp_vectors(HighsLp &lp, int XfromCol, int XtoCol) {
  assert(XfromCol >= 0);
  assert(XtoCol < lp.numCol_);
  assert(XfromCol <= XtoCol);

  int numDeleteCol = XtoCol - XfromCol + 1;
  if (numDeleteCol == 0 || numDeleteCol == lp.numCol_) return;
  //
  // Trivial case is XtoCol = lp.numCol_-1, in which case no columns
  // need be shifted. However, this implies lp.numCol_-numDeleteCol =
  // XfromCol, in which case the loop is vacuous
  for (int col = XfromCol; col < lp.numCol_ - numDeleteCol; col++) {
    lp.colCost_[col] = lp.colCost_[col + numDeleteCol];
    lp.colLower_[col] = lp.colLower_[col + numDeleteCol];
    lp.colUpper_[col] = lp.colUpper_[col + numDeleteCol];
  }
}

void del_cols_from_lp_matrix(HighsLp &lp, int XfromCol, int XtoCol) {
  assert(XfromCol >= 0);
  assert(XtoCol < lp.numCol_);
  assert(XfromCol <= XtoCol);

  int numDeleteCol = XtoCol - XfromCol + 1;
  if (numDeleteCol == 0 || numDeleteCol == lp.numCol_) return;
  //
  // Trivial case is XtoCol = lp.numCol_-1, in which case no columns need be shifted
  // and the loops are vacuous
  int elOs = lp.Astart_[XfromCol];
  int numDeleteEl = lp.Astart_[XtoCol + 1] - elOs;
  for (int el = lp.Astart_[XtoCol + 1]; el < lp.Astart_[lp.numCol_]; el++) {
    lp.Aindex_[el - numDeleteEl] = lp.Aindex_[el];
    lp.Avalue_[el - numDeleteEl] = lp.Avalue_[el];
  }
  for (int col = XfromCol; col <= lp.numCol_ - numDeleteCol; col++) {
    lp.Astart_[col] = lp.Astart_[col + numDeleteCol] - numDeleteEl;
  }

}

void del_rows_from_lp_vectors(HighsLp &lp, int XfromRow, int XtoRow) {
  assert(XfromRow >= 0);
  assert(XtoRow < lp.numRow_);
  assert(XfromRow <= XtoRow);

  int numDeleteRow = XtoRow - XfromRow + 1;
  if (numDeleteRow == 0 || numDeleteRow == lp.numRow_) return;
  //
  // Trivial case is XtoRow = lp.numRow_-1, in which case no rows
  // need be shifted. However, this implies lp.numRow_-numDeleteRow =
  // XfromRow, in which case the loop is vacuous
  for (int row = XfromRow; row < lp.numRow_ - numDeleteRow; row++) {
    lp.rowLower_[row] = lp.rowLower_[row + numDeleteRow];
    lp.rowUpper_[row] = lp.rowUpper_[row + numDeleteRow];
  }
}

void del_rows_from_lp_matrix(HighsLp &lp, int XfromRow, int XtoRow) {
  assert(XfromRow >= 0);
  assert(XtoRow < lp.numRow_);
  assert(XfromRow <= XtoRow);

  int numDeleteRow = XtoRow - XfromRow + 1;
  if (numDeleteRow == 0 || numDeleteRow == lp.numRow_) return;

  int nnz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    int fmEl = lp.Astart_[col];
    lp.Astart_[col] = nnz;
    for (int el = fmEl; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      if (row < XfromRow || row > XtoRow) {
	if (row < XfromRow) {
	  lp.Aindex_[nnz] = row;
	} else {
	  lp.Aindex_[nnz] = row - numDeleteRow;
	}
	lp.Avalue_[nnz] = lp.Avalue_[el];
	nnz++;
      }
    }
  }
  lp.Astart_[lp.numCol_] = nnz;
}

void change_lp_matrix_coefficient(HighsLp &lp, int Xrow, int Xcol, const double XnewValue) {
  int changeElement = -1;
  for (int el = lp.Astart_[Xcol]; el < lp.Astart_[Xcol + 1]; el++) {
    // printf("Column %d: Element %d is row %d. Is it %d?\n", Xcol, el, lp.Aindex_[el], Xrow);
    if (lp.Aindex_[el] == Xrow) {
      changeElement = el;
      break;
    }
  }
  if (changeElement < 0) {
    //    printf("util_changeCoeff: Cannot find row %d in column %d\n", Xrow, Xcol);
    changeElement = lp.Astart_[Xcol + 1];
    int newNumNZ = lp.Astart_[lp.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    lp.Astart_[lp.numCol_], newNumNZ);
    lp.Aindex_.resize(newNumNZ);
    lp.Avalue_.resize(newNumNZ);
    for (int i = Xcol + 1; i <= lp.numCol_; i++) lp.Astart_[i]++;
    for (int el = newNumNZ - 1; el > changeElement; el--) {
      lp.Aindex_[el] = lp.Aindex_[el - 1];
      lp.Avalue_[el] = lp.Avalue_[el - 1];
    }
  }
  lp.Aindex_[changeElement] = Xrow;
  lp.Avalue_[changeElement] = XnewValue;
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
