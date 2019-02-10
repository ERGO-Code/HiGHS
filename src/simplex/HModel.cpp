/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HModel.h"
#include "HConst.h"
#include "HMPSIO.h"
#include "HighsIO.h"
#include "Presolve.h"
#include "HToyIO.h"
#include "HVector.h"

#include "SimplexTimer.h" // For timer
#include "HighsLpUtils.h" // For util_anMl
#include "HighsUtils.h" // For highs_isInfinity

// For compute dual objective alt value
#include "HighsModelObject.h"
#include "HSimplex.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::swap;
using std::fabs;
using std::ofstream;
using std::setprecision;
using std::setw;

// Methods which load whole models, initialise the basis then
// allocate and populate (where possible) work* arrays and
// allocate basis* arrays
HModel::HModel() {
//  clear_solver_lp(highs_model_object);
}

void HModel::rp_basis() {
  printf("\nReporting current basis: solver_lp_->numCol_ = %d; solver_lp_->numRow_ = %d\n", solver_lp_->numCol_,
         solver_lp_->numRow_);
  if (solver_lp_->numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < solver_lp_->numCol_; col++) {
    int var = col;
    if (basis_->nonbasicFlag_[var])
      printf("%6d %6d        %6d %6d\n", var, col, basis_->nonbasicFlag_[var],
             basis_->nonbasicMove_[var]);
    else
      printf("%6d %6d %6d\n", var, col, basis_->nonbasicFlag_[var]);
  }
  if (solver_lp_->numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = solver_lp_->numCol_ + row;
    if (basis_->nonbasicFlag_[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basis_->basicIndex_[row],
             basis_->nonbasicFlag_[var], basis_->nonbasicMove_[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basis_->basicIndex_[row], basis_->nonbasicFlag_[var]);
  }
}

int HModel::get_nonbasicMove(int var) {
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  //  printf("Calling get_nonbasicMove with var = %2d; numTot = %2d\n", var,
  //  numTot); cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  if (!highs_isInfinity(-simplex_info_->workLower_[var])) {
    if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info_->workLower_[var] == simplex_info_->workUpper_[var])
        // Fixed variable so nonbasic move is zero
        return NONBASIC_MOVE_ZE;
      // Boxed variable so nonbasic move is up (from lower bound)
      return NONBASIC_MOVE_UP;
    } else
      // Finite lower bound and infinite upper bound so nonbasic move is up
      // (from lower bound)
      return NONBASIC_MOVE_UP;
  } else
      // Infinite lower bound so nonbasic move depends on whether the upper
      // bound is finite
      if (!highs_isInfinity(simplex_info_->workUpper_[var]))
    // Finite upper bound so nonbasic move is down (from upper bound)
    return NONBASIC_MOVE_DN;
  // Infinite upper bound so free variable: nonbasic move is zero
  return NONBASIC_MOVE_ZE;
}

// ???? Housekeeping done from here down ????

#ifdef HiGHSDEV
void HModel::changeUpdate(int updateMethod) { factor_->change(updateMethod); }
#endif

int HModel::writeToMPS(const char *filename) {
  vector<int> integerColumn;
  int numInt = 0;
  int rtCd =
      writeMPS(filename, solver_lp_->numRow_, solver_lp_->numCol_, numInt, solver_lp_->sense_, solver_lp_->offset_, solver_lp_->Astart_,
               solver_lp_->Aindex_, solver_lp_->Avalue_, solver_lp_->colCost_, solver_lp_->colLower_, solver_lp_->colUpper_,
               solver_lp_->rowLower_, solver_lp_->rowUpper_, integerColumn);
  return rtCd;
}

//>->->->->->->->->->->->->->->->->->->->->->-
// Esoterica!
void HModel::shiftObjectiveValue(double shift) {
  simplex_info_->dualObjectiveValue += shift;
}
//<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-

// Scale a pair of row bound values
void HModel::util_scaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue) {
  (*XrowLowerValue) = (highs_isInfinity(-(*XrowLowerValue)) ? (*XrowLowerValue): (*XrowLowerValue)*scale_->row_[iRow]);
  (*XrowUpperValue) = (highs_isInfinity( (*XrowUpperValue)) ? (*XrowUpperValue): (*XrowUpperValue)*scale_->row_[iRow]);
}

// Scale a pair of column bound values
void HModel::util_scaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue) {
  (*XcolLowerValue) = (highs_isInfinity(-(*XcolLowerValue)) ? (*XcolLowerValue): (*XcolLowerValue)/scale_->col_[iCol]);
  (*XcolUpperValue) = (highs_isInfinity( (*XcolUpperValue)) ? (*XcolUpperValue): (*XcolUpperValue)/scale_->col_[iCol]);
}

// Scale a column cost
void HModel::util_scaleColCostValue(int iCol, double* XcolCostValue) {
  (*XcolCostValue) = (*XcolCostValue)*scale_->col_[iCol];
}

// Unscale a pair of row bound values
void HModel::util_unscaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue) {
  (*XrowLowerValue) = (highs_isInfinity(-(*XrowLowerValue)) ? (*XrowLowerValue): (*XrowLowerValue)/scale_->row_[iRow]);
  (*XrowUpperValue) = (highs_isInfinity( (*XrowUpperValue)) ? (*XrowUpperValue): (*XrowUpperValue)/scale_->row_[iRow]);
}

// Unscale a pair of column bound values
void HModel::util_unscaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue) {
  (*XcolLowerValue) = (highs_isInfinity(-(*XcolLowerValue)) ? (*XcolLowerValue): (*XcolLowerValue)*scale_->col_[iCol]);
  (*XcolUpperValue) = (highs_isInfinity( (*XcolUpperValue)) ? (*XcolUpperValue): (*XcolUpperValue)*scale_->col_[iCol]);
}

// Unscale a column cost
void HModel::util_unscaleColCostValue(int iCol, double* XcolCostValue) {
  (*XcolCostValue) = (*XcolCostValue)/scale_->col_[iCol];
}


// Get the column and row (primal) values and dual (values)
void HModel::util_getPrimalDualValues(vector<double> &XcolValue,
                                      vector<double> &XcolDual,
                                      vector<double> &XrowValue,
                                      vector<double> &XrowDual
				      ) {
  // Take primal solution
  vector<double> value = simplex_info_->workValue_;
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++)
    value[basis_->basicIndex_[iRow]] = simplex_info_->baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info_->workDual_;
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++) dual[basis_->basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < solver_lp_->numCol_; iCol++) {
    value[iCol] *= scale_->col_[iCol];
    dual[iCol] /= (scale_->col_[iCol] / scale_->cost_);
  }
  for (int iRow = 0, iTot = solver_lp_->numCol_; iRow < solver_lp_->numRow_; iRow++, iTot++) {
    value[iTot] /= scale_->row_[iRow];
    dual[iTot] *= (scale_->row_[iRow] * scale_->cost_);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  XcolValue.resize(solver_lp_->numCol_);
  XcolDual.resize(solver_lp_->numCol_);
  XrowValue.resize(solver_lp_->numRow_);
  XrowDual.resize(solver_lp_->numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < solver_lp_->numRow_; i++) XrowValue[i] = -valuePtr[i + solver_lp_->numCol_];
  for (int i = 0; i < solver_lp_->numCol_; i++) XcolValue[i] = valuePtr[i];
  for (int i = 0; i < solver_lp_->numRow_; i++) XrowDual[i] = solver_lp_->sense_ * dual[i + solver_lp_->numCol_];
  for (int i = 0; i < solver_lp_->numCol_; i++) XcolDual[i] = solver_lp_->sense_ * dual[i];
}

void HModel::util_getNonbasicMove(vector<int> &XnonbasicMove) {
  XnonbasicMove = basis_->nonbasicMove_;
}

void HModel::util_getBasicIndexNonbasicFlag(vector<int> &XbasicIndex,
                                            vector<int> &XnonbasicFlag
					    ) {
  XbasicIndex.resize(solver_lp_->numRow_);
  XnonbasicFlag.resize(basis_->nonbasicFlag_.size());
  int basicIndexSz = basis_->basicIndex_.size();
  for (int i = 0; i < basicIndexSz; i++) XbasicIndex[i] = basis_->basicIndex_[i];
  int nonbasicFlagSz = basis_->nonbasicFlag_.size();
  for (int i = 0; i < nonbasicFlagSz; i++) XnonbasicFlag[i] = basis_->nonbasicFlag_[i];
}

// Utilities to get/change costs and bounds
// Get the costs for a contiguous set of columns
void HModel::util_getCosts(HighsLp& lp, int firstcol, int lastcol, double *XcolCost) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) XcolCost[col - firstcol] = lp.colCost_[col];
}

// Get the bounds for a contiguous set of columns
void HModel::util_getColBounds(HighsLp& lp, int firstcol, int lastcol, double *XcolLower,
                               double *XcolUpper) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) {
    if (XcolLower != NULL) XcolLower[col - firstcol] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - firstcol] = lp.colUpper_[col];
  }
}

// Get the bounds for a contiguous set of rows
void HModel::util_getRowBounds(HighsLp& lp, int firstrow, int lastrow, double *XrowLower,
                               double *XrowUpper) {
  assert(0 <= firstrow);
  assert(firstrow <= lastrow);
  assert(lastrow < lp.numRow_);
  for (int row = firstrow; row <= lastrow; ++row) {
    if (XrowLower != NULL) XrowLower[row - firstrow] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - firstrow] = lp.rowUpper_[row];
  }
}

// Possibly change the objective sense
int HModel::util_chgObjSense(const int Xsense) {
  if ((Xsense == OBJSENSE_MINIMIZE) != (solver_lp_->sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the objective sense
    solver_lp_->sense_ = Xsense;
    const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
    for (int var = 0; var < numTot; var++) {
      simplex_info_->workDual_[var] = -simplex_info_->workDual_[var];
      simplex_info_->workCost_[var] = -simplex_info_->workCost_[var];
    }
    simplex_info_->solution_status = SimplexSolutionStatus::UNSET;
  }
  return 0;
}

// Change the costs for all columns
int HModel::util_chgCostsAll(const double *XcolCost) {
  assert(XcolCost != NULL);
  for (int col = 0; col < solver_lp_->numCol_; ++col) {
    solver_lp_->colCost_[col] = XcolCost[col] * scale_->col_[col];
  }
  // Deduce the consequences of new costs
  //  update_solver_lp_status_flags(highs_model, LpAction::NEW_COSTS);
  return 0;
}

// Change the costs for a set of columns
int HModel::util_chgCostsSet(int ncols, const int *XcolCostIndex,
                             const double *XcolCostValues) {
  assert(XcolCostIndex != NULL);
  assert(XcolCostValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolCostIndex[ix];
    assert(0 <= col);
    assert(col < solver_lp_->numCol_);
    solver_lp_->colCost_[col] = XcolCostValues[ix] * scale_->col_[col];
  }
  // Deduce the consequences of new costs
  //  update_solver_lp_status_flags(highs_model, LpAction::NEW_COSTS);
  return 0;
}

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
int HModel::util_chgColBoundsAll(const double *XcolLower,
                                 const double *XcolUpper) {
  assert(XcolLower != NULL);
  assert(XcolUpper != NULL);
  for (int col = 0; col < solver_lp_->numCol_; ++col) {
    double lower = XcolLower[col];
    double upper = XcolUpper[col];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    solver_lp_->colLower_[col] = (highs_isInfinity(-lower) ? lower : lower / scale_->col_[col]);
    solver_lp_->colUpper_[col] = (highs_isInfinity(upper) ? upper : upper / scale_->col_[col]);
    //    printf("[LB; Pr; UB] for column %2d are now [%11g, %11g, %11g] Dual =
    //    %g\n", col, solver_lp_->colLower_[col], simplex_info_->workValue_[col], solver_lp_->colUpper_[col],
    //    simplex_info_->workDual_[col]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}

// Change the bounds for a set of columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
int HModel::util_chgColBoundsSet(int ncols, const int *XcolBoundIndex,
                                 const double *XcolLowerValues,
                                 const double *XcolUpperValues) {
  assert(XcolBoundIndex != NULL);
  assert(XcolLowerValues != NULL);
  assert(XcolUpperValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolBoundIndex[ix];
    assert(0 <= col);
    assert(col < solver_lp_->numCol_);
    double lower = XcolLowerValues[ix];
    double upper = XcolUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    solver_lp_->colLower_[col] = (highs_isInfinity(-lower) ? lower : lower / scale_->col_[col]);
    solver_lp_->colUpper_[col] = (highs_isInfinity(upper) ? upper : upper / scale_->col_[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col,
    //    solver_lp_->colLower_[col], solver_lp_->colUpper_[col], scale_->col_[col]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}

// Change the bounds for all rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
int HModel::util_chgRowBoundsAll(const double *XrowLower,
                                 const double *XrowUpper) {
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < solver_lp_->numRow_; ++row) {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(row + 1);
    solver_lp_->rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale_->row_[row]);
    solver_lp_->rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale_->row_[row]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}

// Change the bounds for a set of rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
int HModel::util_chgRowBoundsSet(int nrows, const int *XrowBoundIndex,
                                 const double *XrowLowerValues,
                                 const double *XrowUpperValues) {
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < nrows; ++ix) {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < solver_lp_->numRow_);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(row + 1);
    solver_lp_->rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale_->row_[row]);
    solver_lp_->rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale_->row_[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row,
    //    solver_lp_->rowLower_[row], solver_lp_->rowUpper_[row]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}

// Utilities to convert model basic/nonbasic status to/from SCIP-like status
// Convert model basic/nonbasic status from SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertBaseStatToWorking(const int *cstat, const int *rstat) {
  int numBasic = 0;
  for (int col = 0; col < solver_lp_->numCol_; col++) {
    int var = col;
    if (cstat[col] == HIGHS_BASESTAT_BASIC) {
      basis_->nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis_->basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == HIGHS_BASESTAT_LOWER) {
      // HIGHS_BASESTAT_LOWER includes fixed variables
      if (solver_lp_->colLower_[col] == solver_lp_->colUpper_[col]) {
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
        continue;
      }
    } else if (cstat[col] == HIGHS_BASESTAT_UPPER) {
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_DN;
      continue;
    } else if (cstat[col] == HIGHS_BASESTAT_ZERO) {
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], solver_lp_->colLower_[col], solver_lp_->colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = solver_lp_->numCol_ + row;
    if (rstat[row] == HIGHS_BASESTAT_BASIC) {
      basis_->nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis_->basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == HIGHS_BASESTAT_LOWER) {
      // HIGHS_BASESTAT_LOWER includes fixed variables
      if (solver_lp_->rowLower_[row] == solver_lp_->rowUpper_[row]) {
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_DN;
        continue;
      }
    } else if (rstat[row] == HIGHS_BASESTAT_UPPER) {
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
      continue;
    } else if (rstat[row] == HIGHS_BASESTAT_ZERO) {
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], solver_lp_->rowLower_[row], solver_lp_->rowUpper_[row]);
#endif
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], solver_lp_->rowLower_[row], solver_lp_->rowUpper_[row]);
      return -(row + 1);
    }
    printf(
        "convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, "
        "nonbasicMove=%d\n",
        row, rstat[row], solver_lp_->rowLower_[row], solver_lp_->rowUpper_[row], basis_->nonbasicMove_[var]);
  }
  assert(numBasic = solver_lp_->numRow_);
  printf("Call simplex_method_.populate_work_arrays();\n");
  // simplex_method.update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  return 0;
}

// Convert model basic/nonbasic status to SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertWorkingToBaseStat(int *cstat, int *rstat) {
  if (cstat != NULL) {
    for (int col = 0; col < solver_lp_->numCol_; col++) {
      int var = col;
      if (!basis_->nonbasicFlag_[var]) {
        cstat[col] = HIGHS_BASESTAT_BASIC;
        continue;
      } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-solver_lp_->colLower_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(solver_lp_->colUpper_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        //	printf("Var %d Move = %d [%g, %g]\n", var, basis_->nonbasicMove_[var],
        // solver_lp_->colLower_[col], solver_lp_->colUpper_[col]);
        if (solver_lp_->colLower_[col] == solver_lp_->colUpper_[col]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(solver_lp_->colUpper_[col]))
#endif
          {
            cstat[col] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-solver_lp_->colLower_[col]) && highs_isInfinity(solver_lp_->colUpper_[col]))
#endif
          {
            cstat[col] = HIGHS_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          col, basis_->nonbasicFlag_[var], basis_->nonbasicMove_[var], solver_lp_->colLower_[col],
          solver_lp_->colUpper_[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL) {
    for (int row = 0; row < solver_lp_->numRow_; row++) {
      int var = solver_lp_->numCol_ + row;
      if (!basis_->nonbasicFlag_[var]) {
        rstat[row] = HIGHS_BASESTAT_BASIC;
        continue;
      }
      // NB nonbasicMove for rows refers to the solver's view where the bounds
      // are switched and negated
      else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN)
      // Free to move only down from -solver_lp_->rowLower_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-solver_lp_->rowLower_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -solver_lp_->rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(solver_lp_->rowUpper_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        if (solver_lp_->rowLower_[row] == solver_lp_->rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(solver_lp_->rowUpper_[row]))
#endif
          {
            rstat[row] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-solver_lp_->rowLower_[row]) && highs_isInfinity(solver_lp_->rowUpper_[row]))
#endif
          {
            rstat[row] = HIGHS_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          row, basis_->nonbasicFlag_[var], basis_->nonbasicMove_[var], solver_lp_->rowLower_[row],
          solver_lp_->rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}

// Utility to get the indices of the basic variables for SCIP
int HModel::util_getBasicIndices(int *bind) {
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = basis_->basicIndex_[row];
    if (var >= solver_lp_->numCol_)
      bind[row] = -(1 + var - solver_lp_->numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

// Methods for brief reports
// is false
void HModel::util_reportNumberIterationObjectiveValue(int i_v) {
  HighsPrintMessage(ML_MINIMAL, "%10d  %20.10e  %2d\n", simplex_info_->iteration_count, simplex_info_->dualObjectiveValue, i_v);
}

void HModel::util_reportSolverOutcome(const char *message) {
  if (simplex_info_->solution_status == SimplexSolutionStatus::OPTIMAL)
    HighsPrintMessage(ML_ALWAYS, "%s: OPTIMAL", message);
  else
    HighsPrintMessage(ML_ALWAYS, "%s: NOT-OPT", message);
  double dualObjectiveValue = simplex_info_->dualObjectiveValue;
#ifdef SCIP_DEV
  double prObjVal = 0; printf("Call simplex_method.compute_primal_objective_function_value\n");
  double dlObjVal =
      abs(prObjVal - dualObjectiveValue) / max(abs(dualObjectiveValue), max(abs(prObjVal), 1.0));
  HighsPrintMessage(ML_MINIMAL, "%32s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f",
		    solver_lp_->model_name_.c_str(),
		    prObjVal,
		    dualObjectiveValue,
		    dlObjVal,
		    simplex_info_->iteration_count,
		    currentRunHighsTime);
#else
  double currentRunHighsTime = timer_->readRunHighsClock();
  HighsPrintMessage(ML_ALWAYS, "%32s %20.10e %10d %10.3f", solver_lp_->model_name_.c_str(), dualObjectiveValue,
         simplex_info_->iteration_count, currentRunHighsTime);
#endif
  HighsPrintMessage(ML_ALWAYS, " [Ph1 = %d; Ph2 = %d; Pr = %d]",
		    simplex_info_->dual_phase1_iteration_count,
		    simplex_info_->dual_phase2_iteration_count,
		    simplex_info_->primal_phase2_iteration_count
		    );
  if (simplex_info_->solution_status == SimplexSolutionStatus::OPTIMAL) {
    HighsPrintMessage(ML_ALWAYS, "\n");
  } else {
    util_reportModelStatus();
    HighsPrintMessage(ML_ALWAYS, "\n");
  }
  // Greppable report line added
  HighsPrintMessage(ML_ALWAYS, "grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s,%d,%d,%d\n",
		    dualObjectiveValue,
		    simplex_info_->iteration_count,
		    currentRunHighsTime,
		    simplex_info_->solution_status,
		    solver_lp_->model_name_.c_str(),
		    simplex_info_->dual_phase1_iteration_count,
		    simplex_info_->dual_phase2_iteration_count,
		    simplex_info_->primal_phase2_iteration_count
		    );
}

// Methods for reporting the model, its solution, row and column data and matrix
// Report the model status
void HModel::util_reportModelStatus() {
  HighsPrintMessage(ML_ALWAYS, "LP status is %2d: ", simplex_info_->solution_status);
  if (simplex_info_->solution_status == SimplexSolutionStatus::UNSET)
    HighsPrintMessage(ML_ALWAYS, "Unset\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::OPTIMAL)
    HighsPrintMessage(ML_ALWAYS, "Optimal\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::INFEASIBLE)
    HighsPrintMessage(ML_ALWAYS, "Infeasible\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::UNBOUNDED)
    HighsPrintMessage(ML_ALWAYS, "Primal unbounded\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::SINGULAR)
    HighsPrintMessage(ML_ALWAYS, "Singular basis\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::FAILED)
    HighsPrintMessage(ML_ALWAYS, "Failed\n");
  else if (simplex_info_->solution_status == SimplexSolutionStatus::OUT_OF_TIME)
    HighsPrintMessage(ML_ALWAYS, "Time limit exceeded\n");
  else
    HighsPrintMessage(ML_ALWAYS, "Unrecognised\n");
}

// The remaining routines are wholly independent of any classes, merely
// printing what's passed in the parameter lists.
//
// ToDo: Consider putting them in a separate class.
void HModel::util_reportRowVecSol(int nrow, vector<double> &XrowLower,
                                  vector<double> &XrowUpper,
                                  vector<double> &XrowPrimal,
                                  vector<double> &XrowDual,
                                  vector<int> &XrowStatus) {
  // Report the LP row data and solution passed to the method, where
  // XrowStatus is the SCIP-like basis status
  if (nrow <= 0) return;
  printf("Row    St      Primal       Lower       Upper        Dual\n");
  for (int row = 0; row < nrow; row++) {
    if (XrowStatus[row] == HIGHS_BASESTAT_BASIC)
      printf("%6d BC", row);
    else if (XrowStatus[row] == HIGHS_BASESTAT_ZERO)
      printf("%6d FR", row);
    else if (XrowStatus[row] == HIGHS_BASESTAT_LOWER) {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6d FX", row);
      else
        printf("%6d LB", row);
    } else if (XrowStatus[row] == HIGHS_BASESTAT_UPPER)
      printf("%6d UB", row);
    else
      printf("%6d ??", row);
    printf(" %11g %11g %11g %11g\n", XrowPrimal[row], XrowLower[row],
           XrowUpper[row], XrowDual[row]);
  }
}

void HModel::util_reportRowMtx(int nrow, vector<int> &XARstart,
                               vector<int> &XARindex,
                               vector<double> &XARvalue) {
  // Report the row-wise matrix passed to the method
  if (nrow <= 0) return;
  printf("Row    Index       Value\n");
  for (int row = 0; row < nrow; row++) {
    printf("%6d Start %8d\n", row, XARstart[row]);
    for (int el = XARstart[row]; el < XARstart[row + 1]; el++) {
      printf("      %6d %11g\n", XARindex[el], XARvalue[el]);
    }
  }
  printf("       Start %8d\n", XARstart[nrow]);
}

void HModel::util_reportColVecSol(int ncol, vector<double> &XcolCost,
                                  vector<double> &colLower,
                                  vector<double> &XcolUpper,
                                  vector<double> &XcolPrimal,
                                  vector<double> &XcolDual,
                                  vector<int> &XcolStatus) {
  // Report the LP column data and solution passed to the method,
  // where XcolStatus is the SCIP-like basis status
  if (ncol <= 0) return;
  printf(
      "Col    St      Primal       Lower       Upper        Dual        "
      "Cost\n");
  for (int col = 0; col < ncol; col++) {
    if (XcolStatus[col] == HIGHS_BASESTAT_BASIC)
      printf("%6d BC", col);
    else if (XcolStatus[col] == HIGHS_BASESTAT_ZERO)
      printf("%6d FR", col);
    else if (XcolStatus[col] == HIGHS_BASESTAT_LOWER) {
      if (colLower[col] == XcolUpper[col])
        printf("%6d FX", col);
      else
        printf("%6d LB", col);
    } else if (XcolStatus[col] == HIGHS_BASESTAT_UPPER)
      printf("%6d UB", col);
    else
      printf("%6d ??", col);
    printf(" %11g %11g %11g %11g %11g\n", XcolPrimal[col], colLower[col],
           XcolUpper[col], XcolDual[col], XcolCost[col]);
  }
}

void HModel::util_reportBasicIndex(const char *message, int numRow, vector<int> &basicIndex) {
  printf("%s: Model has %d basic indices: ", message, numRow);
  for (int i=0; i<numRow; i++){
    printf(" %d", basicIndex[i]);
  }
  printf("\n");
}

void HModel::util_reportModelDa(HighsLp lp, const char *filename) {
  vector<double> wkDseCol;
  wkDseCol.resize(lp.numRow_);
  for (int r_n = 0; r_n < lp.numRow_; r_n++) wkDseCol[r_n] = 0;
  FILE *file = fopen(filename, "w");
  if (file == 0) {
#ifdef HiGHSDEV
    printf("util_reportModelDa: Not opened file OK\n");
#endif
    //      return 1;
  }
  fprintf(file, "%d %d %d\n", lp.numRow_, lp.numCol_, lp.Astart_[lp.numCol_]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColCost, %d, %20g\n", c_n, lp.colCost_[c_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    for (int el_n = lp.Astart_[c_n]; el_n < lp.Astart_[c_n + 1]; el_n++) {
      int r_n = lp.Aindex_[el_n];
      wkDseCol[r_n] = lp.Avalue_[el_n];
    }
    int nnz = 0;
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      if (wkDseCol[r_n] != 0) nnz++;
    }
    fprintf(file, "Mtx_Col_nnz, %d, %d\n", c_n, nnz);
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      if (wkDseCol[r_n] != 0) {
        fprintf(file, "MtxVRow, %d, %20g\n", r_n, wkDseCol[r_n]);
      }
      wkDseCol[r_n] = 0;
    }
  }
  for (int r_n = 0; r_n < lp.numRow_; r_n++)
    fprintf(file, "RowLB, %d, %20g\n", r_n, lp.rowLower_[r_n]);
  for (int r_n = 0; r_n < lp.numRow_; r_n++)
    fprintf(file, "RowUB, %d, %20g\n", r_n, lp.rowUpper_[r_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColLB, %d, %20g\n", c_n, lp.colLower_[c_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColUB, %d, %20g\n", c_n, lp.colUpper_[c_n]);
  fclose(file);
}

#ifdef HiGHSDEV
/*
// Rename appropriately when HModel is split.
void getSolutionFromHModel(const HModel& model, HighsSolution& solution) {
  model.util_getPrimalDualValues(
				 solution.colValue_,
				 solution.colDual_,
				 solution.rowValue_,
				 solution.rowDual_);
}
*/
#endif
