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

void HModel::load_fromArrays(int XnumCol, int Xsense, const double *XcolCost,
                             const double *XcolLower, const double *XcolUpper,
                             int XnumRow, const double *XrowLower,
                             const double *XrowUpper, int XnumNz,
                             const int *XAstart, const int *XAindex,
                             const double *XAvalue) {
  //  printf("load_fromArrays: XnumCol = %d; XnumRow = %d; XnumNz = %d\n",
  //  XnumCol, XnumRow, XnumNz);
  assert(XnumCol > 0);
  assert(XnumRow > 0);

  /*
  solver_lp_->numCol_ = XnumCol;
  solver_lp_->numRow_ = XnumRow;
  solver_lp_->sense_ = Xsense;
  int numNz = XnumNz;
  solver_lp_->colCost_.assign(&XcolCost[0], &XcolCost[0] + solver_lp_->numCol_);
  solver_lp_->colLower_.assign(&XcolLower[0], &XcolLower[0] + solver_lp_->numCol_);
  solver_lp_->colUpper_.resize(XnumCol);
  printf("XnumCol = %d\n", XnumCol);
  for (int iCol=0;iCol<XnumCol;iCol++) {
    printf("ColUpper[%2d]=%g\n",iCol,XcolUpper[iCol]);
  }
  solver_lp_->colUpper_.assign(&XcolUpper[0], &XcolUpper[0] + solver_lp_->numCol_);
  solver_lp_->rowLower_.assign(&XrowLower[0], &XrowLower[0] + solver_lp_->numRow_);
  solver_lp_->rowUpper_.assign(&XrowUpper[0], &XrowUpper[0] + solver_lp_->numRow_);
  solver_lp_->Astart_.assign(&XAstart[0], &XAstart[0] + solver_lp_->numCol_ + 1);
  solver_lp_->Aindex_.assign(&XAindex[0], &XAindex[0] + numNz);
  solver_lp_->Avalue_.assign(&XAvalue[0], &XAvalue[0] + numNz);
  // Assign and initialise the scaling factors

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  */
}

void HModel::replaceWithLogicalBasis() {
  // Replace basis with a logical basis then populate (where possible)
  // work* arrays
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = solver_lp_->numCol_ + row;
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    basis_->basicIndex_[row] = var;
  }
  for (int col = 0; col < solver_lp_->numCol_; col++) {
    basis_->nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
  }
  simplex_info_->num_basic_logicals = solver_lp_->numRow_;

  printf("Call simplex_method_.populate_work_arrays();\n");

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::replaceWithNewBasis(const int *XbasicIndex) {
  // Replace basis with a new basis then populate (where possible)
  // work* arrays

  //  printf("replaceWithNewBasis: \n");
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; var++) {
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
  }
  simplex_info_->num_basic_logicals = 0;
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = XbasicIndex[row];
    if (var >= solver_lp_->numCol_) simplex_info_->num_basic_logicals++;
    basis_->basicIndex_[row] = var;
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
  }

  printf("Call simplex_method_.populate_work_arrays();\n");

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::initFromNonbasic() {
  // Initialise basicIndex from nonbasic* then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  initBasicIndex();
  allocate_WorkAndBaseArrays();
  printf("Call simplex_method_.populate_work_arrays();\n");

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::replaceFromNonbasic() {
  // Initialise basicIndex using nonbasic* then populate (where possible)
  // work* arrays
  initBasicIndex();
  printf("Call simplex_method_.populate_work_arrays();\n");

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::initWithLogicalBasis() {
  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays

  for (int row = 0; row < solver_lp_->numRow_; row++) basis_->basicIndex_[row] = solver_lp_->numCol_ + row;
  for (int col = 0; col < solver_lp_->numCol_; col++) basis_->nonbasicFlag_[col] = 1;
  simplex_info_->num_basic_logicals = solver_lp_->numRow_;

  allocate_WorkAndBaseArrays();
  printf("Call simplex_method_.populate_work_arrays();\n");

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::extendWithLogicalBasis(int firstcol, int lastcol, int firstrow,
                                    int lastrow) {
  // Add nonbasic structurals and basic slacks according to model bounds.
  //
  // NB Assumes that the basis data structures and work vectors on
  // entry are assigned for columns 0..firstcol-1 and rows
  // 0..firstrow-1 and that they constitute a valid basis. Thus they
  // correspond to "firstcol" number of columns and "firstrow" number
  // of rows. Also assumes that solver_lp_->numCol_ and solver_lp_->numRow_ have already been
  // updated to correspond to any additional columns and rows. This is
  // necessary so that generic methods can be used to assign model
  // data to arrays dimensioned 0..numTot
  //
  // Null intervals firstcol...lastcol and firstrow...lastrow are
  // permitted, but this is achieved by setting the "last" to be less
  // than "first" since the latter is used to indicate what's
  // currently in the data structure.

  assert(firstcol >= 0);
  assert(firstrow >= 0);

  // printf("Called extendWithLogicalBasis:\n   solver_lp_->numCol_ =   %d\n   firstcol =
  // %d\n   lastcol =  %d\n   solver_lp_->numRow_ =   %d\n   firstrow = %d\n   lastrow =
  // %d\n", solver_lp_->numCol_, firstcol, lastcol, solver_lp_->numRow_, firstrow, lastrow);
  // Determine the numbers of columns and rows to be added

  int numAddCol = max(lastcol - firstcol + 1, 0);
  int numAddRow = max(lastrow - firstrow + 1, 0);
  int numAddTot = numAddCol + numAddRow;
  if (numAddTot == 0) return;

  // Determine the numbers of columns and rows before and after this method

  int local_oldNumCol = firstcol;
  int local_oldNumRow = firstrow;
  int local_oldNumTot = local_oldNumCol + local_oldNumRow;

  int local_newNumCol = max(local_oldNumCol, lastcol + 1);
  int local_newNumRow = max(local_oldNumRow, lastrow + 1);
  int local_newNumTot = local_newNumCol + local_newNumRow;

  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
#ifdef SCIPDEV
  printf("extendWithLogicalBasis\n");
  printf("solver_lp_->numCol_/Row/Tot = %d/%d/%d\n", solver_lp_->numCol_, solver_lp_->numRow_, numTot);
  printf("local_newNumCol/Row/Tot = %d/%d/%d\n", local_newNumCol,
         local_newNumRow, local_newNumTot);
  cout << flush;
#endif
  // ToDo: Replace references to local_newNum* by references to num* from here
  // on
  assert(local_newNumCol == solver_lp_->numCol_);
  assert(local_newNumRow == solver_lp_->numRow_);
  assert(local_newNumTot == numTot);

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  bool basisOK = nonbasicFlagBasicIndex_OK(local_oldNumCol, local_oldNumRow);
  if (!basisOK)
    printf("HModel::extendWithLogicalBasis: basisOK = %d\n", basisOK);
  assert(basisOK);
#endif

  //  Resize if necessary

  if (solver_lp_->numRow_ > local_oldNumRow) {
    basis_->basicIndex_.resize(solver_lp_->numRow_);

    simplex_info_->baseLower_.resize(solver_lp_->numRow_);
    simplex_info_->baseUpper_.resize(solver_lp_->numRow_);
    simplex_info_->baseValue_.resize(solver_lp_->numRow_);
  }
  if (numTot > local_oldNumTot) {
    basis_->nonbasicFlag_.resize(numTot);
    basis_->nonbasicMove_.resize(numTot);

    simplex_info_->workCost_.resize(numTot);
    simplex_info_->workDual_.resize(numTot);
    simplex_info_->workShift_.resize(numTot);

    simplex_info_->workLower_.resize(numTot);
    simplex_info_->workUpper_.resize(numTot);
    simplex_info_->workRange_.resize(numTot);
    simplex_info_->workValue_.resize(numTot);
  }

  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary

  int rowShift = solver_lp_->numCol_ - local_oldNumCol;
  if (rowShift > 0) {
    // printf("Shifting row data by %d using row=%d..0\n", rowShift,
    // local_oldNumRow-1);cout << flush;
    for (int row = local_oldNumRow - 1; row >= 0; row--) {
      basis_->basicIndex_[row] += rowShift;
      basis_->nonbasicFlag_[solver_lp_->numCol_ + row] = basis_->nonbasicFlag_[local_oldNumCol + row];
      basis_->nonbasicMove_[solver_lp_->numCol_ + row] = basis_->nonbasicMove_[local_oldNumCol + row];

      simplex_info_->workCost_[solver_lp_->numCol_ + row] = simplex_info_->workCost_[local_oldNumCol + row];
      simplex_info_->workDual_[solver_lp_->numCol_ + row] = simplex_info_->workDual_[local_oldNumCol + row];
      simplex_info_->workShift_[solver_lp_->numCol_ + row] = simplex_info_->workShift_[local_oldNumCol + row];

      simplex_info_->workLower_[solver_lp_->numCol_ + row] = simplex_info_->workLower_[local_oldNumCol + row];
      simplex_info_->workUpper_[solver_lp_->numCol_ + row] = simplex_info_->workUpper_[local_oldNumCol + row];
      simplex_info_->workRange_[solver_lp_->numCol_ + row] = simplex_info_->workRange_[local_oldNumCol + row];
      simplex_info_->workValue_[solver_lp_->numCol_ + row] = simplex_info_->workValue_[local_oldNumCol + row];

      // printf("Setting basicIndex[%2d] = %2d; basis_->nonbasicFlag_[%2d] = %2d;
      // basis_->nonbasicMove_[%2d] = %2d\n",
      //      row, basicIndex[row],
      //      solver_lp_->numCol_+row, basis_->nonbasicFlag_[local_oldNumCol+row],
      //      solver_lp_->numCol_+row, basis_->nonbasicMove_[local_oldNumCol+row]);cout << flush;
    }
  }
  // rp_basis();
  // printf("After possibly shifting row data\n");
  // Make any new columns nonbasic
  //  printf("Make any new cols nonbasic: %d %d %d\n", solver_lp_->numCol_, firstcol,
  //  lastcol);
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    //    printf("Setting basis_->nonbasicFlag_[%2d] = NONBASIC_FLAG_TRUE; Setting
    //    basis_->nonbasicMove_[%2d] = %2d\n", var, var, get_nonbasicMoveCol(var));
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    //    printf("Calling get_nonbasicMoveCol(%2d)\n", var);
    //    basis_->nonbasicMove_[var] = get_nonbasicMoveCol(var);
  }
  // Make any new rows basic
  //  printf("Make any new rows basic: %d %d %d\n", solver_lp_->numRow_, firstrow, lastrow);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = solver_lp_->numCol_ + row;
    //    printf("Setting basis_->nonbasicFlag_[%2d] = NONBASIC_FLAG_FALSE; Setting
    //    basicIndex[%2d] = %2d\n", var, row, var);
    basis_->nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    basis_->basicIndex_[row] = var;
  }

  // Initialise costs for the new columns and rows
  printf("init_Phase2_col_cost(firstcol, lastcol);\n");
  printf("init_Phase2_row_cost(firstrow, lastrow);\n");

  // Initialise bounds for the new columns and rows
  printf("init_Phase2_col_bound(firstcol, lastcol);\n");
  printf("init_Phase2_row_bound(firstrow, lastrow);\n");

  // Initialise values (and nonbasicMove) for the new columns
  printf("Call init_value_from_nonbasic(firstcol, lastcol);\n");

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  basisOK = nonbasicFlagBasicIndex_OK(solver_lp_->numCol_, solver_lp_->numRow_);
  assert(basisOK);
#endif

  simplex_info_->num_basic_logicals += numAddRow;

  //  rp_basis();

  // Deduce the consequences of adding new columns and/or rows
  //  if (numAddCol) update_solver_lp_status_flags(highs_model, LpAction::NEW_COLS);
  //  if (numAddRow) update_solver_lp_status_flags(highs_model, LpAction::NEW_ROWS);
}

bool HModel::nonbasicFlagBasicIndex_OK(int XnumCol, int XnumRow) {
  assert(XnumCol >= 0);
  assert(XnumRow >= 0);
  //  printf("Called nonbasicFlagBasicIndex_OK(%d, %d)\n", XnumCol,
  //  XnumRow);cout << flush;
  int XnumTot = XnumCol + XnumRow;
  int numBasic = 0;
  if (XnumTot > 0) {
    for (int var = 0; var < XnumTot; var++)
      if (!basis_->nonbasicFlag_[var]) numBasic++;
  }
  assert(numBasic == XnumRow);
  if (numBasic != XnumRow) return false;
  if (XnumRow > 0) {
    for (int row = 0; row < XnumRow; row++) {
      assert(!basis_->nonbasicFlag_[basis_->basicIndex_[row]]);
      if (basis_->nonbasicFlag_[basis_->basicIndex_[row]]) return false;
    }
  }
  return true;
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

void HModel::initBasicIndex() {
  int numBasic = 0;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; var++) {
    if (!basis_->nonbasicFlag_[var]) {
      assert(numBasic < solver_lp_->numRow_);
      basis_->basicIndex_[numBasic] = var;
      numBasic++;
    }
  }
  assert(numBasic = solver_lp_->numRow_ - 1);
}

void HModel::allocate_WorkAndBaseArrays() {
  // Allocate bounds and solution spaces
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  simplex_info_->workCost_.resize(numTot);
  simplex_info_->workDual_.resize(numTot);
  simplex_info_->workShift_.resize(numTot);

  simplex_info_->workLower_.resize(numTot);
  simplex_info_->workUpper_.resize(numTot);
  simplex_info_->workRange_.resize(numTot);
  simplex_info_->workValue_.resize(numTot);

  simplex_info_->baseLower_.resize(solver_lp_->numRow_);
  simplex_info_->baseUpper_.resize(solver_lp_->numRow_);
  simplex_info_->baseValue_.resize(solver_lp_->numRow_);
}

// ???? Housekeeping done from here down ????
// For the solver: methods to call INVERT and form dual and primal activities
// Call INVERT

#ifdef HiGHSDEV
void HModel::changeUpdate(int updateMethod) { factor_->change(updateMethod); }
#endif

#ifdef HiGHSDEV
// Checking methods Check loading of a model from arrays of data -
// just loads using arrays from an MPS read so optimality is check
void HModel::check_load_fromArrays() {
  // Use the arrays read from an MPS file to test the routine to
  // read a model passed by arrays. First copy the data.
  int XnumCol = solver_lp_->numCol_;
  int XnumRow = solver_lp_->numRow_;
  int XnumNz = solver_lp_->Astart_[solver_lp_->numCol_];
  int Xsense = solver_lp_->sense_;
  vector<double> XcolCost;
  vector<double> colLower;
  vector<double> XcolUpper;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;

  XcolCost.assign(&solver_lp_->colCost_[0], &solver_lp_->colCost_[0] + XnumCol);
  colLower.assign(&solver_lp_->colLower_[0], &solver_lp_->colLower_[0] + XnumCol);
  XcolUpper.assign(&solver_lp_->colUpper_[0], &solver_lp_->colUpper_[0] + XnumCol);
  XrowLower.assign(&solver_lp_->rowLower_[0], &solver_lp_->rowLower_[0] + XnumRow);
  XrowUpper.assign(&solver_lp_->rowUpper_[0], &solver_lp_->rowUpper_[0] + XnumRow);
  XAstart.assign(&solver_lp_->Astart_[0], &solver_lp_->Astart_[0] + XnumCol + 1);
  XAindex.assign(&solver_lp_->Aindex_[0], &solver_lp_->Aindex_[0] + XnumNz);
  XAvalue.assign(&solver_lp_->Avalue_[0], &solver_lp_->Avalue_[0] + XnumNz);

  //  clear_solver_lp(highs_model_object);
  load_fromArrays(XnumCol, Xsense, &XcolCost[0], &colLower[0],
                  &XcolUpper[0], XnumRow, &XrowLower[0], &XrowUpper[0], XnumNz,
                  &XAstart[0], &XAindex[0], &XAvalue[0]);
}

// Check that what's loaded from postsolve is correct
void HModel::check_load_fromPostsolve() {
  //  printf("Checking load_fromPostsolve\n");
  bool ok;
  ok = nonbasicFlagBasicIndex_OK(solver_lp_->numCol_, solver_lp_->numRow_);
  printf(
      "HModel::check_load_fromPostsolve: return from nonbasicFlagBasicIndex_OK "
      "= %d\n",
      ok);
  assert(ok);
  ok = allNonbasicMoveVsWorkArrays_OK();
  printf(
      "HModel::check_load_fromPostsolve: return from "
      "allNonbasicMoveVsWorkArrays_OK = %d\n",
      ok);
  assert(ok);
}
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

// Utilities to add, extract and delete columns and rows
// Add a contiguous set of columns to the model data---making them nonbasic
void HModel::util_addCols(int ncols, const double *XcolCost,
                          const double *colLower, const double *XcolUpper,
                          int nnonz, const int *XAstart, const int *XAindex,
                          const double *XAvalue) {
  assert(ncols >= 0);
  assert(nnonz >= 0);
  // ToDo How to check that solver_lp_->Astart_[solver_lp_->numCol_] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
  cout << flush;
#endif

  if (ncols == 0) return;

  int nwNumCol = solver_lp_->numCol_ + ncols;
  solver_lp_->colCost_.resize(nwNumCol);
  solver_lp_->colLower_.resize(nwNumCol);
  solver_lp_->colUpper_.resize(nwNumCol);
  scale_->col_.resize(nwNumCol);
  solver_lp_->Astart_.resize(nwNumCol + 1);

  // Note that the new columns must have starts, even if they have no entries
  // (yet)
  for (int col = 0; col < ncols; col++) {
    solver_lp_->colCost_[solver_lp_->numCol_ + col] = XcolCost[col];
    solver_lp_->colLower_[solver_lp_->numCol_ + col] = colLower[col];
    solver_lp_->colUpper_[solver_lp_->numCol_ + col] = XcolUpper[col];
    scale_->col_[solver_lp_->numCol_ + col] = 1.0;
    //    printf("In HModel::util_addCols: column %d: setting
    //    solver_lp_->Astart_[solver_lp_->numCol_+col+1] = %d \n", col, solver_lp_->Astart_[solver_lp_->numCol_]); cout << flush;
    solver_lp_->Astart_[solver_lp_->numCol_ + col + 1] = solver_lp_->Astart_[solver_lp_->numCol_];
  }

  //  printf("In HModel::util_addCols: nnonz = %d; cuNnonz = %d\n", nnonz,
  //  solver_lp_->Astart_[solver_lp_->numCol_]); cout << flush;
  if (nnonz > 0) {
    // Determine the current number of nonzeros
    int cuNnonz = solver_lp_->Astart_[solver_lp_->numCol_];

    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    // solver_lp_->Astart_.resize(nwNumCol+1);
    solver_lp_->Aindex_.resize(nwNnonz);
    solver_lp_->Avalue_.resize(nwNnonz);

    // Add the new columns
    for (int col = 0; col < ncols; col++) {
      //      printf("In HModel::util_addCols: column %d: setting
      //      solver_lp_->Astart_[solver_lp_->numCol_+col] = %d = %d + %d\n",
      //             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz); cout
      //             << flush;
      solver_lp_->Astart_[solver_lp_->numCol_ + col] = XAstart[col] + cuNnonz;
    }
    //    printf("In HModel::util_addCols: setting solver_lp_->Astart_[solver_lp_->numCol_+ncols] = %d\n",
    //    nwNnonz);
    cout << flush;
    solver_lp_->Astart_[solver_lp_->numCol_ + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++) {
      int row = XAindex[el];
      assert(row >= 0);
      assert(row < solver_lp_->numRow_);
      solver_lp_->Aindex_[cuNnonz + el] = row;
      solver_lp_->Avalue_[cuNnonz + el] = XAvalue[el];
    }
  }
  // Increase the number of columns and total number of variables in the model
  solver_lp_->numCol_ += ncols;
  //  numTot += ncols;

  //  printf("In HModel::util_addCols: Model now has solver_lp_->Astart_[%d] = %d
  //  nonzeros\n", solver_lp_->numCol_, solver_lp_->Astart_[solver_lp_->numCol_]);
  cout << flush;

  // Update the basis and work vectors correponding to new nonbasic columns
  extendWithLogicalBasis(solver_lp_->numCol_ - ncols, solver_lp_->numCol_ - 1, solver_lp_->numRow_, -1);
}

// Delete the model data for a contiguous set of columns
void HModel::util_deleteCols(int firstcol, int lastcol) {
  assert(firstcol >= 0);
  assert(lastcol < solver_lp_->numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_deleteCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Trivial cases are
  //
  // colStep = 0, in which case no columns are removed
  //
  // lastcol = solver_lp_->numCol_-1, in which case no columns need be
  // shifted. However, this implies solver_lp_->numCol_-colStep=firstcol, in which
  // case the loop is vacuous
  int colStep = lastcol - firstcol + 1;
  if (colStep) {
    for (int col = firstcol; col < solver_lp_->numCol_ - colStep; col++) {
      solver_lp_->colCost_[col] = solver_lp_->colCost_[col + colStep];
      solver_lp_->colLower_[col] = solver_lp_->colLower_[col + colStep];
      solver_lp_->colUpper_[col] = solver_lp_->colUpper_[col + colStep];
      scale_->col_[col] = scale_->col_[col + colStep];
    }
  }
  // Trivial cases are
  //
  // colstep = 0, in which case no columns are removed so elStep = 0
  //
  // lastcol = solver_lp_->numCol_-1, in which case no columns need be
  // shifted and the loops are vacuous
  if (colStep) {
    int elOs = solver_lp_->Astart_[firstcol];
    int elStep = solver_lp_->Astart_[lastcol + 1] - elOs;
    //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1,
    //    solver_lp_->Astart_[lastcol+1], solver_lp_->numCol_+1, solver_lp_->Astart_[solver_lp_->numCol_]-1);
    for (int el = solver_lp_->Astart_[lastcol + 1]; el < solver_lp_->Astart_[solver_lp_->numCol_]; el++) {
      //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n",
      //        el-elStep, solver_lp_->Aindex_[el-elStep], el, solver_lp_->Aindex_[el]);
      solver_lp_->Aindex_[el - elStep] = solver_lp_->Aindex_[el];
      solver_lp_->Avalue_[el - elStep] = solver_lp_->Avalue_[el];
    }
    for (int col = firstcol; col <= solver_lp_->numCol_ - colStep; col++) {
      //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col,
      //    solver_lp_->Astart_[col], col+colStep,  solver_lp_->Astart_[col+colStep]-elStep);
      solver_lp_->Astart_[col] = solver_lp_->Astart_[col + colStep] - elStep;
    }
  }

  // Reduce the number of columns and total number of variables in the model
  solver_lp_->numCol_ -= colStep;
  //  numTot -= colStep;

  // ToDo Determine consequences for basis when deleting columns
  // Invalidate matrix copies
  simplex_info_->solver_lp_has_matrix_col_wise = false;
  simplex_info_->solver_lp_has_matrix_row_wise = false;
}

// Delete the model data for a set of columns
void HModel::util_deleteColset(vector<int> &dstat) {
  printf("util_deleteColset is not implemented");
  assert(1 == 0);
}

// Extract the model data for a contiguous set of columns
void HModel::util_extractCols(int firstcol, int lastcol, double *colLower,
                              double *XcolUpper, int *nnonz, int *XAstart,
                              int *XAindex, double *XAvalue) {
  assert(firstcol >= 0);
  assert(lastcol < solver_lp_->numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Determine the number of columns to be extracted
  // int numExtractCols = lastcol-firstcol+1;
  // printf("Extracting %d columns\n", numExtractCols);cout << flush;
  int elOs = solver_lp_->Astart_[firstcol];
  for (int col = firstcol; col <= lastcol; col++) {
    //    printf("Extracting column %d\n", col);cout << flush;
    colLower[col - firstcol] = solver_lp_->colLower_[col];
    XcolUpper[col - firstcol] = solver_lp_->colUpper_[col];
    XAstart[col - firstcol] = solver_lp_->Astart_[col] - elOs;
  }
  for (int el = solver_lp_->Astart_[firstcol]; el < solver_lp_->Astart_[lastcol + 1]; el++) {
    XAindex[el - elOs] = solver_lp_->Aindex_[el];
    XAvalue[el - elOs] = solver_lp_->Avalue_[el];
  }
  *nnonz = solver_lp_->Astart_[lastcol + 1] - elOs;
}

// Add a contiguous set of rows to the model data---making them basic
void HModel::util_addRows(int nrows, const double *XrowLower,
                          const double *XrowUpper, int nnonz,
                          const int *XARstart, const int *XARindex,
                          const double *XARvalue) {
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || solver_lp_->numCol_ > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

  if (nrows == 0) return;

  int nwNumRow = solver_lp_->numRow_ + nrows;
  solver_lp_->rowLower_.resize(nwNumRow);
  solver_lp_->rowUpper_.resize(nwNumRow);
  scale_->row_.resize(nwNumRow);

  for (int row = 0; row < nrows; row++) {
    solver_lp_->rowLower_[solver_lp_->numRow_ + row] = XrowLower[row];
    solver_lp_->rowUpper_[solver_lp_->numRow_ + row] = XrowUpper[row];
    scale_->row_[solver_lp_->numRow_ + row] = 1.0;
  }
  // NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last
  // entry
  if (nnonz > 0) {
    int cuNnonz = solver_lp_->Astart_[solver_lp_->numCol_];
    vector<int> Alength;
    Alength.assign(solver_lp_->numCol_, 0);
    for (int el = 0; el < nnonz; el++) {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); cout <<
      //      flush;
      assert(col >= 0);
      assert(col < solver_lp_->numCol_);
      Alength[col]++;
    }
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    solver_lp_->Aindex_.resize(nwNnonz);
    solver_lp_->Avalue_.resize(nwNnonz);

    // Add the new rows
    // Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = solver_lp_->numCol_ - 1; col >= 0; col--) {
      // printf("Column %2d has additional length %2d\n", col,
      // Alength[col]);cout << flush;
      int Astart_Colp1 = nwEl;
      nwEl -= Alength[col];
      // printf("Shift: nwEl = %2d\n", nwEl); cout << flush;
      for (int el = solver_lp_->Astart_[col + 1] - 1; el >= solver_lp_->Astart_[col]; el--) {
        nwEl--;
        // printf("Shift: Over-writing solver_lp_->Aindex_[%2d] with solver_lp_->Aindex_[%2d]=%2d\n",
        // nwEl, el, solver_lp_->Aindex_[el]); cout << flush;
        solver_lp_->Aindex_[nwEl] = solver_lp_->Aindex_[el];
        solver_lp_->Avalue_[nwEl] = solver_lp_->Avalue_[el];
      }
      solver_lp_->Astart_[col + 1] = Astart_Colp1;
    }
    // printf("After shift: nwEl = %2d\n", nwEl); cout << flush;
    assert(nwEl == 0);
    // util_reportColMtx(solver_lp_->numCol_, solver_lp_->Astart_, solver_lp_->Aindex_, solver_lp_->Avalue_);

    // Insert the new entries
    for (int row = 0; row < nrows; row++) {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++) {
        int col = XARindex[el];
        nwEl = solver_lp_->Astart_[col + 1] - Alength[col];
        Alength[col]--;
        // printf("Insert: row = %2d; col = %2d; solver_lp_->Astart_[col+1]-Alength[col] =
        // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
        // solver_lp_->Astart_[col+1]-Alength[col], Alength[col], nwEl); cout << flush;
        assert(nwEl >= 0);
        assert(el >= 0);
        // printf("Insert: Over-writing solver_lp_->Aindex_[%2d] with solver_lp_->Aindex_[%2d]=%2d\n",
        // nwEl, el, solver_lp_->Aindex_[el]); cout << flush;
        solver_lp_->Aindex_[nwEl] = solver_lp_->numRow_ + row;
        solver_lp_->Avalue_[nwEl] = XARvalue[el];
      }
    }
  }
  // Increase the number of rows and total number of variables in the model
  solver_lp_->numRow_ += nrows;
  //  numTot += nrows;

  // Update the basis and work vectors correponding to new basic rows
  extendWithLogicalBasis(solver_lp_->numCol_, -1, solver_lp_->numRow_ - nrows, solver_lp_->numRow_ - 1);
}

// Delete the model data for a contiguous set of rows
void HModel::util_deleteRows(int firstrow, int lastrow) {
  assert(firstrow >= 0);
  assert(lastrow < solver_lp_->numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
  cout << flush;
#endif
  // Trivial cases are
  //
  // rowStep = 0, in which case no rows are removed
  //
  // lastrow = solver_lp_->numRow_-1, in which case no rows need be
  // shifted. However, this implies solver_lp_->numRow_-rowStep=firstrow, in which
  // case the loop is vacuous. However, they still have to be removed
  // from the matrix unless all rows are to be removed
  int rowStep = lastrow - firstrow + 1;
  bool allRows = rowStep == solver_lp_->numRow_;
#ifdef HiGHSDEV
  if (allRows) printf("In model.util_deleteRows, aa rows are being removed)\n");
#endif
  if (rowStep) {
    // Was: for (int row = firstrow; row < lastrow; row++) - surely wrong!
    for (int row = firstrow; row < solver_lp_->numRow_ - rowStep; row++) {
      solver_lp_->rowLower_[row] = solver_lp_->rowLower_[row + rowStep];
      solver_lp_->rowUpper_[row] = solver_lp_->rowUpper_[row + rowStep];
      //    scale_->row_[row] = scale_->row_[row + rowStep];
    }
    if (!allRows) {
      int nnz = 0;
      for (int col = 0; col < solver_lp_->numCol_; col++) {
        int fmEl = solver_lp_->Astart_[col];
        solver_lp_->Astart_[col] = nnz;
        for (int el = fmEl; el < solver_lp_->Astart_[col + 1]; el++) {
          int row = solver_lp_->Aindex_[el];
          if (row < firstrow || row > lastrow) {
            if (row < firstrow) {
              solver_lp_->Aindex_[nnz] = row;
            } else {
              solver_lp_->Aindex_[nnz] = row - rowStep;
            }
            solver_lp_->Avalue_[nnz] = solver_lp_->Avalue_[el];
            nnz++;
          }
        }
      }
      solver_lp_->Astart_[solver_lp_->numCol_] = nnz;
    }
  }

  // Reduce the number of rows and total number of variables in the model
  solver_lp_->numRow_ -= rowStep;
  //  numTot -= rowStep;

  // Determine consequences for basis when deleting rows
  //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS);
}

// Delete the model data for a set of rows
void HModel::util_deleteRowset(int *dstat) {
  bool rp = false;
  if (rp) {
    printf("Called model.util_deleteRowSet\n");
    cout << flush;
    printf("Before\n");
  }
  //  solver_lp_->reportLp();

  int newRow = 0;
  // Look through the rows removing any being deleted and shifting data
  // for the rest
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = solver_lp_->numCol_ + row;
      int newVar = solver_lp_->numCol_ + newRow;
      dstat[row] = newRow;
      solver_lp_->rowLower_[newRow] = solver_lp_->rowLower_[row];
      solver_lp_->rowUpper_[newRow] = solver_lp_->rowUpper_[row];
      //    scale_->row_[row] = scale_->row_[rowStep+row];
      basis_->nonbasicFlag_[newVar] = basis_->nonbasicFlag_[var];
      basis_->nonbasicMove_[newVar] = basis_->nonbasicMove_[var];
      simplex_info_->workCost_[newVar] = simplex_info_->workCost_[var];
      simplex_info_->workShift_[newVar] = simplex_info_->workShift_[var];
      simplex_info_->workLower_[newVar] = simplex_info_->workLower_[var];
      simplex_info_->workUpper_[newVar] = simplex_info_->workUpper_[var];
      simplex_info_->workRange_[newVar] = simplex_info_->workRange_[var];
      simplex_info_->workValue_[newVar] = simplex_info_->workValue_[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, solver_lp_->rowLower_[newRow], solver_lp_->rowUpper_[newRow],
            basis_->nonbasicFlag_[newVar], basis_->nonbasicMove_[newVar]);
      newRow++;
    } else {
      // Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row,
               dstat[row], solver_lp_->numCol_ + row);
    }
  }

  if (rp) {
    printf("After\n");
    for (int row = 0; row < solver_lp_->numRow_; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  // Look through the column-wise matrix, removing entries
  // corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < solver_lp_->numCol_; col++) {
    int fmEl = solver_lp_->Astart_[col];
    solver_lp_->Astart_[col] = nnz;
    for (int el = fmEl; el < solver_lp_->Astart_[col + 1]; el++) {
      int row = solver_lp_->Aindex_[el];
      if (dstat[row] >= 0) {
        solver_lp_->Aindex_[nnz] = dstat[row];
        solver_lp_->Avalue_[nnz] = solver_lp_->Avalue_[el];
        nnz++;
      }
    }
  }
  solver_lp_->Astart_[solver_lp_->numCol_] = nnz;

  // Reduce the number of rows and total number of variables in the model
  int dlNumRow = solver_lp_->numRow_ - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", solver_lp_->numRow_, dlNumRow,
           newRow);
#endif
  solver_lp_->numRow_ -= dlNumRow;
  //  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; var++) {
    if (!basis_->nonbasicFlag_[var]) {
      basis_->basicIndex_[numBasic] = var;
      numBasic++;
      if (numBasic > newRow) {
        basisOK = false;
        break;
      }
    }
  }

  if (rp) {
    printf("Now have %d cols; %d rows and %d total\n", solver_lp_->numCol_, solver_lp_->numRow_, numTot);
    for (int row = 0; row < solver_lp_->numRow_; row++)
      printf("Basic variable in row %2d is %2d\n", row, basis_->basicIndex_[row]);
    for (int col = 0; col < solver_lp_->numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, basis_->nonbasicFlag_[col]);
    for (int row = 0; row < solver_lp_->numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             solver_lp_->numCol_ + row, basis_->nonbasicFlag_[solver_lp_->numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(solver_lp_->numCol_, solver_lp_->numRow_);
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so
    //    basisOK\n"); cout<<flush;
#endif
    // Determine consequences for basis when deleting rows to leave an OK basis
    //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS_BASIS_OK);
  } else {
    assert(basisOK);
#ifdef SCIP_DEV
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
#endif
    // Determine consequences for basis when deleting rows to leave no basis
  //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS);
  }
}

// Extract the model data for a contiguous set of rows
void HModel::util_extractRows(int firstrow, int lastrow, double *XrowLower,
                              double *XrowUpper, int *nnonz, int *XARstart,
                              int *XARindex, double *XARvalue) {
  assert(firstrow >= 0);
  assert(lastrow < solver_lp_->numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
  cout << flush;
#endif
  // Determine the number of rows to be extracted
  int numExtractRows = lastrow - firstrow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);cout << flush;
  for (int row = firstrow; row <= lastrow; row++) {
    // printf("Extracting row %d\n", row);cout << flush;
    XrowLower[row - firstrow] = solver_lp_->rowLower_[row];
    XrowUpper[row - firstrow] = solver_lp_->rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow],
    // XrowUpper[row-firstrow]);cout << flush;
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = solver_lp_->Astart_[0]; el < solver_lp_->Astart_[solver_lp_->numCol_]; el++) {
    int row = solver_lp_->Aindex_[el];
    if (row >= firstrow && row <= lastrow) XARlength[row - firstrow] += 1;
  }
  XARstart[0] = 0;
  // printf("Start of row %2d is %d\n", 0, XARstart[0]);cout << flush;
  // printf("Length of row %2d is %d\n", 0, XARlength[0]);cout << flush;
  for (int row = 0; row < numExtractRows - 1; row++) {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
    // printf("Start of row %2d is %d\n", row+1, XARstart[row+1]);cout << flush;
    // printf("Length of row %2d is %d\n", row+1, XARlength[row+1]);cout <<
    // flush;
  }
  XARlength[numExtractRows - 1] = 0;

  for (int col = 0; col < solver_lp_->numCol_; col++) {
    for (int el = solver_lp_->Astart_[col]; el < solver_lp_->Astart_[col + 1]; el++) {
      int row = solver_lp_->Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);cout <<
      // flush;
      if (row >= firstrow && row <= lastrow) {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, solver_lp_->Avalue_[el]);cout << flush;
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = solver_lp_->Avalue_[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);cout << flush;
}

// Change a single coefficient in the matrix
void HModel::util_changeCoeff(int row, int col, const double newval) {
  assert(row >= 0 && row < solver_lp_->numRow_);
  assert(col >= 0 && col < solver_lp_->numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col,
         newval);
  cout << flush;
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  row, col, newval);cout << flush;

  //  solver_lp_->reportLp();
  int cg_el = -1;
  for (int el = solver_lp_->Astart_[col]; el < solver_lp_->Astart_[col + 1]; el++) {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //    solver_lp_->Aindex_[el], row);
    if (solver_lp_->Aindex_[el] == row) {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0) {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n",
    //    row, col);
    cg_el = solver_lp_->Astart_[col + 1];
    int nwNnonz = solver_lp_->Astart_[solver_lp_->numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    solver_lp_->Astart_[solver_lp_->numCol_], nwNnonz);
    solver_lp_->Aindex_.resize(nwNnonz);
    solver_lp_->Avalue_.resize(nwNnonz);
    for (int i = col + 1; i <= solver_lp_->numCol_; i++) solver_lp_->Astart_[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--) {
      solver_lp_->Aindex_[el] = solver_lp_->Aindex_[el - 1];
      solver_lp_->Avalue_[el] = solver_lp_->Avalue_[el - 1];
    }
  }
  solver_lp_->Avalue_[cg_el] = newval;

  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  //  update_solver_lp_status_flags(highs_model, LpAction::NEW_ROWS);
  //  solver_lp_->reportLp();
}

// Get a single coefficient from the matrix
void HModel::util_getCoeff(HighsLp lp, int row, int col, double *val) {
  assert(row >= 0 && row < lp.numRow_);
  assert(col >= 0 && col < lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);
  cout << flush;
#endif
  //  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);cout <<
  //  flush;

  cout << val << endl;

  int get_el = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    //  printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //  lp.Aindex_[el], row);cout << flush;
    if (lp.Aindex_[el] == row) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    //  printf("model.util_getCoeff: Cannot find row %d in column %d\n", row,
    //  col);cout << flush;
    *val = 0;
  } else {
    //  printf("model.util_getCoeff: Found row %d in column %d as element %d:
    //  value %g\n", row, col, get_el, lp.Avalue_[get_el]);cout << flush;
    *val = lp.Avalue_[get_el];
  }
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
