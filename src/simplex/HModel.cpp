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

  populate_WorkArrays();

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

  populate_WorkArrays();

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::initFromNonbasic() {
  // Initialise basicIndex from nonbasic* then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  initBasicIndex();
  allocate_WorkAndBaseArrays();
  populate_WorkArrays();

  // Deduce the consequences of a new basis
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
}

void HModel::replaceFromNonbasic() {
  // Initialise basicIndex using nonbasic* then populate (where possible)
  // work* arrays
  initBasicIndex();
  populate_WorkArrays();

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
  populate_WorkArrays();

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
  initPh2ColCost(firstcol, lastcol);
  initPh2RowCost(firstrow, lastrow);

  // Initialise bounds for the new columns and rows
  initPh2ColBound(firstcol, lastcol);
  initPh2RowBound(firstrow, lastrow);

  // Initialise values (and nonbasicMove) for the new columns
  initValueFromNonbasic(firstcol, lastcol);

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

void HModel::setup_for_solve() {
  //  timer.reset();
  if (solver_lp_->numRow_ == 0) return;

  //  report_solver_lp_status_flags(highs_mode_object);cout<<flush;
  //  printf("In setup_for_solve: basis_->valid_ = %d \n", basis_->valid_);cout<<flush;
  if (basis_->valid_) {
    // Model has a basis so just count the number of basic logicals
    setup_numBasicLogicals();
  } else {
    // Model has no basis: set up a logical basis then populate (where
    // possible) work* arrays
    replaceWithLogicalBasis();
    //    printf("Called replaceWithLogicalBasis\n");cout<<flush;
  }

  if (!(simplex_info_->solver_lp_has_matrix_col_wise && simplex_info_->solver_lp_has_matrix_row_wise)) {
    // Make a copy of col-wise matrix for HMatrix and create its row-wise matrix
    if (simplex_info_->num_basic_logicals == solver_lp_->numRow_) {
      matrix_->setup_lgBs(solver_lp_->numCol_, solver_lp_->numRow_, &solver_lp_->Astart_[0], &solver_lp_->Aindex_[0], &solver_lp_->Avalue_[0]);
      //      printf("Called matrix_->setup_lgBs\n");cout<<flush;
    } else {
      matrix_->setup(solver_lp_->numCol_, solver_lp_->numRow_, &solver_lp_->Astart_[0], &solver_lp_->Aindex_[0], &solver_lp_->Avalue_[0],
                   &basis_->nonbasicFlag_[0]);
      //      printf("Called matrix_->setup\n");cout<<flush;
    }
    // Indicate that there is a colum-wise and row-wise copy of the
    // matrix: can't be done in matrix_->setup_lgBs
    //    simplex_info_->solver_lp_has_matrix_col_wise = true;
    //    simplex_info_->solver_lp_has_matrix_row_wise = true;
  }

    // TODO Put something in to skip factor_->setup
    // Initialise factor arrays, passing &basis_->basicIndex_[0] so that its
    // address can be copied to the internal Factor pointer
    factor_->setup(solver_lp_->numCol_, solver_lp_->numRow_, &solver_lp_->Astart_[0], &solver_lp_->Aindex_[0], &solver_lp_->Avalue_[0],
                 &basis_->basicIndex_[0]);
    // Indicate that the model has factor arrays: can't be done in factor.setup
    //simplex_info_->solver_lp_has_factor_arrays = true;
}

bool HModel::OKtoSolve(int level, int phase) {
  //  printf("Called OKtoSolve(%1d, %1d)\n", level, phase);
  bool ok;
  // Level 0: Minimal check - just look at flags. This means we trust them!
  ok =
    basis_->valid_ &&
    simplex_info_->solver_lp_has_matrix_col_wise &&
    simplex_info_->solver_lp_has_matrix_row_wise &&
    //    simplex_info_->solver_lp_has_factor_arrays &&
    simplex_info_->solver_lp_has_dual_steepest_edge_weights &&
    simplex_info_->solver_lp_has_invert;
  // TODO: Eliminate the following line ASAP!!!
  ok = true;
  if (!ok) {
    if (!basis_->valid_)
      printf("Not OK to solve since basis_->valid_ = %d\n",
	     basis_->valid_);
    if (!simplex_info_->solver_lp_has_matrix_col_wise)
      printf("Not OK to solve since simplex_info_->solver_lp_has_matrix_col_wise = %d\n",
             simplex_info_->solver_lp_has_matrix_col_wise);
    if (!simplex_info_->solver_lp_has_matrix_row_wise)
      printf("Not OK to solve since simplex_info_->solver_lp_has_matrix_row_wise = %d\n",
             simplex_info_->solver_lp_has_matrix_row_wise);
    //    if (!simplex_info_->solver_lp_has_factor_arrays)
    //      printf("Not OK to solve since simplex_info_->solver_lp_has_factor_arrays = %d\n",
    //             simplex_info_->solver_lp_has_factor_arrays);
    if (!simplex_info_->solver_lp_has_dual_steepest_edge_weights)
      printf("Not OK to solve since simplex_info_->solver_lp_has_dual_steepest_edge_weights = %d\n",
	     simplex_info_->solver_lp_has_dual_steepest_edge_weights);
    if (!simplex_info_->solver_lp_has_invert)
      printf("Not OK to solve since simplex_info_->solver_lp_has_invert = %d\n",
	     simplex_info_->solver_lp_has_invert);
    cout << flush;
  }
#ifdef HiGHSDEV
  assert(ok);
#endif
  if (level <= 0) return ok;
  // Level 1: Basis and data check
  ok = nonbasicFlagBasicIndex_OK(solver_lp_->numCol_, solver_lp_->numRow_);
  if (!ok) {
    printf("Error in nonbasicFlag and basicIndex\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  ok = workArrays_OK(phase);
  if (!ok) {
    printf("Error in workArrays\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; ++var) {
    if (basis_->nonbasicFlag_[var]) {
      // Nonbasic variable
      ok = oneNonbasicMoveVsWorkArrays_OK(var);
      if (!ok) {
        printf("Error in nonbasicMoveVsWorkArrays for variable %d of %d\n", var,
               numTot);
        cout << flush;
#ifdef HiGHSDEV
        assert(ok);
#endif
        return ok;
      }
    }
  }
  if (level <= 1) return ok;
  printf("OKtoSolve(%1d) not implemented\n", level);
  cout << flush;
  return ok;
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

bool HModel::workArrays_OK(int phase) {
  //  printf("Called workArrays_OK(%d)\n", phase);cout << flush;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < solver_lp_->numCol_; ++col) {
      int var = col;
      if (!highs_isInfinity(-simplex_info_->workLower_[var])) {
        ok = simplex_info_->workLower_[var] == solver_lp_->colLower_[col];
        if (!ok) {
          printf("For col %d, simplex_info_->workLower_ should be %g but is %g\n", col,
                 solver_lp_->colLower_[col], simplex_info_->workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
        ok = simplex_info_->workUpper_[var] == solver_lp_->colUpper_[col];
        if (!ok) {
          printf("For col %d, simplex_info_->workUpper_ should be %g but is %g\n", col,
                 solver_lp_->colUpper_[col], simplex_info_->workUpper_[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < solver_lp_->numRow_; ++row) {
      int var = solver_lp_->numCol_ + row;
      if (!highs_isInfinity(-simplex_info_->workLower_[var])) {
        ok = simplex_info_->workLower_[var] == -solver_lp_->rowUpper_[row];
        if (!ok) {
          printf("For row %d, simplex_info_->workLower_ should be %g but is %g\n", row,
                 -solver_lp_->rowUpper_[row], simplex_info_->workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
        ok = simplex_info_->workUpper_[var] == -solver_lp_->rowLower_[row];
        if (!ok) {
          printf("For row %d, simplex_info_->workUpper_ should be %g but is %g\n", row,
                 -solver_lp_->rowLower_[row], simplex_info_->workUpper_[var]);
          return ok;
        }
      }
    }
  }
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; ++var) {
    ok = simplex_info_->workRange_[var] == (simplex_info_->workUpper_[var] - simplex_info_->workLower_[var]);
    if (!ok) {
      printf("For variable %d, simplex_info_->workRange_ should be %g = %g - %g but is %g\n",
             var, simplex_info_->workUpper_[var] - simplex_info_->workLower_[var], simplex_info_->workUpper_[var],
             simplex_info_->workLower_[var], simplex_info_->workRange_[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!simplex_info_->costs_perturbed) {
    for (int col = 0; col < solver_lp_->numCol_; ++col) {
      int var = col;
      ok = simplex_info_->workCost_[var] == solver_lp_->sense_ * solver_lp_->colCost_[col];
      if (!ok) {
        printf("For col %d, simplex_info_->workLower_ should be %g but is %g\n", col,
               solver_lp_->colLower_[col], simplex_info_->workCost_[var]);
        return ok;
      }
    }
    for (int row = 0; row < solver_lp_->numRow_; ++row) {
      int var = solver_lp_->numCol_ + row;
      ok = simplex_info_->workCost_[var] == 0.;
      if (!ok) {
        printf("For row %d, simplex_info_->workCost_ should be zero but is %g\n", row,
               simplex_info_->workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool HModel::allNonbasicMoveVsWorkArrays_OK() {
  bool ok;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int var = 0; var < numTot; ++var) {
    printf("NonbasicMoveVsWorkArrays: var = %2d; basis_->nonbasicFlag_[var] = %2d\n",
           var, basis_->nonbasicFlag_[var]);
    if (!basis_->nonbasicFlag_[var]) continue;
    ok = oneNonbasicMoveVsWorkArrays_OK(var);
    if (!ok) {
      printf("Error in NonbasicMoveVsWorkArrays for nonbasic variable %d\n",
             var);
      cout << flush;
#ifdef HiGHSDEV
      assert(ok);
#endif
      return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool HModel::oneNonbasicMoveVsWorkArrays_OK(int var) {
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  //  printf("Calling oneNonbasicMoveVsWorkArrays_ok with var = %2d; numTot =
  //  %2d\n Bounds [%11g, %11g] nonbasicMove = %d\n",
  //	 var, numTot, simplex_info_->workLower_[var], simplex_info_->workUpper_[var], basis_->nonbasicMove_[var]);
  // cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  // Make sure we're not checking a basic variable
  if (!basis_->nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-simplex_info_->workLower_[var])) {
    if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info_->workLower_[var] == simplex_info_->workUpper_[var]) {
        // Fixed variable
        ok = basis_->nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          printf(
              "Fixed variable %d (solver_lp_->numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
              "move should be zero but is %d\n",
              var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var], simplex_info_->workUpper_[var],
              basis_->nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info_->workValue_[var] == simplex_info_->workLower_[var];
        if (!ok) {
          printf(
              "Fixed variable %d (solver_lp_->numCol_ = %d) so simplex_info_->work value should be %g but "
              "is %g\n",
              var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          printf(
              "Boxed variable %d (solver_lp_->numCol_ = %d) [%11g, %11g, %11g] range %g so "
              "nonbasic move should be up/down but is  %d\n",
              var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var], simplex_info_->workUpper_[var],
              simplex_info_->workUpper_[var] - simplex_info_->workLower_[var], basis_->nonbasicMove_[var]);
          return ok;
        }
        if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info_->workValue_[var] == simplex_info_->workLower_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (solver_lp_->numCol_ = %d) with NONBASIC_MOVE_UP so work "
                "value should be %g but is %g\n",
                var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info_->workValue_[var] == simplex_info_->workUpper_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (solver_lp_->numCol_ = %d) with NONBASIC_MOVE_DN so work "
                "value should be %g but is %g\n",
                var, solver_lp_->numCol_, simplex_info_->workUpper_[var], simplex_info_->workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (solver_lp_->numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d\n",
            var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var], simplex_info_->workUpper_[var],
            NONBASIC_MOVE_UP, basis_->nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info_->workValue_[var] == simplex_info_->workLower_[var];
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (solver_lp_->numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
      ok = basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (solver_lp_->numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d\n",
            var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var], simplex_info_->workUpper_[var],
            basis_->nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info_->workValue_[var] == simplex_info_->workUpper_[var];
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (solver_lp_->numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, solver_lp_->numCol_, simplex_info_->workUpper_[var], simplex_info_->workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = basis_->nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        printf(
            "Free variable %d (solver_lp_->numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
            "move should be zero but is  %d\n",
            var, solver_lp_->numCol_, simplex_info_->workLower_[var], simplex_info_->workValue_[var], simplex_info_->workUpper_[var],
            basis_->nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info_->workValue_[var] == 0.0;
      if (!ok) {
        printf(
            "Free variable %d (solver_lp_->numCol_ = %d) so work value should be zero but "
            "is %g\n",
            var, solver_lp_->numCol_, simplex_info_->workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void HModel::setup_numBasicLogicals() {
  simplex_info_->num_basic_logicals = 0;
  for (int i = 0; i < solver_lp_->numRow_; i++)
    if (basis_->basicIndex_[i] >= solver_lp_->numCol_) simplex_info_->num_basic_logicals += 1;
  //  printf("Determined simplex_info_->num_basic_logicals = %d of %d\n", simplex_info_->num_basic_logicals,
  //  solver_lp_->numRow_);
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
  // Was workShift.assign(numTot, 0); but shift is populated by call to
  // initCost()
  simplex_info_->workShift_.resize(numTot);

  simplex_info_->workLower_.resize(numTot);
  simplex_info_->workUpper_.resize(numTot);
  simplex_info_->workRange_.resize(numTot);
  simplex_info_->workValue_.resize(numTot);

  simplex_info_->baseLower_.resize(solver_lp_->numRow_);
  simplex_info_->baseUpper_.resize(solver_lp_->numRow_);
  simplex_info_->baseValue_.resize(solver_lp_->numRow_);
}

void HModel::populate_WorkArrays() {
  // Initialize the values
  initCost();
  initBound();
  initValue();
}

void HModel::initCost(int perturb) {
  // Copy the cost
  initPh2ColCost(0, solver_lp_->numCol_ - 1);
  initPh2RowCost(0, solver_lp_->numRow_ - 1);
  // See if we want to skip perturbation
  simplex_info_->costs_perturbed = 0;
  if (perturb == 0 || simplex_info_->perturb_costs == 0) return;
  simplex_info_->costs_perturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < solver_lp_->numCol_; i++) bigc = max(bigc, fabs(simplex_info_->workCost_[i]));
  if (bigc > 100) bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++) boxedRate += (simplex_info_->workRange_[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) bigc = min(bigc, 1.0);
  if (bigc < 1) {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < solver_lp_->numCol_; i++) {
    double lower = solver_lp_->colLower_[i];
    double upper = solver_lp_->colUpper_[i];
    double xpert = (fabs(simplex_info_->workCost_[i]) + 1) * base * (1 + simplex_info_->numTotRandomValue_[i]);
    if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper == HIGHS_CONST_INF) {  // Lower
      simplex_info_->workCost_[i] += xpert;
    } else if (lower == -HIGHS_CONST_INF) {  // Upper
      simplex_info_->workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex_info_->workCost_[i] += (simplex_info_->workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }

  for (int i = solver_lp_->numCol_; i < numTot; i++) {
    simplex_info_->workCost_[i] += (0.5 - simplex_info_->numTotRandomValue_[i]) * 1e-12;
  }
}

void HModel::initBound(int phase) {
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initPh2ColBound(0, solver_lp_->numCol_ - 1);
  initPh2RowBound(0, solver_lp_->numRow_ - 1);
  if (phase == 2) return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HIGHS_CONST_INF;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++) {
    if (simplex_info_->workLower_[i] == -inf && simplex_info_->workUpper_[i] == inf) {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= solver_lp_->numCol_) continue;
      simplex_info_->workLower_[i] = -1000, simplex_info_->workUpper_[i] = 1000;  // FREE
    } else if (simplex_info_->workLower_[i] == -inf) {
      simplex_info_->workLower_[i] = -1, simplex_info_->workUpper_[i] = 0;  // UPPER
    } else if (simplex_info_->workUpper_[i] == inf) {
      simplex_info_->workLower_[i] = 0, simplex_info_->workUpper_[i] = 1;  // LOWER
    } else {
      simplex_info_->workLower_[i] = 0, simplex_info_->workUpper_[i] = 0;  // BOXED or FIXED
    }
    simplex_info_->workRange_[i] = simplex_info_->workUpper_[i] - simplex_info_->workLower_[i];
  }
}

void HModel::initValue() {
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  initValueFromNonbasic(0, numTot - 1);
}

void HModel::initPh2ColCost(int firstcol, int lastcol) {
  // Copy the Phase 2 cost and zero the shift
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    simplex_info_->workCost_[var] = solver_lp_->sense_ * solver_lp_->colCost_[col];
    simplex_info_->workShift_[var] = 0.;
  }
}

void HModel::initPh2RowCost(int firstrow, int lastrow) {
  // Zero the cost and shift
  for (int row = firstrow; row <= lastrow; row++) {
    int var = solver_lp_->numCol_ + row;
    simplex_info_->workCost_[var] = 0;
    simplex_info_->workShift_[var] = 0.;
  }
}

void HModel::initPh2ColBound(int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  assert(firstcol >= 0);
  assert(lastcol < solver_lp_->numCol_);
  for (int col = firstcol; col <= lastcol; col++) {
    simplex_info_->workLower_[col] = solver_lp_->colLower_[col];
    simplex_info_->workUpper_[col] = solver_lp_->colUpper_[col];
    simplex_info_->workRange_[col] = simplex_info_->workUpper_[col] - simplex_info_->workLower_[col];
  }
}

void HModel::initPh2RowBound(int firstrow, int lastrow) {
  // Copy bounds and compute ranges
  assert(firstrow >= 0);
  assert(lastrow < solver_lp_->numRow_);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = solver_lp_->numCol_ + row;
    simplex_info_->workLower_[var] = -solver_lp_->rowUpper_[row];
    simplex_info_->workUpper_[var] = -solver_lp_->rowLower_[row];
    simplex_info_->workRange_[var] = simplex_info_->workUpper_[var] - simplex_info_->workLower_[var];
  }
}

void HModel::initValueFromNonbasic(int firstvar, int lastvar) {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  assert(firstvar >= 0);
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  assert(lastvar < numTot);
  // double dl_pr_act, norm_dl_pr_act;
  // norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++) {
    if (basis_->nonbasicFlag_[var]) {
      // Nonbasic variable
      // double prev_pr_act = simplex_info_->workValue_[var];
      if (simplex_info_->workLower_[var] == simplex_info_->workUpper_[var]) {
        // Fixed
        simplex_info_->workValue_[var] = simplex_info_->workLower_[var];
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      } else if (!highs_isInfinity(-simplex_info_->workLower_[var])) {
        // Finite lower bound so boxed or lower
        if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
          // Finite upper bound so boxed
          if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_UP) {
            // Set at lower
            simplex_info_->workValue_[var] = simplex_info_->workLower_[var];
          } else if (basis_->nonbasicMove_[var] == NONBASIC_MOVE_DN) {
            // Set at upper
            simplex_info_->workValue_[var] = simplex_info_->workUpper_[var];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
            simplex_info_->workValue_[var] = simplex_info_->workLower_[var];
          }
        } else {
          // Lower
          simplex_info_->workValue_[var] = simplex_info_->workLower_[var];
          basis_->nonbasicMove_[var] = NONBASIC_MOVE_UP;
        }
      } else if (!highs_isInfinity(simplex_info_->workUpper_[var])) {
        // Upper
        simplex_info_->workValue_[var] = simplex_info_->workUpper_[var];
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        simplex_info_->workValue_[var] = 0;
        basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      }
      // dl_pr_act = simplex_info_->workValue_[var] - prev_pr_act;
      // norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
      //      %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var, simplex_info_->workLower_[var],
      // simplex_info_->workValue_[var], simplex_info_->workUpper_[var], simplex_info_->workDual_[var], dl_pr_act);
    } else {
      // Basic variable
      basis_->nonbasicMove_[var] = NONBASIC_MOVE_ZE;
    }
  }
  //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
  //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is
  //  %g\n", norm_dl_pr_act);
}

// ???? Housekeeping done from here down ????
// For the solver: methods to call INVERT and form dual and primal activities
// Call INVERT
int HModel::computeFactor() {
#ifdef HiGHSDEV
  double tt0 = 0;
  int iClock = simplex_info_->clock_[InvertClock];
  if (simplex_info_->analyse_invert_time) tt0 = timer_->clock_time[iClock];
#endif
  // TODO Understand why handling noPvC and noPvR in what seem to be
  // different ways ends up equivalent.
  int rankDeficiency = factor_->build();
  if (rankDeficiency) {
    handleRankDeficiency();
    //    simplex_info_->solution_status = SimplexSolutionStatus::SINGULAR;
#ifdef HiGHSDEV
    //    writePivots("failed");
#endif
    //      return rankDeficiency;
  }
  //    printf("INVERT: After %d iterations and %d updates\n", simplex_info_->iteration_count,
  //    simplex_info_->update_count);
  simplex_info_->update_count = 0;

#ifdef HiGHSDEV
  if (simplex_info_->analyse_invert_time) {
    int iClock = simplex_info_->clock_[InvertClock];
    simplex_info_->total_inverts = timer_->clock_num_call[iClock];
    simplex_info_->total_invert_time = timer_->clock_time[iClock];
    double invertTime = simplex_info_->total_invert_time - tt0;
    printf(
        "           INVERT  %4d     on iteration %9d: INVERT  time = %11.4g; "
        "Total INVERT  time = %11.4g\n",
        simplex_info_->total_inverts,
	simplex_info_->iteration_count, invertTime, simplex_info_->total_invert_time);
  }
#endif

  // Now have a representation of B^{-1}, and it is fresh!
  simplex_info_->solver_lp_has_invert = true;
  simplex_info_->solver_lp_has_fresh_invert = true;
  return 0;
}

// Compute the dual activities
void HModel::computeDual() {
  bool an_computeDual_norm2 = false;
  double btranRHS_norm2;
  double btranSol_norm2;
  double workDual_norm2;

  // Create a local buffer for the pi vector
  HVector buffer;
  buffer.setup(solver_lp_->numRow_);
  buffer.clear(); 
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++) {
    buffer.index[iRow] = iRow;
    buffer.array[iRow] =
        simplex_info_->workCost_[basis_->basicIndex_[iRow]] + simplex_info_->workShift_[basis_->basicIndex_[iRow]];
  }
  buffer.count = solver_lp_->numRow_;
  if (an_computeDual_norm2) {
    btranRHS_norm2 = buffer.norm2();
    btranRHS_norm2 = sqrt(btranRHS_norm2);
  }
  //  printf("computeDual: Before BTRAN\n");cout<<flush;
  factor_->btran(buffer, 1);
  //  printf("computeDual: After  BTRAN\n");cout<<flush;
  if (an_computeDual_norm2) {
    btranSol_norm2 = buffer.norm2();
    btranSol_norm2 = sqrt(btranSol_norm2);
  }

  // Create a local buffer for the values of reduced costs
  HVector bufferLong;
  bufferLong.setup(solver_lp_->numCol_);
  bufferLong.clear();
  matrix_->price_by_col(bufferLong, buffer);
  for (int i = 0; i < solver_lp_->numCol_; i++) {
    simplex_info_->workDual_[i] = simplex_info_->workCost_[i] - bufferLong.array[i];
  }
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = solver_lp_->numCol_; i < numTot; i++) {
    simplex_info_->workDual_[i] = simplex_info_->workCost_[i] - buffer.array[i - solver_lp_->numCol_];
  }

  if (an_computeDual_norm2) {
    workDual_norm2 = 0;
    for (int i = 0; i < numTot; i++)
      workDual_norm2 += simplex_info_->workDual_[i] * simplex_info_->workDual_[i];
    workDual_norm2 = sqrt(workDual_norm2);
    //  printf("computeDual: B.pi=c_B has ||c_B||=%11.4g; ||pi||=%11.4g;
    //  ||pi^TA-c||=%11.4g\n", btranRHS_norm2, btranSol_norm2, workDual_norm2);
    double current_dual_feasibility_tolerance = simplex_info_->dual_feasibility_tolerance;
    double new_dual_feasibility_tolerance = workDual_norm2 / 1e16;
    if (new_dual_feasibility_tolerance > 1e-1) {
      printf(
          "Seriously: do you expect to solve an LP with ||pi^TA-c||=%11.4g?\n",
          workDual_norm2);
    } else if (new_dual_feasibility_tolerance > 10 * current_dual_feasibility_tolerance) {
      printf(
	     "||pi^TA-c|| = %12g so solving with dual_feasibility_tolerance = %12g\n",
	     workDual_norm2, new_dual_feasibility_tolerance);
      simplex_info_->dual_feasibility_tolerance = new_dual_feasibility_tolerance;
    }
  }

  // Now have nonbasic duals
  simplex_info_->solver_lp_has_nonbasic_dual_values = true;
}

// Compute the number of dual infeasibilities for the dual algorithm
void HModel::computeDualInfeasInDual(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HIGHS_CONST_INF;
  const double tau_d = simplex_info_->dual_feasibility_tolerance;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!basis_->nonbasicFlag_[i]) continue;
    // Free
    if (simplex_info_->workLower_[i] == -inf && simplex_info_->workUpper_[i] == inf)
      workCount += (fabs(simplex_info_->workDual_[i]) >= tau_d);
    // In dual, assuming that boxed variables will be flipped
    if (simplex_info_->workLower_[i] == -inf || simplex_info_->workUpper_[i] == inf)
      workCount += (basis_->nonbasicMove_[i] * simplex_info_->workDual_[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Compute the number of dual infeasibilities for the primal?? algorithm
void HModel::computeDualInfeasInPrimal(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HIGHS_CONST_INF;
  const double tau_d = simplex_info_->dual_feasibility_tolerance;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!basis_->nonbasicFlag_[i]) continue;
    // Free
    if (simplex_info_->workLower_[i] == -inf && simplex_info_->workUpper_[i] == inf)
      workCount += (fabs(simplex_info_->workDual_[i]) >= tau_d);
    // In primal don't assume flip
    workCount += (basis_->nonbasicMove_[i] * simplex_info_->workDual_[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Correct dual values
void HModel::correctDual(int *freeInfeasCount) {
  const double tau_d = simplex_info_->dual_feasibility_tolerance;
  const double inf = HIGHS_CONST_INF;
  int workCount = 0;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++) {
    if (basis_->nonbasicFlag_[i]) {
      if (simplex_info_->workLower_[i] == -inf && simplex_info_->workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(simplex_info_->workDual_[i]) >= tau_d);
      } else if (basis_->nonbasicMove_[i] * simplex_info_->workDual_[i] <= -tau_d) {
        if (simplex_info_->workLower_[i] != -inf && simplex_info_->workUpper_[i] != inf) {
          // Boxed variable = flip
          flipBound(i);
        } else {
          // Other variable = shift
          simplex_info_->costs_perturbed = 1;
          if (basis_->nonbasicMove_[i] == 1) {
            double random_v = random_->fraction();
            double dual = (1 + random_v) * tau_d;
            double shift = dual - simplex_info_->workDual_[i];
            simplex_info_->workDual_[i] = dual;
            simplex_info_->workCost_[i] = simplex_info_->workCost_[i] + shift;
          } else {
            double dual = -(1 + random_->fraction()) * tau_d;
            double shift = dual - simplex_info_->workDual_[i];
            simplex_info_->workDual_[i] = dual;
            simplex_info_->workCost_[i] = simplex_info_->workCost_[i] + shift;
          }
        }
      }
    }
  }

  *freeInfeasCount = workCount;
}

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void HModel::computePrimal() {
  // Setup a local buffer for the values of basic variables
  HVector buffer;
  buffer.setup(solver_lp_->numRow_);
  buffer.clear();
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  for (int i = 0; i < numTot; i++)
    if (basis_->nonbasicFlag_[i] && simplex_info_->workValue_[i] != 0)
      matrix_->collect_aj(buffer, i, simplex_info_->workValue_[i]);
  factor_->ftran(buffer, 1);

  for (int i = 0; i < solver_lp_->numRow_; i++) {
    int iCol = basis_->basicIndex_[i];
    simplex_info_->baseValue_[i] = -buffer.array[i];
    simplex_info_->baseLower_[i] = simplex_info_->workLower_[iCol];
    simplex_info_->baseUpper_[i] = simplex_info_->workUpper_[iCol];
  }
  // Now have basic primals
  simplex_info_->solver_lp_has_basic_primal_values = true;
}

// Compute the (primal) objective via primal values and costs
double HModel::computePh2Objective(vector<double> &colPrAct) {
  double Ph2Objective = 0;
  for (int i = 0; i < solver_lp_->numCol_; i++) Ph2Objective += colPrAct[i] * solver_lp_->colCost_[i];
  //  printf("Ph2Objective Ph2Objective = %g\n", Ph2Objective);
  Ph2Objective *= scale_->cost_;
  return Ph2Objective;
}

// Compute the (primal) objective via primal values of basic and nonbasic
// columns and their costs
double HModel::computePrObj() {
  double prObj = 0;
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = basis_->basicIndex_[row];
    if (var < solver_lp_->numCol_) prObj += simplex_info_->baseValue_[row] * solver_lp_->colCost_[var];
  }
  for (int col = 0; col < solver_lp_->numCol_; col++)
    if (basis_->nonbasicFlag_[col]) prObj += simplex_info_->workValue_[col] * solver_lp_->colCost_[col];
  prObj *= scale_->cost_;
  return prObj;
}

int HModel::handleRankDeficiency() {
  int rankDeficiency = factor_->rankDeficiency;
  const int *noPvC = factor_->getNoPvC();
  printf("Returned %d = factor_->build();\n", rankDeficiency);
  fflush(stdout);
  vector<int> basicRows;
  const int numTot = solver_lp_->numCol_ + solver_lp_->numRow_;
  basicRows.resize(numTot);
  //    printf("Before - basis_->basicIndex_:"); for (int iRow=0; iRow<solver_lp_->numRow_; iRow++)
  //    printf(" %2d", basis_->basicIndex_[iRow]); printf("\n");
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++) basicRows[basis_->basicIndex_[iRow]] = iRow;
  for (int k = 0; k < rankDeficiency; k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor_->noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = solver_lp_->numCol_ + factor_->noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut, simplex_info_->workLower_[columnOut],
    //      simplex_info_->workUpper_[columnOut]);
    if (basis_->basicIndex_[rowOut] != columnOut) {
      printf("%d = basis_->basicIndex_[rowOut] != noPvC[k] = %d\n", basis_->basicIndex_[rowOut],
             columnOut);
      fflush(stdout);
    }
    int sourceOut = setSourceOutFmBd(columnOut);
    updatePivots(columnIn, rowOut, sourceOut);
    updateMatrix(columnIn, columnOut);
  }
  //    printf("After  - basis_->basicIndex_:"); for (int iRow=0; iRow<solver_lp_->numRow_; iRow++)
  //    printf(" %2d", basis_->basicIndex_[iRow]); printf("\n");
#ifdef HiGHSDEV
  factor_->checkInvert();
#endif
  return 0;
}

int HModel::setSourceOutFmBd(const int columnOut) {
  int sourceOut = 0;
  if (simplex_info_->workLower_[columnOut] != simplex_info_->workUpper_[columnOut]) {
    if (!highs_isInfinity(-simplex_info_->workLower_[columnOut])) {
      // Finite LB so sourceOut = -1 ensures value set to LB if LB < UB
      sourceOut = -1;
      //      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g]
      //      so setting sourceOut = -1\n", columnOut, simplex_info_->workLower_[columnOut],
      //      simplex_info_->workUpper_[columnOut]);
    } else {
      // Infinite LB so sourceOut = 1 ensures value set to UB
      sourceOut = 1;
      if (!highs_isInfinity(simplex_info_->workUpper_[columnOut])) {
        // Free variable => trouble!
        printf("TROUBLE: variable %d leaving the basis is free!\n", columnOut);
      }
    }
  }
  return sourceOut;
}

// Utilities for shifting costs and flipping bounds
// Record the shift in the cost of a particular column
void HModel::shiftCost(int iCol, double amount) {
  simplex_info_->costs_perturbed = 1;
  assert(simplex_info_->workShift_[iCol] == 0);
  simplex_info_->workShift_[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void HModel::shiftBack(int iCol) {
  simplex_info_->workDual_[iCol] -= simplex_info_->workShift_[iCol];
  simplex_info_->workShift_[iCol] = 0;
}

// Flip a primal bound
void HModel::flipBound(int iCol) {
  const int move = basis_->nonbasicMove_[iCol] = -basis_->nonbasicMove_[iCol];
  simplex_info_->workValue_[iCol] = move == 1 ? simplex_info_->workLower_[iCol] : simplex_info_->workUpper_[iCol];
}

// The major model updates. Factor calls factor_->update; Matrix
// calls matrix_->update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HModel::updateFactor(HVector *column, HVector *row_ep, int *iRow,
                          int *hint) {
  //  HighsTimer &timer = highs_model_object->timer_;
  timer_->start(simplex_info_->clock_[UpdateFactorClock]);
  
  factor_->update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  simplex_info_->solver_lp_has_invert = true;
  if (simplex_info_->update_count >= simplex_info_->update_limit) *hint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  timer_->stop(simplex_info_->clock_[UpdateFactorClock]);
}

void HModel::updateMatrix(int columnIn, int columnOut) {
  timer_->start(simplex_info_->clock_[UpdateMatrixClock]);
  matrix_->update(columnIn, columnOut);
  timer_->stop(simplex_info_->clock_[UpdateMatrixClock]);
}

void HModel::updatePivots(int columnIn, int rowOut, int sourceOut) {
  timer_->start(simplex_info_->clock_[UpdatePivotsClock]);
  int columnOut = basis_->basicIndex_[rowOut];

  // Incoming variable
  basis_->basicIndex_[rowOut] = columnIn;
  basis_->nonbasicFlag_[columnIn] = 0;
  basis_->nonbasicMove_[columnIn] = 0;
  simplex_info_->baseLower_[rowOut] = simplex_info_->workLower_[columnIn];
  simplex_info_->baseUpper_[rowOut] = simplex_info_->workUpper_[columnIn];

  // Outgoing variable
  basis_->nonbasicFlag_[columnOut] = 1;
  //  double dlValue;
  //  double vrLb = simplex_info_->workLower_[columnOut];
  //  double vrV = simplex_info_->workValue_[columnOut];
  //  double vrUb = simplex_info_->workUpper_[columnOut];
  if (simplex_info_->workLower_[columnOut] == simplex_info_->workUpper_[columnOut]) {
    //    dlValue = simplex_info_->workLower_[columnOut]-simplex_info_->workValue_[columnOut];
    simplex_info_->workValue_[columnOut] = simplex_info_->workLower_[columnOut];
    basis_->nonbasicMove_[columnOut] = 0;
  } else if (sourceOut == -1) {
    //    dlValue = simplex_info_->workLower_[columnOut]-simplex_info_->workValue_[columnOut];
    simplex_info_->workValue_[columnOut] = simplex_info_->workLower_[columnOut];
    basis_->nonbasicMove_[columnOut] = 1;
  } else {
    //    dlValue = simplex_info_->workUpper_[columnOut]-simplex_info_->workValue_[columnOut];
    simplex_info_->workValue_[columnOut] = simplex_info_->workUpper_[columnOut];
    basis_->nonbasicMove_[columnOut] = -1;
  }
  double nwValue = simplex_info_->workValue_[columnOut];
  double vrDual = simplex_info_->workDual_[columnOut];
  double dlDualObjectiveValue = nwValue*vrDual;
  //  if (abs(nwValue))
  //    printf("HModel::updatePivots columnOut = %6d (%2d): [%11.4g, %11.4g, %11.4g], nwValue = %11.4g, dual = %11.4g, dlObj = %11.4g\n",
  //			   columnOut, basis_->nonbasicMove_[columnOut], vrLb, vrV, vrUb, nwValue, vrDual, dlDualObjectiveValue);
  simplex_info_->updatedDualObjectiveValue += dlDualObjectiveValue;
  simplex_info_->update_count++;
  // Update the number of basic logicals
  if (columnOut < solver_lp_->numCol_) simplex_info_->num_basic_logicals -= 1;
  if (columnIn < solver_lp_->numCol_) simplex_info_->num_basic_logicals += 1;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  simplex_info_->solver_lp_has_invert = false;
  simplex_info_->solver_lp_has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  simplex_info_->solver_lp_has_fresh_rebuild = false;
  timer_->stop(simplex_info_->clock_[UpdatePivotsClock]);
}

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
  populate_WorkArrays();
  // simplex.method_update_solver_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
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
  double prObjVal = computePrObj();
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
  HighsPrintMessage(ML_MINIMAL, "grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s,%d,%d,%d\n",
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
void HModel::util_analyseLpSolution() {
  if (simplex_info_->solution_status != SimplexSolutionStatus::OPTIMAL) return;
  printf("\nAnalysing the model solution\n");
  fflush(stdout);
  const double inf = HIGHS_CONST_INF;
  const double tlValueEr = 1e-8;
  const double tlPrRsduEr = 1e-8;
  const double tlDuRsduEr = 1e-8;
  const double tlPrIfs = simplex_info_->primal_feasibility_tolerance;
  const double tlDuIfs = simplex_info_->dual_feasibility_tolerance;

  // Copy the values of (nonbasic) primal variables and scatter values of primal
  // variables which are basic
  vector<double> value = simplex_info_->workValue_;
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++)
    value[basis_->basicIndex_[iRow]] = simplex_info_->baseValue_[iRow];

  // Copy the values of (nonbasic) dual variables and zero values of dual
  // variables which are basic
  vector<double> dual = simplex_info_->workDual_;
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++) dual[basis_->basicIndex_[iRow]] = 0;

  // Allocate and zero values of row primal activites and column dual activities
  // to check the residuals
  vector<double> sclRowPrAct;
  vector<double> rowPrAct;
  sclRowPrAct.assign(solver_lp_->numRow_, 0);
  rowPrAct.assign(solver_lp_->numRow_, 0);
  vector<double> sclColDuAct;
  vector<double> colDuAct;
  sclColDuAct.assign(solver_lp_->numCol_, 0);
  colDuAct.assign(solver_lp_->numCol_, 0);

  // Determine row primal activites and column dual activities
  for (int iCol = 0; iCol < solver_lp_->numCol_; iCol++) {
    //    printf("\nCol %2d\n", iCol);
    double lcSclColDuAct = -solver_lp_->colCost_[iCol];
    double lcColDuAct = -(solver_lp_->colCost_[iCol] * scale_->cost_) / scale_->col_[iCol];
    for (int en = solver_lp_->Astart_[iCol]; en < solver_lp_->Astart_[iCol + 1]; en++) {
      int iRow = solver_lp_->Aindex_[en];
      double Avalue_En = solver_lp_->Avalue_[en];
      double unsclAvalue_En = Avalue_En / (scale_->col_[iCol] * scale_->row_[iRow]);
      sclRowPrAct[iRow] += Avalue_En * value[iCol];
      rowPrAct[iRow] += unsclAvalue_En * value[iCol] * scale_->col_[iCol];
      //      double lcSum = lcSclColDuAct - Avalue_En*dual[solver_lp_->numCol_+iRow];
      //      printf("Row %2d: %11.4g - (%11.4g*%11.4g=%11.4g) = %11.4g\n",
      //      iRow, lcSclColDuAct, Avalue_En, dual[solver_lp_->numCol_+iRow],
      //      Avalue_En*dual[solver_lp_->numCol_+iRow], lcSum);
      lcSclColDuAct -= Avalue_En * dual[solver_lp_->numCol_ + iRow];
      lcColDuAct -=
          unsclAvalue_En * dual[solver_lp_->numCol_ + iRow] * scale_->cost_ * scale_->row_[iRow];
    }
    sclColDuAct[iCol] = lcSclColDuAct;
    colDuAct[iCol] = lcColDuAct;
  }

  // Look for column residual errors and infeasibilities - primal and dual
  if (solver_lp_->offset_) printf("Primal objective offset is %11.4g\n", solver_lp_->offset_);
  double lcPrObjV = 0;
  double lcValue = 0;

  int numRpFreeRowEr = 0;
  int maxRpFreeRowEr = 100;
  int numRpFreeColEr = 0;
  int maxRpFreeColEr = 100;

  bool rpAllCol = false;
  int numRpCol = 0;
  int mxRpCol = 100;
  bool rpNoCol = false;
  int numColPrIfs = 0;
  double maxColPrIfs = 0;
  double sumColPrIfs = 0;
  int numSclColPrIfs = 0;
  double maxSclColPrIfs = 0;
  double sumSclColPrIfs = 0;
  int numColDuIfs = 0;
  double maxColDuIfs = 0;
  double sumColDuIfs = 0;
  int numSclColDuIfs = 0;
  double maxSclColDuIfs = 0;
  double sumSclColDuIfs = 0;
  int numColDuRsduEr = 0;
  double sumColDuRsduEr = 0;
  double maxColDuRsduEr = 0;
  int numSclColDuRsduEr = 0;
  double sumSclColDuRsduEr = 0;
  double maxSclColDuRsduEr = 0;
  for (int iCol = 0; iCol < solver_lp_->numCol_; iCol++) {
    double sclColValue;
    double sclColDuIfs;
    // Get the unscaled column bounds
    double unsclColLower = solver_lp_->colLower_[iCol];
    double unsclColUpper = solver_lp_->colUpper_[iCol];
    unsclColLower *= unsclColLower == -inf ? 1 : scale_->col_[iCol];
    unsclColUpper *= unsclColUpper == +inf ? 1 : scale_->col_[iCol];
    // Determine the column primal values given nonbasicMove and the bounds -
    // and check the dual residual errors and infeasibilities
    if (basis_->nonbasicFlag_[iCol]) {
      // Nonbasic variable - check that the value array is correct given
      // nonbasicMove and the bounds
      if (basis_->nonbasicMove_[iCol] == NONBASIC_MOVE_UP) {
        // At lower bound
        sclColValue = solver_lp_->colLower_[iCol];
        sclColDuIfs = max(-dual[iCol], 0.);
      } else if (basis_->nonbasicMove_[iCol] == NONBASIC_MOVE_DN) {
        // At upper bound
        sclColValue = solver_lp_->colUpper_[iCol];
        sclColDuIfs = max(dual[iCol], 0.);
      } else {
        // Fixed or free
        if (solver_lp_->colLower_[iCol] == solver_lp_->colUpper_[iCol]) {
          sclColValue = solver_lp_->colUpper_[iCol];
          sclColDuIfs = 0;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!highs_isInfinity(-solver_lp_->colLower_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has lower bound of %g\n",
                  iCol, solver_lp_->colLower_[iCol]);
            }
          }
          if (!highs_isInfinity(solver_lp_->colUpper_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has upper bound of %g\n",
                  iCol, solver_lp_->colUpper_[iCol]);
            }
          }
          sclColValue = value[iCol];
          sclColDuIfs = abs(dual[iCol]);
          //	  if (!freeEr) {printf("Column %7d is free with value %g\n",
          // iCol ,sclColValue);}
        }
      }
      double valueEr = abs(sclColValue - value[iCol]);
      if (valueEr > tlValueEr) {
        printf(
            "Column %7d has value error of %11.4g for sclColValue = %11.4g and "
            "value[iCol] = %11.4g\n",
            iCol, valueEr, sclColValue, value[iCol]);
        sclColValue = value[iCol];
      }

    } else {
      // Basic variable
      sclColValue = value[iCol];
      sclColDuIfs = abs(dual[iCol]);
    }

    lcPrObjV += sclColValue * solver_lp_->colCost_[iCol];

    double unsclColValue = sclColValue * scale_->col_[iCol];
    //      assert(highs_isInfinity(-sclColValue));
    //      assert(highs_isInfinity(sclColValue));
    // Assess primal infeasibility
    // For scaled values
    double sclColPrIfs = max(
        max(solver_lp_->colLower_[iCol] - sclColValue, sclColValue - solver_lp_->colUpper_[iCol]), 0.0);
    if (sclColPrIfs > tlPrIfs) {
      numSclColPrIfs++;
      sumSclColPrIfs += sclColPrIfs;
    }
    maxSclColPrIfs = max(sclColPrIfs, maxSclColPrIfs);
    // For unscaled values
    double colPrIfs = max(
        max(unsclColLower - unsclColValue, unsclColValue - unsclColUpper), 0.0);
    if (colPrIfs > tlPrIfs) {
      numColPrIfs++;
      sumColPrIfs += colPrIfs;
    }
    maxColPrIfs = max(colPrIfs, maxColPrIfs);

    // Assess dual infeasibility
    // In scaled values
    if (sclColDuIfs > tlDuIfs) {
      numSclColDuIfs++;
      sumSclColDuIfs += sclColDuIfs;
    }
    maxSclColDuIfs = max(sclColDuIfs, maxSclColDuIfs);
    // In unscaled values
    double colDuIfs = sclColDuIfs * scale_->cost_ / scale_->col_[iCol];
    if (colDuIfs > tlDuIfs) {
      numColDuIfs++;
      sumColDuIfs += colDuIfs;
    }
    maxColDuIfs = max(colDuIfs, maxColDuIfs);

    // Check column residual errors
    // Using scaled column activities
    double sclColDual = dual[iCol];
    double sclColDuRsduEr = abs(sclColDuAct[iCol] + sclColDual);
    if (sclColDuRsduEr > tlDuRsduEr) {
      /*
      bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
      if (rpCol) {
        numRpCol++;
        printf("Col    %7d has a   dual residual error of %11.4g for
      sclColDuAct[iCol] = %11.4g and -sclColDual = %11.4g\n", iCol,
      sclColDuRsduEr, sclColDuAct[iCol], -sclColDual);
      }
      */
      numSclColDuRsduEr++;
      sumSclColDuRsduEr += sclColDuRsduEr;
    }
    maxSclColDuRsduEr = max(sclColDuRsduEr, maxSclColDuRsduEr);
    // Using unscaled column activities
    double colDual = sclColDual * scale_->cost_ / scale_->col_[iCol];
    double colDuRsduEr = abs(colDuAct[iCol] + colDual);
    if (colDuRsduEr > tlDuRsduEr) {
      /*
      bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
      if (rpCol) {
        numRpCol++;
        printf("Col    %7d has a   dual residual error of %11.4g for
      colDuAct[iCol] = %11.4g and -colDual = %11.4g\n", iCol, colDuRsduEr,
      colDuAct[iCol], -colDual);
      }
      */
      numColDuRsduEr++;
      sumColDuRsduEr += colDuRsduEr;
    }
    maxColDuRsduEr = max(colDuRsduEr, maxColDuRsduEr);

    bool erFd = sclColPrIfs > tlPrIfs || colPrIfs > tlPrIfs ||
                sclColDuIfs > tlDuIfs || colDuIfs > tlDuIfs ||
                sclColDuRsduEr > tlDuRsduEr || colDuRsduEr > tlDuRsduEr;
    bool rpCol = (rpAllCol || (numRpCol < mxRpCol && erFd)) && !rpNoCol;
    if (rpCol) {
      numRpCol++;
      printf("\nCol %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iCol,
             basis_->nonbasicFlag_[iCol], basis_->nonbasicMove_[iCol], scale_->col_[iCol]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          solver_lp_->colLower_[iCol], sclColValue, solver_lp_->colUpper_[iCol], sclColPrIfs, sclColDuIfs,
          sclColDuRsduEr);
      printf(
          "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: %11.4g) "
          "\n",
          unsclColLower, unsclColValue, unsclColUpper, colPrIfs, colDuIfs,
          colDuRsduEr);
    }
  }

  printf(
      "Found %6d   scaled column primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs);
  printf(
      "Found %6d unscaled column primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numColPrIfs, sumColPrIfs, maxColPrIfs);
  printf(
      "Found %6d   scaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs);
  printf(
      "Found %6d unscaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numColDuIfs, sumColDuIfs, maxColDuIfs);
  printf(
      "Found %6d   scaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuRsduEr, sumSclColDuRsduEr, maxSclColDuRsduEr);
  printf(
      "Found %6d unscaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Col,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs, numColPrIfs, sumColPrIfs,
      maxColPrIfs, numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs, numColDuIfs,
      sumColDuIfs, maxColDuIfs, numSclColDuRsduEr, sumSclColDuRsduEr,
      maxSclColDuRsduEr, numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  bool rpAllRow = false;
  int numRpRow = 0;
  int mxRpRow = 100;
  bool rpNoRow = false;
  int numRowPrIfs = 0;
  double sumRowPrIfs = 0;
  double maxRowPrIfs = 0;
  int numSclRowPrIfs = 0;
  double sumSclRowPrIfs = 0;
  double maxSclRowPrIfs = 0;
  int numRowDuIfs = 0;
  double maxRowDuIfs = 0;
  double sumRowDuIfs = 0;
  int numSclRowDuIfs = 0;
  double maxSclRowDuIfs = 0;
  double sumSclRowDuIfs = 0;
  int numRowPrRsduEr = 0;
  double sumRowPrRsduEr = 0;
  double maxRowPrRsduEr = 0;
  int numSclRowPrRsduEr = 0;
  double sumSclRowPrRsduEr = 0;
  double maxSclRowPrRsduEr = 0;
  for (int iRow = 0; iRow < solver_lp_->numRow_; iRow++) {
    double sclRowValue;
    double sclRowDuIfs;
    // Get the unscaled row bounds
    double unsclRowLower = solver_lp_->rowLower_[iRow];
    double unsclRowUpper = solver_lp_->rowUpper_[iRow];
    unsclRowLower *= unsclRowLower == -inf ? 1 : scale_->row_[iRow];
    unsclRowUpper *= unsclRowUpper == +inf ? 1 : scale_->row_[iRow];
    // Determine the row primal values given nonbasicMove and the bounds - and
    // check the dual residual errors and infeasibilities
    if (basis_->nonbasicFlag_[solver_lp_->numCol_ + iRow]) {
      // Nonbasic variable
      if (basis_->nonbasicMove_[solver_lp_->numCol_ + iRow] == NONBASIC_MOVE_DN) {
        // At lower bound
        sclRowValue = solver_lp_->rowLower_[iRow];
        sclRowDuIfs = max(dual[solver_lp_->numCol_ + iRow], 0.);
      } else if (basis_->nonbasicMove_[solver_lp_->numCol_ + iRow] == NONBASIC_MOVE_UP) {
        // At upper bound
        sclRowValue = solver_lp_->rowUpper_[iRow];
        sclRowDuIfs = max(-dual[solver_lp_->numCol_ + iRow], 0.);
      } else {
        // Fixed or free
        if (solver_lp_->rowLower_[iRow] == solver_lp_->rowUpper_[iRow]) {
          sclRowValue = solver_lp_->rowUpper_[iRow];
          sclRowDuIfs = 0.;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!highs_isInfinity(-solver_lp_->rowLower_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has lower bound of %g\n",
                  iRow, solver_lp_->rowLower_[iRow]);
            }
          }
          if (!highs_isInfinity(solver_lp_->rowUpper_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has upper bound of %g\n",
                  iRow, solver_lp_->rowUpper_[iRow]);
            }
          }
          sclRowValue = -value[solver_lp_->numCol_ + iRow];
          sclRowDuIfs = abs(dual[solver_lp_->numCol_ + iRow]);
          //	  if (!freeEr) {printf("Row    %7d is free with value %g\n",
          // iRow, sclRowValue);}
        }
      }
      double valueEr = abs(sclRowValue + value[solver_lp_->numCol_ + iRow]);
      if (valueEr > tlValueEr) {
        printf(
            "Row    %7d has value error of %11.4g for sclRowValue = %11.4g and "
            "-value[solver_lp_->numCol_+iRow] = %11.4g\n",
            iRow, valueEr, sclRowValue, -value[solver_lp_->numCol_ + iRow]);
        sclRowValue = -value[solver_lp_->numCol_ + iRow];
      }
    } else {
      // Basic variable
      sclRowValue = -value[solver_lp_->numCol_ + iRow];
      sclRowDuIfs = abs(dual[solver_lp_->numCol_ + iRow]);
    }
    //      assert(highs_isInfinity(-sclRowValue));
    //      assert(highs_isInfinity(sclRowValue));
    double unsclRowValue = sclRowValue * scale_->row_[iRow];

    // Assess primal infeasibility
    // For scaled values
    double sclRowPrIfs = max(
        max(solver_lp_->rowLower_[iRow] - sclRowValue, sclRowValue - solver_lp_->rowUpper_[iRow]), 0.0);
    if (sclRowPrIfs > tlPrIfs) {
      numSclRowPrIfs++;
      sumSclRowPrIfs += sclRowPrIfs;
    }
    maxSclRowPrIfs = max(sclRowPrIfs, maxSclRowPrIfs);
    // For unscaled values
    double rowPrIfs = max(
        max(unsclRowLower - unsclRowValue, unsclRowValue - unsclRowUpper), 0.0);
    if (rowPrIfs > tlPrIfs) {
      numRowPrIfs++;
      sumRowPrIfs += rowPrIfs;
    }
    maxRowPrIfs = max(rowPrIfs, maxRowPrIfs);

    // Assess dual infeasibility
    // In scaled values
    if (sclRowDuIfs > tlDuIfs) {
      numSclRowDuIfs++;
      sumSclRowDuIfs += sclRowDuIfs;
    }
    maxSclRowDuIfs = max(sclRowDuIfs, maxSclRowDuIfs);
    // In unscaled values
    double rowDuIfs = sclRowDuIfs * scale_->cost_ / scale_->row_[iRow];
    if (rowDuIfs > tlDuIfs) {
      numRowDuIfs++;
      sumRowDuIfs += rowDuIfs;
    }
    maxRowDuIfs = max(rowDuIfs, maxRowDuIfs);

    // Check row residual errors
    // Using scaled row activities
    double sclRowPrRsduEr = abs(sclRowPrAct[iRow] - sclRowValue);
    if (sclRowPrRsduEr > tlPrRsduEr) {
      /*
      bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
      if (rpRow) {
        numRpRow++;
        printf("Row    %7d has a primal residual error of %11.4g for
      sclRowPrAct[iRow] = %11.4g and sclRowValue = %11.4g\n", iRow,
      sclRowPrRsduEr, sclRowPrAct[iRow], sclRowValue);
      }
      */
      numSclRowPrRsduEr++;
      sumSclRowPrRsduEr += sclRowPrRsduEr;
    }
    maxSclRowPrRsduEr = max(sclRowPrRsduEr, maxSclRowPrRsduEr);
    // Using unscaled row activities
    double rowValue = sclRowValue / scale_->row_[iRow];
    double rowPrRsduEr = abs(rowPrAct[iRow] - rowValue);
    if (rowPrRsduEr > tlPrRsduEr) {
      /*
      bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
      if (rpRow) {
        numRpRow++;
        printf("Row    %7d has a primal residual error of %11.4g for
      rowPrAct[iRow] = %11.4g and rowValue = %11.4g\n", iRow, rowPrRsduEr,
      rowPrAct[iRow], rowValue);
      }
      */
      numRowPrRsduEr++;
      sumRowPrRsduEr += rowPrRsduEr;
    }
    maxRowPrRsduEr = max(rowPrRsduEr, maxRowPrRsduEr);

    bool erFd = sclRowPrIfs > tlPrIfs || rowPrIfs > tlPrIfs ||
                sclRowDuIfs > tlDuIfs || rowDuIfs > tlDuIfs ||
                sclRowPrRsduEr > tlPrRsduEr || rowPrRsduEr > tlPrRsduEr;
    bool rpRow = (rpAllRow || (numRpRow < mxRpRow && erFd)) && !rpNoRow;
    if (rpRow) {
      numRpRow++;
      printf("Row %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iRow,
             basis_->nonbasicFlag_[solver_lp_->numCol_ + iRow], basis_->nonbasicMove_[solver_lp_->numCol_ + iRow],
             scale_->row_[iRow]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          solver_lp_->rowLower_[iRow], sclRowValue, solver_lp_->rowUpper_[iRow], sclRowPrIfs, sclRowDuIfs,
          sclRowPrRsduEr);
      printf(
          "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          unsclRowLower, unsclRowValue, unsclRowUpper, rowPrIfs, rowDuIfs,
          rowPrRsduEr);
    }
  }
  printf(
      "Found %6d   scaled    row primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs);
  printf(
      "Found %6d unscaled    row primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numRowPrIfs, sumRowPrIfs, maxRowPrIfs);
  printf(
      "Found %6d   scaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs);
  printf(
      "Found %6d unscaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numRowDuIfs, sumRowDuIfs, maxRowDuIfs);
  printf(
      "Found %6d   scaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclRowPrRsduEr, sumSclRowPrRsduEr, maxSclRowPrRsduEr);
  printf(
      "Found %6d unscaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Row,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs, numRowPrIfs, sumRowPrIfs,
      maxRowPrIfs, numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs, numRowDuIfs,
      sumRowDuIfs, maxRowDuIfs, numSclRowPrRsduEr, sumSclRowPrRsduEr,
      maxSclRowPrRsduEr, numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  lcPrObjV *= scale_->cost_;
  lcPrObjV += solver_lp_->offset_;
  double dualObjectiveValue = simplex_info_->dualObjectiveValue;
  double ObjEr = abs(dualObjectiveValue - lcPrObjV) / max(1.0, fabs(dualObjectiveValue));
  printf(
	 "Relative objective error of %11.4g: dualObjectiveValue = %g; lcPrObjV = %g\n",
	 ObjEr, dualObjectiveValue, lcPrObjV);
}

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
