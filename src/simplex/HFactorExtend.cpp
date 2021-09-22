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
/**@file simplex/HFactorExtend.cpp
 * @brief Types of solution classes
 */
#include <cassert>

#include "simplex/HFactor.h"

using std::fabs;

void HFactor::addCols(const HighsInt num_new_col) {
  invalidAMatrixAction();
  numCol += num_new_col;
}

void HFactor::deleteNonbasicCols(const HighsInt num_deleted_col) {
  invalidAMatrixAction();
  numCol -= num_deleted_col;
}

void HFactor::addRows(const HighsSparseMatrix* ar_matrix) {
  invalidAMatrixAction();
  assert(kExtendInvertWhenAddingRows);
  HighsInt num_new_row = ar_matrix->num_row_;
  HighsInt new_num_row = numRow + num_new_row;
  printf("Adding %" HIGHSINT_FORMAT
         " new rows to HFactor instance: increasing dimension from "
         "%" HIGHSINT_FORMAT " to %" HIGHSINT_FORMAT " \n",
         num_new_row, numRow, new_num_row);

  // Need to know where (if) a column is basic
  vector<HighsInt> in_basis;
  in_basis.assign(numCol, -1);
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    HighsInt iVar = baseIndex[iRow];
    if (iVar < numCol) in_basis[iVar] = iRow;
  }
  for (HighsInt iRow = numRow; iRow < new_num_row; iRow++) {
    HighsInt iVar = baseIndex[iRow];
    assert(iVar >= numCol);
  }
  //  reportLu(kReportLuBoth, true);
  //  reportLu(kReportLuJustL);

  // Create a row-wise sparse matrix containing the new rows of the L
  // matrix - so that a column-wise version can be created (after
  // inserting the rows into the LR matrix) allowing the cew column
  // entries to be inserted efficiently into the L matrix
  HighsSparseMatrix new_lr_rows;
  new_lr_rows.format_ = MatrixFormat::kRowwise;
  new_lr_rows.num_col_ = numRow;
  double expected_density = 0.0;
  HVector rhs;
  rhs.setup(numRow);
  this->LRstart.reserve(new_num_row + 1);
  for (HighsInt inewRow = 0; inewRow < num_new_row; inewRow++) {
    //    printf("\nFor new row %" HIGHSINT_FORMAT "\n", inewRow);
    // Prepare RHS for system U^T.v = r
    rhs.clear();
    rhs.packFlag = true;
    for (HighsInt iEl = ar_matrix->start_[inewRow];
         iEl < ar_matrix->start_[inewRow + 1]; iEl++) {
      HighsInt iCol = ar_matrix->index_[iEl];
      HighsInt basis_index = in_basis[iCol];
      if (basis_index >= 0) {
        rhs.array[basis_index] = ar_matrix->value_[iEl];
        rhs.index[rhs.count++] = basis_index;
      }
    }
    // Solve U^T.v = r
    btranU(rhs, expected_density);
    double local_density = (1.0 * rhs.count) / numRow;
    expected_density = kRunningAverageMultiplier * local_density +
                       (1 - kRunningAverageMultiplier) * expected_density;
    //    printf("New row btranU density: local = %11.4g; expected =  %11.4g\n",
    //    local_density, expected_density);
    rhs.tight();
    //
    // Append v to the matrix containing the new rows of L
    HighsInt rhs_num_nz = rhs.count;
    //    printf("New row has entries:\n");
    for (HighsInt iX = 0; iX < rhs_num_nz; iX++) {
      HighsInt iCol = rhs.index[iX];
      new_lr_rows.index_.push_back(iCol);
      new_lr_rows.value_.push_back(rhs.array[iCol]);
      //      printf("   %4d, %11.4g\n", (int)iCol, rhs.array[iCol]);
    }
    new_lr_rows.start_.push_back(new_lr_rows.index_.size());
    new_lr_rows.num_row_++;

    //
    // Append v to the L matrix
    for (HighsInt iX = 0; iX < rhs_num_nz; iX++) {
      HighsInt iCol = rhs.index[iX];
      LRindex.push_back(iCol);
      LRvalue.push_back(rhs.array[iCol]);
    }
    LRstart.push_back(LRindex.size());
    //    reportLu(kReportLuJustL);
  }
  // Now create a column-wise copy of the new rows
  HighsSparseMatrix new_lr_cols = new_lr_rows;
  new_lr_cols.ensureColwise();
  //
  // Insert the column-wise copy into the L matrix
  //
  // Add pivot indices for the new columns
  this->LpivotIndex.resize(new_num_row);
  for (HighsInt iCol = numRow; iCol < new_num_row; iCol++)
    LpivotIndex[iCol] = iCol;
  //
  // Add starts for the identity columns
  HighsInt l_matrix_new_num_nz = LRindex.size();
  assert(l_matrix_new_num_nz == Lindex.size() + new_lr_cols.index_.size());
  Lstart.resize(new_num_row + 1);
  HighsInt to_el = l_matrix_new_num_nz;
  for (HighsInt iCol = numRow + 1; iCol < new_num_row + 1; iCol++)
    Lstart[iCol] = l_matrix_new_num_nz;
  //
  // Insert the new entries, remembering to offset the index values by
  // numRow, since new_lr_cols only has the new rows
  Lindex.resize(l_matrix_new_num_nz);
  Lvalue.resize(l_matrix_new_num_nz);
  for (HighsInt iCol = numRow - 1; iCol >= 0; iCol--) {
    const HighsInt from_el = Lstart[iCol + 1];
    Lstart[iCol + 1] = to_el;
    for (HighsInt iEl = new_lr_cols.start_[iCol + 1] - 1;
         iEl >= new_lr_cols.start_[iCol]; iEl--) {
      to_el--;
      Lindex[to_el] = numRow + new_lr_cols.index_[iEl];
      Lvalue[to_el] = new_lr_cols.value_[iEl];
    }
    for (HighsInt iEl = from_el - 1; iEl >= Lstart[iCol]; iEl--) {
      to_el--;
      Lindex[to_el] = Lindex[iEl];
      Lvalue[to_el] = Lvalue[iEl];
    }
  }
  assert(to_el == 0);
  this->LpivotLookup.resize(new_num_row);
  for (HighsInt iRow = numRow; iRow < new_num_row; iRow++)
    LpivotLookup[LpivotIndex[iRow]] = iRow;
  // Now update the U matrix with identity rows and columns
  // Allocate space for U factor
  //  reportLu(kReportLuJustL);
  //
  // Now add pivots corresponding to identity columns in U. The use of
  // most arrays in U is self-evident, except for U(R)lastp, which is
  // (one more than) the end of the index/value array that's applied
  // in the *tranU. Here it needs to be equal to the start since there
  // are no non-pivotal entries
  //

  HighsInt UcountX = Uindex.size();
  HighsInt u_pivot_lookup_offset = UpivotIndex.size() - numRow;
  for (HighsInt iRow = numRow; iRow < new_num_row; iRow++) {
    UpivotLookup.push_back(u_pivot_lookup_offset + iRow);
    UpivotIndex.push_back(iRow);
    UpivotValue.push_back(1);
    Ustart.push_back(UcountX);
    Ulastp.push_back(UcountX);
  }

  // Now, to extend UR
  //
  // UR space
  //
  // Borrowing names from buildFinish()
  HighsInt URstuffX = updateMethod == kUpdateMethodFt ? 5 : 0;
  HighsInt ur_size = URindex.size();
  HighsInt URcountX = ur_size + URstuffX * num_new_row;
  URindex.resize(URcountX);
  URvalue.resize(URcountX);

  // Need to refer to just the new UR vectors
  HighsInt ur_cur_num_vec = URstart.size();
  HighsInt ur_new_num_vec = ur_cur_num_vec + num_new_row;
  printf("\nUpdating UR vectors %d - %d\n", (int)ur_cur_num_vec,
         (int)ur_new_num_vec - 1);
  // UR pointer
  //
  // Allow space to the start of new rows, including the start for the
  // fictitious ur_new_num_vec'th row
  URstart.resize(ur_new_num_vec + 1);
  for (HighsInt iRow = ur_cur_num_vec + 1; iRow < ur_new_num_vec + 1; iRow++) {
    URstart[iRow] = ur_size;
  }
  // NB ur_temp plays the role of URlastp when it could be used as
  // temporary storage in buildFinish()
  vector<HighsInt> ur_temp;
  ur_temp.assign(ur_new_num_vec, 0);
  //
  // URspace has its new entries assigned (to URstuffX) as in
  // buildFinish()
  URspace.resize(ur_new_num_vec);
  for (HighsInt iRow = ur_cur_num_vec; iRow < ur_new_num_vec; iRow++)
    URspace[iRow] = URstuffX;
  // Compute ur_temp exactly as in buildFinish()
  for (HighsInt k = 0; k < UcountX; k++) ur_temp[UpivotLookup[Uindex[k]]]++;
  HighsInt iStart = ur_size;
  URstart[ur_cur_num_vec] = iStart;
  for (HighsInt iRow = ur_cur_num_vec + 1; iRow < ur_new_num_vec + 1; iRow++) {
    HighsInt gap = ur_temp[iRow - 1] + URstuffX;
    URstart[iRow] = iStart + gap;
    iStart += gap;
    printf("URstart[%d] = %d; gap = %d; iStart = %d\n", (int)iRow,
           (int)URstart[iRow], (int)gap, (int)iStart);
  }
  printf("URcountX = %d; iStart%d\n", (int)URcountX, (int)iStart);
  // Lose the start for the fictitious ur_new_num_vec'th row
  URstart.resize(ur_new_num_vec);
  // Resize URlastp and initialise its new entries to be the URstart
  // values since the rows are empty
  URlastp.resize(ur_new_num_vec);
  for (HighsInt iRow = ur_cur_num_vec; iRow < ur_new_num_vec; iRow++)
    URlastp[iRow] = URstart[iRow];
  //
  // Increase the number of rows in HFactor
  numRow += num_new_row;
  //  reportLu(kReportLuBoth, true);
}
