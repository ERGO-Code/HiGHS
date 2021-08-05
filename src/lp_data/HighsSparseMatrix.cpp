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
/**@file lp_data/HighsSparseMatrix.cpp
 * @brief
 */
#include "lp_data/HighsSparseMatrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "util/HighsSort.h"

using std::fabs;
using std::max;
using std::min;
using std::vector;

bool HighsSparseMatrix::operator==(const HighsSparseMatrix& matrix) const {
  bool equal = true;
  equal = this->format_ == matrix.format_ && equal;
  equal = this->num_col_ == matrix.num_col_ && equal;
  equal = this->num_row_ == matrix.num_row_ && equal;
  equal = this->start_ == matrix.start_ && equal;
  equal = this->index_ == matrix.index_ && equal;
  equal = this->value_ == matrix.value_ && equal;
  return equal;
}

void HighsSparseMatrix::clear() {
  this->num_col_ = 0;
  this->num_row_ = 0;
  this->start_.clear();
  this->p_end_.clear();
  this->index_.clear();
  this->value_.clear();
  this->format_ = MatrixFormat::kColwise;
  this->start_.assign(1, 0);
}

bool HighsSparseMatrix::isRowwise() const {
  return this->format_ == MatrixFormat::kRowwise ||
         this->format_ == MatrixFormat::kRowwisePartitioned;
}

bool HighsSparseMatrix::isColwise() const {
  return this->format_ == MatrixFormat::kColwise;
}

HighsInt HighsSparseMatrix::num_nz() const {
  assert(this->format_ != MatrixFormat::kNone);
  if (this->isColwise()) {
    if ((HighsInt)this->start_.size() < this->num_col_ + 1) {
      printf("%d = this->start_.size() < this->num_col_ + 1 = %d\n",
             (int)this->start_.size(), (int)this->num_col_ + 1);
    }
    assert((HighsInt)this->start_.size() >= this->num_col_ + 1);
    return this->start_[this->num_col_];
  } else {
    assert((HighsInt)this->start_.size() >= this->num_row_ + 1);
    return this->start_[this->num_row_];
  }
}

void HighsSparseMatrix::range(double& min_value, double& max_value) const {
  assert(this->format_ != MatrixFormat::kNone);
  for (HighsInt iEl = 0; iEl < this->start_[this->num_col_]; iEl++) {
    double value = fabs(this->value_[iEl]);
    min_value = min(min_value, value);
    max_value = max(max_value, value);
  }
}

HighsStatus HighsSparseMatrix::setFormat(const MatrixFormat desired_format) {
  if (desired_format == MatrixFormat::kNone) return HighsStatus::kError;
  assert(this->format_ != MatrixFormat::kNone);
  if (this->format_ == desired_format) return HighsStatus::kOk;
  if (desired_format == MatrixFormat::kColwise) {
    this->ensureColWise();
  } else {
    this->ensureRowWise();
  }
  assert(this->format_ == desired_format);
  return HighsStatus::kOk;
}

void HighsSparseMatrix::ensureColWise() {
  assert(this->format_ != MatrixFormat::kNone);
  // Should only call this is orientation is ROWWISE
  assert(this->isRowwise());
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->num_nz();
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
  if (num_nz == 0) {
    // Empty matrix, so just ensure that there are enough zero starts
    // for the new orientation
    this->start_.assign(num_col + 1, 0);
    this->index_.clear();
    this->value_.clear();
  } else {
    // Matrix is non-empty, so transpose it
    //
    // Take a copy of the current matrix - that is rowwise - so that
    // the current matrix is filled colwise
    vector<HighsInt> ARstart = this->start_;
    vector<HighsInt> ARindex = this->index_;
    vector<double> ARvalue = this->value_;
    this->start_.resize(num_col + 1);
    this->index_.resize(num_nz);
    this->value_.resize(num_nz);
    vector<HighsInt> Alength;
    Alength.assign(num_col, 0);
    for (HighsInt iEl = ARstart[0]; iEl < num_nz; iEl++)
      Alength[ARindex[iEl]]++;
    this->start_[0] = 0;
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      this->start_[iCol + 1] = this->start_[iCol] + Alength[iCol];
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      for (HighsInt iEl = ARstart[iRow]; iEl < ARstart[iRow + 1]; iEl++) {
        HighsInt iCol = ARindex[iEl];
        HighsInt iCol_el = this->start_[iCol];
        this->index_[iCol_el] = iRow;
        this->value_[iCol_el] = ARvalue[iEl];
        this->start_[iCol]++;
      }
    }
    this->start_[0] = 0;
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      this->start_[iCol + 1] = this->start_[iCol] + Alength[iCol];
    assert(this->start_[num_col] == num_nz);
  }
  this->format_ = MatrixFormat::kColwise;
  assert((HighsInt)this->start_.size() >= num_col + 1);
  num_nz = this->num_nz();
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
}

void HighsSparseMatrix::ensureRowWise() {
  assert(this->format_ != MatrixFormat::kNone);
  // Should only call this is orientation is COLWISE
  assert(this->isColwise());
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->num_nz();
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
  bool empty_matrix = num_col == 0 || num_row == 0;
  if (num_nz == 0) {
    // Empty matrix, so just ensure that there are enough zero starts
    // for the new orientation
    this->start_.assign(num_row + 1, 0);
    this->index_.clear();
    this->value_.clear();
  } else {
    // Matrix is non-empty, so transpose it
    //
    // Take a copy of the current matrix - that is colwise - so that
    // the current matrix is filled rowwise
    vector<HighsInt> Astart = this->start_;
    vector<HighsInt> Aindex = this->index_;
    vector<double> Avalue = this->value_;
    this->start_.resize(num_row + 1);
    this->index_.resize(num_nz);
    this->value_.resize(num_nz);
    vector<HighsInt> ARlength;
    ARlength.assign(num_row, 0);
    for (HighsInt iEl = Astart[0]; iEl < num_nz; iEl++) ARlength[Aindex[iEl]]++;
    this->start_[0] = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      this->start_[iRow + 1] = this->start_[iRow] + ARlength[iRow];
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      for (HighsInt iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
        HighsInt iRow = Aindex[iEl];
        HighsInt iRow_el = this->start_[iRow];
        this->index_[iRow_el] = iCol;
        this->value_[iRow_el] = Avalue[iEl];
        this->start_[iRow]++;
      }
    }
    this->start_[0] = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      this->start_[iRow + 1] = this->start_[iRow] + ARlength[iRow];
    assert(this->start_[num_row] == num_nz);
  }
  this->format_ = MatrixFormat::kRowwise;
  assert((HighsInt)this->start_.size() >= num_row + 1);
  num_nz = this->num_nz();
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
}

HighsStatus HighsSparseMatrix::addCols(const HighsInt num_new_col,
                                       const HighsInt num_new_nz,
                                       const HighsInt* new_matrix_start,
                                       const HighsInt* new_matrix_index,
                                       const double* new_matrix_value,
                                       const int8_t* in_partition) {
  if (this->format_ == MatrixFormat::kNone) {
    printf("Call with format_ == MatrixFormat::kNone\n");
  }
  assert(this->format_ != MatrixFormat::kNone);
  // Adding columns to a row-wise partitioned matrix needs the
  // partition information
  const bool partitioned = this->format_ == MatrixFormat::kRowwisePartitioned;
  if (partitioned) {
    //    if (in_partition == NULL) { printf("in_partition == NULL\n"); }
    assert(in_partition != NULL);
  }
  // Cannot handle the row-wise case
  assert(!this->isRowwise());
  if (num_new_col < 0) return HighsStatus::kError;
  if (num_new_nz < 0) return HighsStatus::kError;
  if (num_new_col == 0) {
    // No columns are being added, so check that no nonzeros are being
    // added
    assert(num_new_nz == 0);
    if (num_new_nz != 0) return HighsStatus::kError;
    return HighsStatus::kOk;
  }
  // Adding a positive number of columns to a matrix
  if (num_new_nz) {
    // Nonzeros are being added, so ensure that non-null data are
    // being passed
    assert(new_matrix_start != NULL);
    assert(new_matrix_index != NULL);
    assert(new_matrix_value != NULL);
  }
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->num_nz();
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && num_row <= 0) return HighsStatus::kError;

  if (this->format_ == MatrixFormat::kRowwise) {
    // Matrix is currently a standard row-wise matrix, so flip
    // column-wise if there are more new nonzeros than current
    // nonzeros
    if (num_new_nz > num_nz && num_nz == 0) {
      // ToDo remove num_nz == 0 when only using a_matrix
      this->setFormat(MatrixFormat::kColwise);
    }
  }
  // Determine the new number of columns and nonzeros in the matrix
  HighsInt new_num_col = num_col + num_new_col;
  HighsInt new_num_nz = num_nz + num_new_nz;

  if (this->isColwise()) {
    // Matrix is column-wise
    this->start_.resize(new_num_col + 1);
    // Append the starts of the new columns
    if (num_new_nz) {
      // Nontrivial number of nonzeros being added, so use new_matrix_start
      for (HighsInt iNewCol = 0; iNewCol < num_new_col; iNewCol++)
        this->start_[num_col + iNewCol] = num_nz + new_matrix_start[iNewCol];
    } else {
      // No nonzeros being added, so new_matrix_start may be null, but entries
      // of zero are implied.
      for (HighsInt iNewCol = 0; iNewCol < num_new_col; iNewCol++)
        this->start_[num_col + iNewCol] = num_nz;
    }
    this->start_[num_col + num_new_col] = new_num_nz;
    // Update the number of columns
    this->num_col_ += num_new_col;
    // If no nonzeros are being added then there's nothing else to do
    if (num_new_nz <= 0) return HighsStatus::kOk;
    // Adding a non-trivial matrix: resize the column-wise matrix arrays
    // accordingly
    this->index_.resize(new_num_nz);
    this->value_.resize(new_num_nz);
    // Copy in the new indices and values
    for (HighsInt iEl = 0; iEl < num_new_nz; iEl++) {
      this->index_[num_nz + iEl] = new_matrix_index[iEl];
      this->value_[num_nz + iEl] = new_matrix_value[iEl];
    }
  } else {
    // Matrix is row-wise
    assert(1 == 0);
  }
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::addRows(const HighsInt num_new_row,
                                       const HighsInt num_new_nz,
                                       const HighsInt* new_matrix_start,
                                       const HighsInt* new_matrix_index,
                                       const double* new_matrix_value,
                                       const int8_t* in_partition) {
  assert(this->format_ != MatrixFormat::kNone);
  // Adding rows to a row-wise partitioned matrix needs the
  // partition information
  const bool partitioned = this->format_ == MatrixFormat::kRowwisePartitioned;
  if (partitioned) {
    assert(1 == 0);
    assert(in_partition != NULL);
  }
  if (num_new_row < 0) return HighsStatus::kError;
  if (num_new_nz < 0) return HighsStatus::kError;
  if (num_new_row == 0) {
    // No rows are being added, so check that no nonzeros are being
    // added
    assert(num_new_nz == 0);
    if (num_new_nz != 0) return HighsStatus::kError;
    return HighsStatus::kOk;
  }
  // Adding a positive number of rows to a matrix
  if (num_new_nz) {
    // Nonzeros are being added, so ensure that non-null data are
    // being passed
    assert(new_matrix_start != NULL);
    assert(new_matrix_index != NULL);
    assert(new_matrix_value != NULL);
  }
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  HighsInt num_nz = this->num_nz();
  // Check that nonzeros aren't being appended to a matrix with no columns
  if (num_new_nz > 0 && num_col <= 0) return HighsStatus::kError;

  if (this->format_ == MatrixFormat::kColwise) {
    // Matrix is currently a standard col-wise matrix, so flip
    // row-wise if there are more new nonzeros than current nonzeros
    if (num_new_nz > num_nz && num_nz == 0) {
      // ToDo remove num_nz == 0 when only using a_matrix
      this->setFormat(MatrixFormat::kRowwise);
    }
  }
  // Determine the new number of rows and nonzeros in the matrix
  HighsInt new_num_nz = num_nz + num_new_nz;
  HighsInt new_num_row = num_row + num_new_row;

  if (this->isRowwise()) {
    // Matrix is row-wise
    this->start_.resize(new_num_row + 1);
    // Append the starts of the new rows
    if (num_new_nz) {
      // Nontrivial number of nonzeros being added, so use new_matrix_start
      assert(new_matrix_start != NULL);
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
        this->start_[num_row + iNewRow] = num_nz + new_matrix_start[iNewRow];
    } else {
      // No nonzeros being added, so new_matrix_start may be NULL, but entries
      // of zero are implied.
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
        this->start_[num_row + iNewRow] = num_nz;
    }
    this->start_[new_num_row] = new_num_nz;
    if (num_new_nz > 0) {
      // Adding a non-trivial matrix: resize the matrix arrays accordingly
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Copy in the new indices and values
      if (partitioned) {
        // Insert the entries in the partition
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          for (HighsInt iNewEl = new_matrix_start[iNewRow];
               iNewEl < new_matrix_start[iNewRow + 1]; iNewEl++) {
            HighsInt iCol = new_matrix_index[iNewEl];
            if (in_partition[iCol]) {
              HighsInt iEl = this->start_[iRow];
              this->index_[iEl] = new_matrix_index[iNewEl];
              this->value_[iEl] = new_matrix_value[iNewEl];
              this->start_[iRow]++;
            }
          }
        }
        // Use the incremented starts to initialise p_end, save these
        // values and reset the starts
        vector<HighsInt> save_p_end;
        save_p_end.resize(num_new_row);
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          this->start_[iRow] = num_nz + new_matrix_start[iNewRow];
          this->p_end_[iRow] = this->start_[iRow];
          save_p_end[iNewRow] = this->p_end_[iRow];
        }
        // Insert the entries not in the partition
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
          HighsInt iRow = num_row + iNewRow;
          for (HighsInt iNewEl = new_matrix_start[iNewRow];
               iNewEl < new_matrix_start[iNewRow + 1]; iNewEl++) {
            HighsInt iCol = new_matrix_index[iNewEl];
            if (!in_partition[iCol]) {
              HighsInt iEl = this->p_end_[iRow];
              this->index_[iEl] = new_matrix_index[iNewEl];
              this->value_[iEl] = new_matrix_value[iNewEl];
              this->p_end_[iRow]++;
            }
          }
        }
        // Reset p_end using the saved values
        for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++)
          this->p_end_[num_row + iNewRow] = save_p_end[iNewRow];
      } else {
        for (HighsInt iNewEl = 0; iNewEl < num_new_nz; iNewEl++) {
          this->index_[num_nz + iNewEl] = new_matrix_index[iNewEl];
          this->value_[num_nz + iNewEl] = new_matrix_value[iNewEl];
        }
      }
    }
  } else {
    // Storing the matrix column-wise, so have to insert the new rows
    assert(this->format_ == MatrixFormat::kColwise);
    if (num_new_nz) {
      vector<HighsInt> length;
      length.assign(num_col, 0);
      for (HighsInt iEl = 0; iEl < num_new_nz; iEl++)
        length[new_matrix_index[iEl]]++;
      // Determine the new number of nonzeros and resize the column-wise matrix
      // arrays
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Append the new rows
      // Shift the existing columns to make space for the new entries
      HighsInt new_iEl = new_num_nz;
      for (HighsInt iCol = num_col - 1; iCol >= 0; iCol--) {
        HighsInt start_col_plus_1 = new_iEl;
        new_iEl -= length[iCol];
        for (HighsInt iEl = this->start_[iCol + 1] - 1;
             iEl >= this->start_[iCol]; iEl--) {
          new_iEl--;
          this->index_[new_iEl] = this->index_[iEl];
          this->value_[new_iEl] = this->value_[iEl];
        }
        this->start_[iCol + 1] = start_col_plus_1;
      }
      assert(new_iEl == 0);
      // Insert the new entries
      for (HighsInt iNewRow = 0; iNewRow < num_new_row; iNewRow++) {
        HighsInt first_el = new_matrix_start[iNewRow];
        HighsInt last_el =
            (iNewRow < num_new_row - 1 ? new_matrix_start[iNewRow + 1]
                                       : num_new_nz);
        for (HighsInt iEl = first_el; iEl < last_el; iEl++) {
          HighsInt iCol = new_matrix_index[iEl];
          new_iEl = this->start_[iCol + 1] - length[iCol];
          length[iCol]--;
          this->index_[new_iEl] = num_row + iNewRow;
          this->value_[new_iEl] = new_matrix_value[iEl];
        }
      }
    }
  }
  // Update the number of rows
  this->num_row_ += num_new_row;

  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::deleteCols(
    const HighsLogOptions& log_options,
    const HighsIndexCollection& index_collection) {
  assert(this->format_ != MatrixFormat::kNone);
  // Can't handle rowwise matrices yet
  assert(!this->isRowwise());
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0,
                         this->num_col_ - 1, true))
      return HighsStatus::kError;
  }
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = this->num_col_;
  HighsInt new_num_col = 0;
  HighsInt new_num_nz = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      new_num_col = delete_from_col;
      new_num_nz = this->start_[delete_from_col];
    }
    // Ensure that the starts of the deleted columns are zeroed to
    // avoid redundant start information for columns whose indices
    // are't used after the deletion takes place. In particular, if
    // all columns are deleted then something must be done to ensure
    // that the matrix isn't magially recreated by increasing the
    // number of columns from zero when there are no rows in the
    // matrix.
    for (HighsInt col = delete_from_col; col <= delete_to_col; col++)
      this->start_[col] = 0;
    // Shift the starts - both in place and value - to account for the
    // columns and nonzeros removed
    const HighsInt keep_from_el = this->start_[keep_from_col];
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      this->start_[new_num_col] = new_num_nz + this->start_[col] - keep_from_el;
      new_num_col++;
    }
    for (HighsInt el = keep_from_el; el < this->start_[keep_to_col + 1]; el++) {
      this->index_[new_num_nz] = this->index_[el];
      this->value_[new_num_nz] = this->value_[el];
      new_num_nz++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  // Ensure that the start of the spurious last column is zeroed so
  // that it doesn't give a positive number of matrix entries if the
  // number of columns in the matrix is increased when there are no
  // rows in the matrix.
  this->start_[this->num_col_] = 0;
  this->start_[new_num_col] = new_num_nz;
  this->start_.resize(new_num_col + 1);
  this->index_.resize(new_num_nz);
  this->value_.resize(new_num_nz);
  // Update the number of columns
  this->num_col_ = new_num_col;
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::deleteRows(
    const HighsLogOptions& log_options,
    const HighsIndexCollection& index_collection) {
  assert(this->format_ != MatrixFormat::kNone);
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0,
                         this->num_row_ - 1, true))
      return HighsStatus::kError;
  }
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt row_dim = this->num_row_;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  // Set up a row mask to indicate the new row index of kept rows and
  // -1 for deleted rows so that the kept entries in the column-wise
  // matrix can be identified and have their correct row index.
  vector<HighsInt> new_index;
  new_index.resize(this->num_row_);
  HighsInt new_num_row = 0;
  bool mask = index_collection.is_mask_;
  const HighsInt* row_mask = index_collection.mask_;
  if (!mask) {
    keep_to_row = -1;
    current_set_entry = 0;
    for (HighsInt k = from_k; k <= to_k; k++) {
      updateIndexCollectionOutInIndex(index_collection, delete_from_row,
                                      delete_to_row, keep_from_row, keep_to_row,
                                      current_set_entry);
      if (k == from_k) {
        // Account for any initial rows being kept
        for (HighsInt row = 0; row < delete_from_row; row++) {
          new_index[row] = new_num_row;
          new_num_row++;
        }
      }
      for (HighsInt row = delete_from_row; row <= delete_to_row; row++) {
        new_index[row] = -1;
      }
      for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
        new_index[row] = new_num_row;
        new_num_row++;
      }
      if (keep_to_row >= row_dim - 1) break;
    }
  } else {
    for (HighsInt row = 0; row < this->num_row_; row++) {
      if (row_mask[row]) {
        new_index[row] = -1;
      } else {
        new_index[row] = new_num_row;
        new_num_row++;
      }
    }
  }
  HighsInt new_num_nz = 0;
  for (HighsInt col = 0; col < this->num_col_; col++) {
    HighsInt from_el = this->start_[col];
    this->start_[col] = new_num_nz;
    for (HighsInt el = from_el; el < this->start_[col + 1]; el++) {
      HighsInt row = this->index_[el];
      HighsInt new_row = new_index[row];
      if (new_row >= 0) {
        this->index_[new_num_nz] = new_row;
        this->value_[new_num_nz] = this->value_[el];
        new_num_nz++;
      }
    }
  }
  this->start_[this->num_col_] = new_num_nz;
  this->start_.resize(this->num_col_ + 1);
  this->index_.resize(new_num_nz);
  this->value_.resize(new_num_nz);
  // Update the number of rows
  this->num_row_ = new_num_row;
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::assessDimensions(
    const HighsLogOptions& log_options, const std::string matrix_name) {
  if (this->format_ == MatrixFormat::kNone) return HighsStatus::kOk;
  // Use error_found to track whether an error has been found in multiple tests
  bool error_found = false;
  // Identify main dimensions
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);
  HighsInt num_vec;
  if (this->format_ == MatrixFormat::kColwise) {
    num_vec = this->num_col_;
  } else {
    num_vec = this->num_row_;
  }
  vector<HighsInt>& matrix_start = this->start_;
  vector<HighsInt>& matrix_p_end = this->p_end_;
  vector<HighsInt>& matrix_index = this->index_;
  vector<double>& matrix_value = this->value_;
  // Assess main dimensions
  bool legal_num_vec = num_vec >= 0;
  if (!legal_num_vec) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix has illegal number of vectors = %" HIGHSINT_FORMAT
                 "\n",
                 matrix_name.c_str(), num_vec);
    error_found = true;
  }
  HighsInt matrix_start_size = matrix_start.size();
  bool legal_matrix_start_size = matrix_start_size >= num_vec + 1;
  if (!legal_matrix_start_size) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix has illegal start vector size = %" HIGHSINT_FORMAT
                 " < %" HIGHSINT_FORMAT "\n",
                 matrix_name.c_str(), matrix_start_size, num_vec + 1);
    error_found = true;
  }
  // Check whether the first start is zero
  if (matrix_start[0]) {
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s matrix start vector begins with %" HIGHSINT_FORMAT
                 " rather than 0\n",
                 matrix_name.c_str(), matrix_start[0]);
    error_found = true;
  }
  if (this->format_ == MatrixFormat::kRowwisePartitioned) {
    HighsInt matrix_p_end_size = matrix_p_end.size();
    bool legal_matrix_p_end_size = matrix_p_end_size >= num_vec + 1;
    if (!legal_matrix_p_end_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal p_end vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_p_end_size, num_vec + 1);
      error_found = true;
    }
  }
  // Possibly check the sizes of the index and value vectors. Can only
  // do this with the number of nonzeros, and this is only known if
  // the start vector has a legal size. Setting num_nz = 0 otherwise
  // means that all tests pass, as they just check that the sizes of
  // the index and value vectors are non-negative.
  HighsInt num_nz = 0;
  if (legal_matrix_start_size) num_nz = matrix_start[num_vec];
  bool legal_num_nz = num_nz >= 0;
  if (!legal_num_nz) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix has illegal number of nonzeros = %" HIGHSINT_FORMAT
                 "\n",
                 matrix_name.c_str(), num_nz);
    error_found = true;
  } else {
    HighsInt matrix_index_size = matrix_index.size();
    HighsInt matrix_value_size = matrix_value.size();
    bool legal_matrix_index_size = matrix_index_size >= num_nz;
    bool legal_matrix_value_size = matrix_value_size >= num_nz;
    if (!legal_matrix_index_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal index vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_index_size, num_nz);
      error_found = true;
    }
    if (!legal_matrix_value_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal value vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_value_size, num_nz);
      error_found = true;
    }
  }
  HighsStatus return_status = HighsStatus::kOk;
  if (error_found) return_status = HighsStatus::kError;
  return return_status;
}

HighsStatus HighsSparseMatrix::assess(const HighsLogOptions& log_options,
                                      const std::string matrix_name,
                                      const double small_matrix_value,
                                      const double large_matrix_value) {
  if (this->format_ == MatrixFormat::kNone) return HighsStatus::kOk;
  if (this->assessDimensions(log_options, matrix_name) == HighsStatus::kError)
    return HighsStatus::kError;
  // Identify main dimensions
  HighsInt vec_dim;
  HighsInt num_vec;
  if (this->isColwise()) {
    vec_dim = this->num_row_;
    num_vec = this->num_col_;
  } else {
    vec_dim = this->num_col_;
    num_vec = this->num_row_;
  }
  const bool partitioned = this->format_ == MatrixFormat::kRowwisePartitioned;
  vector<HighsInt>& matrix_start = this->start_;
  vector<HighsInt>& matrix_p_end = this->p_end_;
  vector<HighsInt>& matrix_index = this->index_;
  vector<double>& matrix_value = this->value_;

  const HighsInt num_nz = matrix_start[num_vec];

  bool error_found = false;
  bool warning_found = false;

  // Assess the starts
  // Set up previous_start for a fictitious previous empty packed vector
  HighsInt previous_start = matrix_start[0];
  // Set up this_start to be the first start in case num_vec = 0
  HighsInt this_start = matrix_start[0];
  HighsInt this_p_end;
  if (partitioned) this_p_end = matrix_p_end[0];
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    this_start = matrix_start[ix];
    HighsInt next_start = matrix_start[ix + 1];
    bool this_start_too_small = this_start < previous_start;
    if (this_start_too_small) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal start of %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT
                   " = "
                   "previous start\n",
                   matrix_name.c_str(), ix, this_start, previous_start);
      return HighsStatus::kError;
    }
    if (partitioned) {
      this_p_end = matrix_p_end[ix];
      bool this_p_end_too_small = this_p_end < this_start;
      if (this_p_end_too_small) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     " has illegal partition end of %" HIGHSINT_FORMAT
                     " < %" HIGHSINT_FORMAT
                     " = "
                     " start\n",
                     matrix_name.c_str(), ix, this_p_end, this_start);
        return HighsStatus::kError;
      }
    }
    previous_start = this_start;
  }
  bool this_start_too_big = this_start > num_nz;
  if (this_start_too_big) {
    highsLogUser(log_options, HighsLogType::kError,
                 "%s matrix packed vector %" HIGHSINT_FORMAT
                 " has illegal start of %" HIGHSINT_FORMAT
                 " > %" HIGHSINT_FORMAT
                 " = "
                 "number of nonzeros\n",
                 matrix_name.c_str(), num_vec, this_start, num_nz);
    return HighsStatus::kError;
  }
  if (partitioned) {
    bool this_p_end_too_big = this_p_end > num_nz;
    if (this_p_end_too_big) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal partition end of %" HIGHSINT_FORMAT
                   " > %" HIGHSINT_FORMAT
                   " = "
                   "number of nonzeros\n",
                   matrix_name.c_str(), num_vec, this_p_end, num_nz);
      return HighsStatus::kError;
    }
  }
  // Assess the indices and values
  // Count the number of acceptable indices/values
  HighsInt num_new_nz = 0;
  HighsInt num_small_values = 0;
  double max_small_value = 0;
  double min_small_value = kHighsInf;
  // Set up a zeroed vector to detect duplicate indices
  vector<HighsInt> check_vector;
  if (vec_dim > 0) check_vector.assign(vec_dim, 0);
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    HighsInt from_el = matrix_start[ix];
    HighsInt to_el = matrix_start[ix + 1];
    // Account for any index-value pairs removed so far
    matrix_start[ix] = num_new_nz;
    for (HighsInt el = from_el; el < to_el; el++) {
      // Check the index
      HighsInt component = matrix_index[el];
      // Check that the index is non-negative
      bool legal_component = component >= 0;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is illegal index %" HIGHSINT_FORMAT "\n",
                     matrix_name.c_str(), ix, el, component);
        return HighsStatus::kError;
      }
      // Check that the index does not exceed the vector dimension
      legal_component = component < vec_dim;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is illegal index "
                     "%12" HIGHSINT_FORMAT " >= %" HIGHSINT_FORMAT
                     " = vector dimension\n",
                     matrix_name.c_str(), ix, el, component, vec_dim);
        return HighsStatus::kError;
      }
      // Check that the index has not already ocurred
      legal_component = check_vector[component] == 0;
      if (!legal_component) {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s matrix packed vector %" HIGHSINT_FORMAT
                     ", entry %" HIGHSINT_FORMAT
                     ", is duplicate index %" HIGHSINT_FORMAT "\n",
                     matrix_name.c_str(), ix, el, component);
        return HighsStatus::kError;
      }
      // Indicate that the index has occurred
      check_vector[component] = 1;
      // Check the value
      double abs_value = fabs(matrix_value[el]);
      /*
      // Check that the value is not zero
      bool zero_value = abs_value == 0;
      if (zero_value) {
        highsLogUser(log_options, HighsLogType::kError,
                        "%s matrix packed vector %" HIGHSINT_FORMAT ", entry %"
      HIGHSINT_FORMAT ", is zero\n", matrix_name.c_str(),  ix, el); return
      HighsStatus::kError;
      }
      */
      // Check that the value is not too large
      bool large_value = abs_value > large_matrix_value;
      if (large_value) {
        highsLogUser(
            log_options, HighsLogType::kError,
            "%s matrix packed vector %" HIGHSINT_FORMAT
            ", entry %" HIGHSINT_FORMAT ", is large value |%g| >= %g\n",
            matrix_name.c_str(), ix, el, abs_value, large_matrix_value);
        return HighsStatus::kError;
      }
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
        if (max_small_value < abs_value) max_small_value = abs_value;
        if (min_small_value > abs_value) min_small_value = abs_value;
        num_small_values++;
      }
      if (ok_value) {
        // Shift the index and value of the OK entry to the new
        // position in the index and value vectors, and increment
        // the new number of nonzeros
        matrix_index[num_new_nz] = matrix_index[el];
        matrix_value[num_new_nz] = matrix_value[el];
        num_new_nz++;
      } else {
        // Zero the check_vector entry since the small value
        // _hasn't_ occurred
        check_vector[component] = 0;
      }
    }
    // Zero check_vector
    for (HighsInt el = matrix_start[ix]; el < num_new_nz; el++)
      check_vector[matrix_index[el]] = 0;
#ifdef HiGHSDEV
    // NB This is very expensive so shouldn't be true
    const bool check_check_vector = false;
    if (check_check_vector) {
      // Check zeroing of check vector
      for (HighsInt component = 0; component < vec_dim; component++) {
        if (check_vector[component]) error_found = true;
      }
      if (error_found)
        highsLogUser(log_options, HighsLogType::kError,
                     "assessMatrix: check_vector not zeroed\n");
    }
#endif
  }
  if (num_small_values) {
    if (partitioned) {
      // Shouldn't happen with a parttioned row-wise matrix since its
      // values should be OK and the code above doesn't handle p_end
      highsLogUser(
          log_options, HighsLogType::kError,
          "%s matrix packed partitioned vector contains %" HIGHSINT_FORMAT
          " |values| in [%g, %g] "
          "less than %g: ignored\n",
          matrix_name.c_str(), num_small_values, min_small_value,
          max_small_value, small_matrix_value);
      error_found = true;
      assert(num_small_values = 0);
    }
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s matrix packed vector contains %" HIGHSINT_FORMAT
                 " |values| in [%g, %g] "
                 "less than %g: ignored\n",
                 matrix_name.c_str(), num_small_values, min_small_value,
                 max_small_value, small_matrix_value);
    warning_found = true;
  }
  matrix_start[num_vec] = num_new_nz;
  HighsStatus return_status = HighsStatus::kOk;
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  return return_status;
}

void HighsSparseMatrix::scaleCol(const HighsInt col, const double colScale) {
  assert(this->format_ != MatrixFormat::kNone);
  assert(col >= 0);
  assert(col < this->num_col_);
  assert(colScale);

  if (this->isColwise()) {
    for (HighsInt iEl = this->start_[col]; iEl < this->start_[col + 1]; iEl++)
      this->value_[iEl] *= colScale;
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        if (this->index_[iEl] == col) this->value_[iEl] *= colScale;
      }
    }
  }
}

void HighsSparseMatrix::scaleRow(const HighsInt row, const double rowScale) {
  assert(this->format_ != MatrixFormat::kNone);
  assert(row >= 0);
  assert(row < this->num_row_);
  assert(rowScale);

  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        if (this->index_[iEl] == row) this->value_[iEl] *= rowScale;
      }
    }
  } else {
    for (HighsInt iEl = this->start_[row]; iEl < this->start_[row + 1]; iEl++)
      this->value_[iEl] *= rowScale;
  }
}

void HighsSparseMatrix::applyScale(const SimplexScale& scale) {
  assert(this->format_ != MatrixFormat::kNone);
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        HighsInt iRow = this->index_[iEl];
        this->value_[iEl] *= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        HighsInt iCol = this->index_[iEl];
        this->value_[iEl] *= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  }
}

void HighsSparseMatrix::unapplyScale(const SimplexScale& scale) {
  assert(this->format_ != MatrixFormat::kNone);
  if (this->isColwise()) {
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      for (HighsInt iEl = this->start_[iCol]; iEl < this->start_[iCol + 1];
           iEl++) {
        HighsInt iRow = this->index_[iEl];
        this->value_[iEl] /= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      for (HighsInt iEl = this->start_[iRow]; iEl < this->start_[iRow + 1];
           iEl++) {
        HighsInt iCol = this->index_[iEl];
        this->value_[iEl] /= (scale.col[iCol] * scale.row[iRow]);
      }
    }
  }
}

void HighsSparseMatrix::createPartition(const HighsSparseMatrix& matrix,
                                        const int8_t* in_partition) {
  assert(matrix.format_ != MatrixFormat::kNone);
  assert(this->format_ != MatrixFormat::kNone);
}

void HighsSparseMatrix::priceByColumn(HVector& result,
                                      const HVector& vector) const {
  assert(this->isColwise());
}

void HighsSparseMatrix::priceByRow(HVector& result,
                                   const HVector& vector) const {
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
}

void HighsSparseMatrix::priceByRowWithSwitch(
    HVector& result, const HVector& vector, const double expected_density,
    const HighsInt from_row, const double switch_density) const {
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
}

void HighsSparseMatrix::update(const HighsInt var_in, const HighsInt var_out) {
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
}

double HighsSparseMatrix::computeDot(const HVector& vector,
                                     const HighsInt use_col) const {
  assert(this->format_ != MatrixFormat::kNone);
  return 0.0;
}

void HighsSparseMatrix::collectAj(HVector& vector, const HighsInt use_col,
                                  const double multiplier) const {
  assert(this->format_ != MatrixFormat::kNone);
}

void HighsSparseMatrix::priceByRowDenseResult(HVector& result,
                                              const HVector& vector,
                                              const HighsInt from_row) {
  assert(this->format_ == MatrixFormat::kRowwisePartitioned);
}
