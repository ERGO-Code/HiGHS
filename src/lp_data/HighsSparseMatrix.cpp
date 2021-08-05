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

void HighsSparseMatrix::initialise() {
  this->clear();
  //  this->format_ = MatrixFormat::kColwise;
  //  this->start_.assign(1, 0);
}

void HighsSparseMatrix::clear() {
  this->format_ = MatrixFormat::kNone;
  this->num_col_ = 0;
  this->num_row_ = 0;
  this->start_.clear();
  this->p_end_.clear();
  this->index_.clear();
  this->value_.clear();
}

void HighsSparseMatrix::range(double& min_value, double& max_value) const {
  for (HighsInt iEl = 0; iEl < this->start_[this->num_col_]; iEl++) {
    double value = fabs(this->value_[iEl]);
    min_value = min(min_value, value);
    max_value = max(max_value, value);
  }
}

HighsStatus HighsSparseMatrix::setFormat(const MatrixFormat desired_format) {
  if (desired_format == MatrixFormat::kNone) return HighsStatus::kError;
  if (this->format_ == desired_format) return HighsStatus::kOk;
  if (this->num_col_ == 0 && this->num_row_ == 0) {
    // No rows or columns, so either orientation is possible and has
    // identical data: just requires the start of the fictitious
    // row/column 0
    this->start_.assign(1, 0);
    this->format_ = desired_format;
  } else {
    // Any LP with positive numbers of rows or columns must have an orientation
    assert(this->format_ != MatrixFormat::kNone);
    if (desired_format == MatrixFormat::kColwise) {
      this->ensureColWise();
    } else {
      this->ensureRowWise();
    }
  }
  assert(this->format_ == desired_format);
  return HighsStatus::kOk;
}

void HighsSparseMatrix::ensureColWise() {
  // Should only call this is orientation is ROWWISE
  assert(this->format_ == MatrixFormat::kRowwise);
  HighsInt num_nz;
  bool empty_matrix = this->num_col_ == 0 || this->num_row_ == 0;
  if (!empty_matrix) {
    // Matrix is probably non-empty
    assert((HighsInt)this->start_.size() >= this->num_row_ + 1);
    num_nz = this->start_[this->num_row_];
    assert(num_nz >= 0);
    assert((HighsInt)this->index_.size() >= num_nz);
    assert((HighsInt)this->value_.size() >= num_nz);
    empty_matrix = num_nz == 0;
    if (!empty_matrix) {
      // Matrix is non-empty, so transpose it
      //
      // Take a copy of the current matrix - that is rowwise - so that
      // the current matrix is filled colwise
      vector<HighsInt> ARstart = this->start_;
      vector<HighsInt> ARindex = this->index_;
      vector<double> ARvalue = this->value_;
      this->start_.resize(this->num_col_ + 1);
      this->index_.resize(num_nz);
      this->value_.resize(num_nz);
      vector<HighsInt> Alength;
      Alength.assign(this->num_col_, 0);
      for (HighsInt iEl = ARstart[0]; iEl < num_nz; iEl++)
        Alength[ARindex[iEl]]++;
      this->start_[0] = 0;
      for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
        this->start_[iCol + 1] = this->start_[iCol] + Alength[iCol];
      for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
        for (HighsInt iEl = ARstart[iRow]; iEl < ARstart[iRow + 1]; iEl++) {
          HighsInt iCol = ARindex[iEl];
          HighsInt iCol_el = this->start_[iCol];
          this->index_[iCol_el] = iRow;
          this->value_[iCol_el] = ARvalue[iEl];
          this->start_[iCol]++;
        }
      }
      this->start_[0] = 0;
      for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
        this->start_[iCol + 1] = this->start_[iCol] + Alength[iCol];
      assert(this->start_[this->num_col_] == num_nz);
    }
  }
  if (empty_matrix) {
    // Matrix is empty, so set up empty column-wise structure
    this->start_.assign(this->num_col_ + 1, 0);
    this->index_.clear();
    this->value_.clear();
  }
  assert((HighsInt)this->start_.size() >= this->num_col_ + 1);
  num_nz = this->start_[this->num_col_];
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
  this->format_ = MatrixFormat::kColwise;
}

void HighsSparseMatrix::ensureRowWise() {
  // Should only call this is orientation is COLWISE
  assert(this->format_ == MatrixFormat::kColwise);
  HighsInt num_nz;
  bool empty_matrix = this->num_col_ == 0 || this->num_row_ == 0;
  if (!empty_matrix) {
    // Matrix is probably non-empty
    assert((HighsInt)this->start_.size() >= this->num_col_ + 1);
    num_nz = this->start_[this->num_col_];
    assert(num_nz >= 0);
    assert((HighsInt)this->index_.size() >= num_nz);
    assert((HighsInt)this->value_.size() >= num_nz);
    empty_matrix = num_nz == 0;
    if (!empty_matrix) {
      // Matrix is non-empty, so transpose it
      //
      // Take a copy of the current matrix - that is colwise - so that
      // the current matrix is filled rowwise
      vector<HighsInt> Astart = this->start_;
      vector<HighsInt> Aindex = this->index_;
      vector<double> Avalue = this->value_;
      this->start_.resize(this->num_row_ + 1);
      this->index_.resize(num_nz);
      this->value_.resize(num_nz);
      vector<HighsInt> ARlength;
      ARlength.assign(this->num_row_, 0);
      for (HighsInt iEl = Astart[0]; iEl < num_nz; iEl++)
        ARlength[Aindex[iEl]]++;
      this->start_[0] = 0;
      for (HighsInt iRow = 0; iRow < this->num_row_; iRow++)
        this->start_[iRow + 1] = this->start_[iRow] + ARlength[iRow];
      for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
        for (HighsInt iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
          HighsInt iRow = Aindex[iEl];
          HighsInt iRow_el = this->start_[iRow];
          this->index_[iRow_el] = iCol;
          this->value_[iRow_el] = Avalue[iEl];
          this->start_[iRow]++;
        }
      }
      this->start_[0] = 0;
      for (HighsInt iRow = 0; iRow < this->num_row_; iRow++)
        this->start_[iRow + 1] = this->start_[iRow] + ARlength[iRow];
      assert(this->start_[this->num_row_] == num_nz);
    }
  }
  if (empty_matrix) {
    // Matrix is empty, so set up empty row-wise structure
    this->start_.assign(this->num_row_ + 1, 0);
    this->index_.clear();
    this->value_.clear();
  }
  assert((HighsInt)this->start_.size() >= this->num_row_ + 1);
  num_nz = this->start_[this->num_row_];
  assert(num_nz >= 0);
  assert((HighsInt)this->index_.size() >= num_nz);
  assert((HighsInt)this->value_.size() >= num_nz);
  this->format_ = MatrixFormat::kRowwise;
}

HighsStatus HighsSparseMatrix::addCols(const HighsInt num_new_col,
                                       const HighsInt num_new_nz,
                                       const HighsInt* new_matrix_start,
                                       const HighsInt* new_matrix_index,
                                       const double* new_matrix_value) {
  if (num_new_col < 0) return HighsStatus::kError;
  if (num_new_col == 0) return HighsStatus::kOk;
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && num_row <= 0) return HighsStatus::kError;
  // Adding a positive number of columns to a matrix
  if (this->format_ == MatrixFormat::kNone) {
    // LP is currently empty, store the matrix column-wise
    assert(num_col == 0 && num_row == 0);
    this->format_ = MatrixFormat::kColwise;
  } else {
    // Ensure that the matrix is stored column-wise
    this->setFormat();
  }
  // Determine the new number of columns in the matrix and resize the
  // starts accordingly.
  HighsInt new_num_col = num_col + num_new_col;
  this->start_.resize(new_num_col + 1);
  // If adding columns to an empty LP then introduce the start for the
  // fictitious column 0
  if (num_col == 0) this->start_[0] = 0;

  // Determine the current number of nonzeros and the new number of nonzeros
  HighsInt current_num_nz = this->start_[num_col];
  HighsInt new_num_nz = current_num_nz + num_new_nz;

  // Append the starts of the new columns
  if (num_new_nz) {
    // Nontrivial number of nonzeros being added, so use new_matrix_start
    assert(new_matrix_start != NULL);
    for (HighsInt iCol = 0; iCol < num_new_col; iCol++)
      this->start_[num_col + iCol] = current_num_nz + new_matrix_start[iCol];
  } else {
    // No nonzeros being added, so new_matrix_start may be null, but entries of
    // zero are implied.
    for (HighsInt iCol = 0; iCol < num_new_col; iCol++)
      this->start_[num_col + iCol] = current_num_nz;
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
    this->index_[current_num_nz + iEl] = new_matrix_index[iEl];
    this->value_[current_num_nz + iEl] = new_matrix_value[iEl];
  }
  return HighsStatus::kOk;
}

HighsStatus HighsSparseMatrix::addRows(const HighsInt num_new_row,
                                       const HighsInt num_new_nz,
                                       const HighsInt* new_matrix_start,
                                       const HighsInt* new_matrix_index,
                                       const double* new_matrix_value) {
  if (num_new_row < 0) return HighsStatus::kError;
  if (num_new_row == 0) return HighsStatus::kOk;
  HighsInt num_col = this->num_col_;
  HighsInt num_row = this->num_row_;
  // Check that nonzeros aren't being appended to a matrix with no columns
  if (num_new_nz > 0 && num_col <= 0) return HighsStatus::kError;
  // Adding a positive number of rows to a matrix
  HighsInt current_num_nz = 0;
  if (this->format_ == MatrixFormat::kNone) {
    // LP is currently empty, store the matrix row-wise
    assert(num_col == 0 && num_row == 0);
    this->format_ = MatrixFormat::kRowwise;
  } else if (this->format_ == MatrixFormat::kColwise) {
    assert(num_col > 0);
    assert((HighsInt)this->start_.size() >= num_col);
    current_num_nz = this->start_[num_col];
    if (current_num_nz == 0) {
      // Matrix is currently empty and stored column-wise. It can be
      // converted trivially to row-wise storage so that rows can be
      // added easily.
      //
      // It's possible that the model could have columns and (empty)
      // rows - hence the assignment of zero starts for rows
      // 0...num_row.
      //
      // However, this allows efficient handling of the (common) case
      // where a modeller defines variables without constraints, and
      // then constraints one-by-one.
      this->format_ = MatrixFormat::kRowwise;
      this->start_.assign(num_row + 1, 0);
    }
  }
  if (this->format_ == MatrixFormat::kRowwise) {
    HighsInt new_num_row = num_row + num_new_row;
    this->start_.resize(new_num_row + 1);
    // If adding rows to an empty matrix then introduce the start for the
    // fictitious row 0
    if (num_row == 0) this->start_[0] = 0;

    // Determine the current number of nonzeros and the new number of nonzeros
    HighsInt current_num_nz = this->start_[num_row];
    HighsInt new_num_nz = current_num_nz + num_new_nz;

    // Append the starts of the new rows
    if (num_new_nz) {
      // Nontrivial number of nonzeros being added, so use new_matrix_start
      assert(new_matrix_start != NULL);
      for (HighsInt iRow = 0; iRow < num_new_row; iRow++)
        this->start_[num_row + iRow] = current_num_nz + new_matrix_start[iRow];
    } else {
      // No nonzeros being added, so new_matrix_start may be null, but entries
      // of zero are implied.
      for (HighsInt iRow = 0; iRow < num_new_row; iRow++)
        this->start_[num_row + iRow] = current_num_nz;
    }
    this->start_[num_row + num_new_row] = new_num_nz;

    if (num_new_nz > 0) {
      // Adding a non-trivial matrix: resize the matrix arrays accordingly
      this->index_.resize(new_num_nz);
      this->value_.resize(new_num_nz);
      // Copy in the new indices and values
      for (HighsInt el = 0; el < num_new_nz; el++) {
        this->index_[current_num_nz + el] = new_matrix_index[el];
        this->value_[current_num_nz + el] = new_matrix_value[el];
      }
    }
  } else {
    // Storing the matrix column-wise, so have to insert the new rows
    assert(this->format_ == MatrixFormat::kColwise);
    vector<HighsInt> Alength;
    Alength.assign(num_col, 0);
    for (HighsInt el = 0; el < num_new_nz; el++)
      Alength[new_matrix_index[el]]++;
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    HighsInt new_num_nz = current_num_nz + num_new_nz;
    this->index_.resize(new_num_nz);
    this->value_.resize(new_num_nz);
    // Append the new rows
    // Shift the existing columns to make space for the new entries
    HighsInt new_el = new_num_nz;
    for (HighsInt col = num_col - 1; col >= 0; col--) {
      HighsInt start_col_plus_1 = new_el;
      new_el -= Alength[col];
      for (HighsInt el = this->start_[col + 1] - 1; el >= this->start_[col];
           el--) {
        new_el--;
        this->index_[new_el] = this->index_[el];
        this->value_[new_el] = this->value_[el];
      }
      this->start_[col + 1] = start_col_plus_1;
    }
    assert(new_el == 0);
    // Insert the new entries
    for (HighsInt row = 0; row < num_new_row; row++) {
      HighsInt first_el = new_matrix_start[row];
      HighsInt last_el =
          (row < num_new_row - 1 ? new_matrix_start[row + 1] : num_new_nz);
      for (HighsInt el = first_el; el < last_el; el++) {
        HighsInt col = new_matrix_index[el];
        new_el = this->start_[col + 1] - Alength[col];
        Alength[col]--;
        this->index_[new_el] = num_row + row;
        this->value_[new_el] = new_matrix_value[el];
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
  HighsStatus return_status = HighsStatus::kOk;
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
  bool legal_matrix_start_size = false;
  // Don't expect the matrix_start_size to be legal if there are no vectors
  if (num_vec > 0) {
    legal_matrix_start_size = matrix_start_size >= num_vec + 1;
    if (!legal_matrix_start_size) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix has illegal start vector size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   matrix_name.c_str(), matrix_start_size, num_vec + 1);
      error_found = true;
    }
  }
  if (matrix_start_size > 0) {
    // Check whether the first start is zero
    if (matrix_start[0]) {
      highsLogUser(log_options, HighsLogType::kWarning,
                   "%s matrix start vector begins with %" HIGHSINT_FORMAT
                   " rather than 0\n",
                   matrix_name.c_str(), matrix_start[0]);
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
  if (error_found)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;
  return return_status;
}

HighsStatus HighsSparseMatrix::assess(const HighsLogOptions& log_options,
                                      const std::string matrix_name,
                                      const double small_matrix_value,
                                      const double large_matrix_value) {
  // Identify main dimensions
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);
  if (this->assessDimensions(log_options, matrix_name) == HighsStatus::kError)
    return HighsStatus::kError;
  // Identify main dimensions
  HighsInt vec_dim;
  HighsInt num_vec;
  if (this->format_ == MatrixFormat::kColwise) {
    vec_dim = this->num_row_;
    num_vec = this->num_col_;
  } else {
    vec_dim = this->num_col_;
    num_vec = this->num_row_;
  }
  vector<HighsInt>& matrix_start = this->start_;
  vector<HighsInt>& matrix_index = this->index_;
  vector<double>& matrix_value = this->value_;

  const HighsInt num_nz = matrix_start[num_vec];
  if (num_vec <= 0) return HighsStatus::kOk;
  if (num_nz <= 0) return HighsStatus::kOk;

  HighsStatus return_status = HighsStatus::kOk;
  bool error_found = false;
  bool warning_found = false;

  // Assess the starts
  // Set up previous_start for a fictitious previous empty packed vector
  HighsInt previous_start = matrix_start[0];
  for (HighsInt ix = 0; ix < num_vec; ix++) {
    HighsInt this_start = matrix_start[ix];
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
    bool this_start_too_big = this_start > num_nz;
    if (this_start_too_big) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%s matrix packed vector %" HIGHSINT_FORMAT
                   " has illegal start of %" HIGHSINT_FORMAT
                   " > %" HIGHSINT_FORMAT
                   " = "
                   "number of nonzeros\n",
                   matrix_name.c_str(), ix, this_start, num_nz);
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
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s matrix packed vector contains %" HIGHSINT_FORMAT
                 " |values| in [%g, %g] "
                 "less than %g: ignored\n",
                 matrix_name.c_str(), num_small_values, min_small_value,
                 max_small_value, small_matrix_value);
    warning_found = true;
  }
  matrix_start[num_vec] = num_new_nz;
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

void HighsSparseMatrix::scaleCol(const HighsInt col, const double colScale) {
  assert(col >= 0);
  assert(col < this->num_col_);
  assert(colScale);
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);

  if (this->format_ == MatrixFormat::kColwise) {
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
  assert(row >= 0);
  assert(row < this->num_row_);
  assert(rowScale);
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);

  if (this->format_ == MatrixFormat::kColwise) {
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
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);
  if (this->format_ == MatrixFormat::kColwise) {
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
  assert(this->format_ == MatrixFormat::kColwise ||
         this->format_ == MatrixFormat::kRowwise);
  if (this->format_ == MatrixFormat::kColwise) {
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
}

void HighsSparseMatrix::priceByColumn(HVector& result,
				      const HVector& vector) const {
}

void HighsSparseMatrix::priceByRow(HVector& result,
				   const HVector& vector) const {
}

void HighsSparseMatrix::priceByRowWithSwitch(HVector& result,
					     const HVector& vector,
					     const double expected_density,
					     const HighsInt from_row,
					     const double switch_density) const {
}

void HighsSparseMatrix::update(const HighsInt var_in,
			       const HighsInt var_out) {
}

double HighsSparseMatrix::computeDot(const HVector& vector,
				     const HighsInt use_col) const {
  return 0.0;
}

void HighsSparseMatrix::collectAj(HVector& vector,
				  const HighsInt use_col,
				  const double multiplier) const {
}

void HighsSparseMatrix::priceByRowDenseResult(HVector& result, const HVector& vector,
					      const HighsInt from_row) {
}


