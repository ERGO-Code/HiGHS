/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsMiqp.cpp
 * @brief
 */
#include "lp_data/HighsMiqp.h"

#include "Highs.h"

HighsStatus Highs::solveMiqp() {
  HighsStatus status = HighsStatus::kOk;
  // Currently trying to solve a MIQP yields an error return
  status = this->optimizeModel();
  assert(status == HighsStatus::kError);

  // A handy const reference to the incumbent HighsLp
  HighsLp& lp = this->model_.lp_;
  // Extract the original number of columns and rows in the incumbent
  // HighsLp
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  //
  // Note that HighsInt is the HiGHS integer: there is a build flag
  // (HIGHSINT64) that enables HiGHS to use 64-bit integers. Hence
  // printing of a HighsInt requires a cast to int if using %d,
  // otherise the PRId64d value can be used
  //
  // Illustrate the use of highsLogUser to perform logging: no other
  // method of permanent logging must be used
  highsLogUser(
      options_.log_options, HighsLogType::kInfo,
      "Highs::solveMiqp Solving an MIQP with %d columns and %" HIGHSINT_FORMAT
      " rows\n",
      int(num_col), num_row);

  // Move the Hessian into a local instance
  HighsHessian hessian = std::move(this->model_.hessian_);
  // Study the HighsHessian class (in highs/model/HighsHessian.h)!

  // Log information about the Hessian: also gives an illustration of
  // the Hessian properties
  assert(hessian.dim() == num_col);
  logHessian(options_.log_options, hessian);

  // Clear the Hessian in the incumbent HighsModel so that it is just
  // a MIP
  this->model_.hessian_.clear();

  // Lambda to find the maximum value in an std::vector<HighsInt>
  auto maxIndex = [&](const std::vector<HighsInt>& index) {
    HighsInt max_index = 0;
    for (HighsInt k = 0; k < int(index.size()); k++)
      max_index = std::max(index[k], max_index);
    return max_index;
  };

  // Linearising the quadratic terms in the Hessian requires rows and
  // (possibly) columns to be added to the incumbent HighsLp
  //
  // Add a column and row to the problem
  std::vector<HighsInt> index = {0, 1, lp.num_row_};
  std::vector<double> value = {1, 1, 1};
  assert(index.size() == value.size());
  HighsInt num_new_index = index.size();

  // Row indices must not exceed lp.num_row_ - 1, and this call
  // illustrates what happens if they do!
  status = this->addCol(1.0, 0, kHighsInf, num_new_index, index.data(),
                        value.data());
  assert(status == HighsStatus::kError);
  // If addCol/addRow returns an error, the incumbent HighsLp is not
  // changed
  assert(lp.num_col_ == num_col);

  // Now give a legal value to index[2]
  index[2] = lp.num_row_ - 1;

  assert(maxIndex(index) < lp.num_row_);
  status = this->addCol(1.0, 0, kHighsInf, num_new_index, index.data(),
                        value.data());
  assert(status == HighsStatus::kOk);
  assert(lp.num_col_ == num_col + 1);

  // The new column is continuous (HighsVarType::kContinuous) by
  // default. This is how to make it an integer variable.
  status = this->changeColIntegrality(1, HighsVarType::kInteger);
  assert(status == HighsStatus::kOk);

  // Now add a row
  index[2] = lp.num_col_ - 1;
  assert(maxIndex(index) < lp.num_col_);
  status =
      this->addRow(0, kHighsInf, num_new_index, index.data(), value.data());
  assert(status == HighsStatus::kOk);
  assert(lp.num_row_ == num_row + 1);

  // Once the linearisations have been added, the MILP is solved by
  // calling Highs::optimizeModel()
  this->reportModelStats();
  const HighsStatus run_status = this->optimizeModel();

  // Delete the columns and rows that have been added
  status = this->deleteCols(num_col, lp.num_col_ - 1);
  assert(status == HighsStatus::kOk);

  status = this->deleteRows(num_row, lp.num_row_ - 1);
  assert(status == HighsStatus::kOk);

  // Restore the Hessian in the incumbent HighsModel
  this->model_.hessian_ = std::move(hessian);
  this->reportModelStats();

  // Check that the columns and rows added have been removed
  assert(lp.num_col_ == num_col);
  assert(lp.num_row_ == num_row);
  // Currently still have to ensure that trying to solve a MIQP yields
  // an error return
  return HighsStatus::kError;  // run_status;
}

void logHessian(const HighsLogOptions log_options,
                const HighsHessian& hessian) {
  // The Hessians received here will be triangular, with the first
  // nonzero in each column being the diagonal entry (possibly with
  // zero value). [Yes, hessian.numNz() is confusing...]
  assert(hessian.format_ == HessianFormat::kTriangular);
  bool diagonal = true;
  HighsInt num_diagonal_nz = 0;
  for (HighsInt iCol = 0; iCol < hessian.dim(); iCol++) {
    HighsInt iEl = hessian.start_[iCol];
    HighsInt iRow = hessian.index_[iEl];
    assert(iRow == iCol);
    if (hessian.value_[iEl]) num_diagonal_nz++;
    if (hessian.start_[iCol + 1] > iEl + 1) diagonal = false;
  }
  assert(diagonal == hessian.isDiagonal());
  highsLogUser(log_options, HighsLogType::kInfo,
               "Highs::solveMiqp Hessian has %d entries, and is %sdiagonal,"
               " with %d / %d nonzero diagonal entries\n",
               int(hessian.numNz()), diagonal ? "" : "not ",
               int(num_diagonal_nz), int(hessian.dim()));
}
