/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolveInitialSweep.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include "presolve/HPresolveUtils.h"
#include "presolve/HighsPostsolveStack.h"
#include "util/HighsCDouble.h"

#define CHECKED_CALL(call)                                \
  do {                                                    \
    HPresolveInitialSweep::Result __result = call;        \
    if (__result != HPresolveInitialSweep::Result::kOk) { \
      return __result;                                    \
    }                                                     \
  } while (0)

namespace presolve {

HPresolveInitialSweep::HPresolveInitialSweep(HighsLp& model,
                                             const HighsOptions& options,
                                             double primal_feastol)
    : model_(&model),
      options_(&options),
      primal_feastol_(primal_feastol),
      num_deleted_rows_(0),
      num_deleted_cols_(0) {}

HPresolveInitialSweep::Result HPresolveInitialSweep::checkColBounds(
    HighsInt col, bool& isFixed) {
  double boundDiff = model_->col_upper_[col] - model_->col_lower_[col];
  double max_abs_col_value = 0;
  for (HighsInt iEl = model_->a_matrix_.start_[col];
       iEl < model_->a_matrix_.start_[col + 1]; iEl++)
    max_abs_col_value =
        std::max(std::abs(model_->a_matrix_.value_[iEl]), max_abs_col_value);
  isFixed = false;
  assert(boundDiff >= 0);
  if (boundDiff <= primal_feastol_ &&
      (boundDiff <= options_->small_matrix_value ||
       max_abs_col_value * boundDiff <= primal_feastol_)) {
    if (std::abs(model_->col_lower_[col]) == kHighsInf)
      return Result::kDualInfeasible;
    isFixed = true;
  }
  return Result::kOk;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::emptyCol(
    HighsPostsolveStack& postsolve_stack, HighsInt col) {
  const HighsInt col_nnz =
      model_->a_matrix_.start_[col + 1] - model_->a_matrix_.start_[col];
  assert(col_nnz == 0);
  double cost = model_->col_cost_[col];
  const double lower = model_->col_lower_[col];
  const double upper = model_->col_upper_[col];

  if ((cost > 0 && lower == -kHighsInf) || (cost < 0 && upper == kHighsInf)) {
    if (std::abs(cost) <= options_->dual_feasibility_tolerance)
      cost = 0;
    else
      return Result::kDualInfeasible;
  }
  double fixval = kHighsInf;
  if (cost > 0) {
    fixval = lower;
    if (fixval == -kHighsInf) return Result::kDualInfeasible;
  } else if (cost < 0 || std::abs(upper) < std::abs(lower)) {
    fixval = upper;
    if (fixval == kHighsInf) return Result::kDualInfeasible;
  } else if (lower != -kHighsInf) {
    fixval = lower;
    if (fixval == -kHighsInf) return Result::kDualInfeasible;
  } else {
    fixval = 0.0;
  }
  assert(fixval != kHighsInf);
  postsolve_stack.removedModelFixedCol(col, fixval, cost, col_nnz, nullptr,
                                       nullptr);
  num_deleted_cols_++;
  if (col == model_->fme_obj_col_) model_->fme_obj_col_ = -1;
  return Result::kOk;
}

void HPresolveInitialSweep::removeFixedCol(HighsInt col) {
  double fixval = model_->col_lower_[col];
  num_deleted_cols_++;
  if (col == model_->fme_obj_col_) model_->fme_obj_col_ = -1;
  for (HighsInt iEl = model_->a_matrix_.start_[col];
       iEl < model_->a_matrix_.start_[col + 1]; iEl++) {
    HighsInt colrow = model_->a_matrix_.index_[iEl];
    double colval = model_->a_matrix_.value_[iEl];
    if (model_->row_lower_[colrow] != -kHighsInf)
      model_->row_lower_[colrow] -= colval * fixval;
    if (model_->row_upper_[colrow] != kHighsInf)
      model_->row_upper_[colrow] -= colval * fixval;
  }
}

HPresolveInitialSweep::Result HPresolveInitialSweep::emptyRow(
    HighsPostsolveStack& postsolve_stack, HighsInt row) {
  if (model_->row_upper_[row] < -primal_feastol_ ||
      model_->row_lower_[row] > primal_feastol_)
    return Result::kPrimalInfeasible;
  postsolve_stack.redundantRow(row);
  return Result::kOk;
}

double HPresolveInitialSweep::getMaxAbsColVal(HighsInt col) const {
  double maxVal = 0.0;
  for (HighsInt iEl = model_->a_matrix_.start_[col];
       iEl < model_->a_matrix_.start_[col + 1]; iEl++)
    maxVal = std::max(std::abs(model_->a_matrix_.value_[iEl]), maxVal);
  return maxVal;
}

bool HPresolveInitialSweep::isRedundant(HighsInt row, double sumLower,
                                        double sumUpper) const {
  return sumLower >= model_->row_lower_[row] - primal_feastol_ &&
         sumUpper <= model_->row_upper_[row] + primal_feastol_;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::singletonRow(
    HighsPostsolveStack& postsolve_stack, HighsInt row, HighsInt col,
    double val) {
  num_deleted_rows_++;

  assert(std::abs(val) > options_->small_matrix_value);

  double lb, ub;
  bool lowerTightened, upperTightened;
  const bool isIntegral =
      model_->integrality_[col] != HighsVarType::kContinuous;
  SingletonRowResult sr = computeSingletonRowBounds(
      val, model_->row_lower_[row], model_->row_upper_[row],
      model_->col_lower_[col], model_->col_upper_[col], primal_feastol_,
      getMaxAbsColVal(col), isIntegral, lb, ub, lowerTightened, upperTightened);
  if (sr == SingletonRowResult::kRedundant) {
    postsolve_stack.redundantRow(row);
    return Result::kOk;
  }
  if (sr == SingletonRowResult::kPrimalInfeasible)
    return Result::kPrimalInfeasible;

  postsolve_stack.singletonRow(row, col, val, lowerTightened, upperTightened);

  model_->col_lower_[col] = lb;
  model_->col_upper_[col] = ub;
  return Result::kOk;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::run(
    HighsPostsolveStack& postsolve_stack) {
  const bool have_col_names = model_->col_names_.size() > 0;
  const bool have_row_names = model_->row_names_.size() > 0;
  const HighsInt original_num_col = model_->num_col_;
  const HighsInt original_num_row = model_->num_row_;

  HighsInt num_fixed_col = 0;
  HighsInt num_empty_col = 0;
  HighsInt num_empty_row = 0;
  HighsInt num_singleton_row = 0;
  HighsInt num_redundant_row = 0;
  HighsInt num_col = 0;
  HighsInt nnz = 0;
  bool isFixed;

  std::vector<HighsInt> newColIndex(model_->num_col_);
  std::vector<HighsInt> row_count(model_->num_row_, 0);
  std::vector<HighsInt> col_of_row(model_->num_row_, -1);
  std::vector<double> val_of_row(model_->num_row_, 0);
  std::vector<HighsCDouble> implied_row_lower(model_->num_row_, 0);
  std::vector<HighsCDouble> implied_row_upper(model_->num_row_, 0);

  // Column pass: remove empty and fixed columns, compress in place
  for (HighsInt iCol = 0; iCol < model_->num_col_; iCol++) {
    HighsInt col_nnz =
        model_->a_matrix_.start_[iCol + 1] - model_->a_matrix_.start_[iCol];
    CHECKED_CALL(checkColBounds(iCol, isFixed));
    if (col_nnz == 0) {
      newColIndex[iCol] = -1;
      num_empty_col++;
      CHECKED_CALL(emptyCol(postsolve_stack, iCol));
    } else if (isFixed) {
      newColIndex[iCol] = -1;
      num_fixed_col++;
      HighsInt iEl = model_->a_matrix_.start_[iCol];
      postsolve_stack.removedModelFixedCol(
          iCol, model_->col_lower_[iCol], model_->col_cost_[iCol], col_nnz,
          &model_->a_matrix_.index_[iEl], &model_->a_matrix_.value_[iEl]);
      removeFixedCol(iCol);
    } else {
      newColIndex[iCol] = num_col;
      model_->col_cost_[num_col] = model_->col_cost_[iCol];
      model_->col_lower_[num_col] = model_->col_lower_[iCol];
      model_->col_upper_[num_col] = model_->col_upper_[iCol];
      model_->integrality_[num_col] = model_->integrality_[iCol];
      if (have_col_names)
        model_->col_names_[num_col] = std::move(model_->col_names_[iCol]);
      HighsInt from_os = model_->a_matrix_.start_[iCol];
      HighsInt new_col_start = nnz;
      for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
        HighsInt iRow = model_->a_matrix_.index_[from_os + iEl];
        double value = model_->a_matrix_.value_[from_os + iEl];
        row_count[iRow]++;
        col_of_row[iRow] = num_col;
        val_of_row[iRow] = value;
        model_->a_matrix_.index_[nnz] = iRow;
        model_->a_matrix_.value_[nnz] = value;
        nnz++;
        double row_lower_bnd = value > 0 ? model_->col_lower_[num_col]
                                         : model_->col_upper_[num_col];
        double row_upper_bnd = value > 0 ? model_->col_upper_[num_col]
                                         : model_->col_lower_[num_col];
        if (std::abs(row_lower_bnd) == kHighsInf)
          implied_row_lower[iRow] = static_cast<HighsCDouble>(-kHighsInf);
        else if (static_cast<double>(implied_row_lower[iRow]) > -kHighsInf)
          implied_row_lower[iRow] +=
              static_cast<HighsCDouble>(value) * row_lower_bnd;
        if (std::abs(row_upper_bnd) == kHighsInf)
          implied_row_upper[iRow] = static_cast<HighsCDouble>(kHighsInf);
        else if (static_cast<double>(implied_row_upper[iRow]) < kHighsInf)
          implied_row_upper[iRow] +=
              static_cast<HighsCDouble>(value) * row_upper_bnd;
      }
      model_->a_matrix_.start_[num_col] = new_col_start;
      num_col++;
    }
  }
  model_->a_matrix_.start_[num_col] = nnz;
  HighsInt num_removed_cols = num_empty_col + num_fixed_col;
  assert(num_col + num_removed_cols == model_->num_col_);
  model_->col_cost_.resize(num_col);
  model_->col_lower_.resize(num_col);
  model_->col_upper_.resize(num_col);
  model_->integrality_.resize(num_col);
  if (have_col_names) model_->col_names_.resize(num_col);
  model_->num_col_ = num_col;
  model_->a_matrix_.num_col_ = num_col;
  model_->a_matrix_.start_.resize(num_col + 1);
  model_->a_matrix_.index_.resize(nnz);
  model_->a_matrix_.value_.resize(nnz);
  if (model_->fme_obj_col_ >= 0)
    model_->fme_obj_col_ = newColIndex[model_->fme_obj_col_];
  postsolve_stack.compressColIndexMap(newColIndex);

  // Row pass: count empty, singleton, and redundant rows
  for (HighsInt iRow = 0; iRow < model_->num_row_; iRow++) {
    if (row_count[iRow] == 0)
      num_empty_row++;
    else if (row_count[iRow] == 1)
      num_singleton_row++;
    else if (isRedundant(iRow, static_cast<double>(implied_row_lower[iRow]),
                         static_cast<double>(implied_row_upper[iRow])))
      num_redundant_row++;
  }

  HighsInt num_removed_rows =
      num_empty_row + num_singleton_row + num_redundant_row;
  if (num_empty_row || num_singleton_row || num_redundant_row) {
    HighsInt num_row = 0;
    std::vector<HighsBool> has_singleton_row(model_->num_col_, false);
    std::vector<HighsInt> newRowIndex(model_->num_row_);
    for (HighsInt iRow = 0; iRow < model_->num_row_; iRow++) {
      if (row_count[iRow] <= 1) {
        newRowIndex[iRow] = -1;
        if (row_count[iRow] == 0) {
          CHECKED_CALL(emptyRow(postsolve_stack, iRow));
          num_deleted_rows_++;
        } else {
          has_singleton_row[col_of_row[iRow]] = true;
          assert(val_of_row[iRow]);
          CHECKED_CALL(singletonRow(postsolve_stack, iRow, col_of_row[iRow],
                                    val_of_row[iRow]));
        }
      } else {
        if (isRedundant(iRow, static_cast<double>(implied_row_lower[iRow]),
                        static_cast<double>(implied_row_upper[iRow]))) {
          postsolve_stack.redundantRow(iRow);
          newRowIndex[iRow] = -1;
          num_deleted_rows_++;
          continue;
        }
        newRowIndex[iRow] = num_row;
        model_->row_lower_[num_row] = model_->row_lower_[iRow];
        model_->row_upper_[num_row] = model_->row_upper_[iRow];
        if (have_row_names)
          model_->row_names_[num_row] = std::move(model_->row_names_[iRow]);
        num_row++;
      }
    }
    assert(num_row + num_removed_rows == model_->num_row_);

    if (num_redundant_row == 0) {
      // Only removing singleton row entries — compress column-wise
      nnz = 0;
      HighsInt from_col = 0;
      auto shiftCols = [&](const HighsInt to_col) {
        for (HighsInt iCol = from_col; iCol < to_col; iCol++) {
          HighsInt from_os = model_->a_matrix_.start_[iCol];
          HighsInt col_nnz = model_->a_matrix_.start_[iCol + 1] - from_os;
          HighsInt new_col_start = nnz;
          for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
            HighsInt iRow = model_->a_matrix_.index_[from_os + iEl];
            HighsInt newRow = newRowIndex[iRow];
            assert(newRow >= 0);
            model_->a_matrix_.index_[nnz] = newRow;
            model_->a_matrix_.value_[nnz] =
                model_->a_matrix_.value_[from_os + iEl];
            nnz++;
          }
          model_->a_matrix_.start_[iCol] = new_col_start;
        }
      };
      for (HighsInt iCol0 = 0; iCol0 < model_->num_col_; iCol0++) {
        if (!has_singleton_row[iCol0]) continue;
        shiftCols(iCol0);
        HighsInt from_os = model_->a_matrix_.start_[iCol0];
        HighsInt col_nnz = model_->a_matrix_.start_[iCol0 + 1] - from_os;
        HighsInt new_col_start = nnz;
        bool found_row_singleton = false;
        for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
          HighsInt iRow = model_->a_matrix_.index_[from_os + iEl];
          HighsInt newRow = newRowIndex[iRow];
          if (newRow >= 0) {
            model_->a_matrix_.index_[nnz] = newRow;
            model_->a_matrix_.value_[nnz] =
                model_->a_matrix_.value_[from_os + iEl];
            nnz++;
          } else {
            assert(row_count[iRow] == 1);
            assert(col_of_row[iRow] == iCol0);
            assert(val_of_row[iRow] == model_->a_matrix_.value_[from_os + iEl]);
            found_row_singleton = true;
          }
        }
        assert(found_row_singleton);
        model_->a_matrix_.start_[iCol0] = new_col_start;
        from_col = iCol0 + 1;
      }
      shiftCols(model_->num_col_);
      model_->a_matrix_.start_[num_col] = nnz;
    } else {
      // Also removing redundant rows — convert to rowwise and compress
      nnz = 0;
      HighsInt from_row = 0;
      num_row = 0;
      auto shiftRows = [&](const HighsInt to_row) {
        for (HighsInt iRow = from_row; iRow < to_row; iRow++) {
          HighsInt new_row_start = nnz;
          for (HighsInt iEl = model_->a_matrix_.start_[iRow];
               iEl < model_->a_matrix_.start_[iRow + 1]; iEl++) {
            model_->a_matrix_.index_[nnz] = model_->a_matrix_.index_[iEl];
            model_->a_matrix_.value_[nnz] = model_->a_matrix_.value_[iEl];
            nnz++;
          }
          model_->a_matrix_.start_[num_row] = new_row_start;
          num_row++;
        }
      };
      model_->a_matrix_.ensureRowwise();
      for (HighsInt iRow0 = 0; iRow0 < model_->num_row_; iRow0++) {
        if (newRowIndex[iRow0] >= 0) continue;
        shiftRows(iRow0);
        from_row = iRow0 + 1;
      }
      shiftRows(model_->num_row_);
      assert(num_row + num_removed_rows == model_->num_row_);
      model_->a_matrix_.start_[num_row] = nnz;
      model_->a_matrix_.num_row_ = num_row;
      model_->a_matrix_.ensureColwise();
    }
    model_->row_lower_.resize(num_row);
    model_->row_upper_.resize(num_row);
    if (have_row_names) model_->row_names_.resize(num_row);
    model_->num_row_ = num_row;
    model_->a_matrix_.num_row_ = num_row;
    model_->a_matrix_.index_.resize(nnz);
    model_->a_matrix_.value_.resize(nnz);
    postsolve_stack.compressRowIndexMap(newRowIndex);
  }

  if (num_fixed_col || num_empty_col)
    highsLogUser(
        options_->log_options, HighsLogType::kInfo,
        "Initial sweep removes %d + %d = %d / %d empty + fixed columns\n",
        int(num_empty_col), int(num_fixed_col), int(num_removed_cols),
        int(original_num_col));
  if (num_empty_row || num_singleton_row || num_redundant_row)
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "Initial sweep identifies %d + %d + %d = %d / %d empty + "
                 "singleton + redundant rows\n",
                 int(num_empty_row), int(num_singleton_row),
                 int(num_redundant_row), int(num_removed_rows),
                 int(original_num_row));

  return Result::kOk;
}

}  // namespace presolve
