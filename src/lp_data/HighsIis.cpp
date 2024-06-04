/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsIis.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "lp_data/HighsIis.h"

#include "Highs.h"

void HighsIis::invalidate() {
  this->valid = false;
  this->strategy = kIisStrategyMin;
  this->col_index.clear();
  this->row_index.clear();
  this->col_bound.clear();
  this->row_bound.clear();
}

bool iisInconsistentBounds(const HighsLp& lp, const HighsOptions& options,
                           HighsIis& iis) {
  iis.invalidate();
  const bool col_priority =
      options.iis_strategy == kIisStrategyFromRayColPriority ||
      options.iis_strategy == kIisStrategyFromLpColPriority;
  for (HighsInt k = 0; k < 2; k++) {
    if ((col_priority && k == 0) || (!col_priority && k == 1)) {
      // Loop over columns first
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (lp.col_lower_[iCol] - lp.col_upper_[iCol] >
            2 * options.primal_feasibility_tolerance) {
          iis.col_index.push_back(iCol);
          break;
        }
      }
      if (iis.col_index.size() > 0) break;
    } else {
      // Loop over rows first
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        if (lp.row_lower_[iRow] - lp.row_upper_[iRow] >
            2 * options.primal_feasibility_tolerance) {
          iis.row_index.push_back(iRow);
          break;
        }
      }
      if (iis.row_index.size() > 0) break;
    }
  }
  HighsInt num_iis_col = iis.col_index.size();
  HighsInt num_iis_row = iis.row_index.size();
  // If none found then return false
  if (num_iis_col + num_iis_row == 0) return false;
  // Should have found exactly 1
  assert((num_iis_col == 1 || num_iis_row == 1) &&
         num_iis_col + num_iis_row < 2);
  assert(lp.a_matrix_.isColwise());
  if (num_iis_col > 0) {
    // Found inconsistent column
    HighsInt iCol = iis.col_index[0];
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      iis.row_index.push_back(lp.a_matrix_.index_[iEl]);

  } else {
    // Found inconsistent row
    HighsInt iRow = iis.row_index[0];
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
        if (lp.a_matrix_.index_[iEl] == iRow) iis.col_index.push_back(iCol);
    }
  }
  iis.valid = true;
  iis.strategy = options.iis_strategy;
  return true;
}

HighsStatus getIisData(const HighsLp& lp, const HighsOptions& options,
                       const std::vector<double>& dual_ray_value,
                       HighsIis& iis) {
  // Check for inconsistent column and row bounds should have been
  // done earlier
  assert(!iisInconsistentBounds(lp, options, iis));
  std::vector<HighsInt> from_row;
  std::vector<HighsInt> from_col;
  std::vector<HighsInt> to_row;
  to_row.assign(lp.num_row_, -1);
  assert(lp.a_matrix_.isColwise());
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (dual_ray_value[iRow]) {
      to_row[iRow] = from_row.size();
      from_row.push_back(iRow);
    }
    printf("getIisData: dual_ray_value[%2d] = %g; to_row[%2d] = %d\n",
           int(iRow), dual_ray_value[iRow], int(iRow), to_row[iRow]);
  }
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    bool use_col = false;
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      use_col = use_col || to_row[lp.a_matrix_.index_[iEl]] >= 0;
    if (use_col) from_col.push_back(iCol);
  }
  HighsLp to_lp;
  HighsInt to_num_col = from_col.size();
  HighsInt to_num_row = from_row.size();
  to_lp.num_col_ = to_num_col;
  to_lp.num_row_ = to_num_row;
  for (HighsInt iCol = 0; iCol < to_num_col; iCol++) {
    printf("getIisData: from_col[%2d] = %d\n", int(iCol), int(from_col[iCol]));
    to_lp.col_cost_.push_back(0);
    to_lp.col_lower_.push_back(lp.col_lower_[from_col[iCol]]);
    to_lp.col_upper_.push_back(lp.col_upper_[from_col[iCol]]);
    for (HighsInt iEl = lp.a_matrix_.start_[from_col[iCol]];
         iEl < lp.a_matrix_.start_[from_col[iCol] + 1]; iEl++) {
      HighsInt iRow = lp.a_matrix_.index_[iEl];
      if (to_row[iRow] >= 0) {
        to_lp.a_matrix_.index_.push_back(to_row[iRow]);
        to_lp.a_matrix_.value_.push_back(lp.a_matrix_.value_[iEl]);
      }
    }
    to_lp.a_matrix_.start_.push_back(to_lp.a_matrix_.index_.size());
  }
  for (HighsInt iRow = 0; iRow < to_num_row; iRow++) {
    printf("getIisData: from_row[%2d] = %d\n", int(iRow), int(from_row[iRow]));
    to_lp.row_lower_.push_back(lp.row_lower_[from_row[iRow]]);
    to_lp.row_upper_.push_back(lp.row_upper_[from_row[iRow]]);
  }
  Highs highs;
  highs.setOptionValue("presolve", kHighsOffString);
  HighsStatus status = highs.passModel(to_lp);
  assert(status == HighsStatus::kOk);
  status = highs.run();
  if (status != HighsStatus::kOk) return status;
  assert(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  iis.valid = true;
  iis.strategy = options.iis_strategy;
  iis.col_index = from_col;
  iis.row_index = from_row;

  return HighsStatus::kOk;
}
