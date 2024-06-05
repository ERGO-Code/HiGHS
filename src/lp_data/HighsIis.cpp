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
  this->valid_ = false;
  this->strategy_ = kIisStrategyMin;
  this->col_index_.clear();
  this->row_index_.clear();
  this->col_bound_.clear();
  this->row_bound_.clear();
}

void HighsIis::removeCol(const HighsInt col) {
  HighsInt num_col = this->col_index_.size();
  assert(col < num_col);
  this->col_index_[col] = this->col_index_[num_col-1];
  this->col_index_.resize(num_col-1);  
}

void HighsIis::removeRow(const HighsInt row) {
  HighsInt num_row = this->row_index_.size();
  assert(row < num_row);
  this->row_index_[row] = this->row_index_[num_row-1];
  this->row_index_.resize(num_row-1);  
}

bool HighsIis::inconsistentBounds(const HighsLp& lp, const HighsOptions& options) {
  this->invalidate();
  const bool col_priority =
      options.iis_strategy == kIisStrategyFromRayColPriority ||
      options.iis_strategy == kIisStrategyFromLpColPriority;
  for (HighsInt k = 0; k < 2; k++) {
    if ((col_priority && k == 0) || (!col_priority && k == 1)) {
      // Loop over columns first
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (lp.col_lower_[iCol] - lp.col_upper_[iCol] >
            2 * options.primal_feasibility_tolerance) {
          this->col_index_.push_back(iCol);
          break;
        }
      }
      if (this->col_index_.size() > 0) break;
    } else {
      // Loop over rows first
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        if (lp.row_lower_[iRow] - lp.row_upper_[iRow] >
            2 * options.primal_feasibility_tolerance) {
          this->row_index_.push_back(iRow);
          break;
        }
      }
      if (this->row_index_.size() > 0) break;
    }
  }
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  // If none found then return false
  if (num_iis_col + num_iis_row == 0) return false;
  // Should have found exactly 1
  assert((num_iis_col == 1 || num_iis_row == 1) &&
         num_iis_col + num_iis_row < 2);
  assert(lp.a_matrix_.isColwise());
  if (num_iis_col > 0) {
    // Found inconsistent column
    HighsInt iCol = this->col_index_[0];
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      this->row_index_.push_back(lp.a_matrix_.index_[iEl]);

  } else {
    // Found inconsistent row
    HighsInt iRow = this->row_index_[0];
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
           iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
        if (lp.a_matrix_.index_[iEl] == iRow) this->col_index_.push_back(iCol);
    }
  }
  this->valid_ = true;
  this->strategy_ = options.iis_strategy;
  return true;
}

HighsStatus HighsIis::getData(const HighsLp& lp, const HighsOptions& options,
			      const std::vector<double>& dual_ray_value) {
  // Check for inconsistent column and row bounds should have been
  // done earlier
  assert(!this->inconsistentBounds(lp, options));

  if (options.iis_strategy == kIisStrategyFromRayRowPriority ||
      options.iis_strategy == kIisStrategyFromRayColPriority) {
    // Identify the LP corresponding to the ray
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
      printf("HighsIis::getData: dual_ray_value[%2d] = %g; to_row[%2d] = %d\n",
	     int(iRow), dual_ray_value[iRow], int(iRow), to_row[iRow]);
    }
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      bool use_col = false;
      for (HighsInt iEl = lp.a_matrix_.start_[iCol];
	   iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
	use_col = use_col || to_row[lp.a_matrix_.index_[iEl]] >= 0;
      if (use_col) from_col.push_back(iCol);
    }
    HighsInt to_num_col = from_col.size();
    HighsInt to_num_row = from_row.size();
    HighsLp to_lp;
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
    if (this->compute(to_lp, options) != HighsStatus::kOk) return HighsStatus::kError;
    // IIS col/row information is for to_lp, so indirect the values
    // into the original LP
    for (HighsInt iCol = 0; iCol < HighsInt(this->col_index_.size()); iCol++)
      this->col_index_[iCol] = from_col[iCol];
    for (HighsInt iRow = 0; iRow < HighsInt(this->row_index_.size()); iRow++)
      this->row_index_[iRow] = from_row[iRow];
  } else {
    // Use the whole LP
    if (this->compute(lp, options) != HighsStatus::kOk) return HighsStatus::kError;
  }
  return HighsStatus::kOk;
}

HighsStatus HighsIis::compute(const HighsLp& lp, const HighsOptions& options) {
  this->invalidate();
  const HighsLogOptions& log_options = options.log_options;
  const bool row_priority =
    options.iis_strategy == kIisStrategyFromRayRowPriority ||
    options.iis_strategy == kIisStrategyFromLpRowPriority;
  // Initially all columns and rows are candidates for the IIS
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    this->col_index_.push_back(iCol);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
    this->row_index_.push_back(iRow);
  Highs highs;
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("output_flag", false);
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);
  highs.setOptionValue("output_flag", true);
  highs.writeModel("");
  highs.setOptionValue("output_flag", false);
  // Zero the objective
  std::vector<double> cost;
  cost.assign(lp.num_col_, 0);
  status = highs.changeColsCost(0, lp.num_col_-1, cost.data());
  assert(status == HighsStatus::kOk);
  status = highs.run();
  if (status != HighsStatus::kOk) return status;
  assert(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  // Pass twice: rows before columns, or columns before rows, according to row_priority
  for (HighsInt k = 0; k < 2; k++) {
    const bool row_deletion = (row_priority && k == 0) || (!row_priority && k == 1);
    std::string type = row_deletion ? "row" : "col";
    // Perform deletion pass
    HighsInt num_cdd = row_deletion ? this->row_index_.size() : this->col_index_.size();
    for (HighsInt iCdd = 0; iCdd < num_cdd; iCdd++) {
      for (;;) {
	const HighsInt iX = row_deletion ? this->row_index_[iCdd] : this->col_index_[iCdd];
	const double lower = row_deletion ? lp.row_lower_[iX] : lp.col_lower_[iX];
	const double upper = row_deletion ? lp.row_upper_[iX] : lp.col_upper_[iX];
	// Record whether a bound can be dropped: by default it's
	// possible, and only not possible if the bound is finite and the
	// LP remains infeasible if the bound is dropped
	bool drop_lower = true;
	bool drop_upper = true;
	if (lower > -kHighsInf) {
	  // Drop the lower bound temporarily
	  status = row_deletion ? highs.changeRowBounds(iX, -kHighsInf, upper) : highs.changeColBounds(iX, -kHighsInf, upper);
	  assert(status == HighsStatus::kOk);
	  status = highs.run();
	  assert(status == HighsStatus::kOk);
	  HighsModelStatus model_status = highs.getModelStatus();
	  if (model_status == HighsModelStatus::kOptimal) {
	    // Now feasible, so restore the lower bound and indicate that
	    // it cannot be dropped permanently
	    status = row_deletion ? highs.changeRowBounds(iX, lower, upper) : highs.changeColBounds(iX, lower, upper);
	    assert(status == HighsStatus::kOk);
	    drop_lower = false;
	  } else {
	    assert(model_status == HighsModelStatus::kInfeasible);
	    // Bound can be dropped permanently
	  }
	}	
	if (upper < kHighsInf) {
	  // Drop the upper bound temporarily
	  status = row_deletion ? highs.changeRowBounds(iX, lower, kHighsInf) : highs.changeColBounds(iX, lower, kHighsInf);
	  assert(status == HighsStatus::kOk);
	  status = highs.run();
	  assert(status == HighsStatus::kOk);
	  HighsModelStatus model_status = highs.getModelStatus();
	  if (model_status == HighsModelStatus::kOptimal) {
	    // Now feasible, so restore the upper bound and indicate that
	    // it cannot be dropped permanently
	    status = row_deletion ? highs.changeRowBounds(iX, lower, upper) : highs.changeColBounds(iX, lower, upper);
	    assert(status == HighsStatus::kOk);
	    drop_upper = false;
	  } else {
	    assert(model_status == HighsModelStatus::kInfeasible);
	    // Bound can be dropped permanently
	  }
	}
	if (drop_lower && drop_upper) {
	  // Both bounds can be dropped, so remove from the set of
	  // candidates
	  status = row_deletion ? highs.changeRowBounds(iX, -kHighsInf, kHighsInf) : highs.changeColBounds(iX, -kHighsInf, kHighsInf);
	  assert(status == HighsStatus::kOk);
	  if (row_deletion) {
	    this->removeRow(iCdd);
	  } else {
	    this->removeCol(iCdd);
	  }	   
	  num_cdd--;
	  highsLogUser(log_options, HighsLogType::kInfo, "Dropped  %s %d from candidate set\n", type.c_str(), int(iX));
	  if (iCdd >= num_cdd) break;
	} else {
	  highsLogUser(log_options, HighsLogType::kInfo, "Retained %s %d in   candidate set\n", type.c_str(), int(iX));
	  break;
	}
      }
    }
  }
  this->valid_ = true;
  this->strategy_ = options.iis_strategy;
  return HighsStatus::kOk;
}
