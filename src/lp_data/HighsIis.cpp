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

std::string HighsIis::iisBoundStatusToString(HighsInt bound_status) const {
  if (bound_status == kIisBoundStatusDropped) return "Dropped";
  if (bound_status == kIisBoundStatusNull) return "   Null";
  if (bound_status == kIisBoundStatusFree) return "   Free";
  if (bound_status == kIisBoundStatusLower) return "  Lower";
  if (bound_status == kIisBoundStatusUpper) return "  Upper";
  if (bound_status == kIisBoundStatusBoxed) return "  Boxed";
  return "*****";
}

void HighsIis::report(const std::string message, const HighsLp& lp) const {
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  if (num_iis_col > 10 || num_iis_row > 10) return;
  printf("\nIIS %s\n===\n", message.c_str());
  printf("Column: ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++) printf("%9d ", iCol);
  printf("\nStatus: ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9s ", iisBoundStatusToString(this->col_bound_[iCol]).c_str());
  printf("\nLower:  ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9.2g ", lp.col_lower_[iCol]);
  printf("\nUpper:  ");
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++)
    printf("%9.2g ", lp.col_upper_[iCol]);
  printf("\n");
  printf("Row:    Status     Lower     Upper\n");
  for (HighsInt iRow = 0; iRow < num_iis_row; iRow++)
    printf("%2d   %9s %9.2g %9.2g\n", int(iRow),
           iisBoundStatusToString(this->row_bound_[iRow]).c_str(),
           lp.row_lower_[iRow], lp.row_upper_[iRow]);
  printf("\n");
}

void HighsIis::addCol(const HighsInt col, const HighsInt status) {
  this->col_index_.push_back(col);
  this->col_bound_.push_back(status);
}

void HighsIis::addRow(const HighsInt row, const HighsInt status) {
  this->row_index_.push_back(row);
  this->row_bound_.push_back(status);
}

void HighsIis::removeCol(const HighsInt col) {
  HighsInt num_col = this->col_index_.size();
  assert(col < num_col);
  this->col_index_[col] = this->col_index_[num_col - 1];
  this->col_index_.resize(num_col - 1);
}

void HighsIis::removeRow(const HighsInt row) {
  HighsInt num_row = this->row_index_.size();
  assert(row < num_row);
  this->row_index_[row] = this->row_index_[num_row - 1];
  this->row_index_.resize(num_row - 1);
}

bool HighsIis::trivial(const HighsLp& lp, const HighsOptions& options) {
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
          this->addCol(iCol, kIisBoundStatusBoxed);
          break;
        }
      }
      if (this->col_index_.size() > 0) break;
    } else {
      // Loop over rows first
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        if (lp.row_lower_[iRow] - lp.row_upper_[iRow] >
            2 * options.primal_feasibility_tolerance) {
          this->addRow(iRow, kIisBoundStatusBoxed);
          break;
        }
      }
      if (this->row_index_.size() > 0) break;
    }
  }
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  // If one is found then we're done
  if (num_iis_col + num_iis_row > 0) {
    // Should have found exactly 1
    assert((num_iis_col == 1 || num_iis_row == 1) &&
           num_iis_col + num_iis_row < 2);
    this->valid_ = true;
    this->strategy_ = options.iis_strategy;
    return true;
  }
  // Now look for empty rows that cannot have zero activity
  std::vector<HighsInt> count;
  count.assign(lp.num_row_, 0);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      count[lp.a_matrix_.index_[iEl]]++;
  }
  assert(this->row_index_.size() == 0);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (count[iRow] > 0) continue;
    if (lp.row_lower_[iRow] > options.primal_feasibility_tolerance) {
      this->addRow(iRow, kIisBoundStatusLower);
    } else if (lp.row_upper_[iRow] < -options.primal_feasibility_tolerance) {
      this->addRow(iRow, kIisBoundStatusUpper);
    }
    if (this->row_index_.size() > 0) {
      this->valid_ = true;
      this->strategy_ = options.iis_strategy;
      return true;
    }
  }
  return false;
}

HighsStatus HighsIis::getData(const HighsLp& lp, const HighsOptions& options,
                              const HighsBasis& basis,
                              const double* dual_ray_value) {
  // Check for trivial IIS should have been done earlier
  assert(!this->trivial(lp, options));

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
      to_lp.col_cost_.push_back(0);
      to_lp.col_lower_.push_back(lp.col_lower_[from_col[iCol]]);
      to_lp.col_upper_.push_back(lp.col_upper_[from_col[iCol]]);
      to_lp.col_names_.push_back(lp.col_names_[from_col[iCol]]);
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
      to_lp.row_lower_.push_back(lp.row_lower_[from_row[iRow]]);
      to_lp.row_upper_.push_back(lp.row_upper_[from_row[iRow]]);
      to_lp.row_names_.push_back(lp.row_names_[from_row[iRow]]);
    }
    if (this->compute(to_lp, options) != HighsStatus::kOk)
      return HighsStatus::kError;
    // IIS col/row information is for to_lp, so indirect the values
    // into the original LP
    for (HighsInt iCol = 0; iCol < HighsInt(this->col_index_.size()); iCol++)
      this->col_index_[iCol] = from_col[this->col_index_[iCol]];
    for (HighsInt iRow = 0; iRow < HighsInt(this->row_index_.size()); iRow++)
      this->row_index_[iRow] = from_row[this->row_index_[iRow]];
  } else {
    // Use the whole LP
    if (this->compute(lp, options, &basis) != HighsStatus::kOk)
      return HighsStatus::kError;
  }
  this->report("On exit", lp);
  return HighsStatus::kOk;
}

HighsStatus HighsIis::compute(const HighsLp& lp, const HighsOptions& options,
                              const HighsBasis* basis) {
  this->invalidate();
  const HighsLogOptions& log_options = options.log_options;
  const bool row_priority =
      options.iis_strategy == kIisStrategyFromRayRowPriority ||
      options.iis_strategy == kIisStrategyFromLpRowPriority;
  // Initially all columns and rows are candidates for the IIS
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) this->addCol(iCol);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) this->addRow(iRow);
  Highs highs;
  //  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("presolve", kHighsOffString);
  const HighsLp& incumbent_lp = highs.getLp();
  HighsStatus run_status = highs.passModel(lp);
  assert(run_status == HighsStatus::kOk);
  if (basis) highs.setBasis(*basis);
  //
  // Elasticity filter
  //
  // Construct the e-LP:
  //
  // Constraints L <= Ax <= U; l <= x <= u
  //
  // Transformed to
  //
  // L <= Ax + e_L - e_U <= U,
  //
  // l <=  x + e_l - e_u <= u,
  //
  // where the elastic variables are not used if the corresponding bound is infinite.
  //
  // x is free, and the objective is the sum of the elastic variables.
  //
  // Determine the number of lower and upper elastic variables for
  // columns and rows
  //
  // col_of_ecol lists the column indices corresponding to the entries in
  // bound_of_col_of_ecol so that the results can be interpreted
  //
  // row_of_ecol lists the row indices corresponding to the entries in
  // bound_of_row_of_ecol so that the results can be interpreted
  std::vector<HighsInt> col_of_ecol;
  std::vector<HighsInt> row_of_ecol;
  std::vector<double> bound_of_row_of_ecol;
  std::vector<double> bound_of_col_of_ecol;
  std::vector<double> erow_lower;
  std::vector<double> erow_upper;
  std::vector<HighsInt> erow_start;
  std::vector<HighsInt> erow_index;
  std::vector<double> erow_value;
  // Accumulate names for ecols and erows, re-using ecol_name for the
  // names of row ecols after defining the names of col ecols
  std::vector<std::string> ecol_name;
  std::vector<std::string> erow_name;
  // When defining names, need to know the column number 
  HighsInt previous_num_col = highs.getNumCol();
  HighsInt previous_num_row = highs.getNumRow();
  HighsInt col_ecol_offset = previous_num_col;
  erow_start.push_back(0);
  HighsInt evar_ix = lp.num_col_;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
    // Free columns have no erow
    if (lower <= -kHighsInf && upper >= kHighsInf) continue;
    erow_lower.push_back(lower);
    erow_upper.push_back(upper);
    erow_name.push_back("row_"+std::to_string(iCol)+"_"+lp.col_names_[iCol]+"_erow");
    // Define the entry for x[iCol]
    erow_index.push_back(iCol);
    erow_value.push_back(1);
    if (lower > -kHighsInf) {
      // New e_l variable 
      col_of_ecol.push_back(iCol);
      ecol_name.push_back("col_"+std::to_string(iCol)+"_"+lp.col_names_[iCol]+"_lower");
      bound_of_col_of_ecol.push_back(lower);
      erow_index.push_back(evar_ix);
      erow_value.push_back(1);
      evar_ix++;
    }
    if (upper < kHighsInf) {
      // New e_u variable 
      col_of_ecol.push_back(iCol);
      ecol_name.push_back("col_"+std::to_string(iCol)+"_"+lp.col_names_[iCol]+"_upper");
      bound_of_col_of_ecol.push_back(upper);
      erow_index.push_back(evar_ix);
      erow_value.push_back(-1);
      evar_ix++;
    }
    erow_start.push_back(erow_index.size());
    HighsInt row_nz = erow_start[erow_start.size()-1] - erow_start[erow_start.size()-2];
    printf("eRow for column %d has %d nonzeros\n", int(iCol), int(row_nz));
    assert(row_nz == 2 || row_nz == 3);
  }
  HighsInt num_new_col = col_of_ecol.size();
  HighsInt num_new_row = erow_start.size()-1;
  HighsInt num_new_nz = erow_start[num_new_row];
  printf("Elasticity filter: For columns there are %d variables and %d constraints\n", int(num_new_col), int(num_new_row));
  const bool write_model = true;
  // Free the original columns
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  col_lower.assign(lp.num_col_, -kHighsInf);
  col_upper.assign(lp.num_col_, kHighsInf);
  run_status = highs.changeColsBounds(0, lp.num_col_-1, col_lower.data(), col_upper.data());
  assert(run_status == HighsStatus::kOk);
  // Add the new columns
  std::vector<double> ecol_cost;
  std::vector<double> ecol_lower;
  std::vector<double> ecol_upper;
  ecol_cost.assign(num_new_col, 1);
  ecol_lower.assign(num_new_col, 0);
  ecol_upper.assign(num_new_col, kHighsInf);
  run_status = highs.addCols(num_new_col, ecol_cost.data(), ecol_lower.data(), ecol_upper.data(),
			     0, nullptr, nullptr, nullptr);
  assert(run_status == HighsStatus::kOk);
  // Add the new rows
  run_status = highs.addRows(num_new_row, erow_lower.data(), erow_upper.data(),
			     num_new_nz, erow_start.data(), erow_index.data(), erow_value.data());
  assert(run_status == HighsStatus::kOk);
  for (HighsInt iCol = 0; iCol < num_new_col; iCol++) {
    highs.passColName(previous_num_col+iCol, ecol_name[iCol]);
  }
  for (HighsInt iRow = 0; iRow < num_new_row; iRow++) {
    highs.passRowName(previous_num_row+iRow, erow_name[iRow]);
  }
  if (write_model) {
    printf("\nAfter adding e-rows\n=============\n");
    bool output_flag;
    run_status = highs.getOptionValue("output_flag", output_flag);
    highs.setOptionValue("output_flag", true);
    highs.writeModel("");
    highs.setOptionValue("output_flag", output_flag);
  }
  // Add the columns corresponding to the e_L and e_U variables for
  // the constraints
  ecol_name.clear();
  std::vector<HighsInt> ecol_start;
  std::vector<HighsInt> ecol_index;
  std::vector<double> ecol_value;
  ecol_start.push_back(0);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    if (lower > -kHighsInf) {
      // Create an e-var for the row lower bound
      row_of_ecol.push_back(iRow);
      ecol_name.push_back("row_"+std::to_string(iRow)+"_"+lp.row_names_[iRow]+"_lower");
      bound_of_row_of_ecol.push_back(lower);
      // Define the sub-matrix column
      ecol_index.push_back(iRow);
      ecol_value.push_back(1);
      ecol_start.push_back(ecol_index.size());
      evar_ix++;
    }
    if (upper < kHighsInf) {
      // Create an e-var for the row upper bound
      row_of_ecol.push_back(iRow);
      ecol_name.push_back("row_"+std::to_string(iRow)+"_"+lp.row_names_[iRow]+"_upper");
      bound_of_row_of_ecol.push_back(upper);
      // Define the sub-matrix column
      ecol_index.push_back(iRow);
      ecol_value.push_back(-1);
      ecol_start.push_back(ecol_index.size());
      evar_ix++;
    }
  }
  num_new_col = ecol_start.size()-1;
  num_new_nz = ecol_start[num_new_col];
  ecol_cost.assign(num_new_col, 1);
  ecol_lower.assign(num_new_col, 0);
  ecol_upper.assign(num_new_col, kHighsInf);
  previous_num_col = highs.getNumCol();
  HighsInt row_ecol_offset = previous_num_col;
  run_status = highs.addCols(num_new_col, ecol_cost.data(), ecol_lower.data(), ecol_upper.data(),
			     num_new_nz, ecol_start.data(), ecol_index.data(), ecol_value.data());
  assert(run_status == HighsStatus::kOk);
  for (HighsInt iCol = 0; iCol < num_new_col; iCol++) {
    highs.passColName(previous_num_col+iCol, ecol_name[iCol]);
  }

  if (write_model) {
    bool output_flag;
    printf("\nAfter adding e-cols\n=============\n");
    run_status = highs.getOptionValue("output_flag", output_flag);
    highs.setOptionValue("output_flag", true);
    highs.writeModel("");
    highs.setOptionValue("output_flag", output_flag);
  }
  run_status = highs.run();
  assert(run_status == HighsStatus::kOk);
  highs.writeSolution("", kSolutionStylePretty);
  HighsModelStatus model_status = highs.getModelStatus();
  assert(model_status == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  // Now fix e-variables that are positive and re-solve until e-LP is infeasible
  HighsInt loop_k = 0;
  for (;;) {
    printf("\nElasticity filter pass %d\n==============\n", int(loop_k));
    HighsInt num_fixed = 0;
    for (HighsInt eCol = 0; eCol < col_of_ecol.size(); eCol++) {
      HighsInt iCol = col_of_ecol[eCol];
      if (solution.col_value[col_ecol_offset+eCol] > options.primal_feasibility_tolerance) {
      printf("E-col %2d (column %2d) corresponds to column %2d with bound %g and has solution value %g\n",
	     int(eCol), int(col_ecol_offset+eCol), int(iCol), bound_of_col_of_ecol[eCol], solution.col_value[col_ecol_offset+eCol]);
	highs.changeColBounds(col_ecol_offset+eCol, 0, 0);
	num_fixed++;
      }
    }
    for (HighsInt eCol = 0; eCol < row_of_ecol.size(); eCol++) {
      HighsInt iRow = row_of_ecol[eCol];
      if (solution.col_value[row_ecol_offset+eCol] > options.primal_feasibility_tolerance) {
      printf("E-row %2d (column %2d) corresponds to    row %2d with bound %g and has solution value %g\n",
	     int(eCol), int(row_ecol_offset+eCol), int(iRow), bound_of_row_of_ecol[eCol], solution.col_value[row_ecol_offset+eCol]);
	highs.changeColBounds(row_ecol_offset+eCol, 0, 0);
	num_fixed++;
      }
    }
    assert(num_fixed>0);
    run_status = highs.run();
    assert(run_status == HighsStatus::kOk);
    highs.writeSolution("", kSolutionStylePretty);
    HighsModelStatus model_status = highs.getModelStatus();
    if (model_status == HighsModelStatus::kInfeasible) break;
    loop_k++;
    if (loop_k>10) assert(1666==1999);
  }

  HighsInt num_enforced_col_ecol = 0;
  HighsInt num_enforced_row_ecol = 0;
  for (HighsInt eCol = 0; eCol < col_of_ecol.size(); eCol++) {
    HighsInt iCol = col_of_ecol[eCol];
    if (incumbent_lp.col_upper_[col_ecol_offset+eCol] == 0) {
      num_enforced_col_ecol++;
      printf("Col e-col %2d (column %2d) corresponds to column %2d with bound %g and is enforced\n",
	     int(eCol), int(col_ecol_offset+eCol), int(iCol), bound_of_col_of_ecol[eCol]);
    }
  }
  for (HighsInt eCol = 0; eCol < row_of_ecol.size(); eCol++) {
    HighsInt iRow = row_of_ecol[eCol];
    if (incumbent_lp.col_upper_[row_ecol_offset+eCol] == 0) {
      num_enforced_row_ecol++;
      printf("Row e-col %2d (column %2d) corresponds to    row %2d with bound %g and is enforced\n",
	     int(eCol), int(row_ecol_offset+eCol), int(iRow), bound_of_row_of_ecol[eCol]);
    }
  }
  printf("\nElasticity filter after %d passes enforces bounds on %d cols and %d rows\n", int(loop_k), int(num_enforced_col_ecol), int(num_enforced_row_ecol));

  assert(666==999);
  // Zero the objective
  std::vector<double> cost;
  cost.assign(lp.num_col_, 0);
  run_status = highs.changeColsCost(0, lp.num_col_ - 1, cost.data());
  assert(run_status == HighsStatus::kOk);
  // Solve the LP
  if (basis) highs.setBasis(*basis);
  run_status = highs.run();
  assert(run_status == HighsStatus::kOk);
  if (run_status != HighsStatus::kOk) return run_status;
  assert(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  const bool use_sensitivity_filter = false;
  if (use_sensitivity_filter) {
    highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
    //
    highs.setOptionValue("output_flag", true);
    // Solve the LP
    run_status = highs.run();
    highs.writeSolution("", kSolutionStylePretty);
  }

  // Pass twice: rows before columns, or columns before rows, according to
  // row_priority
  for (HighsInt k = 0; k < 2; k++) {
    const bool row_deletion =
        (row_priority && k == 0) || (!row_priority && k == 1);
    std::string type = row_deletion ? "Row" : "Col";
    // Perform deletion pass
    HighsInt num_index = row_deletion ? lp.num_row_ : lp.num_col_;
    for (HighsInt iX = 0; iX < num_index; iX++) {
      const HighsInt ix_status =
          row_deletion ? this->row_bound_[iX] : this->col_bound_[iX];
      if (ix_status == kIisBoundStatusDropped ||
          ix_status == kIisBoundStatusFree)
        continue;
      double lower = row_deletion ? lp.row_lower_[iX] : lp.col_lower_[iX];
      double upper = row_deletion ? lp.row_upper_[iX] : lp.col_upper_[iX];
      // Record whether the upper bound has been dropped due to the lower bound
      // being kept
      bool drop_upper = false;
      if (lower > -kHighsInf) {
        // Drop the lower bound temporarily
        run_status = row_deletion
                         ? highs.changeRowBounds(iX, -kHighsInf, upper)
                         : highs.changeColBounds(iX, -kHighsInf, upper);
        assert(run_status == HighsStatus::kOk);
        run_status = highs.run();
        assert(run_status == HighsStatus::kOk);
        HighsModelStatus model_status = highs.getModelStatus();
        if (model_status == HighsModelStatus::kOptimal) {
          // Now feasible, so restore the lower bound
          run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                                    : highs.changeColBounds(iX, lower, upper);
          assert(run_status == HighsStatus::kOk);
          // If the lower bound must be kept, then any finite upper bound
          // must be dropped
          const bool apply_reciprocal_rule = false;
          if (apply_reciprocal_rule) {
            if (upper < kHighsInf) {
              // Drop the upper bound permanently
              upper = kHighsInf;
              run_status = row_deletion
                               ? highs.changeRowBounds(iX, lower, kHighsInf)
                               : highs.changeColBounds(iX, lower, upper);
              assert(run_status == HighsStatus::kOk);
              drop_upper = true;
            }
            //	    continue;
          }
        } else {
          // Bound can be dropped permanently
          assert(model_status == HighsModelStatus::kInfeasible);
          lower = -kHighsInf;
        }
      }
      if (upper < kHighsInf) {
        // Drop the upper bound temporarily
        run_status = row_deletion ? highs.changeRowBounds(iX, lower, kHighsInf)
                                  : highs.changeColBounds(iX, lower, kHighsInf);
        assert(run_status == HighsStatus::kOk);
        run_status = highs.run();
        assert(run_status == HighsStatus::kOk);
        HighsModelStatus model_status = highs.getModelStatus();
        // If the upper bound has been dropped due to the reciprical
        // rule, the LP must be infeasible
        if (drop_upper) assert(model_status == HighsModelStatus::kInfeasible);
        if (model_status == HighsModelStatus::kOptimal) {
          // Now feasible, so restore the upper bound
          run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                                    : highs.changeColBounds(iX, lower, upper);
          assert(run_status == HighsStatus::kOk);
        } else {
          // Bound can be dropped permanently
          assert(model_status == HighsModelStatus::kInfeasible);
          upper = kHighsInf;
        }
      }
      const bool debug_bound_change = true;
      if (debug_bound_change) {
        // Check bounds have been changed correctly
        double check_lower;
        double check_upper;
        double check_cost;
        HighsInt check_num_ix;
        HighsInt check_num_nz;
        run_status =
            row_deletion
                ? highs.getRows(iX, iX, check_num_ix, &check_lower,
                                &check_upper, check_num_nz, nullptr, nullptr,
                                nullptr)
                : highs.getCols(iX, iX, check_num_ix, &check_cost, &check_lower,
                                &check_upper, check_num_nz, nullptr, nullptr,
                                nullptr);
        assert(run_status == HighsStatus::kOk);
        assert(check_lower == lower);
        assert(check_upper == upper);
      }
      HighsInt iss_bound_status = kIisBoundStatusNull;
      if (lower <= -kHighsInf) {
        if (upper >= kHighsInf) {
          if (row_deletion) {
            // Free rows can be dropped
            iss_bound_status = kIisBoundStatusDropped;
          } else {
            // Free columns can only be dropped if they are empty
            iss_bound_status = kIisBoundStatusFree;
          }
        } else {
          iss_bound_status = kIisBoundStatusUpper;
        }
      } else {
        if (upper >= kHighsInf) {
          iss_bound_status = kIisBoundStatusLower;
        } else {
          // FX or BX: shouldn't happen
          iss_bound_status = kIisBoundStatusBoxed;
        }
      }
      assert(iss_bound_status != kIisBoundStatusNull);
      assert(iss_bound_status != kIisBoundStatusBoxed);
      if (row_deletion) {
        this->row_bound_[iX] = iss_bound_status;
      } else {
        this->col_bound_[iX] = iss_bound_status;
      }
      highsLogUser(log_options, HighsLogType::kInfo, "%s %d has status %s\n",
                   type.c_str(), int(iX),
                   iisBoundStatusToString(iss_bound_status).c_str());
    }
    if (k == 1) continue;
    // End of first pass: look to simplify second pass
    this->report("End of deletion", incumbent_lp);
    if (row_deletion) {
      // Mark empty columns as dropped
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        bool empty_col = true;
        for (HighsInt iEl = lp.a_matrix_.start_[iCol];
             iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
          if (this->row_bound_[lp.a_matrix_.index_[iEl]] !=
              kIisBoundStatusDropped) {
            empty_col = false;
            break;
          }
        }
        if (empty_col) {
          highsLogUser(log_options, HighsLogType::kInfo,
                       "Col %d has status Dropped: Empty\n", int(iCol));
          this->col_bound_[iCol] = kIisBoundStatusDropped;
          run_status = highs.changeColBounds(iCol, -kHighsInf, kHighsInf);
          assert(run_status == HighsStatus::kOk);
        }
      }
    }
    this->report("End of pass 1", incumbent_lp);
  }
  this->report("End of pass 2", incumbent_lp);
  HighsInt iss_num_col = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (this->col_bound_[iCol] != kIisBoundStatusDropped) {
      this->col_index_[iss_num_col] = this->col_index_[iCol];
      this->col_bound_[iss_num_col] = this->col_bound_[iCol];
      iss_num_col++;
    }
  }
  HighsInt iss_num_row = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (this->row_bound_[iRow] != kIisBoundStatusDropped) {
      this->row_index_[iss_num_row] = this->row_index_[iRow];
      this->row_bound_[iss_num_row] = this->row_bound_[iRow];
      iss_num_row++;
    }
  }
  this->col_index_.resize(iss_num_col);
  this->col_bound_.resize(iss_num_col);
  this->row_index_.resize(iss_num_row);
  this->row_bound_.resize(iss_num_row);
  this->valid_ = true;
  this->strategy_ = options.iis_strategy;
  return HighsStatus::kOk;
}
