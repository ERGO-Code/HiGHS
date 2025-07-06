/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsIis.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "Highs.h"
#include <tuple>

void HighsIis::invalidate() {
  this->valid_ = false;
  this->strategy_ = kIisStrategyMin;
  this->col_index_.clear();
  this->row_index_.clear();
  this->col_bound_.clear();
  this->row_bound_.clear();
  this->info_.clear();
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
  for (HighsInt iCol = 0; iCol < num_iis_col; iCol++) printf("%9d ", int(iCol));
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
      //      options.iis_strategy == kIisStrategyFromRayColPriority ||
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
                              const std::vector<HighsInt>& infeasible_row) {
  // Check for trivial IIS should have been done earlier
  assert(!this->trivial(lp, options));
  // The number of infeasible rows must be positive
  assert(infeasible_row.size() > 0);
  // Identify the LP corresponding to the set of infeasible rows
  std::vector<HighsInt> from_row = infeasible_row;
  std::vector<HighsInt> from_col;
  std::vector<HighsInt> to_row;
  to_row.assign(lp.num_row_, -1);
  assert(lp.a_matrix_.isColwise());
  // Determine how to detect whether a row is in infeasible_row and
  // (then) gather information about it
  for (HighsInt iX = 0; iX < HighsInt(infeasible_row.size()); iX++)
    to_row[infeasible_row[iX]] = iX;
  // Identify the columns (from_col) with nonzeros in the infeasible
  // rows
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
  to_lp.a_matrix_.num_col_ = to_lp.num_col_;
  to_lp.a_matrix_.num_row_ = to_lp.num_row_;
  const bool has_col_names = lp.col_names_.size() > 0;
  for (HighsInt iCol = 0; iCol < to_num_col; iCol++) {
    to_lp.col_cost_.push_back(0);
    to_lp.col_lower_.push_back(lp.col_lower_[from_col[iCol]]);
    to_lp.col_upper_.push_back(lp.col_upper_[from_col[iCol]]);
    if (has_col_names)
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
  const bool has_row_names = lp.row_names_.size() > 0;
  for (HighsInt iRow = 0; iRow < to_num_row; iRow++) {
    to_lp.row_lower_.push_back(lp.row_lower_[from_row[iRow]]);
    to_lp.row_upper_.push_back(lp.row_upper_[from_row[iRow]]);
    if (has_row_names)
      to_lp.row_names_.push_back(lp.row_names_[from_row[iRow]]);
  }
  if (this->compute(to_lp, options) != HighsStatus::kOk)
    return HighsStatus::kError;
  // Indirect the values into the original LP
  for (HighsInt& colindex : this->col_index_) colindex = from_col[colindex];
  for (HighsInt& rowindex : this->row_index_) rowindex = from_row[rowindex];
  if (kIisDevReport) this->report("On exit", lp);
  return HighsStatus::kOk;
}

HighsStatus HighsIis::compute(const HighsLp& lp, const HighsOptions& options,
                              const HighsBasis* basis) {
  const HighsLogOptions& log_options = options.log_options;
  const bool row_priority =
      //      options.iis_strategy == kIisStrategyFromRayRowPriority ||
      options.iis_strategy == kIisStrategyFromLpRowPriority;
  // Initially all columns and rows are candidates for the IIS
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) this->addCol(iCol);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) this->addRow(iRow);
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", kIisDevReport);
  highs.setOptionValue("presolve", kHighsOffString);
  const HighsLp& incumbent_lp = highs.getLp();
  const HighsBasis& incumbent_basis = highs.getBasis();
  const HighsSolution& solution = highs.getSolution();
  HighsStatus run_status = highs.passModel(lp);
  assert(run_status == HighsStatus::kOk);
  if (basis) highs.setBasis(*basis);

  // Zero the objective
  std::vector<double> cost;
  cost.assign(lp.num_col_, 0);
  run_status = highs.changeColsCost(0, lp.num_col_ - 1, cost.data());
  assert(run_status == HighsStatus::kOk);
  // Solve the LP
  if (basis) highs.setBasis(*basis);
  const bool use_sensitivity_filter = false;
  std::vector<double> primal_phase1_dual;
  bool row_deletion = false;
  HighsInt iX = -1;
  bool drop_lower = false;

  // Lambda for gathering data when solving an LP
  auto solveLp = [&]() -> HighsStatus {
    HighsIisInfo iis_info;
    iis_info.simplex_time = -highs.getRunTime();
    iis_info.simplex_iterations = -info.simplex_iteration_count;
    run_status = highs.run();
    assert(run_status == HighsStatus::kOk);
    if (run_status != HighsStatus::kOk) return run_status;
    HighsModelStatus model_status = highs.getModelStatus();
    if (use_sensitivity_filter &&
        model_status == HighsModelStatus::kInfeasible) {
      printf("\nHighsIis::compute %s deletion for %d and %s bound\n",
             row_deletion ? "Row" : "Col", int(iX),
             drop_lower ? "Lower" : "Upper");
      bool output_flag;
      highs.getOptionValue("output_flag", output_flag);
      highs.setOptionValue("output_flag", true);
      HighsInt simplex_strategy;
      highs.getOptionValue("simplex_strategy", simplex_strategy);
      highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
      // Solve the LP
      run_status = highs.run();
      if (run_status != HighsStatus::kOk) return run_status;
      highs.writeSolution("", kSolutionStylePretty);
      const HighsInt* basic_index = highs.getBasicVariablesArray();
      std::vector<double> rhs;
      rhs.assign(lp.num_row_, 0);
      // Get duals for nonbasic rows, and initialise duals so that basic duals
      // are zero
      assert(101 == 202);

      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        HighsInt iVar = basic_index[iRow];
        const double lower = iVar < lp.num_col_
                                 ? lp.col_lower_[iVar]
                                 : lp.row_lower_[iVar - lp.num_col_];
        const double upper = iVar < lp.num_col_
                                 ? lp.col_upper_[iVar]
                                 : lp.row_upper_[iVar - lp.num_col_];
        const double value = iVar < lp.num_col_
                                 ? solution.col_value[iVar]
                                 : solution.row_value[iVar - lp.num_col_];
        if (value < lower - options.primal_feasibility_tolerance) {
          rhs[iRow] = -1;
        } else if (value > upper + options.primal_feasibility_tolerance) {
          rhs[iRow] = 1;
        }
      }
      HVector pi;
      pi.setup(lp.num_row_);
      highs.getBasisTransposeSolve(rhs.data(), &pi.array[0], NULL, NULL);
      pi.count = lp.num_row_;
      std::vector<double> reduced_costs_value;
      std::vector<HighsInt> reduced_costs_index;
      lp.a_matrix_.productTransposeQuad(reduced_costs_value,
                                        reduced_costs_index, pi);

      primal_phase1_dual = highs.getPrimalPhase1Dual();
      HighsInt num_zero_dual = 0;
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        const HighsBasisStatus status = incumbent_basis.col_status[iCol];
        const double dual = primal_phase1_dual[iCol];
        const double lower = lp.col_lower_[iCol];
        const double upper = lp.col_upper_[iCol];
        const double value = solution.col_value[iCol];
        if (status != HighsBasisStatus::kBasic &&
            std::fabs(dual) < options.dual_feasibility_tolerance) {
          num_zero_dual++;
          // Small dual for nonbasic variable
          printf(
              "HighsIis::compute Column %d [%g, %g, %g] with status %s has "
              "dual %g\n",
              int(iCol), lower, value, upper,
              highs.basisStatusToString(status).c_str(), dual);
          //	  assert(123 == 456);
        }
      }
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        const HighsBasisStatus status = incumbent_basis.row_status[iRow];
        const double dual = primal_phase1_dual[lp.num_col_ + iRow];
        const double lower = lp.row_lower_[iRow];
        const double upper = lp.row_upper_[iRow];
        const double value = solution.row_value[iRow];
        if (status != HighsBasisStatus::kBasic &&
            std::fabs(dual) < options.dual_feasibility_tolerance) {
          num_zero_dual++;
          // Small dual for nonbasic variable
          printf(
              "HighsIis::compute Row    %d [%g, %g, %g] with status %s has "
              "dual %g\n",
              int(iRow), lower, value, upper,
              highs.basisStatusToString(status).c_str(), dual);
          //	  assert(123 == 456);
        }
      }
      highs.setOptionValue("output_flag", output_flag);
      highs.setOptionValue("simplex_strategy", simplex_strategy);
      assert(!num_zero_dual);
    }
    iis_info.simplex_time += highs.getRunTime();
    iis_info.simplex_iterations += info.simplex_iteration_count;
    this->info_.push_back(iis_info);
    return run_status;
  };

  run_status = solveLp();
  if (run_status != HighsStatus::kOk) return run_status;

  assert(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  // Pass twice: rows before columns, or columns before rows, according to
  // row_priority
  for (HighsInt k = 0; k < 2; k++) {
    row_deletion = (row_priority && k == 0) || (!row_priority && k == 1);
    std::string type = row_deletion ? "Row" : "Col";
    // Perform deletion pass
    HighsInt num_index = row_deletion ? lp.num_row_ : lp.num_col_;
    for (iX = 0; iX < num_index; iX++) {
      const HighsInt ix_status =
          row_deletion ? this->row_bound_[iX] : this->col_bound_[iX];
      if (ix_status == kIisBoundStatusDropped ||
          ix_status == kIisBoundStatusFree)
        continue;
      double lower = row_deletion ? lp.row_lower_[iX] : lp.col_lower_[iX];
      double upper = row_deletion ? lp.row_upper_[iX] : lp.col_upper_[iX];
      // Record whether the upper bound has been dropped due to the lower bound
      // being kept
      if (lower > -kHighsInf) {
        // Drop the lower bound temporarily
        bool drop_lower = true;
        run_status = row_deletion
                         ? highs.changeRowBounds(iX, -kHighsInf, upper)
                         : highs.changeColBounds(iX, -kHighsInf, upper);
        assert(run_status == HighsStatus::kOk);
        // Solve the LP
        run_status = solveLp();
        if (run_status != HighsStatus::kOk) return run_status;
        HighsModelStatus model_status = highs.getModelStatus();
        if (model_status == HighsModelStatus::kOptimal) {
          // Now feasible, so restore the lower bound
          run_status = row_deletion ? highs.changeRowBounds(iX, lower, upper)
                                    : highs.changeColBounds(iX, lower, upper);
          assert(run_status == HighsStatus::kOk);
          // If the lower bound must be kept, then any finite upper bound
          // must be dropped
          const bool apply_reciprocal_rule = true;
          if (apply_reciprocal_rule) {
            if (upper < kHighsInf) {
              // Drop the upper bound permanently
              upper = kHighsInf;
              run_status = row_deletion
                               ? highs.changeRowBounds(iX, lower, upper)
                               : highs.changeColBounds(iX, lower, upper);
              assert(run_status == HighsStatus::kOk);
            }
            assert(upper >= kHighsInf);
            // Since upper = kHighsInf, allow the loop to run so that
            // bound status is set as if upper were set to kHighsInf
            // by relaxing it and finding that the LP was still
            // infeasible
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
        // Solve the LP
        run_status = solveLp();
        if (run_status != HighsStatus::kOk) return run_status;
        HighsModelStatus model_status = highs.getModelStatus();
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
    if (kIisDevReport) this->report("End of deletion", incumbent_lp);
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
    if (kIisDevReport) this->report("End of pass 1", incumbent_lp);
  }
  if (kIisDevReport) this->report("End of pass 2", incumbent_lp);
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

bool HighsIis::rowValueBounds(const HighsLp& lp, const HighsOptions& options) {
  // Look for infeasible rows based on row value bounds
  this->invalidate();
  std::vector<double> lower_value;
  std::vector<double> upper_value;
  if (lp.a_matrix_.isColwise()) {
    lower_value.assign(lp.num_row_, 0);
    upper_value.assign(lp.num_row_, 0);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      const double lower = lp.col_lower_[iCol];
      const double upper = lp.col_upper_[iCol];
      for (HighsInt iEl = lp.a_matrix_.start_[iCol]; iEl < lp.a_matrix_.start_[iCol+1]; iEl++) {
	HighsInt iRow = lp.a_matrix_.index_[iEl];
	double value = lp.a_matrix_.value_[iEl];
	if (value > 0) {
	  lower_value[iRow] += value * lower;
	  upper_value[iRow] += value * upper;
	} else {
	  lower_value[iRow] += value * upper;
	  upper_value[iRow] += value * lower;
	}
      }
    }
  } else {
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      double lower_row_value = 0;
      double upper_row_value = 0;
      for (HighsInt iEl = lp.a_matrix_.start_[iRow]; iEl < lp.a_matrix_.start_[iRow+1]; iEl++) {
	HighsInt iCol = lp.a_matrix_.index_[iEl];
	const double lower = lp.col_lower_[iCol];
	const double upper = lp.col_upper_[iCol];
	double value = lp.a_matrix_.value_[iEl];
	if (value > 0) {
	  lower_row_value += value * lower;
	  upper_row_value += value * upper;
	} else {
	  lower_row_value += value * upper;
	  upper_row_value += value * lower;
	}
      }
      lower_value[iRow] = lower_row_value;
      upper_value[iRow] = upper_row_value;
    }
  }
  bool below_lower;
  bool above_upper;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    below_lower = upper_value[iRow] < lp.row_lower_[iRow] - options.primal_feasibility_tolerance;
    above_upper = lower_value[iRow] > lp.row_upper_[iRow] + options.primal_feasibility_tolerance;
    if (below_lower || above_upper) {
      this->row_index_.push_back(iRow);
      if (below_lower) {
	this->row_bound_.push_back(kIisBoundStatusLower);
      } else {
	this->row_bound_.push_back(kIisBoundStatusUpper);
      }
      break;
    }     
  }
  if (this->row_index_.size() == 0) return false;
  assert(below_lower || above_upper);
  assert(!(below_lower && above_upper));
  // Found an infeasible row
  HighsInt iRow = this->row_index_[0];
  if (below_lower) {
    printf("Row %d has maximum row value of %g, below lower bound of %g\n",
	   int(iRow), upper_value[iRow], lp.row_lower_[iRow]);
  } else {
    printf("Row %d has minimum row value of %g, above upper bound of %g\n",
	   int(iRow), lower_value[iRow], lp.row_upper_[iRow]);
  }
  std::vector<std::tuple<double, HighsInt, double>> row_data;
  double activity = 0;
  // Lambda for adding row data tuples
  auto addRowData = [&](HighsInt iCol, double value) {
    double bound;
    if (below_lower) {
      if (value > 0) {
	bound = lp.col_upper_[iCol];
      } else {
	bound = lp.col_lower_[iCol];
      }
    } else {
      if (value > 0) {
	bound = lp.col_lower_[iCol];
      } else {
	bound = lp.col_upper_[iCol];
      }
    }
    double term = bound * value;
    activity += term;
    printf("addRowData: tuple(%g, %d, %g)\n", term, int(iCol), value);
    row_data.push_back(std::make_tuple(term, iCol, value));
  };
  if (lp.a_matrix_.isColwise()) {
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = lp.a_matrix_.start_[iCol]; iEl < lp.a_matrix_.start_[iCol+1]; iEl++) {
	if (lp.a_matrix_.index_[iEl] == iRow) 
	  addRowData(iCol, lp.a_matrix_.value_[iEl]);
      }
    }
  } else {
    for (HighsInt iEl = lp.a_matrix_.start_[iRow]; iEl < lp.a_matrix_.start_[iRow+1]; iEl++)
      addRowData(lp.a_matrix_.index_[iEl], lp.a_matrix_.value_[iEl]);
  }
  if (below_lower) {
    assert(std::fabs(activity-upper_value[iRow]) < 1e-5);
  } else {
    assert(std::fabs(activity-lower_value[iRow]) < 1e-5);
  }
  std::sort(row_data.begin(), row_data.end());
  HighsInt row_data_size = row_data.size();
  for (HighsInt iisCol = 0; iisCol < row_data_size; iisCol++) {
    
    printf("row_data: tuple(%g, %d, %g)\n",
	   std::get<0>(row_data[iisCol]),
	   int(std::get<1>(row_data[iisCol])),
	   std::get<2>(row_data[iisCol]));
  }

  // Determine whether any entries can be removed
  if (below_lower) {
    const double lower = lp.row_lower_[iRow] - options.primal_feasibility_tolerance;
    // Loop through the terms from most positive to most negative
    // removing them from the activity, until it is within the bound,
    // then use that column and any others in the IIS
    for (HighsInt iisCol = row_data_size-1; iisCol >= 0; iisCol--) {
      const double term = std::get<0>(row_data[iisCol]);
      const HighsInt iCol = std::get<1>(row_data[iisCol]);
      const double value = std::get<2>(row_data[iisCol]);
      double activity_m_term = activity-term;
      printf("Deletion loop: iisCol = %2d; iCol = %3d; term = %11.4g; value = %11.4g; lower = %11.4g; activity = %11.4g; activity - term = %11.4g\n",
	     int(iisCol), int(iCol), term, value, lower, activity, activity_m_term);
      assert(activity < lower);
      if (activity - term <= lower) {
	// Column not in IIS
	activity -= term;
      } else {
	// Column in IIS
	this->col_index_.push_back(iCol);
	double bound;
	if (value > 0) {
	  this->col_bound_.push_back(kIisBoundStatusUpper);
	  bound = lp.col_upper_[iCol];
	} else {
	  this->col_bound_.push_back(kIisBoundStatusLower);
	  bound = lp.col_lower_[iCol];
	}
	if (bound * value != term) {
	  printf("Deletion loop %d: bound(%g) * value(%g) = %g != %g = term\n",
		 int(iisCol), bound, value, bound * value, term);
	}
	assert(bound * value == term);
      }
    }
  } else {
    const double upper = lp.row_upper_[iRow] + options.primal_feasibility_tolerance;
    // Loop through the terms from most negative to most positive,
    // removing them from the activity, until it is within the bound,
    // then use that column and any others in the IIS
    for (HighsInt iisCol = 0; iisCol < row_data_size; iisCol++) {
      const double term = std::get<0>(row_data[iisCol]);
      const HighsInt iCol = std::get<1>(row_data[iisCol]);
      const double value = std::get<2>(row_data[iisCol]);
      printf("Deletion loop: iisCol = %2d; iCol = %3d; term = %11.4g; value = %11.4g; upper = %11.4g; activity = %11.4g; activity - term = %11.4g\n",
	     int(iisCol), int(iCol), term, value, upper, activity, activity-term);
      assert(activity > upper);
      if (activity - term >= upper) {
	// Column not in IIS
	activity -= term;
      } else {
	// Column in IIS
	this->col_index_.push_back(iCol);
	double bound;
	if (value > 0) {
	  this->col_bound_.push_back(kIisBoundStatusLower);
	  bound = lp.col_lower_[iCol];
	} else {
	  this->col_bound_.push_back(kIisBoundStatusUpper);
	  bound = lp.col_upper_[iCol];
	}
	assert(bound * value == term);
      }
    }
  }
  // There must be at least one column in the IIS
  assert(this->col_index_.size() > 0);
  this->valid_ = true;
  return this->valid_;
}

bool HighsIis::ok(const HighsLp& lp, const HighsOptions& options) const {
  if (!this->valid_) return true;
  const HighsLogOptions& log_options = options.log_options;
  Highs h;
  h.passOptions(options);
  h.passModel(lp);
  h.run();
  if (h.getModelStatus() != HighsModelStatus::kInfeasible) {
    highsLogUser(log_options, HighsLogType::kError, "HighsIis: Given LP is not infeasible\n");
    return false;
  }
  auto optimalOrUnbounded = [&]() -> bool {
    h.writeModel("");
    h.run();
    return
      h.getModelStatus() == HighsModelStatus::kOptimal ||
      h.getModelStatus() == HighsModelStatus::kUnbounded;
  };
  HighsInt num_iis_col = this->col_index_.size();
  HighsInt num_iis_row = this->row_index_.size();
  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = this->col_index_[iisCol];
    if (this->col_bound_[iisCol] == kIisBoundStatusLower ||
	this->col_bound_[iisCol] == kIisBoundStatusBoxed) {
      h.changeColBounds(iCol, -kHighsInf, lp.col_upper_[iCol]);
      if (!optimalOrUnbounded()) {
	highsLogUser(log_options, HighsLogType::kError,
		     "HighsIis: IIS column %d (LP column %d): relaxing lower bound of %g does not yield LP with optimal or unbounded status\n",
		     int(iisCol), int(iCol), lp.col_lower_[iCol]);
	return false;
      }
      h.changeColBounds(iCol, lp.col_lower_[iCol], lp.col_upper_[iCol]);
    }
    if (this->col_bound_[iisCol] == kIisBoundStatusUpper ||
	this->col_bound_[iisCol] == kIisBoundStatusBoxed) {
      h.changeColBounds(iCol, lp.col_lower_[iCol], kHighsInf);
      if (!optimalOrUnbounded()) {
	highsLogUser(log_options, HighsLogType::kError,
		     "HighsIis: IIS column %d (LP column %d): relaxing upper bound of %g does not yield LP with optimal or unbounded status\n",
		     int(iisCol), int(iCol), lp.col_upper_[iCol]);
	return false;
      }
      h.changeColBounds(iCol, lp.col_lower_[iCol], lp.col_upper_[iCol]);
    }
  }
  return true;
  
  
}
