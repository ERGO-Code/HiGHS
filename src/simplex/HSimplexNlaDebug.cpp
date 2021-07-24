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
/**@file simplex/HSimplexNlaDebug.cpp
 *
 * @brief Debugging code for simplex NLA
 */
#include "simplex/HSimplexNlaDebug.h"

#include "util/HighsRandom.h"

//#include <stdio.h>

const double kSolveLargeError = 1e-12;
const double kSolveExcessiveError = sqrt(kSolveLargeError);

const double kInverseLargeError = 1e-12;
const double kInverseExcessiveError = sqrt(kInverseLargeError);

HighsDebugStatus debugCheckInvert(const HSimplexNla& simplex_nla,
                                  const bool force) {
  if (simplex_nla.options_->highs_debug_level < kHighsDebugLevelCostly &&
      !force)
    return HighsDebugStatus::kNotChecked;
  if (force)
    highsLogDev(simplex_nla.options_->log_options, HighsLogType::kInfo,
                "CheckNlaINVERT:   Forcing debug\n");

  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;
  return_status = HighsDebugStatus::kOk;

  const HighsInt num_row = simplex_nla.lp_->num_row_;
  const HighsInt num_col = simplex_nla.lp_->num_col_;
  const vector<HighsInt>& a_start = simplex_nla.lp_->a_start_;
  const vector<HighsInt>& a_index = simplex_nla.lp_->a_index_;
  const vector<double>& a_value = simplex_nla.lp_->a_value_;
  const HighsInt* base_index = simplex_nla.base_index_;
  const HighsOptions* options = simplex_nla.options_;
  HVector column;
  HVector rhs;
  column.setup(num_row);
  rhs.setup(num_row);
  double expected_density = 1;

  // Solve for a random solution
  HighsRandom random;
  column.clear();
  rhs.clear();
  column.count = -1;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    rhs.index[rhs.count++] = iRow;
    double value = random.fraction();
    column.array[iRow] = value;
    HighsInt iCol = base_index[iRow];
    if (iCol < num_col) {
      for (HighsInt iEl = a_start[iCol]; iEl < a_start[iCol + 1]; iEl++) {
        HighsInt index = a_index[iEl];
        rhs.array[index] += value * a_value[iEl];
      }
    } else {
      HighsInt index = iCol - num_col;
      rhs.array[index] += value;
    }
  }
  simplex_nla.ftran(rhs, expected_density);

  double solve_error_norm = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double solve_error = fabs(rhs.array[iRow] - column.array[iRow]);
    solve_error_norm = std::max(solve_error, solve_error_norm);
  }
  std::string value_adjective;
  HighsLogType report_level;
  return_status = HighsDebugStatus::kOk;

  if (solve_error_norm) {
    if (solve_error_norm > kSolveExcessiveError) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (solve_error_norm > kSolveLargeError) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }

    if (force) report_level = HighsLogType::kInfo;

    highsLogDev(
        options->log_options, report_level,
        "CheckINVERT:   %-9s (%9.4g) norm for random solution solve error\n",
        value_adjective.c_str(), solve_error_norm);
  }

  if (options->highs_debug_level < kHighsDebugLevelExpensive)
    return return_status;

  expected_density = 0;
  double inverse_error_norm = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iCol = base_index[iRow];
    column.clear();
    column.packFlag = true;
    if (iCol < num_col) {
      for (HighsInt k = a_start[iCol]; k < a_start[iCol + 1]; k++) {
        HighsInt index = a_index[k];
        column.array[index] = a_value[k];
        column.index[column.count++] = index;
      }
    } else {
      HighsInt index = iCol - num_col;
      column.array[index] = 1.0;
      column.index[column.count++] = index;
    }

    simplex_nla.ftran(column, expected_density);

    double inverse_column_error_norm = 0;
    for (HighsInt lc_iRow = 0; lc_iRow < num_row; lc_iRow++) {
      double value = column.array[lc_iRow];
      double ckValue;
      if (lc_iRow == iRow) {
        ckValue = 1;
      } else {
        ckValue = 0;
      }
      double inverse_error = fabs(value - ckValue);
      inverse_column_error_norm =
          std::max(inverse_error, inverse_column_error_norm);
    }
    inverse_error_norm =
        std::max(inverse_column_error_norm, inverse_error_norm);
  }
  if (inverse_error_norm) {
    if (inverse_error_norm > kInverseExcessiveError) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (inverse_error_norm > kInverseLargeError) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }
    highsLogDev(options->log_options, report_level,
                "CheckINVERT:   %-9s (%9.4g) norm for inverse error\n",
                value_adjective.c_str(), inverse_error_norm);
  }

  return return_status;
}
