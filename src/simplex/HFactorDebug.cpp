/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HFactorDebug.cpp
 * @brief
 */

#include "simplex/HFactorDebug.h"

#include "simplex/HVector.h"
#include "util/HighsRandom.h"

const double solve_large_error = 1e-12;
const double solve_excessive_error = sqrt(solve_large_error);

const double inverse_large_error = 1e-12;
const double inverse_excessive_error = sqrt(inverse_large_error);

HighsDebugStatus debugCheckInvert(const HighsOptions& options,
                                  const HFactor& factor, const bool force) {
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;
  if (force)
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "CheckINVERT:   Forcing debug\n");

  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;
  return_status = HighsDebugStatus::kOk;
  const HighsInt numRow = factor.numRow;
  const HighsInt numCol = factor.numCol;
  const HighsInt* Astart = factor.getAstart();
  const HighsInt* Aindex = factor.getAindex();
  const double* Avalue = factor.getAvalue();
  const HighsInt* baseIndex = factor.getBaseIndex();

  HVector column;
  HVector rhs;
  column.setup(numRow);
  rhs.setup(numRow);
  double rhsDensity = 1;

  // Solve for a random solution
  HighsRandom random;
  column.clear();
  rhs.clear();
  column.count = -1;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    rhs.index[rhs.count++] = iRow;
    double value = random.fraction();
    column.array[iRow] = value;
    HighsInt iCol = baseIndex[iRow];
    if (iCol < numCol) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt index = Aindex[k];
        rhs.array[index] += value * Avalue[k];
      }
    } else {
      HighsInt index = iCol - numCol;
      rhs.array[index] += value;
    }
  }
  factor.ftran(rhs, rhsDensity);
  double solve_error_norm = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    double solve_error = fabs(rhs.array[iRow] - column.array[iRow]);
    solve_error_norm = std::max(solve_error, solve_error_norm);
  }
  std::string value_adjective;
  HighsLogType report_level;
  return_status = HighsDebugStatus::kOk;

  if (solve_error_norm) {
    if (solve_error_norm > solve_excessive_error) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (solve_error_norm > solve_large_error) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }

    if (force) report_level = HighsLogType::kInfo;

    highsLogDev(
        options.log_options, report_level,
        "CheckINVERT:   %-9s (%9.4g) norm for random solution solve error\n",
        value_adjective.c_str(), solve_error_norm);
  }

  if (options.highs_debug_level < kHighsDebugLevelExpensive)
    return return_status;

  double columnDensity = 0;
  double inverse_error_norm = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    HighsInt iCol = baseIndex[iRow];
    column.clear();
    column.packFlag = true;
    if (iCol < numCol) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt index = Aindex[k];
        column.array[index] = Avalue[k];
        column.index[column.count++] = index;
      }
    } else {
      HighsInt index = iCol - numCol;
      column.array[index] = 1.0;
      column.index[column.count++] = index;
    }
    factor.ftran(column, columnDensity);
    double inverse_column_error_norm = 0;
    for (HighsInt lc_iRow = 0; lc_iRow < numRow; lc_iRow++) {
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
    if (inverse_error_norm > inverse_excessive_error) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (inverse_error_norm > inverse_large_error) {
      value_adjective = "Large";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "Small";
      report_level = HighsLogType::kInfo;
    }
    highsLogDev(options.log_options, report_level,
                "CheckINVERT:   %-9s (%9.4g) norm for inverse error\n",
                value_adjective.c_str(), inverse_error_norm);
  }

  return return_status;
}

void debugReportRankDeficiency(
    const HighsInt call_id, const HighsInt highs_debug_level,
    const HighsLogOptions& log_options, const HighsInt numRow,
    const vector<HighsInt>& permute, const vector<HighsInt>& iwork,
    const HighsInt* baseIndex, const HighsInt rank_deficiency,
    const vector<HighsInt>& noPvR, const vector<HighsInt>& noPvC) {
  if (highs_debug_level == kHighsDebugLevelNone) return;
  if (call_id == 0) {
    if (numRow > 123) return;
    highsLogDev(log_options, HighsLogType::kWarning, "buildRankDeficiency0:");
    highsLogDev(log_options, HighsLogType::kWarning, "\nIndex  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\nPerm   ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  permute[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\nIwork  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  iwork[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\nBaseI  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  baseIndex[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
  } else if (call_id == 1) {
    if (rank_deficiency > 100) return;
    highsLogDev(log_options, HighsLogType::kWarning, "buildRankDeficiency1:");
    highsLogDev(log_options, HighsLogType::kWarning, "\nIndex  ");
    for (HighsInt i = 0; i < rank_deficiency; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\nnoPvR  ");
    for (HighsInt i = 0; i < rank_deficiency; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  noPvR[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\nnoPvC  ");
    for (HighsInt i = 0; i < rank_deficiency; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  noPvC[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
    if (numRow > 123) return;
    highsLogDev(log_options, HighsLogType::kWarning, "Index  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\nIwork  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  iwork[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
  } else if (call_id == 2) {
    if (numRow > 123) return;
    highsLogDev(log_options, HighsLogType::kWarning, "buildRankDeficiency2:");
    highsLogDev(log_options, HighsLogType::kWarning, "\nIndex  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\nPerm   ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  permute[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
  }
}

void debugReportRankDeficientASM(
    const HighsInt highs_debug_level, const HighsLogOptions& log_options,
    const HighsInt numRow, const vector<HighsInt>& MCstart,
    const vector<HighsInt>& MCcountA, const vector<HighsInt>& MCindex,
    const vector<double>& MCvalue, const vector<HighsInt>& iwork,
    const HighsInt rank_deficiency, const vector<HighsInt>& noPvC,
    const vector<HighsInt>& noPvR) {
  if (highs_debug_level == kHighsDebugLevelNone) return;
  if (rank_deficiency > 10) return;
  double* ASM;
  ASM = (double*)malloc(sizeof(double) * rank_deficiency * rank_deficiency);
  for (HighsInt i = 0; i < rank_deficiency; i++) {
    for (HighsInt j = 0; j < rank_deficiency; j++) {
      ASM[i + j * rank_deficiency] = 0;
    }
  }
  for (HighsInt j = 0; j < rank_deficiency; j++) {
    HighsInt ASMcol = noPvC[j];
    HighsInt start = MCstart[ASMcol];
    HighsInt end = start + MCcountA[ASMcol];
    for (HighsInt en = start; en < end; en++) {
      HighsInt ASMrow = MCindex[en];
      HighsInt i = -iwork[ASMrow] - 1;
      if (i < 0 || i >= rank_deficiency) {
        highsLogDev(log_options, HighsLogType::kWarning,
                    "STRANGE: 0 > i = %" HIGHSINT_FORMAT " || %" HIGHSINT_FORMAT
                    " = i >= rank_deficiency = %" HIGHSINT_FORMAT "\n",
                    i, i, rank_deficiency);
      } else {
        if (noPvR[i] != ASMrow) {
          highsLogDev(log_options, HighsLogType::kWarning,
                      "STRANGE: %" HIGHSINT_FORMAT
                      " = noPvR[i] != ASMrow = %" HIGHSINT_FORMAT "\n",
                      noPvR[i], ASMrow);
        }
        highsLogDev(log_options, HighsLogType::kWarning,
                    "Setting ASM(%2" HIGHSINT_FORMAT ", %2" HIGHSINT_FORMAT
                    ") = %11.4g\n",
                    i, j, MCvalue[en]);
        ASM[i + j * rank_deficiency] = MCvalue[en];
      }
    }
  }
  highsLogDev(log_options, HighsLogType::kWarning, "ASM:                    ");
  for (HighsInt j = 0; j < rank_deficiency; j++)
    highsLogDev(log_options, HighsLogType::kWarning, " %11" HIGHSINT_FORMAT "",
                j);
  highsLogDev(log_options, HighsLogType::kWarning,
              "\n                        ");
  for (HighsInt j = 0; j < rank_deficiency; j++)
    highsLogDev(log_options, HighsLogType::kWarning, " %11" HIGHSINT_FORMAT "",
                noPvC[j]);
  highsLogDev(log_options, HighsLogType::kWarning,
              "\n                        ");
  for (HighsInt j = 0; j < rank_deficiency; j++)
    highsLogDev(log_options, HighsLogType::kWarning, "------------");
  highsLogDev(log_options, HighsLogType::kWarning, "\n");
  for (HighsInt i = 0; i < rank_deficiency; i++) {
    highsLogDev(log_options, HighsLogType::kWarning,
                "%11" HIGHSINT_FORMAT " %11" HIGHSINT_FORMAT "|", i, noPvR[i]);
    for (HighsInt j = 0; j < rank_deficiency; j++) {
      highsLogDev(log_options, HighsLogType::kWarning, " %11.4g",
                  ASM[i + j * rank_deficiency]);
    }
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
  }
  free(ASM);
}

void debugReportMarkSingC(const HighsInt call_id,
                          const HighsInt highs_debug_level,
                          const HighsLogOptions& log_options,
                          const HighsInt numRow, const vector<HighsInt>& iwork,
                          const HighsInt* baseIndex) {
  if (highs_debug_level == kHighsDebugLevelNone) return;
  if (numRow > 123) return;
  if (call_id == 0) {
    highsLogDev(log_options, HighsLogType::kWarning, "\nMarkSingC1");
    highsLogDev(log_options, HighsLogType::kWarning, "\nIndex  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\niwork  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  iwork[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\nBaseI  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  baseIndex[i]);
  } else if (call_id == 1) {
    highsLogDev(log_options, HighsLogType::kWarning, "\nMarkSingC2");
    highsLogDev(log_options, HighsLogType::kWarning, "\nIndex  ");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  i);
    highsLogDev(log_options, HighsLogType::kWarning, "\nNwBaseI");
    for (HighsInt i = 0; i < numRow; i++)
      highsLogDev(log_options, HighsLogType::kWarning, " %2" HIGHSINT_FORMAT "",
                  baseIndex[i]);
    highsLogDev(log_options, HighsLogType::kWarning, "\n");
  }
}

void debugLogRankDeficiency(
    const HighsInt highs_debug_level, const HighsLogOptions& log_options,
    const HighsInt rank_deficiency, const HighsInt basis_matrix_num_el,
    const HighsInt invert_num_el, const HighsInt& kernel_dim,
    const HighsInt kernel_num_el, const HighsInt nwork) {
  if (highs_debug_level == kHighsDebugLevelNone) return;
  if (!rank_deficiency) return;
  highsLogDev(
      log_options, HighsLogType::kWarning,
      "Rank deficiency %1" HIGHSINT_FORMAT ": basis_matrix (%" HIGHSINT_FORMAT
      " el); INVERT (%" HIGHSINT_FORMAT " el); kernel (%" HIGHSINT_FORMAT
      " "
      "dim; %" HIGHSINT_FORMAT " el): nwork = %" HIGHSINT_FORMAT "\n",
      rank_deficiency, basis_matrix_num_el, invert_num_el, kernel_dim,
      kernel_num_el, nwork);
}

void debugPivotValueAnalysis(const HighsInt highs_debug_level,
                             const HighsLogOptions& log_options,
                             const HighsInt numRow,
                             const vector<double>& UpivotValue) {
  if (highs_debug_level < kHighsDebugLevelCheap) return;
  double min_pivot = kHighsInf;
  double mean_pivot = 0;
  double max_pivot = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    double abs_pivot = fabs(UpivotValue[iRow]);
    min_pivot = min(abs_pivot, min_pivot);
    max_pivot = max(abs_pivot, max_pivot);
    mean_pivot += log(abs_pivot);
  }
  mean_pivot = exp(mean_pivot / numRow);
  if (highs_debug_level > kHighsDebugLevelCheap || min_pivot < 1e-8)
    highsLogDev(log_options, HighsLogType::kError,
                "InvertPivotAnalysis: %" HIGHSINT_FORMAT
                " pivots: Min %g; Mean "
                "%g; Max %g\n",
                numRow, min_pivot, mean_pivot, max_pivot);
}
