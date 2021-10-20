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
/**@file lp_data/HFactorDebug.cpp
 * @brief
 */

#include "simplex/HFactorDebug.h"

#include <algorithm>
#include <cmath>

#include "simplex/HVector.h"
#include "util/HighsRandom.h"

using std::fabs;
using std::max;
using std::min;

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
