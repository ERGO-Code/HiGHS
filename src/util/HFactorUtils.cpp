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
/**@file util/HFactorUtils.cpp
 * @brief Types of solution classes
 */
#include "util/HFactor.h"

void HFactor::invalidAMatrixAction() {
  this->a_matrix_valid = false;
  refactor_info_.clear();
}

void HFactor::reportLu(const HighsInt l_u_or_both, const bool full) const {
  if (l_u_or_both < kReportLuJustL || l_u_or_both > kReportLuBoth) return;
  if (l_u_or_both & 1) {
    printf("L");
    if (full) printf(" - full");
    printf(":\n");

    if (full) reportIntVector("l_pivot_lookup", l_pivot_lookup);
    if (full) reportIntVector("l_pivot_index", l_pivot_index);
    reportIntVector("l_start", l_start);
    reportIntVector("l_index", l_index);
    reportDoubleVector("l_value", l_value);
    if (full) {
      reportIntVector("lr_start", lr_start);
      reportIntVector("lr_index", lr_index);
      reportDoubleVector("lr_value", lr_value);
    }
  }
  if (l_u_or_both & 2) {
    printf("U");
    if (full) printf(" - full");
    printf(":\n");
    if (full) reportIntVector("u_pivot_lookup", u_pivot_lookup);
    reportIntVector("u_pivot_index", u_pivot_index);
    reportDoubleVector("u_pivot_value", u_pivot_value);
    reportIntVector("u_start", u_start);
    if (full) reportIntVector("u_last_p", u_last_p);
    reportIntVector("u_index", u_index);
    reportDoubleVector("u_value", u_value);
    if (full) {
      reportIntVector("ur_start", ur_start);
      reportIntVector("ur_lastp", ur_lastp);
      reportIntVector("ur_space", ur_space);
      for (HighsInt iRow = 0; iRow < ur_start.size(); iRow++) {
        const HighsInt start = ur_start[iRow];
        const HighsInt end = ur_lastp[iRow];
        if (start >= end) continue;
        printf("UR    Row %2d: ", (int)iRow);
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11d ", (int)ur_index[iEl]);
        printf("\n              ");
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11.4g ", ur_value[iEl]);
        printf("\n");
      }
      //      reportIntVector("ur_index", ur_index);
      //      reportDoubleVector("ur_value", ur_value);
    }
  }
  if (l_u_or_both == 3 && full) {
    reportDoubleVector("PFpivotValue", PFpivotValue);
    reportIntVector("PFpivotIndex", PFpivotIndex);
    reportIntVector("PFstart", PFstart);
    reportIntVector("PFindex", PFindex);
    reportDoubleVector("PFvalue", PFvalue);
  }
}

void HFactor::reportIntVector(const std::string name,
                              const vector<HighsInt> entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11d ", (int)entry[iEn]);
  }
  printf("\n");
}
void HFactor::reportDoubleVector(const std::string name,
                                 const vector<double> entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11.4g ", entry[iEn]);
  }
  printf("\n");
}
