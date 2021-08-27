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
/**@file simplex/HFactorUtils.cpp
 * @brief Types of solution classes
 */
#include "simplex/HFactor.h"

void HFactor::invalidAMatrixAction() {
  this->a_matrix_valid = false;
  refactor_info_.clear();
}

void HFactor::reportLu(const HighsInt l_u_or_both, const bool full) const {
  if (l_u_or_both<kReportLuJustL || l_u_or_both>kReportLuBoth) return;
  if (l_u_or_both & 1) {
    printf("L");
    if (full) printf(" - full");
    printf(":\n");

    if (full) reportIntVector("LpivotLookup", LpivotLookup);
    if (full) reportIntVector("LpivotIndex", LpivotIndex);
    reportIntVector("Lstart", Lstart);
    reportIntVector("Lindex", Lindex);
    reportDoubleVector("Lvalue", Lvalue);
    if (full) {
      reportIntVector("LRstart", LRstart);
      reportIntVector("LRindex", LRindex);
      reportDoubleVector("LRvalue", LRvalue);
    }
  }
  if (l_u_or_both & 2) {
    printf("U");
    if (full) printf(" - full");
    printf(":\n");
    if (full) reportIntVector("UpivotLookup", UpivotLookup);
    reportIntVector("UpivotIndex", UpivotIndex);
    reportDoubleVector("UpivotValue", UpivotValue);
    reportIntVector("Ustart", Ustart);
    if (full) reportIntVector("Ulastp", Ulastp);
    reportIntVector("Uindex", Uindex);
    reportDoubleVector("Uvalue", Uvalue);
    if (full) {
      reportIntVector("URstart", URstart);
      reportIntVector("URlastp", URlastp);
      reportIntVector("URspace", URspace);
      reportIntVector("URindex", URindex);
      reportDoubleVector("URvalue", URvalue);
    }
  }
}

void HFactor::reportIntVector(const std::string name,
                              const vector<HighsInt> entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0) printf("\n                                  ");
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
    if (iEn > 0 && iEn % 10 == 0) printf("\n                                  ");
    printf("%11.4g ", entry[iEn]);
  }
  printf("\n");
}
