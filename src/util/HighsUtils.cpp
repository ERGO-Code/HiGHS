/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "util/HighsUtils.h"

#include <cmath>
#include <cstdio>
#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h"

/*
int getOmpNumThreads() {
  return omp_get_num_threads()
}
*/

double getNorm2(const std::vector<double> values) {
  double sum = 0;
  int values_size = values.size();
  for (int i = 0; i < values_size; i++) sum += values[i] * values[i];
  return sum;
}

bool highs_isInfinity(double val) {
  if (val >= HIGHS_CONST_INF) return true;
  return false;
}

double highs_relative_difference(const double v0, const double v1) {
  return fabs(v0 - v1) / std::max(v0, std::max(v1, 1.0));
}

#ifdef HiGHSDEV
void analyseVectorValues(const char* message, int vecDim,
                         const std::vector<double>& vec, bool analyseValueList,
                         std::string model_name) {
  if (vecDim == 0) return;
  double log10 = log(10.0);
  const int nVK = 20;
  int nNz = 0;
  int nPosInfV = 0;
  int nNegInfV = 0;
  std::vector<int> posVK;
  std::vector<int> negVK;
  posVK.resize(nVK + 1, 0);
  negVK.resize(nVK + 1, 0);

  const int VLsMxZ = 10;
  std::vector<int> VLsK;
  std::vector<double> VLsV;
  VLsK.resize(VLsMxZ, 0);
  VLsV.resize(VLsMxZ, 0);
  // Ensure that 1.0 and -1.0 are counted
  const int PlusOneIx = 0;
  const int MinusOneIx = 1;
  bool excessVLsV = false;
  int VLsZ = 2;
  VLsV[PlusOneIx] = 1.0;
  VLsV[MinusOneIx] = -1.0;

  for (int ix = 0; ix < vecDim; ix++) {
    double v = vec[ix];
    double absV = std::fabs(v);
    int log10V;
    if (absV > 0) {
      // Nonzero value
      nNz++;
      if (highs_isInfinity(-v)) {
        //-Inf value
        nNegInfV++;
      } else if (highs_isInfinity(v)) {
        //+Inf value
        nPosInfV++;
      } else {
        // Finite nonzero value
        if (absV == 1) {
          log10V = 0;
        } else if (absV == 10) {
          log10V = 1;
        } else if (absV == 100) {
          log10V = 2;
        } else if (absV == 1000) {
          log10V = 3;
        } else {
          log10V = log(absV) / log10;
        }
        if (log10V >= 0) {
          int k = std::min(log10V, nVK);
          posVK[k]++;
        } else {
          int k = std::min(-log10V, nVK);
          negVK[k]++;
        }
      }
    }
    if (analyseValueList) {
      if (v == 1.0) {
        VLsK[PlusOneIx]++;
      } else if (v == -1.0) {
        VLsK[MinusOneIx]++;
      } else {
        int fdIx = -1;
        for (int ix = 2; ix < VLsZ; ix++) {
          if (v == VLsV[ix]) {
            fdIx = ix;
            break;
          }
        }
        if (fdIx == -1) {
          // New value
          if (VLsZ < VLsMxZ) {
            fdIx = VLsZ;
            VLsV[fdIx] = v;
            VLsK[fdIx]++;
            VLsZ++;
          } else {
            excessVLsV = true;
          }
        } else {
          // Existing value
          VLsK[fdIx]++;
        }
      }
    }
  }
  printf("%s of dimension %d with %d nonzeros (%3d%%): Analysis\n", message,
         vecDim, nNz, 100 * nNz / vecDim);
  if (nNegInfV > 0) printf("%12d values are -Inf\n", nNegInfV);
  if (nPosInfV > 0) printf("%12d values are +Inf\n", nPosInfV);
  int k = nVK;
  int vK = posVK[k];
  if (vK > 0) printf("%12d values satisfy 10^(%3d) <= v < Inf\n", vK, k);
  for (int k = nVK - 1; k >= 0; k--) {
    int vK = posVK[k];
    if (vK > 0)
      printf("%12d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, k, k + 1);
  }
  for (int k = 1; k <= nVK; k++) {
    int vK = negVK[k];
    if (vK > 0)
      printf("%12d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, -k, 1 - k);
  }
  vK = vecDim - nNz;
  if (vK > 0) printf("%12d values are zero\n", vK);
  if (analyseValueList) {
    printf("           Value distribution:");
    if (excessVLsV) printf(" More than %d different values", VLsZ);
    printf("\n            Value        Count\n");
    for (int ix = 0; ix < VLsZ; ix++) {
      int pct = ((100.0 * VLsK[ix]) / vecDim) + 0.5;
      printf("     %12g %12d (%3d%%)\n", VLsV[ix], VLsK[ix], pct);
    }
    printf("grep_value_distrib,%s,%d", model_name.c_str(), VLsZ);
    printf(",");
    if (excessVLsV) printf("!");
    for (int ix = 0; ix < VLsZ; ix++) printf(",%g", VLsV[ix]);
    printf("\n");
  }
}

void analyseMatrixSparsity(const char* message, int numCol, int numRow,
                           const std::vector<int>& Astart,
                           const std::vector<int>& Aindex) {
  if (numCol == 0) return;
  std::vector<int> rowCount;
  std::vector<int> colCount;

  rowCount.assign(numRow, 0);
  colCount.resize(numCol);

  for (int col = 0; col < numCol; col++) {
    colCount[col] = Astart[col + 1] - Astart[col];
    for (int el = Astart[col]; el < Astart[col + 1]; el++)
      rowCount[Aindex[el]]++;
  }
  const int maxCat = 10;
  std::vector<int> CatV;
  std::vector<int> rowCatK;
  std::vector<int> colCatK;
  CatV.resize(maxCat + 1);
  rowCatK.assign(maxCat + 1, 0);
  colCatK.assign(maxCat + 1, 0);

  CatV[1] = 1;
  for (int cat = 2; cat < maxCat + 1; cat++) {
    CatV[cat] = 2 * CatV[cat - 1];
  }

  int maxRowCount = 0;
  int maxColCount = 0;
  for (int col = 0; col < numCol; col++) {
    maxColCount = std::max(colCount[col], maxColCount);
    int fdCat = maxCat;
    for (int cat = 0; cat < maxCat - 1; cat++) {
      if (colCount[col] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    colCatK[fdCat]++;
  }

  for (int row = 0; row < numRow; row++) {
    maxRowCount = std::max(rowCount[row], maxRowCount);
    int fdCat = maxCat;
    for (int cat = 0; cat < maxCat - 1; cat++) {
      if (rowCount[row] < CatV[cat + 1]) {
        fdCat = cat;
        break;
      }
    }
    rowCatK[fdCat]++;
  }

  printf("\n%s\n\n", message);
  int lastRpCat;
  for (int cat = 0; cat < maxCat + 1; cat++) {
    if (colCatK[cat]) lastRpCat = cat;
  }
  int cat = maxCat;
  if (colCatK[cat]) lastRpCat = cat;
  int sumK = 0;
  int pct;
  double v;
  int sumPct = 0;
  for (int cat = 0; cat < lastRpCat; cat++) {
    sumK += colCatK[cat];
    v = 100 * colCatK[cat];
    v = v / numCol + 0.5;
    pct = v;
    sumPct += pct;
    printf("%12d (%3d%%) columns of count in [%3d, %3d]\n", colCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  sumK += colCatK[cat];
  v = 100 * colCatK[cat];
  v = v / numCol + 0.5;
  pct = v;
  sumPct += pct;
  if (cat == maxCat) {
    printf("%12d (%3d%%) columns of count in [%3d, inf]\n", colCatK[cat], pct,
           CatV[cat]);
  } else {
    printf("%12d (%3d%%) columns of count in [%3d, %3d]\n", colCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }
  printf("Max count is %d / %d\n\n", maxColCount, numRow);

  lastRpCat = -1;
  for (int cat = 0; cat < maxCat + 1; cat++) {
    if (rowCatK[cat]) lastRpCat = cat;
  }
  cat = maxCat;
  if (rowCatK[cat]) lastRpCat = cat;
  sumK = 0;
  pct = 0;
  v = 0;
  sumPct = 0;
  for (int cat = 0; cat < lastRpCat; cat++) {
    sumK += rowCatK[cat];
    v = 100 * rowCatK[cat];
    v = v / numRow + 0.5;
    pct = v;
    sumPct += pct;
    printf("%12d (%3d%%)    rows of count in [%3d, %3d]\n", rowCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }

  cat = lastRpCat;
  sumK += rowCatK[cat];
  v = 100 * rowCatK[cat];
  v = v / numRow + 0.5;
  pct = v;
  sumPct += pct;
  if (cat == maxCat) {
    printf("%12d (%3d%%)    rows of count in [%3d, inf]\n", rowCatK[cat], pct,
           CatV[cat]);
  } else {
    printf("%12d (%3d%%)    rows of count in [%3d, %3d]\n", rowCatK[cat], pct,
           CatV[cat], CatV[cat + 1] - 1);
  }
  printf("Max count is %d / %d\n", maxRowCount, numCol);
}

bool initialiseValueDistribution(const double min_value_limit,
                                 const double max_value_limit,
                                 const double base_value_limit,
                                 HighsValueDistribution& value_distribution) {
  if (min_value_limit <= 0) return false;
  if (max_value_limit < min_value_limit) return false;
  int num_count;
  if (min_value_limit == max_value_limit) {
    // For counting values below and above a value
    num_count = 1;
  } else {
    if (base_value_limit <= 0) return false;
    const double log_ratio = log(max_value_limit / min_value_limit);
    const double log_base_value_limit = log(base_value_limit);
    //    printf("initialiseValueDistribution: log_ratio = %g;
    //    log_base_value_limit = %g; log_ratio/log_base_value_limit = %g\n",
    //	   log_ratio, log_base_value_limit, log_ratio/log_base_value_limit);
    num_count = log_ratio / log_base_value_limit + 1;
  }
  //  printf("initialiseValueDistribution: num_count = %d\n", num_count);
  value_distribution.count_.assign(num_count + 1, 0);
  value_distribution.limit_.assign(num_count, 0);
  value_distribution.limit_[0] = min_value_limit;
  //  printf("Interval  0 is [%10.4g, %10.4g)\n", 0.0,
  //  value_distribution.limit_[0]);
  for (int i = 1; i < num_count; i++) {
    value_distribution.limit_[i] =
        base_value_limit * value_distribution.limit_[i - 1];
    //    printf("Interval %2d is [%10.4g, %10.4g)\n", i,
    //    value_distribution.limit_[i-1], value_distribution.limit_[i]);
  }
  //  printf("Interval %2d is [%10.4g, inf)\n", num_count,
  //  value_distribution.limit_[num_count-1]);
  value_distribution.num_count_ = num_count;
  value_distribution.num_zero_ = 0;
  value_distribution.num_one_ = 0;
  value_distribution.min_value_ = HIGHS_CONST_INF;
  value_distribution.max_value_ = 0;
  return true;
}

bool updateValueDistribution(const double value,
                             HighsValueDistribution& value_distribution) {
  if (value_distribution.num_count_ < 0) return false;
  const double abs_value = fabs(value);
  value_distribution.min_value_ =
      std::min(abs_value, value_distribution.min_value_);
  value_distribution.max_value_ =
      std::max(abs_value, value_distribution.max_value_);
  if (!abs_value) {
    value_distribution.num_zero_++;
    return true;
  }
  if (abs_value == 1.0) {
    value_distribution.num_one_++;
    return true;
  }
  for (int i = 0; i < value_distribution.num_count_; i++) {
    if (abs_value < value_distribution.limit_[i]) {
      value_distribution.count_[i]++;
      return true;
    }
  }
  value_distribution.count_[value_distribution.num_count_]++;
  return true;
}

double doublePercentage(const int of, const int in) {
  return ((100.0 * of) / in);
}

int integerPercentage(const int of, const int in) {
  const double double_percentage = ((100.0 * of) / in) + 0.4999;
  return (int)double_percentage;
}

bool printValueDistribution(std::string value_name,
                            const HighsValueDistribution& value_distribution,
                            const int mu) {
  const int num_count = value_distribution.num_count_;
  if (num_count < 0) return false;
  bool not_reported_ones = true;
  int sum_count = value_distribution.num_zero_ + value_distribution.num_one_;
  double sum_percentage = 0;
  const double min_value = value_distribution.min_value_;
  for (int i = 0; i < num_count + 1; i++)
    sum_count += value_distribution.count_[i];
  if (!sum_count) return false;
  printf("     Minimum %svalue is %10.4g", value_name.c_str(), min_value);
  if (mu > 0) {
    printf("  corresponding to  %10d / %10d\n", (int)(min_value * mu), mu);
  } else {
    printf("\n");
  }
  printf("     Maximum %svalue is %10.4g", value_name.c_str(),
         value_distribution.max_value_);
  if (mu > 0) {
    printf("  corresponding to  %10d / %10d\n",
           (int)(value_distribution.max_value_ * mu), mu);
  } else {
    printf("\n");
  }
  int sum_report_count = 0;
  double percentage;
  int int_percentage;
  int count = value_distribution.num_zero_;
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) are %10.4g\n", count, value_name.c_str(),
           int_percentage, 0.0);
    sum_report_count += count;
  }
  count = value_distribution.count_[0];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) in (%10.4g, %10.4g)", count,
           value_name.c_str(), int_percentage, 0.0,
           value_distribution.limit_[0]);
    sum_report_count += count;
    if (mu > 0) {
      printf(" corresponding to (%10d, %10d)\n", 0,
             (int)(value_distribution.limit_[0] * mu));
    } else {
      printf("\n");
    }
  }
  for (int i = 1; i < num_count; i++) {
    if (not_reported_ones && value_distribution.limit_[i - 1] >= 1.0) {
      count = value_distribution.num_one_;
      if (count) {
        percentage = doublePercentage(count, sum_count);
        sum_percentage += percentage;
        int_percentage = percentage;
        printf("%12d %svalues (%3d%%) are             %10.4g", count,
               value_name.c_str(), int_percentage, 1.0);
        sum_report_count += count;
        if (mu > 0) {
          printf(" corresponding to %10d\n", mu);
        } else {
          printf("\n");
        }
      }
      not_reported_ones = false;
    }
    count = value_distribution.count_[i];
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) in [%10.4g, %10.4g)", count,
             value_name.c_str(), int_percentage,
             value_distribution.limit_[i - 1], value_distribution.limit_[i]);
      sum_report_count += count;
      if (mu > 0) {
        printf(" corresponding to [%10d, %10d)\n",
               (int)(value_distribution.limit_[i - 1] * mu),
               (int)(value_distribution.limit_[i] * mu));
      } else {
        printf("\n");
      }
    }
  }
  if (not_reported_ones && value_distribution.limit_[num_count - 1] >= 1.0) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) are             %10.4g", count,
             value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        printf("  corresponding to  %10d\n", mu);
      } else {
        printf("\n");
      }
    }
    not_reported_ones = false;
  }
  count = value_distribution.count_[num_count];
  if (count) {
    percentage = doublePercentage(count, sum_count);
    sum_percentage += percentage;
    int_percentage = percentage;
    printf("%12d %svalues (%3d%%) in [%10.4g,        inf)", count,
           value_name.c_str(), int_percentage,
           value_distribution.limit_[num_count - 1]);
    sum_report_count += count;
    if (mu > 0) {
      printf(" corresponding to [%10d,        inf)\n",
             (int)(value_distribution.limit_[num_count - 1] * mu));
    } else {
      printf("\n");
    }
  }
  if (not_reported_ones) {
    count = value_distribution.num_one_;
    if (count) {
      percentage = doublePercentage(count, sum_count);
      sum_percentage += percentage;
      int_percentage = percentage;
      printf("%12d %svalues (%3d%%) are             %10.4g", count,
             value_name.c_str(), int_percentage, 1.0);
      sum_report_count += count;
      if (mu > 0) {
        printf("  corresponding to  %10d\n", mu);
      } else {
        printf("\n");
      }
    }
  }
  printf("%12d %svalues\n", sum_count, value_name.c_str());
  if (sum_report_count != sum_count)
    printf("ERROR: %d = sum_report_count != sum_count = %d\n", sum_report_count,
           sum_count);
  return true;
}
#endif
