/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "util/HighsUtils.h"
#include "HConfig.h"
#include "lp_data/HConst.h"

#include <stdio.h>
#include <cmath>
#include <vector>

bool highs_isInfinity(double val) {
  if (val >= HIGHS_CONST_INF) return true;
  return false;
}

#ifdef HiGHSDEV
void analyseVectorValues(const char* message, int vecDim,
                         const std::vector<double>& vec,
                         bool analyseValueList) {
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

  lastRpCat;
  for (int cat = 0; cat < maxCat + 1; cat++) {
    if (rowCatK[cat]) lastRpCat = cat;
  }
  cat = maxCat;
  if (rowCatK[cat]) lastRpCat = cat;
  sumK = 0;
  pct;
  v;
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

#endif
