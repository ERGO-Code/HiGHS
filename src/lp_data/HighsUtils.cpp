/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "HighsLp.h"
#include "HighsIO.h"
#include "HighsUtils.h"
#include "HConst.h"

#include <cmath>

bool highs_isInfinity(double val) {
  if (val >= HIGHS_CONST_INF) return true;
  return false;
}

#ifdef HiGHSDEV
void HighsUtils::util_anVecV(const char *message, int vecDim, std::vector<double> &vec,
                         bool anVLs) {
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
    double absV = abs(v);
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
    if (anVLs) {
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
  if (nNegInfV > 0) printf("   %7d values are -Inf\n", nNegInfV);
  if (nPosInfV > 0) printf("   %7d values are +Inf\n", nPosInfV);
  int k = nVK;
  int vK = posVK[k];
  if (vK > 0) printf("   %7d values satisfy 10^(%3d) <= v < Inf\n", vK, k);
  for (int k = nVK - 1; k >= 0; k--) {
    int vK = posVK[k];
    if (vK > 0)
      printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, k, k + 1);
  }
  for (int k = 1; k <= nVK; k++) {
    int vK = negVK[k];
    if (vK > 0)
      printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, -k, 1 - k);
  }
  vK = vecDim - nNz;
  if (vK > 0) printf("   %7d values are zero\n", vK);
  if (anVLs) {
    printf("           Value distribution:");
    if (excessVLsV) printf(" More than %d different values", VLsZ);
    printf("\n           Value    Count\n");
    for (int ix = 0; ix < VLsZ; ix++) {
      int pct = ((100.0 * VLsK[ix]) / vecDim) + 0.5;
      printf("     %11.4g %8d (%3d%%)\n", VLsV[ix], VLsK[ix], pct);
    }
  }
}
#endif

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void HighsUtils::reportLp(HighsLp lp) {
  reportLpBrief(lp);
  reportLpColVec(lp);
  reportLpRowVec(lp);
  reportLpColMtx(lp);
}

// Report the LP briefly
void HighsUtils::reportLpBrief(HighsLp lp) {
  reportLpDimensions(lp);
  reportLpObjSense(lp);
}

// Report the LP dimensions
void HighsUtils::reportLpDimensions(HighsLp lp) {
  HighsPrintMessage(HighsMessageType::INFO,
                    "LP has %d columns, %d rows and %d nonzeros\n",
                    lp.numCol_, lp.numRow_, lp.Astart_[lp.numCol_]);
}

// Report the LP objective sense
void HighsUtils::reportLpObjSense(HighsLp lp) {
  if (lp.sense_ == OBJSENSE_MINIMIZE)
    HighsPrintMessage(HighsMessageType::INFO, "Objective sense is minimize\n");
  else if (lp.sense_ == OBJSENSE_MAXIMIZE)
    HighsPrintMessage(HighsMessageType::INFO, "Objective sense is maximize\n");
  else
    HighsPrintMessage(HighsMessageType::INFO,
                      "Objective sense is ill-defined as %d\n", lp.sense_);
}

// Report the vectors of LP column data
void HighsUtils::reportLpColVec(HighsLp lp) {
  if (lp.numCol_ <= 0) return;
  HighsPrintMessage(HighsMessageType::INFO,
                    "  Column        Lower        Upper         Cost\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(HighsMessageType::INFO, "%8d %12g %12g %12g\n", iCol,
                      lp.colLower_[iCol], lp.colUpper_[iCol], lp.colCost_[iCol]);
  }
}

// Report the vectors of LP row data
void HighsUtils::reportLpRowVec(HighsLp lp) {
  if (lp.numRow_ <= 0) return;
  HighsPrintMessage(HighsMessageType::INFO,
                    "     Row        Lower        Upper\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsPrintMessage(HighsMessageType::INFO, "%8d %12g %12g\n", iRow,
                      lp.rowLower_[iRow], lp.rowUpper_[iRow]);
  }
}

// Report the LP column-wise matrix
void HighsUtils::reportLpColMtx(HighsLp lp) {
  if (lp.numCol_ <= 0) return;
  HighsPrintMessage(HighsMessageType::INFO,
                    "Column Index              Value\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(HighsMessageType::INFO, "    %8d Start   %10d\n", iCol,
                      lp.Astart_[iCol]);
    for (int el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
      HighsPrintMessage(HighsMessageType::INFO, "          %8d %12g\n",
                        lp.Aindex_[el], lp.Avalue_[el]);
    }
  }
  HighsPrintMessage(HighsMessageType::INFO, "             Start   %10d\n",
                    lp.Astart_[lp.numCol_]);
}

/*
void HighsUtils::reportLpSolution(HighsModelObject highs_model) {
  HighsLp lp = highs_model.lp_scaled_;
  reportLpBrief(lp);
  //  model->util_reportModelStatus(lp);
  assert(lp.numCol_ > 0);
  assert(lp.numRow_ > 0);
  vector<double> colPrimal(lp.numCol_);
  vector<double> colDual(lp.numCol_);
  vector<int> colStatus(lp.numCol_);
  vector<double> rowPrimal(lp.numRow_);
  vector<double> rowDual(lp.numRow_);
  vector<int> rowStatus(lp.numRow_);
  //  util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  //  if (util_convertWorkingToBaseStat(&colStatus[0], &rowStatus[0])) return;
  //  util_reportColVecSol(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_, colPrimal, colDual, colStatus);
  //  util_reportRowVecSol(lp.numRow_, lp.rowLower_, lp.rowUpper_, rowPrimal, rowDual, rowStatus);
}
*/

#ifdef HiGHSDEV
void HighsUtils::util_anMl(HighsLp lp, const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_anVecV("Column costs", lp.numCol_, lp.colCost_, false);
  util_anVecV("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_anVecV("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_anVecV("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_anVecV("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_anVecV("Matrix entries", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  util_anMlBd("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_anMlBd("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
}

void HighsUtils::util_anMlBd(const char *message, int numBd, std::vector<double> &lower,
                         std::vector<double> &upper) {
  if (numBd == 0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix = 0; ix < numBd; ix++) {
    if (highs_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Finite lower bound and infinite upper bound: Lb
        numLb++;
      } else {
        // Finite lower bound and   finite upper bound:
        if (lower[ix] < upper[ix]) {
          // Distinct finite bounds: Bx
          numBx++;
        } else {
          // Equal finite bounds: Fx
          numFx++;
        }
      }
    }
  }
  printf("Analysing %d %s bounds\n", numBd, message);
  if (numFr > 0)
    printf("   Free:  %7d (%3d%%)\n", numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    printf("   LB:    %7d (%3d%%)\n", numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    printf("   UB:    %7d (%3d%%)\n", numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    printf("   Boxed: %7d (%3d%%)\n", numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    printf("   Fixed: %7d (%3d%%)\n", numFx, (100 * numFx) / numBd);
  printf("grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  printf("grep_CharMl,%d,%d,%d,%d,%d,%d\n", numBd, numFr, numLb, numUb, numBx,
         numFx);
}

#endif
