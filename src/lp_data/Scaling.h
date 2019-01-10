/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Scaling.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_SCALING_H_
#define LP_DATA_SCALING_H_

#include "HConfig.h"
#include "HConst.h"
#include "HighsModelObject.h"
#include "HighsLpUtils.h"
#include <cassert>
#include <iostream>

// Limits on scaling factors
  const double minAlwScale = 1 / 1024.0;
  const double maxAlwScale = 1024.0;
  const double maxAlwCostScale = maxAlwScale;
  const double minAlwColScale = minAlwScale;
  const double maxAlwColScale = maxAlwScale;
  const double minAlwRowScale = minAlwScale;
  const double maxAlwRowScale = maxAlwScale;

#ifdef HiGHSDEV
  // Information on large costs
  const double tlLargeCo = 1e5;
  int numLargeCo;
  vector<int> largeCostFlag;
  double largeCostScale;
#endif

void scaleHighsModelInit(HighsModelObject &highs_model) {
  highs_model.scale_.col_.assign(highs_model.lp_scaled_.numCol_, 1);
  highs_model.scale_.row_.assign(highs_model.lp_scaled_.numRow_, 1);
  highs_model.scale_.cost_ = 1;
#ifdef HiGHSDEV
  //  largeCostScale = 1;
#endif
}

void scaleCosts(HighsModelObject &highs_model) {
  // Scale the costs by no less than minAlwCostScale
  double costScale = highs_model.scale_.cost_;
  double maxNzCost = 0;
  for (int iCol = 0; iCol < highs_model.lp_scaled_.numCol_; iCol++) {
    if (highs_model.lp_scaled_.colCost_[iCol]) {
      maxNzCost = max(fabs(highs_model.lp_scaled_.colCost_[iCol]), maxNzCost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  costScale = 1;
  const double ln2 = log(2.0);
  // Scale the costs if the max cost is positive and outside the range [1/16,
  // 16]
  if ((maxNzCost > 0) && ((maxNzCost < (1.0 / 16)) || (maxNzCost > 16))) {
    costScale = maxNzCost;
    costScale = pow(2.0, floor(log(costScale) / ln2 + 0.5));
    costScale = min(costScale, maxAlwCostScale);
  }
#ifdef HiGHSDEV
  printf(
      "MaxNzCost = %11.4g: scaling all costs by %11.4g\ngrep_CostScale,%g,%g\n",
      maxNzCost, costScale, maxNzCost, costScale);
#endif
  if (costScale == 1) return;
  // Scale the costs (and record of maxNzCost) by costScale, being at most
  // maxAlwCostScale
  for (int iCol = 0; iCol < highs_model.lp_scaled_.numCol_; iCol++) {
    highs_model.lp_scaled_.colCost_[iCol] /= costScale;
  }
  maxNzCost /= costScale;

#ifdef HiGHSDEV
  bool alwLargeCostScaling = false;
  /*
  if (alwLargeCostScaling && (numLargeCo > 0)) {
    // Scale any large costs by largeCostScale, being at most (a further)
    // maxAlwCostScale
    largeCostScale = maxNzCost;
    largeCostScale = pow(2.0, floor(log(largeCostScale) / ln2 + 0.5));
    largeCostScale = min(largeCostScale, maxAlwCostScale);
    printf(
        "   Scaling all |cost| > %11.4g by %11.4g\ngrep_LargeCostScale,%g,%g\n",
        tlLargeCo, largeCostScale, tlLargeCo, largeCostScale);
    for (int iCol = 0; iCol < highs_model.lp_scaled_.numCol_; iCol++) {
      if (largeCostFlag[iCol]) {
        highs_model.lp_scaled_.colCost_[iCol] /= largeCostScale;
      }
    }
  }
  */
  printf("After cost scaling\n");
  //  utils.util_analyseVectorValues("Column costs", highs_model.lp_scaled_.numCol_, highs_model.lp_scaled_.colCost_, false);
#endif
}

void scaleLp(HighsModelObject &highs_model) {
  // Scale the LP highs_model.lp_scaled_, assuming all data are in place
  // Reset all scaling to 1
  HighsSimplexInfo &simplex_info = highs_model.simplex_info_;
  HighsTimer &timer = highs_model.timer_;
  timer.start(timer.scaleClock);
  scaleHighsModelInit(highs_model);
  double *colScale = &highs_model.scale_.col_[0];
  double *rowScale = &highs_model.scale_.row_[0];
  int numCol = highs_model.lp_scaled_.numCol_;
  int numRow = highs_model.lp_scaled_.numRow_;

  // Allow a switch to/from the original scaling rules
  bool originalScaling = true;
  bool alwCostScaling = true;
  if (originalScaling) alwCostScaling = false;


  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double inf = HIGHS_CONST_INF;
  double min0 = inf, max0 = 0;
  for (int k = 0, AnX = highs_model.lp_scaled_.Astart_[numCol]; k < AnX; k++) {
    double value = fabs(highs_model.lp_scaled_.Avalue_[k]);
    min0 = min(min0, value);
    max0 = max(max0, value);
  }
  bool noScaling = min0 >= 0.2 && max0 <= 5;
  //   noScaling = false; printf("!!!! FORCE SCALING !!!!\n");
  if (noScaling) {
    // No matrix scaling, but possible cost scaling
#ifdef HiGHSDEV
    printf("grep_Scaling,%s,Obj,0,Row,1,1,Col,1,1,0\n", highs_model.modelName.c_str());
#endif
    // Possibly scale the costs
    if (!originalScaling && alwCostScaling) scaleCosts(highs_model);
    timer.stop(timer.scaleClock);
    return;
  }
  // See if we want to include cost include if minimum nonzero cost is less than
  // 0.1
  double minNzCost = inf;
  for (int i = 0; i < numCol; i++) {
    if (highs_model.lp_scaled_.colCost_[i]) minNzCost = min(fabs(highs_model.lp_scaled_.colCost_[i]), minNzCost);
  }
  bool includeCost = false;
  //  if (originalScaling)
  includeCost = minNzCost < 0.1;

  // Search up to 6 times
  vector<double> rowMin(numRow, inf);
  vector<double> rowMax(numRow, 1 / inf);
  for (int search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < numCol; iCol++) {
      // For column scale (find)
      double colMin = inf;
      double colMax = 1 / inf;
      double myCost = fabs(highs_model.lp_scaled_.colCost_[iCol]);
      if (includeCost && myCost != 0)
        colMin = min(colMin, myCost), colMax = max(colMax, myCost);
      for (int k = highs_model.lp_scaled_.Astart_[iCol]; k < highs_model.lp_scaled_.Astart_[iCol + 1]; k++) {
        double value = fabs(highs_model.lp_scaled_.Avalue_[k]) * rowScale[highs_model.lp_scaled_.Aindex_[k]];
        colMin = min(colMin, value), colMax = max(colMax, value);
      }
      colScale[iCol] = 1 / sqrt(colMin * colMax);
      if (!originalScaling) {
        // Ensure that column scale factor is not excessively large or small
        colScale[iCol] =
            min(max(minAlwColScale, colScale[iCol]), maxAlwColScale);
      }
      // For row scale (only collect)
      for (int k = highs_model.lp_scaled_.Astart_[iCol]; k < highs_model.lp_scaled_.Astart_[iCol + 1]; k++) {
        int iRow = highs_model.lp_scaled_.Aindex_[k];
        double value = fabs(highs_model.lp_scaled_.Avalue_[k]) * colScale[iCol];
        rowMin[iRow] = min(rowMin[iRow], value);
        rowMax[iRow] = max(rowMax[iRow], value);
      }
    }

    // For row scale (find)
    for (int iRow = 0; iRow < numRow; iRow++) {
      rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
      if (!originalScaling) {
        // Ensure that row scale factor is not excessively large or small
        rowScale[iRow] =
            min(max(minAlwRowScale, rowScale[iRow]), maxAlwRowScale);
      }
    }
    rowMin.assign(numRow, inf);
    rowMax.assign(numRow, 1 / inf);
  }

  // Make it numerical better
  // Also determine the max and min row and column scaling factors
  double minColScale = inf;
  double maxColScale = 1 / inf;
  double minRowScale = inf;
  double maxRowScale = 1 / inf;
  const double ln2 = log(2.0);
  for (int iCol = 0; iCol < numCol; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
    minColScale = min(colScale[iCol], minColScale);
    maxColScale = max(colScale[iCol], maxColScale);
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));
    minRowScale = min(rowScale[iRow], minRowScale);
    maxRowScale = max(rowScale[iRow], maxRowScale);
  }
#ifdef HiGHSDEV
  bool excessScaling =
      (minColScale < minAlwColScale) || (maxColScale > maxAlwColScale) ||
      (minRowScale < minAlwRowScale) || (maxRowScale > maxAlwRowScale);

  printf("grep_Scaling,%s,%d,%d,Obj,%g,%d,Row,%g,%g,Col,%g,%g,%d\n",
         highs_model.modelName.c_str(), originalScaling, alwCostScaling, minNzCost,
         includeCost, minColScale, maxColScale, minRowScale, maxRowScale,
         excessScaling);
#endif

  // Apply scaling to matrix and bounds
  for (int iCol = 0; iCol < numCol; iCol++)
    for (int k = highs_model.lp_scaled_.Astart_[iCol]; k < highs_model.lp_scaled_.Astart_[iCol + 1]; k++)
      highs_model.lp_scaled_.Avalue_[k] *= (colScale[iCol] * rowScale[highs_model.lp_scaled_.Aindex_[k]]);

  for (int iCol = 0; iCol < numCol; iCol++) {
    highs_model.lp_scaled_.colLower_[iCol] /= highs_model.lp_scaled_.colLower_[iCol] == -inf ? 1 : colScale[iCol];
    highs_model.lp_scaled_.colUpper_[iCol] /= highs_model.lp_scaled_.colUpper_[iCol] == +inf ? 1 : colScale[iCol];
    highs_model.lp_scaled_.colCost_[iCol] *= colScale[iCol];
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    highs_model.lp_scaled_.rowLower_[iRow] *= highs_model.lp_scaled_.rowLower_[iRow] == -inf ? 1 : rowScale[iRow];
    highs_model.lp_scaled_.rowUpper_[iRow] *= highs_model.lp_scaled_.rowUpper_[iRow] == +inf ? 1 : rowScale[iRow];
  }
  // Deduce the consequences of scaling the LP
  //  mlFg_Update(mlFg_action_ScaleLP);
#ifdef HiGHSDEV
  // Analyse the scaled LP
  //  if (simplex_info.analyse_lp) {
  //    util_analyseLp(highs_model.lp_scaled_, "Scaled");
  //  }
  //  if (mlFg_scaledLP) {
  //  utils.util_analyseVectorValues("Column scaling factors", numCol, colScale, false);
  //  utils.util_analyseVectorValues("Row scaling factors", numRow, rowScale, false);
  //  }
#endif
  // Possibly scale the costs
  if (!originalScaling && alwCostScaling) scaleCosts(highs_model);
  timer.stop(timer.scaleClock);
}

#endif
