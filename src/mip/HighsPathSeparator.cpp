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
/**@file mip/HighsPathSeparator.cpp
 */

#include "mip/HighsPathSeparator.h"

#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

enum class RowType : int8_t {
  kUnusuable = -2,
  kGeq = -1,
  kEq = 0,
  kLeq = 1,
};

void HighsPathSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& lpAggregator,
                                            HighsTransformedLp& transLp,
                                            HighsCutPool& cutpool) {
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  const HighsLp& lp = lpRelaxation.getLp();
  const HighsSolution& lpSolution = lpRelaxation.getSolution();

  randgen.initialise(mip.options_mip_->random_seed +
                     lpRelaxation.getNumLpIterations());
  std::vector<RowType> rowtype;
  rowtype.resize(lp.num_row_);
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    if (lp.row_lower_[i] == lp.row_upper_[i]) {
      rowtype[i] = RowType::kEq;
      continue;
    }

    double lowerslack = kHighsInf;
    double upperslack = kHighsInf;

    if (lp.row_lower_[i] != -kHighsInf)
      lowerslack = lpSolution.row_value[i] - lp.row_lower_[i];

    if (lp.row_upper_[i] != kHighsInf)
      upperslack = lp.row_upper_[i] - lpSolution.row_value[i];

    if (lowerslack > mip.mipdata_->feastol &&
        upperslack > mip.mipdata_->feastol)
      rowtype[i] = RowType::kUnusuable;
    else if (lowerslack < upperslack)
      rowtype[i] = RowType::kGeq;
    else
      rowtype[i] = RowType::kLeq;
  }

  std::vector<HighsInt> numContinuous(lp.num_row_);

  size_t maxAggrRowSize = 0;
  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;

    maxAggrRowSize += lp.a_start_[col + 1] - lp.a_start_[col];
    for (HighsInt i = lp.a_start_[col]; i != lp.a_start_[col + 1]; ++i)
      ++numContinuous[lp.a_index_[i]];
  }

  std::vector<std::pair<HighsInt, double>> colSubstitutions(
      lp.num_col_, std::make_pair(-1, 0.0));

  // identify equality rows where only a single continuous variable with nonzero
  // transformed solution value is present. Mark those columns and remember the
  // rows so that we can always substitute such columns away using this equation
  // and block the equation from being used as a start row
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    if (rowtype[i] == RowType::kEq && numContinuous[i] == 1) {
      HighsInt len;
      const HighsInt* rowinds;
      const double* rowvals;

      lpRelaxation.getRow(i, len, rowinds, rowvals);

      HighsInt j;
      for (j = 0; j != len; ++j) {
        if (mip.variableType(rowinds[j]) != HighsVarType::kContinuous) continue;
        if (transLp.boundDistance(rowinds[j]) == 0.0) continue;

        break;
      }

      HighsInt col = rowinds[j];
      double val = rowvals[j];

      assert(mip.variableType(rowinds[j]) == HighsVarType::kContinuous);
      assert(transLp.boundDistance(col) > 0.0);

      if (colSubstitutions[col].first != -1) continue;

      colSubstitutions[col].first = i;
      colSubstitutions[col].second = val;
      rowtype[i] = RowType::kUnusuable;
    }
  }

  // for each continuous variable with nonzero transformed solution value
  // remember the <= and == rows where it is present with a positive coefficient
  // in its set of in-arc rows. Treat >= rows as <= rows with reversed sign
  // The reason to only store one set of rows for one sign of the coefficients
  // is that this directs the selection to be more diverse. Consider
  // aggregations of 2 rows where we allow both directions. When one of the rows
  // is used as start row we can always select the other one. When we only
  // project out variables with negative coefficients we give the aggregation
  // path an orientation and avoid this situation. For each aggregation of rows
  // the cut generation will try the reversed orientation in any case too.

  std::vector<std::pair<HighsInt, double>> inArcRows;
  inArcRows.reserve(maxAggrRowSize);
  std::vector<std::pair<HighsInt, int>> colInArcs(lp.num_col_);

  std::vector<std::pair<HighsInt, double>> outArcRows;
  outArcRows.reserve(maxAggrRowSize);
  std::vector<std::pair<HighsInt, int>> colOutArcs(lp.num_col_);

  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;
    if (colSubstitutions[col].first != -1) continue;

    colInArcs[col].first = inArcRows.size();
    colOutArcs[col].first = outArcRows.size();
    for (HighsInt i = lp.a_start_[col]; i != lp.a_start_[col + 1]; ++i) {
      switch (rowtype[lp.a_index_[i]]) {
        case RowType::kUnusuable:
          continue;
        case RowType::kLeq:
          if (lp.a_value_[i] < 0)
            inArcRows.emplace_back(lp.a_index_[i], lp.a_value_[i]);
          else
            outArcRows.emplace_back(lp.a_index_[i], lp.a_value_[i]);
          break;
        case RowType::kGeq:
        case RowType::kEq:
          if (lp.a_value_[i] > 0)
            inArcRows.emplace_back(lp.a_index_[i], lp.a_value_[i]);
          else
            outArcRows.emplace_back(lp.a_index_[i], lp.a_value_[i]);
          break;
      }
    }

    colInArcs[col].second = inArcRows.size();
    colOutArcs[col].second = outArcRows.size();
  }

  HighsCutGeneration cutGen(lpRelaxation, cutpool);
  std::vector<HighsInt> baseRowInds;
  std::vector<double> baseRowVals;
  const HighsInt maxPathLen = 6;
  std::vector<std::pair<std::vector<HighsInt>, std::vector<double>>>
      aggregatedPath;
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    switch (rowtype[i]) {
      case RowType::kUnusuable:
        continue;
      case RowType::kLeq:
        lpAggregator.addRow(i, -1.0);
        break;
      default:
        lpAggregator.addRow(i, 1.0);
    }

    HighsInt currPathLen = 1;
    const double maxWeight = 1. / mip.mipdata_->feastol;
    const double minWeight = mip.mipdata_->feastol;

    auto checkWeight = [&](double w) {
      w = std::abs(w);
      return w <= maxWeight && w >= minWeight;
    };

    aggregatedPath.clear();

    while (currPathLen != maxPathLen) {
      lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);
      HighsInt baseRowLen = baseRowInds.size();
      bool addedSubstitutionRows = false;

      HighsInt bestOutArcCol = -1;
      double outArcColVal = 0.0;
      double outArcColBoundDist = 0.0;

      HighsInt bestInArcCol = -1;
      double inArcColVal = 0.0;
      double inArcColBoundDist = 0.0;

      for (HighsInt j = 0; j != baseRowLen; ++j) {
        HighsInt col = baseRowInds[j];
        if (col >= lp.num_col_ || transLp.boundDistance(col) == 0.0 ||
            lpRelaxation.isColIntegral(col))
          continue;

        if (colSubstitutions[col].first != -1) {
          addedSubstitutionRows = true;
          lpAggregator.addRow(colSubstitutions[col].first,
                              -baseRowVals[j] / colSubstitutions[col].second);
          continue;
        }

        if (addedSubstitutionRows) continue;

        if (baseRowVals[j] < 0) {
          if (colInArcs[col].first == colInArcs[col].second) continue;
          if (bestOutArcCol == -1 ||
              transLp.boundDistance(col) > outArcColBoundDist) {
            bestOutArcCol = col;
            outArcColVal = baseRowVals[j];
            outArcColBoundDist = transLp.boundDistance(col);
          }
        } else {
          if (colOutArcs[col].first == colOutArcs[col].second) continue;
          if (bestInArcCol == -1 ||
              transLp.boundDistance(col) > inArcColBoundDist) {
            bestInArcCol = col;
            inArcColVal = baseRowVals[j];
            inArcColBoundDist = transLp.boundDistance(col);
          }
        }
      }

      if (addedSubstitutionRows) continue;

      double rhs = 0;

      bool success = cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

      lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
      if (!aggregatedPath.empty() || bestOutArcCol != -1 || bestInArcCol != -1)
        aggregatedPath.emplace_back(baseRowInds, baseRowVals);
      rhs = 0;

      success =
          success || cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

      if (success || (bestOutArcCol == -1 && bestInArcCol == -1)) break;

      ++currPathLen;
      // we prefer to use an out edge if the bound distances are equal in
      // feasibility tolerance otherwise we choose an inArc. This tie breaking
      // is arbitrary, but we should direct the substitution to prefer one
      // direction to increase diversity.
      if (bestInArcCol == -1 ||
          (bestOutArcCol != -1 &&
           outArcColBoundDist >= inArcColBoundDist - mip.mipdata_->feastol)) {
        HighsInt inArcRow = randgen.integer(colInArcs[bestOutArcCol].first,
                                            colInArcs[bestOutArcCol].second);

        HighsInt row = inArcRows[inArcRow].first;
        double weight = -outArcColVal / inArcRows[inArcRow].second;

        if (!checkWeight(weight)) {
          bool success = false;
          for (HighsInt nextRow = inArcRow + 1;
               nextRow < colInArcs[bestOutArcCol].second && !success;
               ++nextRow) {
            row = inArcRows[nextRow].first;
            weight = -outArcColVal / inArcRows[nextRow].second;
            success = checkWeight(weight);
          }

          for (HighsInt nextRow = colInArcs[bestOutArcCol].first;
               nextRow < inArcRow && !success; ++nextRow) {
            row = inArcRows[nextRow].first;
            weight = -outArcColVal / inArcRows[nextRow].second;
            success = checkWeight(weight);
          }

          if (!success) {
            if (bestInArcCol == -1)
              break;
            else
              goto check_out_arc_col;
          }
        }

        lpAggregator.addRow(row, weight);
      } else {
      check_out_arc_col:
        HighsInt outArcRow = randgen.integer(colOutArcs[bestInArcCol].first,
                                             colOutArcs[bestInArcCol].second);

        HighsInt row = outArcRows[outArcRow].first;
        double weight = -inArcColVal / outArcRows[outArcRow].second;

        if (!checkWeight(weight)) {
          bool success = false;
          for (HighsInt nextRow = outArcRow + 1;
               nextRow < colOutArcs[bestInArcCol].second && !success;
               ++nextRow) {
            row = outArcRows[nextRow].first;
            weight = -inArcColVal / outArcRows[nextRow].second;
            success = checkWeight(weight);
          }

          for (HighsInt nextRow = colOutArcs[bestInArcCol].first;
               nextRow < outArcRow && !success; ++nextRow) {
            row = outArcRows[nextRow].first;
            weight = -inArcColVal / outArcRows[nextRow].second;
            success = checkWeight(weight);
          }

          if (!success) break;
        }

        lpAggregator.addRow(row, weight);
      }
    }

    // if the path has length at least 2 try to separate a path mixing cut
    HighsInt pathLen = aggregatedPath.size();
    if (pathLen > 1) {
      // generate path mixing cut
      HighsHashTable<HighsInt, HighsInt> indexPos;

      std::vector<HighsInt> inds;
      std::vector<double> solval;
      std::vector<double> upper;
      std::vector<uint8_t> isIntegral;
      inds.reserve(lp.num_col_ + lp.num_row_);
      solval.reserve(lp.num_col_ + lp.num_row_);
      upper.reserve(lp.num_col_ + lp.num_row_);
      isIntegral.reserve(lp.num_col_ + lp.num_row_);

      std::vector<double> rhs(pathLen);
      std::vector<double> tmpUpper;
      std::vector<double> tmpSolval;

      double delta = 1.0;

      for (HighsInt k = 0; k < pathLen; ++k) {
        bool integralPositive = false;

        if (!transLp.transform(aggregatedPath[k].second, tmpUpper, tmpSolval,
                               aggregatedPath[k].first, rhs[k],
                               integralPositive)) {
          pathLen = k;
          break;
        }

        if (rhs[k] > kHighsTiny ||
            (k > 0 && rhs[k - 1] <= rhs[k] + mip.mipdata_->feastol)) {
          pathLen = k;
          break;
        }
        rhs[k] = std::min(0., rhs[k]);
        delta = std::max(std::abs(rhs[k]), delta);

        HighsInt len = aggregatedPath[k].first.size();
        for (HighsInt j = 0; j < len; ++j) {
          HighsInt index = aggregatedPath[k].first[j];
          HighsInt* pos = &indexPos[index];
          if (*pos == 0) {
            inds.push_back(index);
            solval.push_back(tmpSolval[j]);
            upper.push_back(tmpUpper[j]);
            isIntegral.push_back(lpRelaxation.isColIntegral(index));
            if (isIntegral.back())
              delta = std::max(std::abs(aggregatedPath[k].second[j]), delta);
            *pos = inds.size();
          } else {
            assert(inds[*pos - 1] == index);
            assert(solval[*pos - 1] == tmpSolval[j]);
            assert(upper[*pos - 1] == tmpUpper[j]);
            if (isIntegral[*pos - 1])
              delta = std::max(std::abs(aggregatedPath[k].second[j]), delta);
          }
        }
      }

      if (pathLen > 1) {
        delta = std::exp2(std::ceil(std::log2(delta + 1.0)));

        HighsInt numInds = inds.size();

        std::vector<double> valueMatrix;
        valueMatrix.resize(pathLen * numInds, 0.0);

        HighsCDouble cutRhs = 0.0;
        std::vector<double> cutVals(numInds);
        std::vector<double> maxFrac(numInds);
        std::vector<HighsCDouble> downSum(numInds);
        std::vector<HighsCDouble> fSum(numInds);

        double fLast = 0;
        double scale = -1.0 / delta;

        for (HighsInt k = 0; k < pathLen; ++k) {
          double f = rhs[k] * scale;
          HighsCDouble fDiff = HighsCDouble(f) - fLast;
          HighsInt len = aggregatedPath[k].first.size();
          cutRhs += fDiff;
          for (HighsInt j = 0; j < len; ++j) {
            HighsInt i = indexPos[aggregatedPath[k].first[j]] - 1;
            assert(i >= 0);

            double gj = aggregatedPath[k].second[j] * scale;

            switch (isIntegral[i]) {
              case 0:
                cutVals[i] = std::max(cutVals[i], gj);
                break;
              case 1: {
                double gjdown = std::floor(gj);
                double hj = gj - gjdown;
                maxFrac[i] = std::max(maxFrac[i], hj);
                downSum[i] += fDiff * gjdown;
                fSum[i] += fDiff;

                if (fSum[i] < maxFrac[i]) {
                  cutVals[i] = double(downSum[i] + fSum[i]);
                } else {
                  cutVals[i] = double(downSum[i] + maxFrac[i]);
                }
                break;
              }
            }
          }

          if (k > 0) {
            double viol = double(cutRhs);
            for (HighsInt j = 0; j < numInds; ++j) {
              viol -= solval[j] * cutVals[j];
            }

            viol *= delta;

            if (viol > 10 * mip.mipdata_->feastol) {
              scale = -delta;
              double rhs = double(cutRhs) * scale;
              for (HighsInt j = 0; j < numInds; ++j) {
                cutVals[j] *= scale;
              }

              for (HighsInt j = numInds - 1; j >= 0; --j) {
                if (std::abs(cutVals[j]) <= mip.mipdata_->epsilon) {
                  --numInds;
                  std::swap(cutVals[j], cutVals[numInds]);
                  std::swap(inds[j], inds[numInds]);
                }
              }

              cutVals.resize(numInds);
              inds.resize(numInds);

              if (transLp.untransform(cutVals, inds, rhs)) {
                HighsInt cutLen = inds.size();
                mip.mipdata_->debugSolution.checkCut(
                    inds.data(), cutVals.data(), cutLen, rhs);

                // compute violation in untransformed space again
                double viol = -rhs;
                for (HighsInt j = 0; j < cutLen; ++j)
                  viol += cutVals[j] *
                          lpRelaxation.getSolution().col_value[inds[j]];

                if (viol > 10 * mip.mipdata_->feastol) {
                  mip.mipdata_->domain.tightenCoefficients(
                      inds.data(), cutVals.data(), cutLen, rhs);
                  cutpool.addCut(mip, inds.data(), cutVals.data(), inds.size(),
                                 rhs);
                }
              }
              // printf("cut is violated for k = %d\n", k);
              break;
            }
          }
          fLast = f;
        }
      }
    }

    lpAggregator.clear();
  }
}
