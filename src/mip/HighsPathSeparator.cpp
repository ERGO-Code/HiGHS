/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsPathSeparator.cpp
 * @author Leona Gottwald
 */

#include "mip/HighsPathSeparator.h"

#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

enum class RowType : int8_t {
  Unusuable = -2,
  Geq = -1,
  Eq = 0,
  Leq = 1,
};

void HighsPathSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& lpAggregator,
                                            HighsTransformedLp& transLp,
                                            HighsCutPool& cutpool) {
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  const HighsLp& lp = lpRelaxation.getLp();
  const HighsSolution& lpSolution = lpRelaxation.getSolution();

  randgen.initialise(mip.options_mip_->highs_random_seed +
                     lpRelaxation.getNumLpIterations());
  std::vector<RowType> rowtype;
  rowtype.resize(lp.numRow_);
  for (HighsInt i = 0; i != lp.numRow_; ++i) {
    if (lp.rowLower_[i] == lp.rowUpper_[i]) {
      rowtype[i] = RowType::Eq;
      continue;
    }

    double lowerslack = HIGHS_CONST_INF;
    double upperslack = HIGHS_CONST_INF;

    if (lp.rowLower_[i] != -HIGHS_CONST_INF)
      lowerslack = lpSolution.row_value[i] - lp.rowLower_[i];

    if (lp.rowUpper_[i] != HIGHS_CONST_INF)
      upperslack = lp.rowUpper_[i] - lpSolution.row_value[i];

    if (lowerslack > mip.mipdata_->feastol &&
        upperslack > mip.mipdata_->feastol)
      rowtype[i] = RowType::Unusuable;
    else if (lowerslack < upperslack)
      rowtype[i] = RowType::Geq;
    else
      rowtype[i] = RowType::Leq;
  }

  std::vector<HighsInt> numContinuous(lp.numRow_);

  size_t maxAggrRowSize = 0;
  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;

    maxAggrRowSize += lp.Astart_[col + 1] - lp.Astart_[col];
    for (HighsInt i = lp.Astart_[col]; i != lp.Astart_[col + 1]; ++i)
      ++numContinuous[lp.Aindex_[i]];
  }

  std::vector<std::pair<HighsInt, double>> colSubstitutions(
      lp.numCol_, std::make_pair(-1, 0.0));

  // identify equality rows where only a single continuous variable with nonzero
  // transformed solution value is present. Mark those columns and remember the
  // rows so that we can always substitute such columns away using this equation
  // and block the equation from being used as a start row
  for (HighsInt i = 0; i != lp.numRow_; ++i) {
    if (rowtype[i] == RowType::Eq && numContinuous[i] == 1) {
      HighsInt len;
      const HighsInt* rowinds;
      const double* rowvals;

      lpRelaxation.getRow(i, len, rowinds, rowvals);

      HighsInt j;
      for (j = 0; j != len; ++j) {
        if (mip.variableType(rowinds[j]) != HighsVarType::CONTINUOUS) continue;
        if (transLp.boundDistance(rowinds[j]) == 0.0) continue;

        break;
      }

      HighsInt col = rowinds[j];
      double val = rowvals[j];

      assert(mip.variableType(rowinds[j]) == HighsVarType::CONTINUOUS);
      assert(transLp.boundDistance(col) > 0.0);

      if (colSubstitutions[col].first != -1) continue;

      colSubstitutions[col].first = i;
      colSubstitutions[col].second = val;
      rowtype[i] = RowType::Unusuable;
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
  std::vector<std::pair<HighsInt, int>> colInArcs(lp.numCol_);

  std::vector<std::pair<HighsInt, double>> outArcRows;
  outArcRows.reserve(maxAggrRowSize);
  std::vector<std::pair<HighsInt, int>> colOutArcs(lp.numCol_);

  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;
    if (colSubstitutions[col].first != -1) continue;

    colInArcs[col].first = inArcRows.size();
    colOutArcs[col].first = outArcRows.size();
    for (HighsInt i = lp.Astart_[col]; i != lp.Astart_[col + 1]; ++i) {
      switch (rowtype[lp.Aindex_[i]]) {
        case RowType::Unusuable:
          continue;
        case RowType::Leq:
          if (lp.Avalue_[i] < 0)
            inArcRows.emplace_back(lp.Aindex_[i], lp.Avalue_[i]);
          else
            outArcRows.emplace_back(lp.Aindex_[i], lp.Avalue_[i]);
          break;
        case RowType::Geq:
        case RowType::Eq:
          if (lp.Avalue_[i] > 0)
            inArcRows.emplace_back(lp.Aindex_[i], lp.Avalue_[i]);
          else
            outArcRows.emplace_back(lp.Aindex_[i], lp.Avalue_[i]);
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

  for (HighsInt i = 0; i != lp.numRow_; ++i) {
    switch (rowtype[i]) {
      case RowType::Unusuable:
        continue;
      case RowType::Leq:
        lpAggregator.addRow(i, -1.0);
        break;
      default:
        lpAggregator.addRow(i, 1.0);
    }

    HighsInt currPathLen = 1;

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
        if (transLp.boundDistance(col) == 0.0 ||
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

        lpAggregator.addRow(row, weight);
      } else {
        HighsInt outArcRow = randgen.integer(colOutArcs[bestInArcCol].first,
                                             colOutArcs[bestInArcCol].second);

        HighsInt row = outArcRows[outArcRow].first;
        double weight = -inArcColVal / outArcRows[outArcRow].second;

        lpAggregator.addRow(row, weight);
      }
    }

    lpAggregator.clear();
  }
}
