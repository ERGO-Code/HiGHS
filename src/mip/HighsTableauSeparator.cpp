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
/**@file mip/HighsTableauSeparator.cpp
 */

#include "mip/HighsTableauSeparator.h"

#include <algorithm>

#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "pdqsort/pdqsort.h"

void HighsTableauSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                               HighsLpAggregator& lpAggregator,
                                               HighsTransformedLp& transLp,
                                               HighsCutPool& cutpool) {
  std::vector<HighsInt> basisinds;
  Highs& lpSolver = lpRelaxation.getLpSolver();
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  HighsInt numrow = lpRelaxation.numRows();
  basisinds.resize(numrow);
  lpRelaxation.getLpSolver().getBasicVariables(basisinds.data());

  std::vector<HighsInt> nonzeroWeights;
  std::vector<double> rowWeights;
  nonzeroWeights.resize(numrow);
  rowWeights.resize(numrow);
  HighsInt numNonzeroWeights;

  HighsCutGeneration cutGen(lpRelaxation, cutpool);

  std::vector<HighsInt> baseRowInds;
  std::vector<double> baseRowVals;

  const HighsSolution& lpSolution = lpRelaxation.getSolution();
#if 0
  if (mip.mipdata_->objintscale != 0.0 &&
      std::abs(
          std::round(mip.mipdata_->objintscale * lpRelaxation.getObjective()) /
              mip.mipdata_->objintscale -
          lpRelaxation.getObjective()) > 1000 * mip.mipdata_->feastol) {
    HighsInt numRows = 0;
    for (int j = 0; j != numrow; ++j) {
      double weight = mip.mipdata_->objintscale * lpSolution.row_dual[j];
      if (std::abs(weight) <= 10 * mip.mipdata_->epsilon ||
          std::abs(weight) * lpRelaxation.getMaxAbsRowVal(j) <=
              mip.mipdata_->feastol) {
        continue;
      }

      lpAggregator.addRow(j, weight);
    }

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);

    if (baseRowInds.size() - numRows <= 1000 + 0.1 * mip.numCol()) {
      double rhs = 0;
      cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

      lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
      rhs = 0;
      cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);
    }

    lpAggregator.clear();
  }
#endif
  std::vector<std::pair<double, HighsInt>> fractionalBasisvars;
  fractionalBasisvars.reserve(basisinds.size());
  for (HighsInt i = 0; i != HighsInt(basisinds.size()); ++i) {
    if (cutpool.getNumCuts() > mip.options_mip_->mip_pool_soft_limit) break;
    double fractionality;
    if (basisinds[i] < 0) {
      HighsInt row = -basisinds[i] - 1;

      if (!lpRelaxation.isRowIntegral(row)) continue;

      double solval = lpSolution.row_value[row];
      fractionality = std::abs(std::round(solval) - solval);
      fractionality -= lpRelaxation.getRowLen(row) * mip.mipdata_->feastol;
    } else {
      HighsInt col = basisinds[i];
      if (mip.variableType(col) == HighsVarType::kContinuous) continue;

      double solval = lpSolution.col_value[col];
      fractionality = std::abs(std::round(solval) - solval);
    }

    if (fractionality < 1000 * mip.mipdata_->feastol) continue;

    fractionalBasisvars.emplace_back(fractionality, i);
  }

  pdqsort(fractionalBasisvars.begin(), fractionalBasisvars.end(),
          [&](const std::pair<double, HighsInt>& a,
              const std::pair<double, HighsInt>& b) {
            return std::make_tuple(a.first,
                                   HighsHashHelpers::hash(
                                       std::make_pair(a.second, getNumCalls())),
                                   a.second) >
                   std::make_tuple(b.first,
                                   HighsHashHelpers::hash(
                                       std::make_pair(b.second, getNumCalls())),
                                   b.second);
          });

  fractionalBasisvars.resize(std::min(fractionalBasisvars.size(), size_t{200}));
  HighsInt numCuts = cutpool.getNumCuts();
  for (const auto& fracvar : fractionalBasisvars) {
    HighsInt i = fracvar.second;
    if (lpSolver.getBasisInverseRow(i, rowWeights.data(), &numNonzeroWeights,
                                    nonzeroWeights.data()) != HighsStatus::kOk)
      continue;

    // already handled by other separator
    if (numNonzeroWeights == 1) continue;

    double maxAbsRowWeight = 0.0;
    for (int j = 0; j != numNonzeroWeights; ++j) {
      int row = nonzeroWeights[j];
      maxAbsRowWeight = std::max(std::abs(rowWeights[row]), maxAbsRowWeight);
    }

    int expshift = 0;
    std::frexp(maxAbsRowWeight, &expshift);
    expshift = -expshift;

    HighsInt numRows = 0;
    for (int j = 0; j != numNonzeroWeights; ++j) {
      HighsInt row = nonzeroWeights[j];
      rowWeights[row] = std::ldexp(rowWeights[row], expshift);
      if (std::abs(rowWeights[row]) <= 10 * mip.mipdata_->epsilon ||
          std::abs(rowWeights[row]) * lpRelaxation.getMaxAbsRowVal(row) <=
              mip.mipdata_->feastol) {
        rowWeights[row] = 0;
      } else
        ++numRows;
    }

    for (int j = 0; j != numNonzeroWeights; ++j) {
      int row = nonzeroWeights[j];
      if (rowWeights[row] == 0) continue;
      lpAggregator.addRow(row, rowWeights[row]);
    }

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);

    if (baseRowInds.size() - numRows > 1000 + 0.1 * mip.numCol()) continue;

    double rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);
    if (mip.mipdata_->domain.infeasible()) break;

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
    rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);
    if (mip.mipdata_->domain.infeasible()) break;

    lpAggregator.clear();
  }
}
