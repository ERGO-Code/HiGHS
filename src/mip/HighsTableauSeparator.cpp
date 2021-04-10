/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTableauSeparator.cpp
 * @author Leona Gottwald
 */

#include "mip/HighsTableauSeparator.h"

#include <algorithm>

#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

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
    } else {
      HighsInt col = basisinds[i];
      if (mip.variableType(col) == HighsVarType::CONTINUOUS) continue;

      double solval = lpSolution.col_value[col];
      fractionality = std::abs(std::round(solval) - solval);
    }

    if (fractionality < 1000 * mip.mipdata_->feastol) continue;

    fractionalBasisvars.emplace_back(fractionality, i);
  }

  std::sort(fractionalBasisvars.begin(), fractionalBasisvars.end(),
            [&fractionalBasisvars](const std::pair<double, HighsInt>& a,
                                   const std::pair<double, HighsInt>& b) {
              return std::make_tuple(
                         a.first,
                         HighsHashHelpers::hash((uint64_t(a.second) << 32) +
                                                fractionalBasisvars.size()),
                         a.second) >
                     std::make_tuple(
                         b.first,
                         HighsHashHelpers::hash((uint64_t(b.second) << 32) +
                                                fractionalBasisvars.size()),
                         b.second);
            });
  HighsInt numCuts = cutpool.getNumCuts();
  for (const auto& fracvar : fractionalBasisvars) {
    HighsInt i = fracvar.second;
    if (lpSolver.getBasisInverseRow(i, rowWeights.data(), &numNonzeroWeights,
                                    nonzeroWeights.data()) != HighsStatus::OK)
      continue;

    // do not use aggregations which use too many rows for reasons of
    // performance and numerical stability
    if (numNonzeroWeights > 100 + 0.15 * lpRelaxation.numRows()) continue;

    if (numNonzeroWeights == 1) {
      lpAggregator.addRow(nonzeroWeights[0], 1);
    } else {
      double maxAbsColWeight = 0.0;
      for (int j = 0; j != numNonzeroWeights; ++j) {
        int row = nonzeroWeights[j];
        double maxRowVal = lpRelaxation.getMaxAbsRowVal(row);
        maxAbsColWeight =
            std::max(std::abs(rowWeights[row]) * maxRowVal, maxAbsColWeight);
      }

      int expshift = 0;
      if (maxAbsColWeight > 1e3 || maxAbsColWeight < 0.5) {
        std::frexp(maxAbsColWeight, &expshift);
        expshift = -expshift + 1;
      }

      for (int j = 0; j != numNonzeroWeights; ++j) {
        int row = nonzeroWeights[j];
        double weight = std::ldexp(rowWeights[row], expshift);
        lpAggregator.addRow(row, weight);
      }
    }

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);

    double rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
    rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

    lpAggregator.clear();

    if (cutpool.getNumCuts() - numCuts >=
        0.1 * mip.options_mip_->mip_pool_soft_limit)
      break;
  }
}
