/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsMachineSchedSeparator.cpp
 */

#include "mip/HighsMachineSchedSeparator.h"

#include <unordered_map>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

bool HighsMachineSchedSeparator::findSingleMachineScheduleClique(
    std::vector<std::vector<double>>& vals,
    std::vector<std::vector<HighsInt>>& inds, std::vector<double>& rhss,
    double& releasedate, const HighsMipSolver& mipsolver) {
  struct pair_hash {
    size_t operator()(const std::pair<HighsInt, HighsInt>& p) const {
      return HighsHashHelpers::hash(p);
    }
  };
  std::vector<HighsCliqueTable::CliqueVar> clique(2);
  enum class ArcType {
    kIfBinOne,   // 0
    kIfBinZero,  // 1
  };
  HighsInt largestDegree = 0;
  HighsInt largestDegreeCol = -1;
  std::vector<HighsInt> degrees(mipsolver.numCol());
  std::unordered_map<std::pair<HighsInt, HighsInt>,
                     std::tuple<double, HighsInt, ArcType>, pair_hash>
      adjacency;

  auto addEntry = [&](HighsInt posCol, HighsInt negCol, HighsInt binCol,
                      double p, ArcType t) {
    // My_ji + s_i - s_j >= p_ji
    // y_ji = val -> s_i >= s_j + p_ji - M * val
    // Make an arc from negCol (j) to posCol (i)
    if (p < 0) return;
    auto it = adjacency.find({negCol, posCol});
    if (it != adjacency.end()) {
      // If there's two overlapping arcs y_ji and y'_ji then we can conclude
      // if y_ji -> y'_ji, i.e., y_ji + ~y'_ji <= 1
      clique[0] = HighsCliqueTable::CliqueVar(
          std::get<1>(it->second),
          std::get<2>(it->second) == ArcType::kIfBinOne);
      clique[1] =
          HighsCliqueTable::CliqueVar(binCol, t == ArcType::kIfBinOne ? 0 : 1);
      mipsolver.mipdata_->cliquetable.addClique(mipsolver, clique.data(), 2);
      if (std::get<0>(it->second) < p) {
        std::get<0>(it->second) = p;
        std::get<1>(it->second) = binCol;
        std::get<2>(it->second) = t;
      }
    } else {
      degrees[posCol]++;
      if (degrees[posCol] > largestDegree ||
          (degrees[posCol] == largestDegree && posCol < largestDegreeCol)) {
        largestDegreeCol = posCol;
        largestDegree = degrees[posCol];
      }
      adjacency.emplace(std::make_pair(negCol, posCol),
                        std::make_tuple(p, binCol, t));
    }
  };

  HighsInt numRows = 0;
  const HighsInt maxRows = std::min(HighsInt{5000}, 2 * mipsolver.numRow());
  adjacency.reserve(maxRows + 2);
  for (HighsInt row = 0; row != mipsolver.numRow(); row++) {
    double rowLower = mipsolver.model_->row_lower_[row];
    double rowUpper = mipsolver.model_->row_upper_[row];
    if (rowLower == rowUpper) continue;
    HighsInt start = mipsolver.mipdata_->ARstart_[row];
    HighsInt end = mipsolver.mipdata_->ARstart_[row + 1];
    if (end - start != 3) continue;
    bool machineSchedRow = true;
    HighsInt posContCol = -1;
    HighsInt negContCol = -1;
    HighsInt binCol = -1;
    double binCoef = 0;
    for (HighsInt i = start; i != end; i++) {
      HighsInt col = mipsolver.mipdata_->ARindex_[i];
      if (mipsolver.mipdata_->domain.col_lower_[col] == -kHighsInf) {
        // TODO: Could some jobs be modelled in reverse?
        machineSchedRow = false;
        break;
      }
      if (mipsolver.mipdata_->domain.isBinary(col)) {
        if (binCol != -1) {
          machineSchedRow = false;
          break;
        }
        binCol = col;
        binCoef = mipsolver.mipdata_->ARvalue_[i];
      } else if (mipsolver.mipdata_->ARvalue_[i] == -1) {
        negContCol = col;
      } else if (mipsolver.mipdata_->ARvalue_[i] == 1) {
        posContCol = col;
      } else {
        machineSchedRow = false;
        break;
      }
    }
    if (!machineSchedRow || binCol == -1 || negContCol == -1 ||
        posContCol == -1 || posContCol == negContCol)
      continue;
    // We want to put the row into form:
    // Mx_ji + s_i - s_j >= p_ji, p_ji >= 0
    // x_ji = 1 -> s_i >= s_j + p_ij - M
    if (rowUpper != kHighsInf) {
      // Given My_ij + s_i - s_j <= d
      // The row becomes (after multiplying by -1):
      // -My_ij + s_j - s_i >= -d
      // Add implication s_j >= s_i + p_ij + M, p_ij + M > 0 when binCol = 1
      // Add implication s_j >= s_i + p_ij, p, p_ij > 0 when binCol = 0
      double rhs_0 = -rowUpper;
      double rhs_1 = -rowUpper + binCoef;
      // TODO: If both rhs_0 > 0 && rhs_1 > 0 then strengthen the claim, i.e.,
      // TODO: min{rhs_0, rhs_1} + y_ji * (max{rhs_0, rhs_1} - min{rhs_0, rhs_1}
      if (rhs_0 > 0 || rhs_1 > 0) {
        if (rhs_0 > 0) {
          addEntry(negContCol, posContCol, binCol, rhs_0, ArcType::kIfBinZero);
        } else {
          addEntry(negContCol, posContCol, binCol, rhs_1, ArcType::kIfBinOne);
        }
        numRows++;
      }
    }
    if (rowLower != -kHighsInf) {
      // Given Mx_ij + s_i - s_j >= d
      // Add implication s_i >= s_j + pji - M, p_ji - M > 0 when binCol = 1
      // Add implication s_i >= s_j + pji, p_ji > 0 when binCol = 0
      double rhs_0 = rowLower;
      double rhs_1 = rowLower - binCoef;
      if (rhs_0 > 0 || rhs_1 > 0) {
        if (rhs_0 > 0) {
          addEntry(posContCol, negContCol, binCol, rhs_0, ArcType::kIfBinZero);
        } else if (rhs_1 > 0) {
          addEntry(posContCol, negContCol, binCol, rhs_1, ArcType::kIfBinOne);
        }
        numRows++;
      }
    }
    if (numRows >= maxRows) break;
  }

  // A clique of size 3 needs at least 6 arcs
  if (numRows <= 5) return false;

  // Greedily search neighbours of largest degree column for a double-sided
  // clique (corresponds to a single machine schedule)
  std::vector<HighsInt> potentialNeighbours;
  potentialNeighbours.reserve(largestDegree);
  for (const auto& arc : adjacency) {
    HighsInt col = std::get<1>(arc.first);
    if (col == largestDegreeCol) {
      potentialNeighbours.emplace_back(std::get<0>(arc.first));
    }
  }
  pdqsort(potentialNeighbours.begin(), potentialNeighbours.end(),
          [&](const HighsInt c1, const HighsInt c2) {
            return degrees[c1] > degrees[c2];
          });

  std::vector<HighsInt> neighbours;
  neighbours.reserve(largestDegree + 1);
  neighbours.emplace_back(largestDegreeCol);
  releasedate = mipsolver.mipdata_->domain.col_lower_[largestDegreeCol];
  std::vector<double> processingTimes;
  processingTimes.resize(largestDegree + 1, kHighsInf);
  // Iterate over potential neighbours and check validity
  for (HighsInt col : potentialNeighbours) {
    bool valid_neighbour = true;
    for (HighsInt neighbour : neighbours) {
      auto fromArc = adjacency.find({col, neighbour});
      auto toArc = adjacency.find({neighbour, col});
      // Need to verify that the to and fromArc exist and are opposites
      if (toArc == adjacency.end() || fromArc == adjacency.end() ||
          std::get<1>(fromArc->second) != std::get<1>(toArc->second) ||
          std::get<2>(fromArc->second) == std::get<2>(toArc->second)) {
        valid_neighbour = false;
        break;
      }
    }
    if (!valid_neighbour) continue;
    size_t newNeighbourIndex = neighbours.size();
    // Extract the processing times from the arcs
    for (size_t i = 0; i != neighbours.size(); ++i) {
      HighsInt neighbour = neighbours[i];
      const auto fromArc = adjacency.find({col, neighbour});
      const auto toArc = adjacency.find({neighbour, col});
      processingTimes[newNeighbourIndex] = std::min(
          processingTimes[newNeighbourIndex], std::get<0>(fromArc->second));
      processingTimes[i] =
          std::min(processingTimes[i], std::get<0>(toArc->second));
    }
    releasedate =
        std::min(mipsolver.mipdata_->domain.col_lower_[col], releasedate);
    neighbours.emplace_back(col);
  }
  if (neighbours.size() < 3) return false;

  // Now populate the actual inequalities
  vals.resize(neighbours.size(), std::vector<double>(neighbours.size()));
  inds.resize(neighbours.size(), std::vector<HighsInt>(neighbours.size()));
  rhss.resize(neighbours.size());
  for (size_t i = 0; i != neighbours.size(); ++i) {
    rhss[i] -= releasedate;
    HighsInt col = neighbours[i];
    for (size_t j = 0; j != neighbours.size(); ++j) {
      size_t jj = j >= i ? j + 1 : j;
      if (jj >= neighbours.size()) continue;
      HighsInt neighbour = neighbours[jj];
      const auto toArc = adjacency.find({neighbour, col});
      assert(toArc != adjacency.end());
      inds[i][j] = std::get<1>(toArc->second);
      vals[i][j] = processingTimes[jj];
      if (std::get<2>(toArc->second) == ArcType::kIfBinZero) {
        rhss[i] -= vals[i][j];
        vals[i][j] *= -1;
      }
    }
    // Put the job start time on the LHS
    inds[i].back() = col;
    vals[i].back() = -1;
  }

  return true;
}

void HighsMachineSchedSeparator::separateLpSolution(
    HighsLpRelaxation& lpRelaxation, HighsLpAggregator& lpAggregator,
    HighsTransformedLp& transLp, HighsCutPool& cutpool) {
  if (separated) return;
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  std::vector<std::vector<double>> vals;
  std::vector<std::vector<HighsInt>> inds;
  std::vector<double> rhss;
  double releasedate;
  has_single_machine_schedule =
      findSingleMachineScheduleClique(vals, inds, rhss, releasedate, mip);
  if (!has_single_machine_schedule) {
    separated = true;
    return;
  }

  // Load the cuts!!!
  const std::vector<double>& lpSolution = lpRelaxation.getSolution().col_value;
  for (size_t i = 0; i != inds.size(); ++i) {
    double viol = -rhss[i];
    for (size_t j = 0; j != inds[i].size(); j++) {
      viol += lpSolution[inds[i][j]] * vals[i][j];
    }
    if (viol >= 10 * mip.mipdata_->feastol)
      cutpool.addCut(mip, inds[i].data(), vals[i].data(), inds[i].size(),
                     rhss[i]);
  }
  separated = true;
}
