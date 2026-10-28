/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "HighsIntegers.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

namespace {

struct KnapsackTerm {
  HighsInt literal;
  int64_t weight;
};

struct KnapsackRow {
  HighsInt row;
  HighsInt direction;
  int64_t capacity;
  int64_t weightsum;
  std::vector<KnapsackTerm> terms;
};

}  // namespace

HPresolve::Result HPresolve::strengthenKnapsacks(
    HighsPostsolveStack& postsolve_stack, bool& strengthenedKnapsacks) {
  strengthenedKnapsacks = false;

  constexpr HighsInt maxRowSize = 1000;
  const int64_t maxValue = std::min<int64_t>(kHighsIInf, int64_t{1} << 50);
  const HighsInt numRows = model->num_row_;
  const HighsInt numCols = model->num_col_;
  HighsInt numNormalized = 0;
  HighsInt numFixed = 0;
  HighsInt numDisaggregated = 0;
  HighsInt numLifted = 0;

  std::vector<KnapsackRow> knapsacks;
  std::vector<std::pair<HighsInt, double>> oldRow;

  auto changeRow = [&](const KnapsackRow& knapsack) {
    strengthenedKnapsacks = true;
    oldRow.clear();
    for (const auto& nz : getRowVector(knapsack.row)) {
      oldRow.emplace_back(nz.index(), nz.value());
    }
    HighsCDouble rhs = knapsack.capacity;
    size_t i = 0;
    size_t j = 0;
    while (i != knapsack.terms.size()) {
      const HighsInt& col = knapsack.terms[i].literal / 2;
      const bool positive = knapsack.terms[i].literal % 2 == 1;
      const double weight = static_cast<double>(knapsack.terms[i].weight);
      // Some term may have been fixed by clique table propagation
      // and still have non-zero knapsack weight
      if (model->col_lower_[col] == model->col_upper_[col])
        rhs -= (positive ? weight : -weight) * model->col_lower_[col];
      if (!positive) rhs -= weight;
      while (oldRow[j].first != col) {
        assert(j != oldRow.size());
        ++j;
      }
      const double delta =
          knapsack.direction *
              static_cast<double>(positive ? weight : -weight) -
          oldRow[j].second;
      if (delta != 0) {
        addToMatrix(knapsack.row, col, delta);
      }
      ++i;
    }
    for (const auto& nz : getRowVector(knapsack.row)) {
      if (model->col_lower_[nz.index()] == model->col_upper_[nz.index()]) {
        rhs += knapsack.direction * nz.value() * model->col_lower_[nz.index()];
      }
    }
    if (knapsack.direction == 1)
      changeRowUpper(knapsack.row, static_cast<double>(rhs));
    else
      changeRowLower(knapsack.row, -static_cast<double>(rhs));
  };

  // Iterate over all rows. Only interested in knapsack rows, i.e.
  // row of pure binary columns with integer coefficients where at least one
  // coefficient > 1
  for (HighsInt row = 0; row < numRows; ++row) {
    if (rowDeleted[row] || rowsize[row] <= 2 || rowsize[row] > maxRowSize ||
        isRanged(row))
      continue;
    const HighsInt direction = model->row_upper_[row] == kHighsInf ? -1 : 1;
    const double rowBound =
        direction == 1 ? model->row_upper_[row] : -model->row_lower_[row];
    if (std::abs(rowBound) > maxValue) continue;
    int64_t overflowBudget =
        maxValue - static_cast<int64_t>(std::ceil(std::abs(rowBound)));
    HighsCDouble capacity = rowBound;
    KnapsackRow knapsack{row, direction, 0, 0, {}};
    bool eligible = true;
    bool hasGeneralIntCoefficients = false;
    bool changed = false;
    for (const auto& nz : getRowVector(row)) {
      const HighsInt col = nz.index();
      const double value = direction * nz.value();
      if (value != std::floor(value) || std::abs(value) > maxValue) {
        eligible = false;
        break;
      }
      // Columns may have been fixed by previous row changes
      if (model->col_lower_[col] == model->col_upper_[col]) {
        const double fixedValue = model->col_lower_[col];
        const double contribution = value * fixedValue;
        if (std::abs(contribution) > overflowBudget ||
            fixedValue != std::floor(fixedValue) ||
            contribution != std::floor(contribution)) {
          eligible = false;
          break;
        }
        overflowBudget -= static_cast<int64_t>(std::abs(contribution));
        capacity -= contribution;
        changed = true;
      } else {
        if (model->integrality_[col] == HighsVarType::kContinuous ||
            model->col_lower_[col] != 0 || model->col_upper_[col] != 1) {
          eligible = false;
          break;
        }
        const int64_t weight = static_cast<int64_t>(std::abs(value));
        if (weight > overflowBudget) {
          eligible = false;
          break;
        }
        overflowBudget -= weight;
        knapsack.weightsum += weight;
        knapsack.terms.push_back({2 * col + (value > 0), weight});
        if (value < 0) capacity += weight;
        if (weight != 1) hasGeneralIntCoefficients = true;
      }
    }
    if (!eligible || !hasGeneralIntCoefficients) continue;
    knapsack.capacity = static_cast<int64_t>(
        static_cast<double>(floor(capacity + primal_feastol)));
    if (knapsack.capacity < 0) return Result::kPrimalInfeasible;
    if (knapsack.weightsum <= knapsack.capacity) {
      postsolve_stack.redundantRow(row);
      removeRow(row);
      continue;
    }

    // Filter out setppc rows
    if (knapsack.capacity <= 1) continue;

    // Reduce coefficient values by common factor if possible
    int64_t gcd = 0;
    for (const KnapsackTerm& term : knapsack.terms) {
      gcd = HighsIntegers::gcd(gcd, term.weight);
      if (gcd == 1) break;
    }
    if (gcd > 1) {
      knapsack.weightsum /= gcd;
      knapsack.capacity /= gcd;
      if (knapsack.capacity <= 1) continue;
      for (KnapsackTerm& term : knapsack.terms) term.weight /= gcd;
      changed = true;
    }

    // Fix any columns with coefficients that are too large
    for (auto& term : knapsack.terms) {
      if (term.weight <= knapsack.capacity) continue;
      const HighsInt col = term.literal / 2;
      if (term.literal % 2 != 0) {
        HPRESOLVE_CHECKED_CALL(changeColUpper(col, 0.0));
      } else {
        HPRESOLVE_CHECKED_CALL(changeColLower(col, 1.0));
      }
      knapsack.weightsum -= term.weight;
      term.weight = 0;
      ++numFixed;
      changed = true;
    }
    knapsack.terms.erase(
        std::remove_if(
            knapsack.terms.begin(), knapsack.terms.end(),
            [](const KnapsackTerm& term) { return term.weight == 0; }),
        knapsack.terms.end());
    if (knapsack.weightsum <= knapsack.capacity) {
      postsolve_stack.redundantRow(row);
      removeRow(row);
      continue;
    }

    if (knapsack.terms.size() <= 2) continue;

    // If a column's absence makes a row redundant then its weight
    // can be reduced
    const int64_t excess = knapsack.weightsum - knapsack.capacity;
    for (auto& term : knapsack.terms) {
      if (term.weight <= excess) continue;
      const int64_t delta = term.weight - excess;
      term.weight = excess;
      knapsack.capacity -= delta;
      knapsack.weightsum -= delta;
      changed = true;
    }

    // (1) If the two smallest items can't fit into the knapsack then
    // the row is a clique.
    // (2) If an item cannot fit in a knapsack along with the smallest item,
    // then its weight can be raised to knapsack capacity.
    // (3) If an item cannot fit in a knapsack with the smallest value
    // larger than minWeight (either some single item larger than minWeight or
    // the two smallest items together), then its weight can be increased to
    // capacity - minWeight.
    int64_t minWeight = maxValue;
    int64_t secondWeight = maxValue;
    int64_t secondUniqueMinWeight = maxValue;
    for (const auto& term : knapsack.terms) {
      if (term.weight < minWeight) {
        secondWeight = minWeight;
        secondUniqueMinWeight = minWeight;
        minWeight = term.weight;
      } else {
        if (term.weight < secondWeight) {
          secondWeight = term.weight;
        }
        if (term.weight > minWeight)
          secondUniqueMinWeight = std::min(term.weight, secondUniqueMinWeight);
      }
    }
    if (minWeight + secondWeight > knapsack.capacity) {
      if (knapsack.capacity != 1) {
        for (auto& term : knapsack.terms) term.weight = 1;
        knapsack.capacity = 1;
        knapsack.weightsum = knapsack.terms.size();
        changed = true;
      }
    } else {
      secondUniqueMinWeight =
          std::min(secondUniqueMinWeight, minWeight + secondWeight);
      for (auto& term : knapsack.terms) {
        if (term.weight < knapsack.capacity &&
            term.weight + minWeight > knapsack.capacity) {
          knapsack.weightsum += knapsack.capacity - term.weight;
          term.weight = knapsack.capacity;
          changed = true;
        } else if (term.weight < knapsack.capacity - minWeight &&
                   term.weight > knapsack.capacity / 2 &&
                   knapsack.capacity - term.weight < secondUniqueMinWeight) {
          knapsack.weightsum += knapsack.capacity - minWeight - term.weight;
          term.weight = knapsack.capacity - minWeight;
          changed = true;
        }
      }
    }
    if (changed) {
      changeRow(knapsack);
      ++numNormalized;
    }
    knapsacks.push_back(std::move(knapsack));
  }

  // Check if knapsack can be disaggregated.
  // Given \sum a_i x_i + \sum b x_j <= b, a_i <= b,
  // The set of rows
  // x_i + \sum x_j <= 1 \forall i
  // has a stronger LP relaxation
  std::vector<std::vector<HighsInt>> newIndices;
  std::vector<std::vector<double>> newValues;
  std::vector<double> newUpper;
  std::vector<HighsInt> disaggregatedRows;
  for (KnapsackRow& knapsack : knapsacks) {
    std::vector<KnapsackTerm> maxWeightItems;
    std::vector<KnapsackTerm> smallItems;
    int64_t smallItemSum = 0;
    for (const KnapsackTerm& term : knapsack.terms) {
      if (term.weight == knapsack.capacity)
        maxWeightItems.push_back({term.literal, 1});
      else {
        smallItems.push_back({term.literal, 1});
        smallItemSum += term.weight;
      }
    }
    if (maxWeightItems.empty() || smallItems.empty() ||
        smallItemSum > knapsack.capacity ||
        maxWeightItems.size() * smallItems.size() >
            10 * knapsack.terms.size()) {
      continue;
    }
    knapsack.capacity = 1;  // Spoof capacity so we don't try lifting afterwards
    for (const KnapsackTerm& term : smallItems) {
      std::vector<KnapsackTerm> clique = maxWeightItems;
      clique.push_back(term);
      std::vector<HighsInt> indices;
      std::vector<double> values;
      double rhs = 1;
      for (const auto& literal : clique) {
        indices.push_back(literal.literal / 2);
        const bool positive = literal.literal % 2 == 1;
        values.push_back(positive ? 1.0 : -1.0);
        if (!positive) rhs -= 1;
      }
      newIndices.push_back(std::move(indices));
      newValues.push_back(std::move(values));
      newUpper.push_back(rhs);
    }
    // TODO: Can we add a clique to the table here or does that cause issues
    // with fixings when no domain is present?
    disaggregatedRows.push_back(knapsack.row);
  }

  if (!addToMatrix(postsolve_stack,
                   std::vector<double>(newUpper.size(), -kHighsInf), newUpper,
                   newIndices, newValues))
    return Result::kOutOfMemory;
  for (HighsInt row : disaggregatedRows) {
    postsolve_stack.redundantRow(row);
    removeRow(row);
    ++numDisaggregated;
    strengthenedKnapsacks = true;
  }

  knapsacks.erase(std::remove_if(knapsacks.begin(), knapsacks.end(),
                                 [](const KnapsackRow& knapsack) {
                                   return knapsack.capacity <= 1;
                                 }),
                  knapsacks.end());

  // Search for a literal u, s.t.,
  // u = 1 -> x_i = 0 for some set of I
  // If \alpha = \sum_{i \in I} a_i - (weightsum - capacity) > 0
  // then we can lift \alpha u into the row

  // The sum of a_i for x_i = 0 when literal = 1, i.e.,
  // the weight from the knapsack that disappears by fixing the literal.
  std::vector<int64_t> removedWeight(2 * numCols, 0);
  // Keeps track of source from where the literal was last visited
  std::vector<HighsInt> lastLiteralSource(2 * numCols, -1);
  // Keep track of literals with non-zero removed weight
  std::vector<HighsInt> candidates;
  for (const KnapsackRow& knapsack : knapsacks) {
    candidates.clear();
    for (const KnapsackTerm& term : knapsack.terms) {
      HighsInt negatedTermLiteral = term.literal ^ 1;
      lastLiteralSource[negatedTermLiteral] = term.literal;
      if (removedWeight[negatedTermLiteral] == 0)
        candidates.push_back(negatedTermLiteral);
      removedWeight[negatedTermLiteral] += term.weight;
      presolveCliqueTable.recordLiteralNeighbourhood(
          term.literal, term.weight, removedWeight, lastLiteralSource,
          candidates, maxRowSize);
    }
    HighsInt bestLiteral = -1;
    int64_t bestWeight = 0;
    for (HighsInt literal : candidates) {
      const HighsInt col = literal / 2;
      const int64_t weight =
          removedWeight[literal] - (knapsack.weightsum - knapsack.capacity);
      if (!colDeleted[col] &&
          model->col_lower_[col] != model->col_upper_[col] &&
          weight > bestWeight) {
        bestLiteral = literal;
        bestWeight = weight;
      }
      removedWeight[literal] = 0;
      lastLiteralSource[literal] = -1;
    }
    if (bestLiteral == -1) continue;
    const double weight = static_cast<double>(bestWeight);
    const bool positive = bestLiteral % 2 == 1;
    addToMatrix(knapsack.row, bestLiteral / 2,
                knapsack.direction * (positive ? weight : -weight));
    if (!positive) {
      if (knapsack.direction == 1) {
        addToRowUpper(knapsack.row, -static_cast<HighsCDouble>(weight));
      } else {
        addToRowLower(knapsack.row, static_cast<HighsCDouble>(weight));
      }
    }
    ++numLifted;
    strengthenedKnapsacks = true;
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Knapsack presolve: %d normalized, %d fixed, %d disaggregated, "
              "%d lifted\n",
              static_cast<int>(numNormalized), static_cast<int>(numFixed),
              static_cast<int>(numDisaggregated), static_cast<int>(numLifted));
  return Result::kOk;
}

}  // namespace presolve