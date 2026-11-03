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
  double weight;
};

struct KnapsackRow {
  HighsInt row;
  HighsInt direction;
  HighsCDouble capacity;
  HighsCDouble weightsum;
  std::vector<KnapsackTerm> terms;
};

bool gcdAllButOne(std::vector<int64_t>& weights, int64_t& b,
                  bool& hasZeroWeight, const HighsInt numTerms) {
  HighsInt numOnes = 0;
  HighsInt numZeros = 0;
  HighsInt numMaxWeight = 0;
  for (const int64_t& weight : weights) {
    if (weight == 1) numOnes++;
    if (weight == 0) numZeros++;
    if (weight == b) numMaxWeight++;
  }
  if (numOnes >= 2 || numZeros + numMaxWeight + numOnes >= numTerms - 2) {
    return false;
  }
  int64_t gcd = 0;
  HighsInt candidate = -1;
  // TODO: Extend logic s.t. search can swap candidate once, see SCIP
  for (HighsInt i = 0; i != numTerms; ++i) {
    if (weights[i] == b || weights[i] == 0) continue;
    const int64_t oldGcd = gcd;
    gcd = HighsIntegers::gcd(gcd, weights[i]);
    if (gcd != 1) continue;
    if (candidate == -1) {
      candidate = i;
      gcd = oldGcd;
    } else {
      return false;
    }
  }
  assert(gcd > 1);
  const int64_t remainderCapacity = b % gcd;
  if (candidate == -1) {
    // Found a GCD for the entire knapsack
    if (remainderCapacity >= 1) {
      const int64_t oldB = b;
      b -= remainderCapacity;
      if (numMaxWeight >= 1) {
        for (HighsInt i = 0; i != numTerms; ++i) {
          if (weights[i] == oldB) weights[i] = b;
        }
      }
    }
  } else {
    // Found a GCD for all terms but one
    int64_t newWeight = 0;
    const int64_t remainderWeight = weights[candidate] % gcd;
    if (remainderWeight > remainderCapacity) {
      newWeight = weights[candidate] - remainderWeight + gcd;
    } else {
      newWeight = weights[candidate] - remainderWeight;
    }
    if (newWeight == 0 && numMaxWeight >= 1) return false;
    if (newWeight == 0) hasZeroWeight = true;
    weights[candidate] = newWeight;
    if (remainderCapacity > 0) {
      const int64_t oldB = b;
      b -= remainderCapacity;
      if (numMaxWeight >= 1) {
        for (HighsInt i = 0; i != numTerms; ++i) {
          if (weights[i] == oldB) weights[i] = b;
        }
      }
    }
  }
  b /= gcd;
  for (HighsInt i = 0; i != numTerms; ++i) {
    weights[i] /= gcd;
  }
  return true;
}

void gcdSimplify(KnapsackRow& knapsack, bool& changed, bool& hasZeroWeight) {
  // Heuristically try recursive all-but-one GCD calls
  const HighsInt numTerms = static_cast<HighsInt>(knapsack.terms.size());
  std::vector<int64_t> weights;
  weights.reserve(numTerms);
  for (const KnapsackTerm& term : knapsack.terms) {
    weights.push_back(static_cast<int64_t>(term.weight));
  }
  int64_t b = static_cast<int64_t>(static_cast<double>(knapsack.capacity));
  while (gcdAllButOne(weights, b, hasZeroWeight, numTerms)) {
    changed = true;
  }
  if (b != static_cast<int64_t>(static_cast<double>(knapsack.capacity))) {
    HighsCDouble newWeightSum = 0;
    for (HighsInt i = 0; i != numTerms; ++i) {
      knapsack.terms[i].weight = static_cast<double>(weights[i]);
      newWeightSum += knapsack.terms[i].weight;
    }
    knapsack.weightsum = newWeightSum;
    knapsack.capacity = static_cast<double>(b);
  }
}

bool recogniseCardinality(KnapsackRow& knapsack, const double feastol,
                          bool& changed) {
  HighsInt numTerms = static_cast<HighsInt>(knapsack.terms.size());
  HighsInt numZeros = 0;
  if (numTerms <= 2) return false;
  std::vector<HighsInt> order;
  order.reserve(numTerms);
  for (HighsInt i = 0; i != numTerms; ++i) {
    order.push_back(i);
    if (knapsack.terms[i].weight == 0) numZeros++;
  }
  HighsInt numNonZeroTerms = numTerms - numZeros;
  pdqsort(order.begin(), order.end(), [&](const HighsInt i, const HighsInt j) {
    return knapsack.terms[i].weight > knapsack.terms[j].weight;
  });
  auto weight = [&](const HighsInt i) -> double& {
    return knapsack.terms[order[i]].weight;
  };
  HighsInt v = 0;
  while (v != numNonZeroTerms && weight(v) == knapsack.capacity) {
    ++v;
  }
  const HighsInt numBig = v;
  if (numBig == numNonZeroTerms) return true;
  if (weight(numBig) == 1 && weight(numNonZeroTerms - 1) == 1) {
    if (numNonZeroTerms - numBig - knapsack.capacity > 0) return true;
    return false;
  }

  auto setCardinality = [&](const double newWeight) {
    changed = true;
    for (HighsInt i = 0; i != numBig; ++i) {
      knapsack.weightsum -= weight(i) - newWeight;
      weight(i) = newWeight;
    }
    for (HighsInt i = numBig; i != numNonZeroTerms; ++i) {
      knapsack.weightsum -= weight(i) - 1;
      weight(i) = 1;
    }
    knapsack.capacity = newWeight;
  };

  HighsCDouble exceedSum = 0;
  while (v != numNonZeroTerms && exceedSum - knapsack.capacity <= 0) {
    exceedSum += weight(v);
    v++;
  }
  const HighsCDouble exceedGap = exceedSum - knapsack.capacity;
  const HighsInt numExceed = v - numBig;

  if (exceedGap > 0) {
    if (exceedGap <= feastol) return false;
    // There is some subset of small items that together exceed capacity
    if (v == numNonZeroTerms) {
      // All small items needed to exceed capacity
      setCardinality(static_cast<double>(numNonZeroTerms - numBig - 1));
      return true;
    }
    HighsCDouble exceedSumBack = 0;
    for (HighsInt i = numNonZeroTerms - 1; i >= numNonZeroTerms - numExceed;
         --i) {
      exceedSumBack += weight(i);
    }
    if (exceedSumBack - knapsack.capacity > feastol) {
      // If the largest n small weights AND the smallest n small weights
      // both exceed capacity, then the row holds for all n combinations
      setCardinality(static_cast<double>(numExceed - 1));
      return true;
    }
  } else {
    // If all small terms don't exceed capacity, then the row can be rewritten
    // as \sum m x_i + \sum x_j <= m, where x_i are large items and
    // x_j are small items (where there's m small items)
    if (static_cast<HighsCDouble>(weight(numBig)) * (numNonZeroTerms - numBig) -
            knapsack.capacity <=
        feastol) {
      // If the largest small term can fit m times into the knapsack
      // then the cardinality upgrade has a stronger LP relaxation
      setCardinality(numNonZeroTerms - numBig);
      return false;
    }
  }
  return false;
}

}  // namespace

HPresolve::Result HPresolve::strengthenKnapsacks(
    HighsPostsolveStack& postsolve_stack, bool& strengthenedKnapsacks) {
  strengthenedKnapsacks = false;

  const HighsInt maxRowSize =
      std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);
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
      const double weight = knapsack.terms[i].weight;
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
          knapsack.direction * (positive ? weight : -weight) - oldRow[j].second;
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

  // Iterate over all rows. Only interested in rows of pure binaries
  // where at least one coefficient differs.
  // Will transform to make all coefficients positive, i.e., a knapsack.
  for (HighsInt row = 0; row < numRows; ++row) {
    if (rowDeleted[row] || rowsize[row] <= 2 || rowsize[row] > maxRowSize ||
        rowsizeInteger[row] + rowsizeImplInt[row] != rowsize[row] ||
        isRanged(row)) {
      continue;
    }
    const HighsInt direction = model->row_upper_[row] == kHighsInf ? -1 : 1;
    const double rowBound =
        direction == 1 ? model->row_upper_[row] : -model->row_lower_[row];
    if (std::abs(rowBound) == kHighsInf) continue;
    HighsCDouble capacity = rowBound;
    KnapsackRow knapsack{row, direction, 0, 0, {}};
    bool eligible = true;
    bool allIntCoefficients = true;
    double firstObservedCoefficient = 0;
    bool allSameCoefficients = true;
    bool changed = false;
    for (const auto& nz : getRowVector(row)) {
      const HighsInt col = nz.index();
      const double value = direction * nz.value();
      const double weight = std::abs(value);
      if (weight == kHighsInf) {
        eligible = false;
        break;
      }
      if (allIntCoefficients) {
        allIntCoefficients = value == std::floor(value);
      }
      // Columns may have been fixed by previous row changes
      if (model->col_lower_[col] == model->col_upper_[col]) {
        const double fixedValue = model->col_lower_[col];
        const double contribution = value * fixedValue;
        if (std::abs(contribution) == kHighsInf) {
          eligible = false;
          break;
        }
        capacity -= contribution;
        changed = true;
      } else {
        if (model->col_lower_[col] != 0 || model->col_upper_[col] != 1) {
          eligible = false;
          break;
        }
        assert(model->integrality_[col] != HighsVarType::kContinuous);
        if (firstObservedCoefficient == 0) firstObservedCoefficient = weight;
        if (allSameCoefficients) {
          allSameCoefficients = weight == firstObservedCoefficient;
        }
        knapsack.weightsum += weight;
        knapsack.terms.push_back({2 * col + (value > 0), weight});
        if (value < 0) capacity += weight;
      }
    }
    if (!eligible || allSameCoefficients) continue;
    knapsack.capacity = capacity;
    if (knapsack.capacity < 0) return Result::kPrimalInfeasible;
    if (knapsack.weightsum - knapsack.capacity <= primal_feastol) {
      postsolve_stack.redundantRow(row);
      removeRow(row);
      continue;
    }

    if (allIntCoefficients) {
      // Round capacity to be integer if all coefficients are integer
      HighsCDouble roundedCapacity = floor(knapsack.capacity);
      if (roundedCapacity < knapsack.capacity + primal_feastol) changed = true;
      knapsack.capacity = roundedCapacity;
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

    // Knapsack may now be redundant
    if (knapsack.weightsum - knapsack.capacity <= primal_feastol) {
      postsolve_stack.redundantRow(row);
      removeRow(row);
      continue;
    }

    if (knapsack.terms.size() <= 2) continue;

    // If a column's absence makes a row redundant then its weight
    // can be reduced
    const HighsCDouble excess = knapsack.weightsum - knapsack.capacity;
    for (auto& term : knapsack.terms) {
      if (term.weight - excess <= primal_feastol) continue;
      const HighsCDouble delta = term.weight - excess;
      term.weight = static_cast<double>(excess);
      knapsack.capacity -= delta;
      knapsack.weightsum -= delta;
      changed = true;
    }

    bool hasZeroWeight = false;
    if (allIntCoefficients && knapsack.capacity < kHighsIInf) {
      // Simplify constraint (or strengthen) using GCD arguments
      gcdSimplify(knapsack, changed, hasZeroWeight);
      // Knapsack may now be redundant
      if (knapsack.weightsum - knapsack.capacity <= primal_feastol) {
        postsolve_stack.redundantRow(row);
        removeRow(row);
        continue;
      }
    }

    // (1) If the two smallest items can't fit into the knapsack then
    // the row is a clique.
    // (2) If an item cannot fit in a knapsack along with the smallest item,
    // then its weight can be raised to knapsack capacity.
    // (3) If an item cannot fit in a knapsack with the smallest value
    // larger than minWeight (either some single item larger than minWeight or
    // the two smallest items together), then its weight can be increased to
    // capacity - minWeight.
    double minWeight = kHighsInf;
    double secondWeight = kHighsInf;
    double secondUniqueMinWeight = kHighsInf;
    for (const auto& term : knapsack.terms) {
      if (term.weight == 0) continue;
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
    HighsCDouble minWeightsCapacityDiff =
        static_cast<HighsCDouble>(minWeight) + secondWeight - knapsack.capacity;
    if (minWeightsCapacityDiff > primal_feastol) {
      if (knapsack.capacity != 1) {
        double newWeightSum = 0;
        for (auto& term : knapsack.terms) {
          if (term.weight != 0) {
            term.weight = 1;
            newWeightSum += 1;
          }
        }
        knapsack.capacity = 1;
        knapsack.weightsum = newWeightSum;
        changed = true;
      }
    } else {
      secondUniqueMinWeight =
          std::min(secondUniqueMinWeight, minWeight + secondWeight);
      for (auto& term : knapsack.terms) {
        const HighsCDouble weightPlusMinWeightExcess =
            static_cast<HighsCDouble>(term.weight) + minWeight -
            knapsack.capacity;
        const HighsCDouble itemPlusSecondUniqueExcess =
            static_cast<HighsCDouble>(secondUniqueMinWeight) + term.weight -
            knapsack.capacity;
        if (term.weight < knapsack.capacity &&
            weightPlusMinWeightExcess > primal_feastol) {
          knapsack.weightsum += knapsack.capacity - term.weight;
          term.weight = static_cast<double>(knapsack.capacity);
          changed = true;
        } else if (weightPlusMinWeightExcess < -primal_feastol &&
                   term.weight > knapsack.capacity / 2 &&
                   itemPlusSecondUniqueExcess > primal_feastol) {
          knapsack.weightsum -= weightPlusMinWeightExcess;
          term.weight = static_cast<double>(knapsack.capacity - minWeight);
          changed = true;
        }
      }
    }

    // Recognise whether the knapsack can be upgraded to a "cardinality" row
    bool skipDisaggregation =
        recogniseCardinality(knapsack, primal_feastol, changed);

    if (changed) {
      changeRow(knapsack);
      ++numNormalized;
      if (hasZeroWeight) {
        knapsack.terms.erase(
            std::remove_if(
                knapsack.terms.begin(), knapsack.terms.end(),
                [](const KnapsackTerm& term) { return term.weight == 0; }),
            knapsack.terms.end());
      }
    }
    if (!skipDisaggregation) {
      knapsacks.push_back(std::move(knapsack));
    }
  }

  // Check if knapsack can be disaggregated.
  // Given \sum a_i x_i + \sum b x_j <= b, \sum a_i < b
  // The set of rows
  // x_i + \sum x_j <= 1 \forall i
  // has a stronger LP relaxation
  HighsInt maxNewNonZeros = numNonzeros() / 4;
  HighsInt numNewNonZeros = 0;
  std::vector<std::vector<HighsInt>> newIndices;
  std::vector<std::vector<double>> newValues;
  std::vector<double> newUpper;
  std::vector<HighsInt> disaggregatedRows;
  for (KnapsackRow& knapsack : knapsacks) {
    if (numNewNonZeros >= maxNewNonZeros) break;
    std::vector<KnapsackTerm> maxWeightItems;
    std::vector<KnapsackTerm> smallItems;
    HighsCDouble smallItemSum = 0;
    for (const KnapsackTerm& term : knapsack.terms) {
      if (term.weight == knapsack.capacity)
        maxWeightItems.push_back({term.literal, 1});
      else {
        smallItems.push_back({term.literal, 1});
        smallItemSum += term.weight;
      }
    }
    if (maxWeightItems.empty() || smallItems.size() <= 1 ||
        smallItems.size() * (maxWeightItems.size() + 1) >=
            10 * (maxWeightItems.size() + smallItems.size()) ||
        smallItemSum - knapsack.capacity > 0) {
      continue;
    }
    numNewNonZeros +=
        static_cast<HighsInt>(smallItems.size() * (maxWeightItems.size() + 1) -
                              maxWeightItems.size() - smallItems.size());
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

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Knapsack presolve: %d normalized, %d fixed, %d disaggregated, "
              "%d lifted\n",
              static_cast<int>(numNormalized), static_cast<int>(numFixed),
              static_cast<int>(numDisaggregated), static_cast<int>(numLifted));
  return Result::kOk;
}

}  // namespace presolve