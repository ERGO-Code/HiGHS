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

struct KnapsackAggregateStats {
  HighsInt numSmall = 0;
  HighsCDouble smallSum = 0;
  double minSmall = kHighsInf;
  double secondSmall = kHighsInf;
  double maxSmall = 0;
  double secondMaxSmall = 0;
  HighsCDouble outsideThirdSum = 0;
  double thirdMax1 = 0;
  double thirdMax2 = 0;
};

double convertWeightToDoubleSafely(const HighsCDouble& weight) {
  double rounded = static_cast<double>(weight);
  const HighsCDouble error = HighsCDouble(rounded) - weight;
  if (error > 0) {
    rounded = std::nextafter(rounded, -kHighsInf);
  }
  return rounded;
}

void tightenKnapsackCapacity(KnapsackRow& knapsack,
                             const HighsCDouble& newCapacity,
                             const double feastol, bool& changed) {
  if (knapsack.capacity - newCapacity <= feastol) return;
  knapsack.capacity = newCapacity;
  changed = true;
}

void increaseKnapsackWeight(KnapsackRow& knapsack, KnapsackTerm& term,
                            const HighsCDouble& newWeight, double feastol,
                            bool& changed) {
  const double weight = convertWeightToDoubleSafely(newWeight);
  if (weight - term.weight <= feastol) return;
  knapsack.weightsum += static_cast<HighsCDouble>(weight) - term.weight;
  term.weight = weight;
  changed = true;
}

bool isLargeKnapsackTerm(const double weight, const HighsCDouble& capacity,
                         const double feastol) {
  return 2 * weight - capacity > feastol;
}

HighsCDouble boundSmallResidual(const KnapsackAggregateStats& stats,
                                const HighsCDouble& residual,
                                const double feastol) {
  if (stats.numSmall == 0 || residual <= feastol) return 0;
  if (stats.smallSum - residual <= 0) return stats.smallSum;

  auto boundDoubleCount = [&](const double count) -> HighsInt {
    if (count >= stats.numSmall) return stats.numSmall;
    if (count <= 0) return 0;
    return static_cast<HighsInt>(count);
  };

  // Get an upper bound on number of small items that can also be selected
  HighsInt maxNumSmallSelected = boundDoubleCount(
      std::floor(static_cast<double>(residual) / stats.minSmall));
  // If sum of small items larger than residual, then we can't select
  // them all. How many would we have to exclude (assuming worst-case)
  // so that we can potentially fit in the residual
  if (stats.smallSum - residual > 0) {
    const HighsInt minNumSmallExcluded = boundDoubleCount(std::ceil(
        static_cast<double>(stats.smallSum - residual) / stats.maxSmall));
    maxNumSmallSelected =
        std::min(maxNumSmallSelected, stats.numSmall - minNumSmallExcluded);
  }

  // Check that selecting one more item is impossible
  // This is just done due to numerics
  if (maxNumSmallSelected < stats.numSmall) {
    const double lower = std::max(
        stats.minSmall * (maxNumSmallSelected + 1),
        static_cast<double>(stats.smallSum) -
            stats.maxSmall * (stats.numSmall - maxNumSmallSelected - 1));
    if (lower - residual <= feastol) {
      maxNumSmallSelected = stats.numSmall;
    }
  }

  // Calculate an upper bound on the remaining residual that can be used
  // by the small items.
  HighsCDouble upper = min(residual, stats.smallSum);
  if (maxNumSmallSelected >= 2) {
    upper = min(upper, static_cast<HighsCDouble>(stats.maxSmall) +
                           static_cast<HighsCDouble>(stats.secondMaxSmall) *
                               (maxNumSmallSelected - 1));
  } else {
    upper = min(
        upper, static_cast<HighsCDouble>(stats.maxSmall) * maxNumSmallSelected);
  }
  upper =
      min(upper, stats.smallSum - static_cast<HighsCDouble>(stats.minSmall) *
                                      (stats.numSmall - maxNumSmallSelected));
  // Check if the two largest items fit in the residual
  if (stats.maxSmall - residual > feastol) {
    HighsCDouble smallSmallSum = stats.smallSum - stats.maxSmall;
    if (stats.numSmall >= 2 && stats.secondMaxSmall - residual > feastol) {
      smallSmallSum -= stats.secondMaxSmall;
    }
    upper = min(upper, smallSmallSum);
  }

  // If all small items together exceed residual, then at least one
  // small item can't be selected with the large item
  if (stats.smallSum - residual > feastol)
    upper = min(upper, stats.smallSum - stats.minSmall);
  // If all small items but the largest can't fit in the residual,
  // then at least two small items must be missing to make them fit.
  if (stats.numSmall >= 2 &&
      stats.smallSum - stats.maxSmall - residual > feastol)
    upper = min(upper, stats.smallSum - stats.minSmall - stats.secondSmall);
  return max(0, upper);
}

bool computeKnapsackAggregateStats(KnapsackRow& knapsack, double feastol,
                                   const HighsCDouble& capacity,
                                   KnapsackAggregateStats& stats) {
  for (const auto& term : knapsack.terms) {
    const double weight = term.weight;
    if (weight == 0) continue;
    // Any weight near capacity means tightening the capacity isn't possible
    if (weight - capacity >= feastol) return false;
    if (!isLargeKnapsackTerm(weight, capacity, feastol)) {
      ++stats.numSmall;
      stats.smallSum += weight;
      if (weight < stats.minSmall) {
        stats.secondSmall = stats.minSmall;
        stats.minSmall = weight;
      } else {
        stats.secondSmall = std::min(stats.secondSmall, weight);
      }
      if (weight > stats.maxSmall) {
        stats.secondMaxSmall = stats.maxSmall;
        stats.maxSmall = weight;
      } else {
        stats.secondMaxSmall = std::max(stats.secondMaxSmall, weight);
      }
    }
    if (3 * static_cast<HighsCDouble>(weight) - capacity > feastol) {
      // These are the two largest entries
      if (weight > stats.thirdMax1) {
        stats.thirdMax2 = stats.thirdMax1;
        stats.thirdMax1 = weight;
      } else {
        stats.thirdMax2 = std::max(stats.thirdMax2, weight);
      }
    } else {
      stats.outsideThirdSum += weight;
    }
  }
  return true;
}

void strengthenKnapsackResiduals(KnapsackRow& knapsack, double feastol,
                                 bool& changed) {
  // Store original capacity. Tightening it changes partition of
  // small and large items
  const HighsCDouble originalCapacity = knapsack.capacity;
  KnapsackAggregateStats stats;
  if (!computeKnapsackAggregateStats(knapsack, feastol, originalCapacity,
                                     stats)) {
    return;
  }

  // Tighten the knapsack capacity if possible
  HighsCDouble capacity = min(originalCapacity, stats.smallSum);
  std::vector<HighsCDouble> strengthenedResiduals;
  strengthenedResiduals.reserve(knapsack.terms.size());
  // At most one large item can be selected. Find worst-case of this
  // happening, where we bound small-item contribution for each large item
  for (const KnapsackTerm& term : knapsack.terms) {
    if (!isLargeKnapsackTerm(term.weight, originalCapacity, feastol)) continue;
    const HighsCDouble boundedResidual =
        boundSmallResidual(stats, originalCapacity - term.weight, feastol);
    strengthenedResiduals.push_back(boundedResidual);
    capacity = max(capacity, term.weight + boundedResidual);
  }

  // At most two items with weight greater than one third of capacity fit.
  capacity =
      min(capacity, stats.outsideThirdSum + stats.thirdMax1 + stats.thirdMax2);
  // Apply the tightened capacity
  tightenKnapsackCapacity(knapsack, capacity, feastol, changed);

  // Try and tighten large terms using arguments on residual capacity
  const bool capacityTightened = originalCapacity - knapsack.capacity > feastol;
  if (capacityTightened) {
    KnapsackAggregateStats updatedStats;
    if (!computeKnapsackAggregateStats(knapsack, feastol, knapsack.capacity,
                                       updatedStats)) {
      return;
    }
    std::swap(stats, updatedStats);
  }

  HighsInt nextResidual = 0;
  for (auto& term : knapsack.terms) {
    if (!isLargeKnapsackTerm(term.weight, knapsack.capacity, feastol)) continue;
    const HighsCDouble smallUpper =
        capacityTightened ? boundSmallResidual(
                                stats, knapsack.capacity - term.weight, feastol)
                          : strengthenedResiduals[nextResidual++];
    increaseKnapsackWeight(knapsack, term, knapsack.capacity - smallUpper,
                           feastol, changed);
  }
}

bool liftKnapsackWeights(KnapsackRow& knapsack, double feastol, bool& changed) {
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
      secondWeight = std::min(secondWeight, term.weight);
      if (term.weight > minWeight)
        secondUniqueMinWeight = std::min(secondUniqueMinWeight, term.weight);
    }
  }
  if (secondWeight == kHighsInf) return false;

  if (HighsCDouble(minWeight) + secondWeight - knapsack.capacity > feastol) {
    // The row is a clique as no two items can be selected
    double newWeightsum = 0;
    for (auto& term : knapsack.terms) {
      if (term.weight == 0) continue;
      if (term.weight != 1) changed = true;
      term.weight = 1;
      newWeightsum += 1;
    }
    if (knapsack.capacity != 1) changed = true;
    knapsack.capacity = 1;
    knapsack.weightsum = newWeightsum;
    return true;
  }

  const HighsCDouble nextSmallSum = min(HighsCDouble(secondUniqueMinWeight),
                                        HighsCDouble(minWeight) + secondWeight);
  for (auto& term : knapsack.terms) {
    if (term.weight == 0) continue;
    const HighsCDouble excess =
        HighsCDouble(term.weight) + minWeight - knapsack.capacity;
    if (excess > feastol) {
      // Large item can't co-exist with the smallest item.
      // Its weight can there be raised to capacity
      increaseKnapsackWeight(knapsack, term, knapsack.capacity, feastol,
                             changed);
    } else if (excess < -feastol &&
               isLargeKnapsackTerm(term.weight, knapsack.capacity, feastol) &&
               term.weight + nextSmallSum - knapsack.capacity > feastol) {
      // Large item can't co-exist with the second unique smallest weight,
      // so can only exist with the smallest item and weight can be increased
      increaseKnapsackWeight(knapsack, term, knapsack.capacity - minWeight,
                             feastol, changed);
    }
  }
  return false;
}

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
  HighsInt firstSeen = -1;
  HighsInt numSeen = 0;
  for (HighsInt i = 0; i != numTerms; ++i) {
    if (weights[i] == b || weights[i] == 0) continue;
    if (++numSeen == 1) firstSeen = i;
    const int64_t oldGcd = gcd;
    gcd = HighsIntegers::gcd(gcd, weights[i]);
    if (gcd != 1) continue;
    if (candidate == -1) {
      candidate = i;
      gcd = oldGcd;
    } else if (numSeen == 3 && candidate != firstSeen) {
      // Excluding the second term failed, so try and exclude the first term
      gcd = HighsIntegers::gcd(weights[candidate], weights[i]);
      if (gcd == 1) return false;
      candidate = firstSeen;
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

void recogniseUnitAndMaxWeightKnapsack(KnapsackRow& knapsack,
                                       const double feastol, bool& changed) {
  // Determine if the knapsack's LP can be tightened and retain its solution
  // set by being represented by a knapsack with only
  // items of weight {1, b}
  HighsInt numTerms = static_cast<HighsInt>(knapsack.terms.size());
  HighsInt numZeros = 0;
  if (numTerms <= 2) return;
  std::vector<HighsInt> order;
  order.reserve(numTerms);
  for (HighsInt i = 0; i != numTerms; ++i) {
    order.push_back(i);
    if (knapsack.terms[i].weight == 0) numZeros++;
  }
  const HighsInt numNonZeroTerms = numTerms - numZeros;
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
  if (numBig == numNonZeroTerms) return;
  if (weight(numBig) == 1 && weight(numNonZeroTerms - 1) == 1) return;

  auto upgradeRow = [&](const double newWeight) {
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
    if (exceedGap <= feastol) {
      return;
    }
    // There is some subset of small items that together exceed capacity
    if (v == numNonZeroTerms) {
      // All small items needed to exceed capacity
      upgradeRow(static_cast<double>(numNonZeroTerms - numBig - 1));
      return;
    }
    HighsCDouble exceedSumBack = 0;
    for (HighsInt i = numNonZeroTerms - 1; i >= numNonZeroTerms - numExceed;
         --i) {
      exceedSumBack += weight(i);
    }
    if (exceedSumBack - knapsack.capacity > feastol) {
      // If the largest n small weights AND the smallest n small weights
      // both exceed capacity, then the row holds for all n combinations
      upgradeRow(static_cast<double>(numExceed - 1));
      return;
    }
  } else {
    // If all small terms don't exceed capacity, then the row can be rewritten
    // as \sum m x_i + \sum x_j <= m, where x_i are large items and
    // x_j are small items (where there's m small items)
    if (static_cast<HighsCDouble>(weight(numBig)) * (numNonZeroTerms - numBig) -
            knapsack.capacity <=
        feastol) {
      // If the largest small term can fit m times into the knapsack
      // then the row upgrade has a strictly stronger LP relaxation
      upgradeRow(numNonZeroTerms - numBig);
    }
  }
}

}  // namespace

HPresolve::Result HPresolve::strengthenKnapsacks(
    HighsPostsolveStack& postsolve_stack) {
  const HighsInt maxRowSize =
      std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);
  const HighsInt numRows = model->num_row_;
  HighsInt numNormalized = 0;
  HighsInt numFixed = 0;

  std::vector<std::pair<HighsInt, double>> oldRow;

  auto changeRow = [&](const KnapsackRow& knapsack) {
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
      if (knapsack.capacity - roundedCapacity > primal_feastol) changed = true;
      knapsack.capacity = roundedCapacity;
    }

    // Fix any columns with coefficients that are too large
    for (auto& term : knapsack.terms) {
      if (term.weight - knapsack.capacity <= primal_feastol) continue;
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

    // Try and tighten RHS and large item weights
    strengthenKnapsackResiduals(knapsack, primal_feastol, changed);
    // Check if knapsack is a clique or some weights can be increased
    const bool cliqueUpgraded =
        liftKnapsackWeights(knapsack, primal_feastol, changed);
    // Recognise whether the knapsack can be upgraded
    if (!cliqueUpgraded) {
      recogniseUnitAndMaxWeightKnapsack(knapsack, primal_feastol, changed);
    }

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
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Knapsack presolve: %d rows changed, %d fixed\n",
              static_cast<int>(numNormalized), static_cast<int>(numFixed));
  return Result::kOk;
}

}  // namespace presolve
