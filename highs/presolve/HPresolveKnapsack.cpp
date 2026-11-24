/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"
#include "util/HighsIntegers.h"

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
  HighsCDouble weightSum;
  std::vector<KnapsackTerm> terms;
};

struct KnapsackStatistics {
  HighsInt numEntries = 0;
  HighsInt numLarge = 0;
  HighsInt numSmall = 0;
  HighsInt numZeros = 0;
  HighsInt numNonZeros = 0;
  HighsInt numOnes = 0;
  double secondSmallestUniqueWeight = 0;
  HighsCDouble smallSum = 0;
  HighsCDouble outsideThirdSum = 0;
};

void clearKnapSackStatistics(KnapsackStatistics& stats) {
  stats.numEntries = 0;
  stats.numLarge = 0;
  stats.numSmall = 0;
  stats.numZeros = 0;
  stats.numNonZeros = 0;
  stats.numOnes = 0;
  stats.secondSmallestUniqueWeight = 0;
  stats.smallSum = 0;
  stats.outsideThirdSum = 0;
}

void tightenKnapsackCapacity(KnapsackRow& knapsack,
                             const HighsCDouble& newCapacity,
                             const double feastol, bool& changed) {
  if (knapsack.capacity - newCapacity <= feastol) return;
  knapsack.capacity = newCapacity;
  changed = true;
}

void increaseKnapsackWeight(KnapsackRow& knapsack, KnapsackTerm& term,
                            const double weight, double feastol,
                            bool& changed) {
  if (weight - term.weight <= feastol) return;
  knapsack.weightSum += static_cast<HighsCDouble>(weight) - term.weight;
  term.weight = weight;
  changed = true;
}

bool isLargeKnapsackTerm(const double weight, const HighsCDouble& capacity,
                         const double feastol) {
  return 2 * weight - capacity > feastol;
}

void computeKnapsackStatistics(KnapsackRow& knapsack,
                               const std::vector<HighsInt>& order,
                               const double feastol,
                               KnapsackStatistics& stats) {
  stats.numEntries = static_cast<HighsInt>(knapsack.terms.size());
  stats.numLarge = stats.numEntries;
  double minWeight = kHighsInf;
  for (HighsInt i = stats.numEntries - 1; i >= 0; --i) {
    const double weight = knapsack.terms[order[i]].weight;
    if (weight == 0) {
      stats.numZeros++;
      stats.numLarge--;
      continue;
    }
    stats.numNonZeros++;
    if (weight == 1) stats.numOnes++;
    if (static_cast<HighsCDouble>(weight) * 3 - knapsack.capacity <= feastol) {
      stats.outsideThirdSum += weight;
    }
    if (!isLargeKnapsackTerm(weight, knapsack.capacity, feastol)) {
      stats.numSmall++;
      stats.numLarge--;
      stats.smallSum += weight;
      if (stats.numSmall == 1) {
        minWeight = weight;
        stats.secondSmallestUniqueWeight = weight;
      } else if (stats.secondSmallestUniqueWeight == minWeight &&
                 weight != minWeight) {
        stats.secondSmallestUniqueWeight = weight;
      }
    }
  }
  if (stats.numSmall >= 2) {
    const double twoSmallestWeights =
        minWeight + knapsack.terms[order[stats.numNonZeros - 2]].weight;
    // If all small weights are equal, the next sum uses two items.
    stats.secondSmallestUniqueWeight =
        stats.secondSmallestUniqueWeight == minWeight
            ? twoSmallestWeights
            : std::min(stats.secondSmallestUniqueWeight, twoSmallestWeights);
  }
}

void tightenRhs(KnapsackRow& knapsack, const std::vector<HighsInt>& order,
                const double feastol, KnapsackStatistics& stats,
                bool& changed) {
  bool changedRhs = false;
  computeKnapsackStatistics(knapsack, order, feastol, stats);
  if (knapsack.terms[order[0]].weight - knapsack.capacity >= -feastol) {
    // Can't tighten capacity to lower than the largest item's weight
    return;
  }
  const HighsCDouble capacity = knapsack.capacity;
  if (isLargeKnapsackTerm(knapsack.terms[order[0]].weight, capacity, feastol)) {
    // Knapsack capacity is min(capacity, largest_item + small_sum)
    // Consider the knapsack
    // 13.5 x_1 + 12 x_2 + 11 x_3 + 3 x_4 + 2 x_5 ≤ 20
    // The first three weights are all large, i.e., can't fit with one another.
    // Assuming the worst-case that we can fit all the small items
    // with the largest large item is then a bound on the capacity.
    // In this example it is tighter than what's provided:
    // 13.5 + 3 + 2 = 18.5 < 20
    tightenKnapsackCapacity(knapsack,
                            knapsack.terms[order[0]].weight + stats.smallSum,
                            feastol, changedRhs);
  }
  if (stats.numNonZeros >= 2 &&
      static_cast<HighsCDouble>(knapsack.terms[order[1]].weight) * 3 -
              capacity >
          feastol) {
    // At most two items above one third of the original capacity can fit.
    // Consider the knapsack
    // 8.5 x_1 + 8 x_2 + 7.5 x_3 + 7 x_4 + 2 x_5 <= 20
    // The first four weights are all larger than 20/3, so no three
    // of them can fit together. Take the largest two and assume that
    // all the remaining smaller items are taken
    // 8.5 + 8 + 2 = 18.5 < 20
    // That lets us tighten the capacity. It's a generalisation of the above
    // case
    tightenKnapsackCapacity(knapsack,
                            stats.outsideThirdSum +
                                knapsack.terms[order[0]].weight +
                                knapsack.terms[order[1]].weight,
                            feastol, changedRhs);
  }
  if (changedRhs) {
    changed = true;
    clearKnapSackStatistics(stats);
    computeKnapsackStatistics(knapsack, order, feastol, stats);
  }
}

bool liftKnapsackWeights(KnapsackRow& knapsack,
                         const std::vector<HighsInt>& order, double feastol,
                         const KnapsackStatistics& stats, bool& changed) {
  if (stats.numNonZeros < 2) return false;
  const double minWeight =
      knapsack.terms[order[stats.numLarge + stats.numSmall - 1]].weight;
  const double secondMinWeight =
      knapsack.terms[order[stats.numLarge + stats.numSmall - 2]].weight;
  if (minWeight + secondMinWeight - knapsack.capacity > feastol) {
    // The row is a clique as no two items can be selected
    // Consider
    // 7 x_1 + 6.5 x_2 + 6 x_3 + 4.5 x_4 <= 10
    // The two smallest items don't fit, i.e., 6 + 4.5 = 10.5 > 10
    // The knapsack can therefore be
    // x_1 + x_2 + x_3 + x_4 <= 1
    double newWeightSum = 0;
    for (auto& term : knapsack.terms) {
      if (term.weight == 0) continue;
      if (term.weight != 1) changed = true;
      term.weight = 1;
      newWeightSum += 1;
    }
    if (knapsack.capacity != 1) changed = true;
    knapsack.capacity = 1;
    knapsack.weightSum = newWeightSum;
    return true;
  }

  for (HighsInt i = 0; i != stats.numLarge; i++) {
    KnapsackTerm& term = knapsack.terms[order[i]];
    if (term.weight + minWeight - knapsack.capacity > feastol) {
      // A large item cannot fit with the smallest item
      // Increase its weight to capacity
      // Consider the knapsack
      // 9.5 x_1 + 6 x_2 + 2 x_3 + 2 x_4 + 1 x_5 <= 10
      // 9.5 + 1 = 10.5 > 10
      // Change to
      // 10 x_1 + 6 x_2 + 2 x_3 + 2 x_4 + 1 x_5 <= 10
      increaseKnapsackWeight(knapsack, term,
                             static_cast<double>(knapsack.capacity), feastol,
                             changed);
    } else if (term.weight + stats.secondSmallestUniqueWeight -
                   knapsack.capacity >
               feastol) {
      // A large item cannot fit with the second-smallest unique value
      // Increase its weight to capacity - minWeight
      // Consider the knapsack
      // 8.5 x_1 + 6 x_2 + 3 x_3 + 1 x_4 + 1 x_5 <= 10
      // 2 is the second-smallest possible value in the knapsack, which
      // comes from selecting the two smallest items
      // 8.5 + 2 = 10.5 > 10
      // Change to
      // 9 x_1 + 6 x_2 + 3 x_3 + 1 x_4 + 1 x_5 <= 10
      increaseKnapsackWeight(knapsack, term,
                             static_cast<double>(knapsack.capacity - minWeight),
                             feastol, changed);
    } else if (stats.numSmall >= 2 &&
               term.weight + stats.smallSum -
                       knapsack.terms[order[stats.numLarge]].weight -
                       knapsack.capacity >
                   feastol) {
      // A large item cannot fit with all small items excluding the largest
      // small item. We therefore know if the item is selected, at least
      // two small items can't be selected.
      // Increase its weight to capacity - smallSum + minWeight +
      // secondMinWeight (provided it's an improvement)
      // Consider the knapsack
      // 6.5 x_1 + 5.5 x_2 + 2.5 x_3 + 2.5 x_4 + 2 x_5 <= 10
      // 6.5 + 2.5 + 2 = 11 > 10 (can't fit with small items excluding the
      // largest small item)
      // 10 - 2.5 = 7.5 > 6.5 (new proposed weight is larger)
      // Change knapsack to
      // 7.5 x_1 + 5.5 x_2 + 2.5 x_3 + 2.5 x_4 + 2 x_5 <= 10
      increaseKnapsackWeight(
          knapsack, term,
          static_cast<double>(knapsack.capacity - stats.smallSum + minWeight +
                              secondMinWeight),
          feastol, changed);
    } else if (term.weight + stats.smallSum - knapsack.capacity > feastol &&
               term.weight + stats.smallSum - minWeight - knapsack.capacity <
                   -feastol) {
      // A large item cannot fit with all small items, but can
      // fit with every small item but the smallest
      // Increase its weight to capacity - smallSum + minWeight (provided
      // it's an improvement)
      // Consider the knapsack
      // 6 x_1 + 5.5 x_2 + 2 x_3 + 2 x_4 + x_5 <= 10
      // 5.5 + 2 + 2 + 1 = 10.5 > 10 (can't fit with all small items)
      // 5.5 + 2 + 2 = 9.5 < 10 (can fit with all small items except
      // the smallest)
      // 10 - 2 - 2 = 6 > 5.5 (new weight is larger)
      // Change knapsack to
      // 6 x_1 + 6 x_2 + 2 x_3 + 2 x_4 + x_5 <= 10
      increaseKnapsackWeight(
          knapsack, term,
          static_cast<double>(knapsack.capacity - stats.smallSum + minWeight),
          feastol, changed);
    } else if (term.weight + stats.smallSum - knapsack.capacity < -feastol) {
      // A large item can fit with all small items
      // Increase its weight to capacity - smallSum (provided it's an
      // improvement)
      // Consider the knapsack
      // 8 x_1 + 7 x_2 + 6.5 x_3 + 2 x_4 + x_5 <= 10
      // 6.5 + 2 + 1 = 9.5 < 10 (large item can co-exist with all small items)
      // Change knapsack to
      // 8 x_1 + 7 x_2 + 7 x_3 + 2 x_4 + x_5 <= 10
      increaseKnapsackWeight(
          knapsack, term,
          static_cast<double>(knapsack.capacity - stats.smallSum), feastol,
          changed);
    }
  }
  return false;
}

bool gcdAllButOne(std::vector<int64_t>& weights, int64_t& b,
                  const HighsInt numTerms) {
  // Consider the knapsack
  // 16 x_1 + 11 x_2 + 8 x_3 + 4 x_4 <= 21
  // The GCD of the four terms is 1, but if we leave out 11 x_2,
  // then the remaining three terms have a GCD of 4.
  // We have that 11 % 4 = 3 and 21 % 4 = 1. We can therefore round up
  // the weight 11 to 12 and round the capacity from 21 to 20.
  // Dividing through by the GCD then gives the strictly stronger knapsack
  // 4 x_1 + 3 x_2 + 2 x_3 + x_4 <= 5
  if (numTerms <= 2) return false;
  HighsInt numOnes = 0;
  HighsInt numZeros = 0;
  HighsInt numMaxWeight = 0;
  for (const int64_t& weight : weights) {
    if (weight == 1) {
      if (++numOnes == 2) return false;
    } else if (weight == 0) {
      numZeros++;
    }
    if (weight == b) {
      numMaxWeight++;
    }
  }
  if (numZeros + numMaxWeight + numOnes >= numTerms - 2) {
    return false;
  }
  std::vector<int64_t> suffixGcd(weights.size() + 1, 1);
  suffixGcd.back() = 0;
  for (HighsInt i = numTerms - 1; i >= 0; --i) {
    suffixGcd[i] = suffixGcd[i + 1];
    if (weights[i] != b && weights[i] != 0) {
      suffixGcd[i] = HighsIntegers::gcd(suffixGcd[i], weights[i]);
    }
    if (suffixGcd[i] == 1) break;
  }

  int64_t gcd = suffixGcd[0];
  HighsInt candidate = -1;
  if (gcd == 1) {
    int64_t prefixGcd = 0;
    for (HighsInt i = 0; i != numTerms && prefixGcd != 1; ++i) {
      if (weights[i] == b || weights[i] == 0) continue;
      gcd = HighsIntegers::gcd(prefixGcd, suffixGcd[i + 1]);
      if (gcd > 1) {
        candidate = i;
        break;
      }
      prefixGcd = HighsIntegers::gcd(prefixGcd, weights[i]);
    }
    if (candidate == -1) return false;
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

void gcdRemoveTail(std::vector<int64_t>& weights, int64_t& b,
                   const HighsInt numTerms, const double feastol) {
  // This function uses the prefix and tail concept
  // Given a sorted knapsack \sum a_i x_i <= b (larger weights first)
  // The prefix are all weights before some item i, and the tail
  // is all weights after and including item i
  // Concrete example:
  // 30 x_1 + 20 x_2 + 20 x_3 + 10 x_4 + 2 x_5 + x_6 <= 43
  // Prefix: 30 x_1 + 20 x_2 + 20 x_3 + 10 x_4
  // Tail: 2 x_5 + x_6
  // We can see that the GCD for the prefix is 10. The
  // tail fits into capacity modulo GCD, and can therefore fit alongside
  // any prefix selection, and are therefore redundant.
  // Rewrite as: 3 x_1 + 2 x_2 + 2 x_3 + x_4 <= 4
  // Second example:
  // 9 x_1 + 8 x_2 + 7 x_3 + 6 x_4 + 2 x_5 + x_6 <= 20
  // Prefix: 9 x_1 + 8 x_2 + 7 x_3 + 6 x_4
  // Tail: 2 x_5 + x_6
  // Every pair of items in the prefix admit the entire tail,
  // but no three items from the prefix can be selected
  // Rewrite as x_1 + x_2 + x_3 + x_4 <= 2
  std::vector<HighsCDouble> tailSum(numTerms + 1, 0);
  for (HighsInt i = numTerms; i > 0; --i) {
    tailSum[i - 1] = tailSum[i] + static_cast<double>(weights[i - 1]);
  }
  int64_t gcd = 0;
  HighsInt maxPrefixItemsWithTail = 0;
  HighsCDouble largeSum = 0;
  for (HighsInt prefixSize = 1; prefixSize < numTerms; ++prefixSize) {
    if (gcd != 1) {
      gcd = HighsIntegers::gcd(gcd, weights[prefixSize - 1]);
      if (gcd > 1 && tailSum[prefixSize] - static_cast<double>(b % gcd) <= 0) {
        b /= gcd;
        for (HighsInt i = 0; i < numTerms; ++i) {
          if (weights[i] == 0) continue;
          weights[i] = i < prefixSize ? weights[i] / gcd : 0;
        }
        return;
      }
    }

    while (maxPrefixItemsWithTail < prefixSize) {
      // Track how many of the largest items can fit in the prefix
      // assuming a full tail
      const HighsCDouble next =
          largeSum + static_cast<double>(weights[maxPrefixItemsWithTail]);
      if (next + tailSum[prefixSize] - static_cast<double>(b) > 0) break;
      largeSum = next;
      ++maxPrefixItemsWithTail;
    }
    // Check if the smallest items in the prefix give a tight bound
    // on the the number of items that can fit with the tail
    if (maxPrefixItemsWithTail == 0 || maxPrefixItemsWithTail == prefixSize ||
        tailSum[prefixSize - maxPrefixItemsWithTail - 1] - tailSum[prefixSize] -
                static_cast<double>(b) <=
            feastol)
      continue;
    // All subsets of maxPrefixItemsWithTail items from the prefix
    // fit with the tail.
    for (HighsInt i = 0; i < numTerms; ++i) {
      if (weights[i] == 0) continue;
      weights[i] = i < prefixSize ? 1 : 0;
    }
    b = maxPrefixItemsWithTail;
    return;
  }
}

bool gcdSimplify(KnapsackRow& knapsack, const std::vector<HighsInt>& order,
                 const double feastol) {
  // Turn all the weights and capacity into actual integer types
  const HighsInt numTerms = static_cast<HighsInt>(knapsack.terms.size());
  std::vector<int64_t> weights;
  weights.reserve(numTerms);
  for (const HighsInt i : order) {
    weights.push_back(static_cast<int64_t>(knapsack.terms[i].weight));
  }
  int64_t b = static_cast<int64_t>(static_cast<double>(knapsack.capacity));
  const int64_t originalB = b;
  // See if some smaller values can be eliminated via GCD arguments
  gcdRemoveTail(weights, b, numTerms, feastol);

  // Try repeated-calls of all-but-one GCD calls until unsuccessful
  while (gcdAllButOne(weights, b, numTerms)) {
  }
  if (b != originalB) {
    HighsCDouble newWeightSum = 0;
    for (HighsInt i = 0; i != numTerms; ++i) {
      KnapsackTerm& term = knapsack.terms[order[i]];
      term.weight = static_cast<double>(weights[i]);
      newWeightSum += term.weight;
    }
    knapsack.weightSum = newWeightSum;
    knapsack.capacity = static_cast<double>(b);
    return true;
  }
  return false;
}

bool tryNormaliseConditionalKnapsack(KnapsackRow& knapsack,
                                     const std::vector<HighsInt>& order,
                                     const KnapsackStatistics& stats,
                                     const double feastol, bool& changed) {
  // Determine if the knapsack's LP can be tightened and retain its solution
  // set by being represented by a knapsack with only weights {1, a}
  // It does this by using the largest or smallest weight as a pivot
  // and seeing if forced selecting or forced non-selection
  // result in a fixed number of maximal items that can be subsequently selected
  // Consider the knapsack
  // 8 x_1 + 5 x_2 + 4 x_3 + 4 x_4 + 3 x_5 <= 13
  // if x_1 = 0 -> any three remaining items fit, but no four items fit
  // if x_1 = 1 -> any one item fits, but no two items fit
  // The knapsack can be rewritten as
  // 2 x_1 + x_2 + x_3 + x_4 + x_5 <= 3
  // We only accept this if it dominates the LP of the original row
  // The check selects a discrete fractional amounts of the
  // pivot item, and then greedily packs in the largest items remaining.
  // Why this is needed: Consider
  // 8 x_1 + 5 x_2 + 5 x_3 + 3 x_4 + 3 x_5 <= 13
  // It can also be transformed like above, but
  // the solution (1/2, 1, 1, 0, 0) would now be valid, despite
  // not being permitted in the original knapsack, i.e., no LP dominance.

  auto weight = [&](const HighsInt i) -> double& {
    return knapsack.terms[order[i]].weight;
  };

  if (stats.numNonZeros < 3) return false;
  std::vector<HighsCDouble> largestSums;
  largestSums.reserve(stats.numNonZeros);

  // See if removing the largest or smallest item changes can result
  // in sets of fixed size being selected
  for (const HighsInt pivot : {HighsInt{0}, stats.numNonZeros - 1}) {
    const HighsInt largestItemStart = pivot == 0 ? 1 : 0;
    const HighsInt smallestItemStart = largestItemStart + stats.numNonZeros - 2;
    const double pivotWeight = weight(pivot);
    largestSums.clear();
    largestSums.emplace_back(0);

    const auto getItemLimit = [&](const HighsCDouble capacity,
                                  const bool computeLargest) -> HighsInt {
      if (capacity < feastol) return -1;
      HighsCDouble smallestSum = 0;
      for (HighsInt i = 0; i != stats.numNonZeros - 1; ++i) {
        if (computeLargest) {
          HighsCDouble largestSum = largestSums.back();
          largestSum += weight(largestItemStart + i);
          largestSums.emplace_back(largestSum);
        }
        smallestSum += weight(smallestItemStart - i);
        if (largestSums[i + 1] - capacity > 0) {
          // Return whether smallestSum also exceeds branch capacity
          return smallestSum - capacity > feastol ? i : -1;
        }
      }
      return stats.numNonZeros - 1;
    };

    const HighsInt countWithoutPivot = getItemLimit(knapsack.capacity, true);
    if (countWithoutPivot < 0) continue;
    const HighsInt countWithPivot =
        getItemLimit(knapsack.capacity - pivotWeight, false);
    if (countWithPivot < 0) continue;

    const HighsInt newPivotWeight = countWithoutPivot - countWithPivot;
    if (pivotWeight == newPivotWeight) continue;

    // Only want to accept normalisation if it will dominate the LP relaxation
    bool dominates = true;
    if (newPivotWeight > 1) {
      for (HighsInt i = 0; i != countWithoutPivot; ++i) {
        if (i >= countWithPivot) {
          const HighsCDouble slack =
              (knapsack.capacity - largestSums[i]) * newPivotWeight -
              static_cast<HighsCDouble>(pivotWeight) * (countWithoutPivot - i);
          if (slack < 0) {
            dominates = false;
            break;
          }
        }
      }
    }
    if (!dominates) continue;

    // New knapsack is \sum x_i + newPivotWeight * x_pivot <= countWithoutPivot
    for (KnapsackTerm& term : knapsack.terms) {
      if (term.weight != 0) term.weight = 1;
    }
    knapsack.terms[order[pivot]].weight = newPivotWeight;
    knapsack.weightSum =
        static_cast<HighsCDouble>(stats.numNonZeros - 1) + newPivotWeight;
    knapsack.capacity = countWithoutPivot;
    changed = true;
    return true;
  }
  return false;
}

bool tryNormaliseKnapsack(KnapsackRow& knapsack,
                          const std::vector<HighsInt>& order,
                          const KnapsackStatistics& stats, const double feastol,
                          bool& changed) {
  // Determine if the knapsack's LP can be tightened and retain its solution
  // set by being represented by a knapsack with only weights {1, a}
  // Consider the knapsack
  // 11 x_1 + 5 x_2 + 5 x_3 + 4 x_4 + 4x_5 <= 11
  // x_1 fills the knapsack by itself, and when it's not selected,
  // the two largest "non-max-weight" items fit, but the three smallest
  // "non-max-weight" items don't fit. A stronger knapsack is then
  // 2 x_1 + x_2 + x_3 + x_4 + x_5 <= 2
  // A cheap guard for checking LP dominance is whether the largest
  // "non-max-weight" item can fit new-rhs many times into the
  // old capacity, e.g., 2 * 5 <= 11

  auto weight = [&](const HighsInt i) -> double& {
    return knapsack.terms[order[i]].weight;
  };

  if (stats.numLarge == 0) return false;
  if (stats.numSmall == 0 ||
      (weight(stats.numLarge - 1) == 1 &&
       weight(stats.numLarge + stats.numSmall - 1) == 1)) {
    return true;
  }

  HighsInt numMaxWeightItems = 0;
  while (numMaxWeightItems != stats.numLarge &&
         weight(numMaxWeightItems) == knapsack.capacity) {
    ++numMaxWeightItems;
  }

  if (numMaxWeightItems == 0 ||
      static_cast<HighsCDouble>(weight(numMaxWeightItems - 1)) +
              weight(stats.numNonZeros - 1) - knapsack.capacity <=
          feastol) {
    // Can't normalise if there's no max weight item
    // or if smallest item is too close to feasibility tolerance
    return false;
  }

  auto upgradeRow = [&](const double newWeight) {
    // The largest non-max-weight item must fit newWeight times for
    // normalised knapsack to dominate original w.r.t. LP relaxation.
    if (static_cast<HighsCDouble>(weight(numMaxWeightItems)) * newWeight -
            knapsack.capacity >
        0)
      return false;
    changed = true;
    for (HighsInt i = 0; i != numMaxWeightItems; ++i) {
      knapsack.weightSum -= weight(i) - newWeight;
      weight(i) = newWeight;
    }
    for (HighsInt i = numMaxWeightItems; i != stats.numNonZeros; ++i) {
      knapsack.weightSum -= weight(i) - 1;
      weight(i) = 1;
    }
    knapsack.capacity = newWeight;
    return true;
  };

  const HighsInt numNonZeros = stats.numNonZeros;
  HighsCDouble exceedSum = 0;
  HighsInt v = numMaxWeightItems;
  while (v != numNonZeros && exceedSum - knapsack.capacity <= 0) {
    exceedSum += weight(v);
    v++;
  }
  const HighsCDouble exceedGap = exceedSum - knapsack.capacity;
  const HighsInt numExceed = v - numMaxWeightItems;

  if (exceedGap > 0) {
    // There is some subset of small items that together exceed capacity
    if (exceedGap <= feastol) {
      return false;
    }
    if (v == numNonZeros) {
      // All small items needed to exceed capacity
      return upgradeRow(numExceed - 1);
    }
    HighsCDouble exceedSumBack = 0;
    for (HighsInt i = numNonZeros - 1; i >= numNonZeros - numExceed; --i) {
      exceedSumBack += weight(i);
    }
    if (exceedSumBack - knapsack.capacity > feastol) {
      // If the largest n small weights AND the smallest n small weights
      // both exceed capacity, then the row holds for all n combinations
      return upgradeRow(numExceed - 1);
    }
    return false;
  }
  return upgradeRow(numNonZeros - numMaxWeightItems);
}

void strengthenTwoWeightKnapsack(KnapsackRow& knapsack,
                                 const std::vector<HighsInt>& order,
                                 const KnapsackStatistics& stats,
                                 const double feastol, bool& changed) {
  if (stats.numNonZeros < 2) return;
  const double large = knapsack.terms[order[0]].weight;
  const double small = knapsack.terms[order[stats.numNonZeros - 1]].weight;
  if (large == small) return;
  HighsInt numSameLarge = 0;
  for (HighsInt i = 0; i != stats.numNonZeros; ++i) {
    const double w = knapsack.terms[order[i]].weight;
    if (w == large)
      ++numSameLarge;
    else if (w != small)
      return;
  }
  const HighsInt numSameSmall = stats.numNonZeros - numSameLarge;
  std::vector<HighsInt> maxSmallItems;
  maxSmallItems.reserve(numSameLarge + 1);
  HighsCDouble maxActivity = 0;
  HighsInt smallItemCount = numSameSmall;
  for (HighsInt i = 0; i <= numSameLarge; ++i) {
    const HighsCDouble largeSum = static_cast<HighsCDouble>(large) * i;
    if (largeSum - knapsack.capacity > feastol) break;
    while (smallItemCount > 0 &&
           largeSum + static_cast<HighsCDouble>(small) * smallItemCount -
                   knapsack.capacity >
               feastol) {
      smallItemCount--;
    }
    maxSmallItems.push_back(smallItemCount);
    const HighsCDouble activity =
        largeSum + static_cast<HighsCDouble>(small) * smallItemCount;
    if (activity - maxActivity > 0) maxActivity = activity;
  }

  // May be able to tighten the knapsack
  tightenKnapsackCapacity(knapsack, maxActivity, feastol, changed);

  // Try and lift
  HighsCDouble largeLift = knapsack.capacity;
  HighsCDouble smallLift = knapsack.capacity;
  HighsInt maxSmallItemsSize = static_cast<HighsInt>(maxSmallItems.size());
  for (HighsInt i = 0; i != maxSmallItemsSize; ++i) {
    const HighsInt s = maxSmallItems[i];
    const HighsCDouble slack = knapsack.capacity -
                               static_cast<HighsCDouble>(large) * i -
                               static_cast<HighsCDouble>(small) * s;
    if (i > 0 && slack - largeLift * i < 0) {
      largeLift = slack / i;
    }
    if (s > 0 && slack - smallLift * s < 0) {
      smallLift = slack / s;
    }
  }
  double newLarge = large;
  double newSmall = small;
  if (largeLift > smallLift && largeLift > feastol) {
    newLarge = static_cast<double>(largeLift + large);
  } else if (smallLift > feastol) {
    newSmall = static_cast<double>(smallLift + small);
  }
  if (newLarge == large && newSmall == small) return;

  for (HighsInt i = 0; i != stats.numNonZeros; ++i) {
    double& w = knapsack.terms[order[i]].weight;
    if (w == 0) continue;
    w = w == large ? newLarge : newSmall;
  }
  knapsack.weightSum = static_cast<HighsCDouble>(newLarge) * numSameLarge +
                       static_cast<HighsCDouble>(newSmall) * numSameSmall;
  changed = true;
}

}  // namespace

HPresolve::Result HPresolve::strengthenKnapsacks(
    HighsPostsolveStack& postsolve_stack) {
  const HighsInt maxRowSize =
      std::max(HighsInt{1000}, (model->num_col_ - numDeletedCols) / 20);
  const HighsInt numRows = model->num_row_;
  HighsInt numRedundantRows = 0;
  HighsInt numFixedCols = 0;
  HighsInt numStrengthenedRows = 0;
  std::vector<HighsInt> oldRowPositions;

  auto checkRedundantKnapsack = [&](const KnapsackRow& knapsack,
                                    const HighsInt row) -> bool {
    if (knapsack.weightSum - knapsack.capacity <= primal_feastol) {
      postsolve_stack.redundantRow(row);
      removeRow(row);
      ++numRedundantRows;
      return true;
    }
    return false;
  };

  auto changeRow = [&](const KnapsackRow& knapsack) {
    ++numStrengthenedRows;
    // Replace the knapsack row with a tightened row in the matrix
    oldRowPositions.clear();
    const auto row = getRowVector(knapsack.row);
    for (auto nz = row.begin(); nz != row.end(); ++nz)
      oldRowPositions.push_back(nz.position());

    HighsCDouble rhs = knapsack.capacity;
    size_t nextTerm = 0;
    for (const HighsInt pos : oldRowPositions) {
      const HighsInt col = Acol[pos];
      const double lower = model->col_lower_[col];
      const bool fixed = lower == model->col_upper_[col];
      if (nextTerm == knapsack.terms.size() ||
          knapsack.terms[nextTerm].literal / 2 != col) {
        // Handle columns that were already fixed while the knapsack
        // was created but still exist in the row
        assert(fixed);
        rhs += knapsack.direction * Avalue[pos] * lower;
        continue;
      }

      const KnapsackTerm& term = knapsack.terms[nextTerm++];
      const bool positive = term.literal % 2 == 1;
      if (fixed) {
        // Column is fixed, so eliminate it from the row
        rhs -= static_cast<HighsCDouble>(term.weight) *
               (positive ? lower : 1 - lower);
        addToMatrix(knapsack.row, col, -Avalue[pos]);
        continue;
      }
      const double newCoef = positive ? term.weight : -term.weight;
      if (!positive) rhs -= term.weight;
      const double delta = knapsack.direction * newCoef - Avalue[pos];
      if (delta != 0) addToMatrix(knapsack.row, col, delta);
    }
    assert(nextTerm == knapsack.terms.size());
    if (knapsack.direction == 1)
      changeRowUpper(knapsack.row, static_cast<double>(rhs));
    else
      changeRowLower(knapsack.row, -static_cast<double>(rhs));
  };

  // Iterate over all rows. Only interested in rows of pure binaries
  // where at least one coefficient differs.
  // Will complement binary columns when necessary to make all
  // coefficients positive, i.e., a knapsack.
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
    KnapsackRow knapsack{row, direction, 0, 0, {}};
    knapsack.capacity = rowBound;
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
      // Columns may have been fixed
      if (model->col_lower_[col] == model->col_upper_[col]) {
        const double fixedValue = model->col_lower_[col];
        const double contribution = value * fixedValue;
        if (std::abs(contribution) == kHighsInf) {
          eligible = false;
          break;
        }
        knapsack.capacity -= contribution;
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
        knapsack.weightSum += weight;
        knapsack.terms.push_back({2 * col + (value > 0), weight});
        if (value < 0) knapsack.capacity += weight;
      }
    }
    if (!eligible || allSameCoefficients) continue;
    if (knapsack.capacity < 0) {
      if (knapsack.capacity < -primal_feastol) return Result::kPrimalInfeasible;
      continue;
    }
    if (checkRedundantKnapsack(knapsack, row)) continue;

    if (allIntCoefficients) {
      // Round capacity to be integer if all coefficients are integer
      knapsack.capacity = floor(knapsack.capacity + primal_feastol);
    }

    // Fix any columns with coefficients that are too large
    for (auto& term : knapsack.terms) {
      if (term.weight - knapsack.capacity <= primal_feastol) continue;
      ++numFixedCols;
      const HighsInt col = term.literal / 2;
      if (term.literal % 2 != 0) {
        HPRESOLVE_CHECKED_CALL(changeColUpper(col, 0.0));
      } else {
        HPRESOLVE_CHECKED_CALL(changeColLower(col, 1.0));
      }
      knapsack.weightSum -= term.weight;
      term.weight = 0;
    }

    if (checkRedundantKnapsack(knapsack, row)) continue;

    // If a column's absence makes a row redundant then its weight
    // can be reduced
    const HighsCDouble excess = knapsack.weightSum - knapsack.capacity;
    for (auto& term : knapsack.terms) {
      if (term.weight - excess <= primal_feastol) continue;
      const HighsCDouble delta = term.weight - excess;
      term.weight = static_cast<double>(excess);
      knapsack.capacity -= delta;
      knapsack.weightSum -= delta;
      changed = true;
    }

    // Sort terms in non-increasing order
    std::vector<HighsInt> order(knapsack.terms.size());
    std::iota(order.begin(), order.end(), HighsInt{0});
    const auto isItemLarger = [&](const HighsInt i, const HighsInt j) {
      return knapsack.terms[i].weight > knapsack.terms[j].weight;
    };
    pdqsort(order.begin(), order.end(), isItemLarger);

    if (allIntCoefficients && knapsack.capacity < kHighsIInf) {
      // Simplify constraint (or strengthen) using GCD arguments
      changed |= gcdSimplify(knapsack, order, primal_feastol);
      if (checkRedundantKnapsack(knapsack, row)) continue;
    }

    KnapsackStatistics stats;
    tightenRhs(knapsack, order, primal_feastol, stats, changed);
    bool lifted = false;
    const bool cliqueUpgraded =
        liftKnapsackWeights(knapsack, order, primal_feastol, stats, lifted);
    changed |= lifted;
    // Recognise whether the knapsack can be upgraded
    if (!cliqueUpgraded) {
      if (lifted) {
        pdqsort(order.begin(), order.end(), isItemLarger);
        if (allIntCoefficients && knapsack.capacity < kHighsIInf) {
          // Simplify constraint (or strengthen) using GCD arguments
          // Lifted coefficients may have enabled new reductions
          const bool gcdSimplified =
              gcdSimplify(knapsack, order, primal_feastol);
          if (checkRedundantKnapsack(knapsack, row)) continue;
          if (gcdSimplified) {
            clearKnapSackStatistics(stats);
            computeKnapsackStatistics(knapsack, order, primal_feastol, stats);
          }
        }
      }
      if (!tryNormaliseKnapsack(knapsack, order, stats, primal_feastol,
                                changed) &&
          !tryNormaliseConditionalKnapsack(knapsack, order, stats,
                                           primal_feastol, changed)) {
        strengthenTwoWeightKnapsack(knapsack, order, stats, primal_feastol,
                                    changed);
      }
    }

    if (changed) {
      changeRow(knapsack);
    }
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Removed %" HIGHSINT_FORMAT
              " rows, "
              "fixed %" HIGHSINT_FORMAT
              " cols, "
              "strengthened %" HIGHSINT_FORMAT " rows\n",
              numRedundantRows, numFixedCols, numStrengthenedRows);

  return Result::kOk;
}

}  // namespace presolve
