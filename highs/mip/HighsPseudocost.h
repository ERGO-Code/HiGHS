/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_PSEUDOCOST_H_
#define HIGHS_PSEUDOCOST_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include "util/HighsInt.h"

class HighsMipSolver;
namespace presolve {
class HighsPostsolveStack;
}

class HighsPseudocost;

constexpr double minThreshold = 1e-6;

struct HighsPseudocostInitialization {
  std::vector<double> pseudocostup;
  std::vector<double> pseudocostdown;
  std::vector<HighsInt> nsamplesup;
  std::vector<HighsInt> nsamplesdown;
  std::vector<double> inferencesup;
  std::vector<double> inferencesdown;
  std::vector<HighsInt> ninferencesup;
  std::vector<HighsInt> ninferencesdown;
  std::vector<double> conflictscoreup;
  std::vector<double> conflictscoredown;
  double cost_total;
  double inferences_total;
  double conflict_avg_score;
  int64_t nsamplestotal;
  int64_t ninferencestotal;

  HighsPseudocostInitialization(const HighsPseudocost& pscost,
                                HighsInt maxCount);
  HighsPseudocostInitialization(
      const HighsPseudocost& pscost, HighsInt maxCount,
      const presolve::HighsPostsolveStack& postsolveStack);
};

class HighsPseudocost {
  friend struct HighsPseudocostInitialization;
  std::vector<double> pseudocostup;
  std::vector<double> pseudocostdown;
  std::vector<HighsInt> nsamplesup;
  std::vector<HighsInt> nsamplesdown;
  std::vector<double> inferencesup;
  std::vector<double> inferencesdown;
  std::vector<HighsInt> ninferencesup;
  std::vector<HighsInt> ninferencesdown;
  std::vector<HighsInt> ncutoffsup;
  std::vector<HighsInt> ncutoffsdown;
  std::vector<double> conflictscoreup;
  std::vector<double> conflictscoredown;
  std::vector<bool> changed;
  std::vector<HighsInt> indschanged;

  double conflict_weight;
  double conflict_avg_score;
  double cost_total;
  double inferences_total;
  int64_t nsamplestotal;
  int64_t ninferencestotal;
  int64_t ncutoffstotal;
  HighsInt minreliable;
  double degeneracyFactor;

 public:
  HighsPseudocost() = default;
  HighsPseudocost(const HighsMipSolver& mipsolver);

  void subtractBase(const HighsPseudocost& base) {
    for (size_t i = 0; i != pseudocostup.size(); ++i) {
      pseudocostup[i] -= base.pseudocostup[i];
      pseudocostdown[i] -= base.pseudocostdown[i];
      nsamplesup[i] -= base.nsamplesup[i];
      nsamplesdown[i] -= base.nsamplesdown[i];
    }
  }

  void increaseConflictWeight() {
    conflict_weight *= 1.02;

    if (conflict_weight > 1000.0) {
      double scale = 1.0 / conflict_weight;
      conflict_weight = 1.0;
      conflict_avg_score *= scale;

      for (size_t i = 0; i != conflictscoreup.size(); ++i) {
        conflictscoreup[i] *= scale;
        conflictscoredown[i] *= scale;
      }
    }
  }

  void setDegeneracyFactor(double degeneracyFactor) {
    assert(degeneracyFactor >= 1.0);
    this->degeneracyFactor = degeneracyFactor;
  }

  void increaseConflictScoreUp(HighsInt col) {
    conflictscoreup[col] += conflict_weight;
    conflict_avg_score += conflict_weight;
    if (!changed[col]) {
      changed[col] = true;
      indschanged.push_back(col);
    }
  }

  void increaseConflictScoreDown(HighsInt col) {
    conflictscoredown[col] += conflict_weight;
    conflict_avg_score += conflict_weight;
    if (!changed[col]) {
      changed[col] = true;
      indschanged.push_back(col);
    }
  }

  void setMinReliable(HighsInt minreliable) { this->minreliable = minreliable; }

  HighsInt getMinReliable() const { return minreliable; }

  HighsInt getNumObservations(HighsInt col) const {
    return nsamplesup[col] + nsamplesdown[col];
  }

  HighsInt getNumObservationsUp(HighsInt col) const { return nsamplesup[col]; }

  HighsInt getNumObservationsDown(HighsInt col) const {
    return nsamplesdown[col];
  }

  void addCutoffObservation(HighsInt col, bool upbranch) {
    ++ncutoffstotal;
    if (upbranch)
      ncutoffsup[col] += 1;
    else
      ncutoffsdown[col] += 1;
    if (!changed[col]) {
      changed[col] = true;
      indschanged.push_back(col);
    }
  }

  void addObservation(HighsInt col, double delta, double objdelta) {
    assert(delta != 0.0);
    assert(objdelta >= 0.0);
    if (delta > 0.0) {
      double unit_gain = objdelta / delta;
      double d = unit_gain - pseudocostup[col];
      nsamplesup[col] += 1;
      pseudocostup[col] += d / nsamplesup[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / static_cast<double>(nsamplestotal);
    } else {
      double unit_gain = -objdelta / delta;
      double d = unit_gain - pseudocostdown[col];
      nsamplesdown[col] += 1;
      pseudocostdown[col] += d / nsamplesdown[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / static_cast<double>(nsamplestotal);
    }
    if (!changed[col]) {
      changed[col] = true;
      indschanged.push_back(col);
    }
  }

  void addInferenceObservation(HighsInt col, HighsInt ninferences,
                               bool upbranch) {
    double d = ninferences - inferences_total;
    ++ninferencestotal;
    inferences_total += d / static_cast<double>(ninferencestotal);
    if (upbranch) {
      d = ninferences - inferencesup[col];
      ninferencesup[col] += 1;
      inferencesup[col] += d / ninferencesup[col];
    } else {
      d = ninferences - inferencesdown[col];
      ninferencesdown[col] += 1;
      inferencesdown[col] += d / ninferencesdown[col];
    }
    if (!changed[col]) {
      changed[col] = true;
      indschanged.push_back(col);
    }
  }

  bool isReliable(HighsInt col) const {
    return std::min(nsamplesup[col], nsamplesdown[col]) >= minreliable;
  }

  bool isReliableUp(HighsInt col) const {
    return nsamplesup[col] >= minreliable;
  }

  bool isReliableDown(HighsInt col) const {
    return nsamplesdown[col] >= minreliable;
  }

  double getAvgPseudocost() const { return cost_total; }

  double getPseudocostUp(HighsInt col, double frac, double offset) const {
    double up = std::ceil(frac) - frac;
    double cost;

    if (nsamplesup[col] == 0 || nsamplesup[col] < minreliable) {
      double weightPs =
          nsamplesup[col] == 0
              ? 0
              : 0.9 + 0.1 * nsamplesup[col] / static_cast<double>(minreliable);
      cost = weightPs * pseudocostup[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostup[col];
    return up * (offset + cost);
  }

  double getPseudocostDown(HighsInt col, double frac, double offset) const {
    double down = frac - std::floor(frac);
    double cost;

    if (nsamplesdown[col] == 0 || nsamplesdown[col] < minreliable) {
      double weightPs = nsamplesdown[col] == 0
                            ? 0
                            : 0.9 + 0.1 * nsamplesdown[col] /
                                        static_cast<double>(minreliable);
      cost = weightPs * pseudocostdown[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostdown[col];

    return down * (offset + cost);
  }

  double getPseudocostUp(HighsInt col, double frac) const {
    double up = std::ceil(frac) - frac;
    if (nsamplesup[col] == 0) return up * cost_total;
    return up * pseudocostup[col];
  }

  double getPseudocostDown(HighsInt col, double frac) const {
    double down = frac - std::floor(frac);
    if (nsamplesdown[col] == 0) return down * cost_total;
    return down * pseudocostdown[col];
  }

  double getConflictScoreUp(HighsInt col) const {
    return conflictscoreup[col] / conflict_weight;
  }

  double getConflictScoreDown(HighsInt col) const {
    return conflictscoredown[col] / conflict_weight;
  }

  double getScore(HighsInt col, double upcost, double downcost) const {
    double costScore = std::max(upcost, minThreshold) *
                       std::max(downcost, minThreshold) /
                       std::max(minThreshold, cost_total * cost_total);
    double inferenceScore =
        std::max(inferencesup[col], minThreshold) *
        std::max(inferencesdown[col], minThreshold) /
        std::max(minThreshold, inferences_total * inferences_total);

    double cutOffScoreUp =
        ncutoffsup[col] /
        std::max(1.0, static_cast<double>(ncutoffsup[col]) +
                          static_cast<double>(nsamplesup[col]));
    double cutOffScoreDown =
        ncutoffsdown[col] /
        std::max(1.0, static_cast<double>(ncutoffsdown[col]) +
                          static_cast<double>(nsamplesdown[col]));
    double avgCutoffs = static_cast<double>(ncutoffstotal) /
                        std::max(1.0, static_cast<double>(ncutoffstotal) +
                                          static_cast<double>(nsamplestotal));

    double cutoffScore = std::max(cutOffScoreUp, minThreshold) *
                         std::max(cutOffScoreDown, minThreshold) /
                         std::max(minThreshold, avgCutoffs * avgCutoffs);

    double conflictScoreUp = conflictscoreup[col] / conflict_weight;
    double conflictScoreDown = conflictscoredown[col] / conflict_weight;
    double conflictScoreAvg =
        conflict_avg_score /
        (conflict_weight * static_cast<double>(conflictscoreup.size()));
    double conflictScore =
        std::max(conflictScoreUp, minThreshold) *
        std::max(conflictScoreDown, minThreshold) /
        std::max(minThreshold, conflictScoreAvg * conflictScoreAvg);

    auto mapScore = [](double score) { return 1.0 - 1.0 / (1.0 + score); };
    return mapScore(costScore) / degeneracyFactor +
           degeneracyFactor *
               (1e-2 * mapScore(conflictScore) +
                1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  double getScore(HighsInt col, double frac) const {
    double upcost = getPseudocostUp(col, frac);
    double downcost = getPseudocostDown(col, frac);

    return getScore(col, upcost, downcost);
  }

  double getScoreUp(HighsInt col, double frac) const {
    double costScore =
        getPseudocostUp(col, frac) / std::max(minThreshold, cost_total);
    double inferenceScore =
        inferencesup[col] / std::max(minThreshold, inferences_total);

    double cutOffScoreUp =
        ncutoffsup[col] /
        std::max(1.0, static_cast<double>(ncutoffsup[col]) +
                          static_cast<double>(nsamplesup[col]));
    double avgCutoffs = static_cast<double>(ncutoffstotal) /
                        std::max(1.0, static_cast<double>(ncutoffstotal) +
                                          static_cast<double>(nsamplestotal));

    double cutoffScore = cutOffScoreUp / std::max(minThreshold, avgCutoffs);

    double conflictScoreUp = conflictscoreup[col] / conflict_weight;
    double conflictScoreAvg =
        conflict_avg_score /
        (conflict_weight * static_cast<double>(conflictscoreup.size()));
    double conflictScore =
        conflictScoreUp / std::max(minThreshold, conflictScoreAvg);

    auto mapScore = [](double score) { return 1.0 - 1.0 / (1.0 + score); };

    return mapScore(costScore) +
           (1e-2 * mapScore(conflictScore) +
            1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  double getScoreDown(HighsInt col, double frac) const {
    double costScore =
        getPseudocostDown(col, frac) / std::max(minThreshold, cost_total);
    double inferenceScore =
        inferencesdown[col] / std::max(minThreshold, inferences_total);

    double cutOffScoreDown =
        ncutoffsdown[col] /
        std::max(1.0, static_cast<double>(ncutoffsdown[col]) +
                          static_cast<double>(nsamplesdown[col]));
    double avgCutoffs = static_cast<double>(ncutoffstotal) /
                        std::max(1.0, static_cast<double>(ncutoffstotal) +
                                          static_cast<double>(nsamplestotal));

    double cutoffScore = cutOffScoreDown / std::max(minThreshold, avgCutoffs);

    double conflictScoreDown = conflictscoredown[col] / conflict_weight;
    double conflictScoreAvg =
        conflict_avg_score /
        (conflict_weight * static_cast<double>(conflictscoredown.size()));
    double conflictScore =
        conflictScoreDown / std::max(minThreshold, conflictScoreAvg);

    auto mapScore = [](double score) { return 1.0 - 1.0 / (1.0 + score); };

    return mapScore(costScore) +
           (1e-2 * mapScore(conflictScore) +
            1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  double getAvgInferencesUp(HighsInt col) const { return inferencesup[col]; }

  double getAvgInferencesDown(HighsInt col) const {
    return inferencesdown[col];
  }

  std::vector<HighsInt> getNSamplesUp() const { return nsamplesup; }

  std::vector<HighsInt> getNSamplesDown() const { return nsamplesdown; }

  std::vector<HighsInt> getNInferencesUp() const { return ninferencesup; }

  std::vector<HighsInt> getNInferencesDown() const { return ninferencesdown; }

  std::vector<HighsInt> getNCutoffsUp() const { return ncutoffsup; }

  std::vector<HighsInt> getNCutoffsDown() const { return ncutoffsdown; }

  void flushPseudoCostObservations(double& curr_observation,
                                   const double& new_observation,
                                   const HighsInt n_new, const HighsInt n_prev,
                                   const HighsInt n_curr, bool inference) {
    const HighsInt n = n_new - n_prev;
    if (n > 0) {
      const double r = static_cast<double>(n) / (n_curr + n);
      const double average = (1 - r) * curr_observation + r * new_observation;
      curr_observation = average;
      if (inference) {
        this->ninferencestotal += n;
      } else {
        this->nsamplestotal += n;
      }
    }
  }

  void flushCutoffObservations(HighsInt& curr_observation,
                               const HighsInt& prev_observation,
                               const HighsInt& new_observation) {
    HighsInt delta = new_observation - prev_observation;
    curr_observation += delta;
    this->ncutoffstotal += delta;
  }

  void flushConflictObservations(double& curr_observation,
                                 double new_observation,
                                 double conflict_weight) {
    double s = this->conflict_weight *
               std::max(curr_observation / this->conflict_weight,
                        new_observation / conflict_weight);
    if (s > curr_observation + minThreshold) {
      this->conflict_avg_score += s - curr_observation;
    }
    curr_observation = s;
  }

  void flushPseudoCost(HighsPseudocost& pseudocost,
                       std::vector<HighsInt>& nsamplesup,
                       std::vector<HighsInt>& nsamplesdown,
                       std::vector<HighsInt>& ninferencesup,
                       std::vector<HighsInt>& ninferencesdown,
                       std::vector<HighsInt>& ncutoffsup,
                       std::vector<HighsInt>& ncutoffsdown) {
    int64_t orig_nsamplestotal = this->nsamplestotal;
    int64_t orig_ninferencestotal = this->ninferencestotal;
    assert(pseudocost.ncutoffsup.size() == ncutoffsup.size() &&
           pseudocost.ncutoffsup.size() == this->ncutoffsup.size());
    for (HighsInt col : pseudocost.indschanged) {
      assert(col >= 0 &&
             col < static_cast<HighsInt>(pseudocost.ncutoffsup.size()));
      flushPseudoCostObservations(this->pseudocostup[col],
                                  pseudocost.pseudocostup[col], nsamplesup[col],
                                  pseudocost.nsamplesup[col],
                                  this->nsamplesup[col], false);
      flushPseudoCostObservations(
          this->pseudocostdown[col], pseudocost.pseudocostdown[col],
          nsamplesdown[col], pseudocost.nsamplesdown[col],
          this->nsamplesdown[col], false);
      flushPseudoCostObservations(
          this->inferencesup[col], pseudocost.inferencesup[col],
          ninferencesup[col], pseudocost.ninferencesup[col],
          this->ninferencesup[col], true);
      flushPseudoCostObservations(
          this->inferencesdown[col], pseudocost.inferencesdown[col],
          ninferencesdown[col], pseudocost.ninferencesdown[col],
          this->ninferencesdown[col], true);
      // Take the max conflict score (no way to guess num observations)
      flushConflictObservations(this->conflictscoreup[col],
                                pseudocost.conflictscoreup[col],
                                pseudocost.conflict_weight);
      flushConflictObservations(this->conflictscoredown[col],
                                pseudocost.conflictscoredown[col],
                                pseudocost.conflict_weight);
      flushCutoffObservations(this->ncutoffsup[col], ncutoffsup[col],
                              pseudocost.ncutoffsup[col]);
      flushCutoffObservations(this->ncutoffsdown[col], ncutoffsdown[col],
                              pseudocost.ncutoffsdown[col]);
      pseudocost.changed[col] = false;
    }
    pseudocost.indschanged.clear();
    if (this->ninferencestotal > orig_ninferencestotal) {
      const double r = static_cast<double>(orig_ninferencestotal) /
                       static_cast<double>(this->ninferencestotal);
      this->inferences_total =
          r * this->inferences_total + (1 - r) * pseudocost.inferences_total;
    }
    if (this->nsamplestotal > orig_nsamplestotal) {
      const double r = static_cast<double>(orig_nsamplestotal) /
                       static_cast<double>(this->nsamplestotal);
      this->cost_total = r * this->cost_total + (1 - r) * pseudocost.cost_total;
    }
  }

  void syncPseudoCost(HighsPseudocost& pseudocost) {
    std::copy(pseudocostup.begin(), pseudocostup.end(),
              pseudocost.pseudocostup.begin());
    std::copy(pseudocostdown.begin(), pseudocostdown.end(),
              pseudocost.pseudocostdown.begin());
    std::copy(nsamplesup.begin(), nsamplesup.end(),
              pseudocost.nsamplesup.begin());
    std::copy(nsamplesdown.begin(), nsamplesdown.end(),
              pseudocost.nsamplesdown.begin());
    std::copy(inferencesup.begin(), inferencesup.end(),
              pseudocost.inferencesup.begin());
    std::copy(inferencesdown.begin(), inferencesdown.end(),
              pseudocost.inferencesdown.begin());
    std::copy(ninferencesup.begin(), ninferencesup.end(),
              pseudocost.ninferencesup.begin());
    std::copy(ninferencesdown.begin(), ninferencesdown.end(),
              pseudocost.ninferencesdown.begin());
    std::copy(ncutoffsup.begin(), ncutoffsup.end(),
              pseudocost.ncutoffsup.begin());
    std::copy(ncutoffsdown.begin(), ncutoffsdown.end(),
              pseudocost.ncutoffsdown.begin());
    std::copy(conflictscoreup.begin(), conflictscoreup.end(),
              pseudocost.conflictscoreup.begin());
    std::copy(conflictscoredown.begin(), conflictscoredown.end(),
              pseudocost.conflictscoredown.begin());
    pseudocost.conflict_weight = conflict_weight;
    pseudocost.conflict_avg_score = conflict_avg_score;
    pseudocost.cost_total = cost_total;
    pseudocost.inferences_total = inferences_total;
    pseudocost.nsamplestotal = nsamplestotal;
    pseudocost.ninferencestotal = ninferencestotal;
    pseudocost.ncutoffstotal = ncutoffstotal;
  }

  void removeChanged() {
    for (HighsInt col : indschanged) {
      changed[col] = false;
    }
    indschanged.clear();
  }
};

#endif
