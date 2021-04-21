/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_PSEUDOCOST_H_
#define HIGHS_PSEUDOCOST_H_

#include <cmath>
#include <limits>
#include <vector>

#include "util/HighsRandom.h"

class HighsMipSolver;

class HighsPseudocost {
  std::vector<double> pseudocostup;
  std::vector<double> pseudocostdown;
  std::vector<int> nsamplesup;
  std::vector<int> nsamplesdown;

  double cost_total;
  size_t nsamplestotal;
  int minreliable;
  unsigned seed;

 public:
  HighsPseudocost(int ncols, unsigned int seed = 0x533D)
      : pseudocostup(ncols),
        pseudocostdown(ncols),
        nsamplesup(ncols),
        nsamplesdown(ncols),
        cost_total(0),
        nsamplestotal(0),
        minreliable(8),
        seed(seed) {}

  void setSeed(unsigned int seed) { seed = seed; }

  void subtractBase(const HighsPseudocost& base) {
    int ncols = pseudocostup.size();

    for (int i = 0; i != ncols; ++i) {
      pseudocostup[i] -= base.pseudocostup[i];
      pseudocostdown[i] -= base.pseudocostdown[i];
      nsamplesup[i] -= base.nsamplesup[i];
      nsamplesdown[i] -= base.nsamplesdown[i];
    }
  }

  void setMinReliable(int minreliable) { this->minreliable = minreliable; }

  int getMinReliable() const { return minreliable; }

  int getNumObservations(int col) const {
    return nsamplesup[col] + nsamplesdown[col];
  }

  void addObservation(int col, double delta, double objdelta) {
    assert(delta != 0.0);
    assert(objdelta >= 0.0);
    if (delta > 0.0) {
      pseudocostup[col] += objdelta / delta;
      nsamplesup[col] += 1;
      cost_total += objdelta / delta;
      ++nsamplestotal;
    } else {
      pseudocostdown[col] -= objdelta / delta;
      cost_total -= objdelta / delta;
      ++nsamplestotal;
      nsamplesdown[col] += 1;
    }
  }

  bool isReliable(int col) const {
    return std::min(nsamplesup[col], nsamplesdown[col]) >= minreliable;
  }

  bool isReliableUp(int col) const { return nsamplesup[col] >= minreliable; }

  bool isReliableDown(int col) const {
    return nsamplesdown[col] >= minreliable;
  }

  double getAvgPseudocost() const {
    return nsamplestotal == 0 ? 0 : cost_total / nsamplestotal;
  }

  double getPseudocostUp(int col, double frac, double offset) const {
    double up = std::ceil(frac) - frac;
    double cost;

    if (nsamplesup[col] == 0 || nsamplesup[col] < minreliable) {
      double weightPs = nsamplesup[col] == 0 ? 0
                                             : 0.75 + 0.25 * nsamplesup[col] /
                                                          (double)minreliable;
      cost = nsamplesup[col] == 0
                 ? 0
                 : weightPs * pseudocostup[col] / nsamplesup[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostup[col] / nsamplesup[col];
    return up * (offset + cost);
  }

  double getPseudocostDown(int col, double frac, double offset) const {
    double down = frac - std::floor(frac);
    double cost;

    if (nsamplesdown[col] == 0 || nsamplesdown[col] < minreliable) {
      double weightPs =
          nsamplesdown[col] == 0
              ? 0
              : 0.75 + 0.25 * nsamplesdown[col] / (double)minreliable;
      cost = nsamplesdown[col] == 0
                 ? 0
                 : weightPs * pseudocostdown[col] / nsamplesdown[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostdown[col] / nsamplesdown[col];

    return down * (offset + cost);
  }

  double getPseudocostUp(int col, double frac) const {
    double up = std::ceil(frac) - frac;
    if (nsamplesup[col] == 0)
      return nsamplestotal == 0 ? 0 : up * cost_total / nsamplestotal;
    return up * pseudocostup[col] / nsamplesup[col];
  }

  double getPseudocostDown(int col, double frac) const {
    double down = frac - std::floor(frac);
    if (nsamplesdown[col] == 0)
      return nsamplestotal == 0 ? 0 : down * cost_total / nsamplestotal;
    return down * pseudocostdown[col] / nsamplesdown[col];
  }

  double getScore(int col, double upcost, double downcost) const {
    return upcost * downcost;
  }

  double getScore(int col, double frac) const {
    double upcost = getPseudocostUp(col, frac);
    double downcost = getPseudocostDown(col, frac);

    return upcost * downcost;
  }
};

#endif