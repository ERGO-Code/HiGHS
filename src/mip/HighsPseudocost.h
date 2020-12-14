
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
  int minreliable;
  unsigned seed;

 public:
  HighsPseudocost(int ncols, unsigned int seed = 0x533D)
      : pseudocostup(ncols),
        pseudocostdown(ncols),
        nsamplesup(ncols),
        nsamplesdown(ncols),
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
    } else {
      pseudocostdown[col] -= objdelta / delta;
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

  double getPseudocostUp(int col, double frac) const {
    if (nsamplesup[col] == 0) return 0;

    double up = std::ceil(frac) - frac;
    return up * pseudocostup[col] / nsamplesup[col];
  }

  double getPseudocostDown(int col, double frac) const {
    if (nsamplesdown[col] == 0) return 0;

    double down = frac - std::floor(frac);
    return down * pseudocostdown[col] / nsamplesdown[col];
  }

  double getScore(int col, double upcost, double downcost) const {
    // each column is assigned a random hash based on the seed
    // for which we compute a random number in [0,1e-5]
    // This number is then added to the final score for random tie breaking
    unsigned hash = (unsigned(col) + seed) * 0x9e3779b9u;
    double randomval =
        1e-5 * hash / (double)std::numeric_limits<unsigned>::max();
    return randomval + upcost * downcost;
  }

  double getScore(int col, double frac) const {
    // each column is assigned a random hash based on the seed
    // for which we compute a random number in [0,1e-5]
    // This number is then added to the final score for random tie breaking
    unsigned hash = (unsigned(col) + seed) * 0x9e3779b9u;
    double randomval =
        1e-5 * hash / (double)std::numeric_limits<unsigned>::max();

    if (nsamplesup[col] == 0 || nsamplesdown[col] == 0) return randomval;

    double up = std::ceil(frac) - frac;
    double down = frac - std::floor(frac);
    return randomval + (up * pseudocostup[col] / nsamplesup[col]) *
                           (down * pseudocostdown[col] / nsamplesdown[col]);
  }
};

#endif