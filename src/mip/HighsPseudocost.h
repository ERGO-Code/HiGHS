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

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

class HighsMipSolver;

class HighsPseudocost {
  std::vector<double> pseudocostup;
  std::vector<double> pseudocostdown;
  std::vector<int> nsamplesup;
  std::vector<int> nsamplesdown;
  std::vector<int> ncutoffsup;
  std::vector<int> ncutoffsdown;

  double cost_total;
  size_t nsamplestotal;
  size_t ncutoffstotal;
  int minreliable;

 public:
  HighsPseudocost() = default;
  HighsPseudocost(const HighsMipSolver& mipsolver);

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

  int getNumObservationsUp(int col) const { return nsamplesup[col]; }

  int getNumObservationsDown(int col) const { return nsamplesdown[col]; }

  void addCutoffObservation(int col, bool upbranch) {
    ++ncutoffstotal;
    if (upbranch)
      ncutoffsup[col] += 1;
    else
      ncutoffsdown[col] += 1;
  }

  void addObservation(int col, double delta, double objdelta) {
    assert(delta != 0.0);
    assert(objdelta >= 0.0);
    if (delta > 0.0) {
      double unit_gain = objdelta / delta;
      double d = unit_gain - pseudocostup[col];
      nsamplesup[col] += 1;
      pseudocostup[col] += d / nsamplesup[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / nsamplestotal;
    } else {
      double unit_gain = -objdelta / delta;
      double d = unit_gain - pseudocostdown[col];
      nsamplesdown[col] += 1;
      pseudocostdown[col] += d / nsamplesdown[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / nsamplestotal;
    }
  }

  bool isReliable(int col) const {
    return std::min(nsamplesup[col] + ncutoffsup[col],
                    nsamplesdown[col] + ncutoffsdown[col]) >= minreliable;
  }

  bool isReliableUp(int col) const {
    return nsamplesup[col] + ncutoffsup[col] >= minreliable;
  }

  bool isReliableDown(int col) const {
    return nsamplesdown[col] + ncutoffsdown[col] >= minreliable;
  }

  double getAvgPseudocost() const { return cost_total; }

  double getPseudocostUp(int col, double frac, double offset) const {
    double up = std::ceil(frac) - frac;
    double cost;

    if (nsamplesup[col] == 0 || nsamplesup[col] < minreliable) {
      double weightPs = nsamplesup[col] == 0 ? 0
                                             : 0.75 + 0.25 * nsamplesup[col] /
                                                          (double)minreliable;
      cost = weightPs * pseudocostup[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostup[col];
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
      cost = weightPs * pseudocostdown[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostdown[col];

    return down * (offset + cost);
  }

  double getPseudocostUp(int col, double frac) const {
    double up = std::ceil(frac) - frac;
    if (nsamplesup[col] == 0) return up * cost_total;
    return up * pseudocostup[col];
  }

  double getPseudocostDown(int col, double frac) const {
    double down = frac - std::floor(frac);
    if (nsamplesdown[col] == 0) return down * cost_total;
    return down * pseudocostdown[col];
  }

  double getScore(int col, double upcost, double downcost) const {
    if (cost_total > 1e-12) {
      double avgCost = getAvgPseudocost();

      double cutoffweightup =
          ncutoffsup[col] /
          std::max(1.0, double(ncutoffsup[col] + nsamplesup[col]));
      double costweightup = 1.0 - cutoffweightup;
      upcost = upcost * costweightup + 2 * avgCost * cutoffweightup;

      double cutoffweightdown =
          ncutoffsdown[col] /
          std::max(1.0, double(ncutoffsdown[col] + nsamplesdown[col]));
      double costweightdown = 1.0 - cutoffweightdown;
      downcost = downcost * costweightdown + 2 * avgCost * cutoffweightdown;
    } else {
      upcost +=
          ncutoffsup[col] == 0
              ? 0.0
              : (ncutoffsup[col] / (double(ncutoffsup[col] + nsamplesup[col])));
      downcost += ncutoffsdown[col] == 0
                      ? 0.0
                      : (ncutoffsdown[col] /
                         (double(ncutoffsdown[col] + nsamplesdown[col])));
    }
    return upcost * downcost;
  }

  double getScore(int col, double frac) const {
    double upcost = getPseudocostUp(col, frac);
    double downcost = getPseudocostDown(col, frac);

    if (cost_total > 1e-12) {
      double avgCost = getAvgPseudocost();

      double cutoffweightup =
          ncutoffsup[col] /
          std::max(1.0, double(ncutoffsup[col] + nsamplesup[col]));
      double costweightup = 1.0 - cutoffweightup;
      upcost = upcost * costweightup + 2 * avgCost * cutoffweightup;

      double cutoffweightdown =
          ncutoffsdown[col] /
          std::max(1.0, double(ncutoffsdown[col] + nsamplesdown[col]));
      double costweightdown = 1.0 - cutoffweightdown;
      downcost = downcost * costweightdown + 2 * avgCost * cutoffweightdown;
    } else {
      upcost +=
          ncutoffsup[col] == 0
              ? 0.0
              : (ncutoffsup[col] / (double(ncutoffsup[col] + nsamplesup[col])));
      downcost += ncutoffsdown[col] == 0
                      ? 0.0
                      : (ncutoffsdown[col] /
                         (double(ncutoffsdown[col] + nsamplesdown[col])));
    }

    return upcost * downcost;
  }
};

#endif
