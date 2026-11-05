/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsRedcostFixing.h
 * @brief reduced cost fixing using the current cutoff bound
 */

#ifndef HIGHS_REDCOST_FIXING_H_
#define HIGHS_REDCOST_FIXING_H_

#include <algorithm>
#include <map>
#include <vector>

#include "HighsConflictPool.h"
#include "mip/HighsDomainChange.h"

class HighsDomain;
class HighsMipSolver;
class HighsLpRelaxation;

// 10xb: flat sorted vector replacing std::multimap<double, HighsInt> for the
// per-column lurking-bound staircases (typically a few dozen entries that are
// rebuilt on every root LP resolve; node-based maps made this churn costly).
struct HighsLurkingStair {
  using value_type = std::pair<double, HighsInt>;
  std::vector<value_type> v;
  using iterator = std::vector<value_type>::iterator;
  using const_iterator = std::vector<value_type>::const_iterator;
  iterator begin() { return v.begin(); }
  iterator end() { return v.end(); }
  const_iterator begin() const { return v.begin(); }
  const_iterator end() const { return v.end(); }
  bool empty() const { return v.empty(); }
  size_t size() const { return v.size(); }
  void clear() { v.clear(); }
  iterator lower_bound(double k) {
    return std::lower_bound(
        v.begin(), v.end(), k,
        [](const value_type& a, double b) { return a.first < b; });
  }
  iterator upper_bound(double k) {
    return std::upper_bound(
        v.begin(), v.end(), k,
        [](double a, const value_type& b) { return a < b.first; });
  }
  iterator erase(iterator a, iterator b) { return v.erase(a, b); }
  iterator erase(iterator a) { return v.erase(a); }
  iterator emplace(double k, HighsInt val) {
    return v.insert(upper_bound(k), value_type(k, val));
  }
};

class HighsRedcostFixing {
  std::vector<HighsLurkingStair> lurkingColUpper;
  std::vector<HighsLurkingStair> lurkingColLower;

 public:
  std::vector<std::pair<double, HighsDomainChange>> getLurkingBounds(
      const HighsMipSolver& mipsolver, const HighsDomain& globaldom) const;

  void propagateRootRedcost(const HighsMipSolver& mipsolver);

  static void propagateRedCost(const HighsMipSolver& mipsolver,
                               HighsDomain& localdomain, HighsDomain& globaldom,
                               const HighsLpRelaxation& lp,
                               HighsConflictPool& conflictpool,
                               HighsPseudocost& pseudocost, double upper_limit);

  void addRootRedcost(const HighsMipSolver& mipsolver,
                      const std::vector<double>& lpredcost, double lpobjective);
};

#endif
