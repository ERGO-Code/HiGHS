#ifndef HIGHS_CONFLICT_POOL_H_
#define HIGHS_CONFLICT_POOL_H_
#include <set>
#include <vector>

#include "mip/HighsDomainChange.h"

class HighsDomain;

class HighsConflictPool {
  // head and tail pointers for linked list of conflict changes for each column
  std::vector<int> head_;
  std::vector<int> tail_;

  std::vector<HighsDomainChange> conflictchange_;
  std::vector<int> conflictindex_;

  // next and prev links for column
  std::vector<int> next_;
  std::vector<int> prev_;
  std::vector<std::pair<int, int>> conflictrange_;

  std::vector<int> deletedindices_;
  std::set<std::pair<int, int>> freespaces_;

 private:
  int getStateContribution(int confchgidx, const HighsDomain& domain) const;

  void unlink(int confchgidx);

  void link(int confchgidx);

 public:
  HighsConflictPool(int ncols) {
    head_.resize(ncols, -1);
    tail_.resize(ncols, -1);
  }

  int addConflict(const HighsDomainChange* conflictchgs, int nchgs);

  void removeConflict(int confid);

  bool conflictPropagates(int state) { return state < 2; }

  template <typename MarkPropagate>
  void updateStates(const HighsDomainChange& domchg, double oldbound,
                    std::vector<int>& states, MarkPropagate&& markPropagate) {
    oldchg.boundval = oldbound;
    int it = head_[domchg.column];

    while (it != -1) {
      conflictchange_[it];

      int oldstate = states[conflictindex_[it]];

      if (domchg.boundtype == HighsBoundType::Lower) {
        if (conflictchange_[it].boundtype == HighsBoundType::Lower) {
          if (conflictchange_[it].boundval > oldbound + 1e-6 &&
              conflictchange_[it].boundval <= domchg.boundval + 1e-6)
            states[conflictindex_[it]] += 1;
        } else {
          if (conflictchange_[it].boundval <= oldbound + 1e-6 &&
              conflictchange_[it].boundval > domchg.boundval + 1e-6)
            states[conflictindex_[it]] -= 1;
        }

      } else {
        if (conflictchange_[it].boundtype == HighsBoundType::Lower) {
          if (conflictchange_[it].boundval <= oldbound + 1e-6 &&
              conflictchange_[it].boundval > domchg.boundval + 1e-6)
            states[conflictindex_[it]] -= 1;
        } else {
          if (conflictchange_[it].boundval > oldbound + 1e-6 &&
              conflictchange_[it].boundval <= domchg.boundval + 1e-6)
            states[conflictindex_[it]] += 1;
        }
      }

      if (states[conflictindex_[it]] < 2) markPropagate(conflictindex_[it]);

      it = next_[it];
    }
  }

  int propagateConflict(int conflict, const HighsDomain& domain,
                        std::vector<HighsDomainChange>& boundchgs,
                        bool& infeasible);
};

#endif