#include "mip/HighsConflictPool.h"

#include "mip/HighsDomain.h"

int HighsConflictPool::addConflict(const HighsDomainChange* conflictchgs,
                                   int nchgs) {
  auto it = freespaces_.lower_bound(std::make_pair(nchgs, 0));

  int confid;

  if (deletedindices_.empty()) {
    confid = conflictrange_.size();
    conflictrange_.resize(confid + 1);
  } else {
    confid = deletedindices_.back();
    deletedindices_.pop_back();
  }

  int& start = conflictrange_[confid].first;
  int& end = conflictrange_[confid].second;
  if (it != freespaces_.end()) {
    start = it->second;
    end = start + nchgs;

    if (nchgs < it->first)
      freespaces_.emplace(it->first - nchgs, start + nchgs);
    else
      freespaces_.erase(it);

    for (int i = 0; i != nchgs; ++i) {
      conflictchange_[start + i] = conflictchgs[i];
      conflictindex_[start + i] = confid;
    }
  } else {
    start = conflictchange_.size();
    end = start + nchgs;

    next_.resize(end, -1);
    prev_.resize(end, -1);
    conflictchange_.insert(conflictchange_.end(), conflictchgs,
                           conflictchgs + nchgs);

    conflictindex_.insert(conflictindex_.end(), confid);
  }
}

void HighsConflictPool::removeConflict(int confid) {
  int start = conflictrange_[confid].first;
  int end = conflictrange_[confid].second;

  for (int i = start; i != end; ++i)
    if (prev_[i] != -1) unlink(i);

  deletedindices_.push_back(confid);
  freespaces_.emplace(end - start, start);
}

int HighsConflictPool::getStateContribution(int confchgidx,
                                            const HighsDomain& domain) const {
  int state;

  if (conflictchange_[confchgidx].boundtype == HighsBoundType::Lower) {
    if (conflictchange_[confchgidx].boundval <
        domain.colLower_[conflictchange_[confchgidx].column] + 1e-6)
      state = 2;
    else if (conflictchange_[confchgidx].boundval >
             domain.colUpper_[conflictchange_[confchgidx].column] + 1e-6)
      state = 0;
    else
      state = 1;
  } else {
    if (conflictchange_[confchgidx].boundval <
        domain.colUpper_[conflictchange_[confchgidx].column] + 1e-6)
      state = 2;
    else if (conflictchange_[confchgidx].boundval >
             domain.colUpper_[conflictchange_[confchgidx].column] + 1e-6)
      state = 0;
    else
      state = 1;
  }

  return state;
}

void HighsConflictPool::unlink(int confchgidx) {
  int next = next_[confchgidx];
  int prev = prev_[confchgidx];

  if (next != -1)
    prev_[next] = prev;
  else if (tail_[conflictchange_[confchgidx].column] == confchgidx)
    tail_[conflictchange_[confchgidx].column] = prev;

  if (prev != confchgidx)
    next_[prev] = next;
  else if (head_[conflictchange_[confchgidx].column] == confchgidx) {
    head_[conflictchange_[confchgidx].column] = next;
    prev_[next] = next;
  }
}

void HighsConflictPool::link(int confchgidx) {
  int col = conflictchange_[confchgidx].column;
  if (tail_[col] != -1) {
    prev_[confchgidx] = tail_[col];
    next_[tail_[col]] = confchgidx;
    tail_[col] = confchgidx;
  } else {
    head_[col] = confchgidx;
    tail_[col] = confchgidx;
    // we link the prev pointer to itself for the head
    // so that we can check whether a conflict change is linked by checking
    // prev_[confchgidx] != -1;
    prev_[confchgidx] = confchgidx;
  }
}

int HighsConflictPool::propagateConflict(
    int conflict, const HighsDomain& domain,
    std::vector<HighsDomainChange>& boundchgs, bool& infeasible) {
  int start = conflictrange_[conflict].first;
  int end = conflictrange_[conflict].second;

  int watchedliterals[2] = {-1, -1};

  int nwatchedliterals = 0;
  for (int i = start; i != end; ++i) {
    if (next_[i] != -1 || prev_[i] != -1)
      watchedliterals[nwatchedliterals++] = i;
  }

  // try to find 2 new literals or replacedments for the current watched
  // literals if those are fixed to zero
  nwatchedliterals = 2;
  int statesum = 0;
  for (int i = 0; i != 2; ++i) {
    int k = watchedliterals[i];

    int state = k == -1 ? 0 : getStateContribution(k, domain);

    if (state == 0) {
      int newliteral = -1;

      if (k != -1) unlink(k);

      // we need to find a new watched literal that has state != 0
      for (int j = start; j != end; ++j) {
        if (j == watchedliterals[0] || j == watchedliterals[1]) continue;

        int statej = getStateContribution(j, domain);
        if (statej != 0) {
          newliteral = j;
          state = statej;
          break;
        }
      }

      watchedliterals[i] = newliteral;

      if (watchedliterals[i] == -1)
        --nwatchedliterals;
      else
        link(newliteral);
    }

    statesum += state;
  }

  // if we found 0 literals to watch, the conflict is infeasible
  // if we found 1 literal to watch, it may be fixed to one if it is not
  // already and if we found 2 literals to watch we can fix nothing for now

  switch (nwatchedliterals) {
    case 0:
      infeasible = true;
      break;
    case 1: {
      int literal =
          watchedliterals[0] != -1 ? watchedliterals[0] : watchedliterals[1];

      if ((conflictchange_[literal].boundtype == HighsBoundType::Lower &&
           conflictchange_[literal].boundval >
               domain.colLower_[conflictchange_[literal].column]) ||
          (conflictchange_[literal].boundtype == HighsBoundType::Upper &&
           conflictchange_[literal].boundval <
               domain.colUpper_[conflictchange_[literal].column]))
        boundchgs.push_back(conflictchange_[literal]);
    }
  }

  return statesum;
}