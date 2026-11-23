#include "Recovery.h"

#include <algorithm>
#include <cassert>
#include <limits>

#include "Parameters.h"

namespace hipo {

// A rung: an action, and how many times it can be chosen. The count is of the
// action, whichever ladder chose it: switch is on two ladders, and once used up
// on one it is used up on the other.
struct Rung {
  RecoveryAction action;
  Int max_times;
};

static const std::vector<Rung> kNumericalRungs = {
    {RecoveryAction::kDropBadPivots, 1},
    {RecoveryAction::kIncreaseReg, kMaxRegChanges},
    {RecoveryAction::kSwitchFactorisation, 1},
};

static const std::vector<Rung> kStagnationRungs = {
    {RecoveryAction::kRecentre, 1},
    {RecoveryAction::kSwitchFactorisation, 1},
    {RecoveryAction::kTerminate, 1},
};

// each row has a floor on its regularisation, so no limit is needed
static const std::vector<Rung> kAbsorbedRungs = {
    {RecoveryAction::kDecreaseReg, std::numeric_limits<Int>::max()},
};

static const std::vector<Rung>& ladderRungs(Recovery::Ladder ladder) {
  switch (ladder) {
    case Recovery::kNumericalLadder:
      return kNumericalRungs;
    case Recovery::kStagnationLadder:
      return kStagnationRungs;
    case Recovery::kAbsorbedLadder:
    default:
      return kAbsorbedRungs;
  }
}

static Int rungIndex(const std::vector<Rung>& rungs, RecoveryAction action) {
  Int index = 0;
  while (index < (Int)rungs.size() && rungs[index].action != action) ++index;
  assert(index < (Int)rungs.size());
  return index;
}

static bool allowedBy(RecoveryAction action,
                      const SolverCapabilities& capabilities) {
  if (action == RecoveryAction::kSwitchFactorisation)
    return capabilities.can_switch;
  if (action == RecoveryAction::kDropBadPivots)
    return capabilities.can_drop_pivots;
  return true;
}

Int Recovery::timesProposed(RecoveryAction action) const {
  Int count = 0;
  for (const ProposedAction& proposed : proposed_actions)
    if (proposed.action == action) ++count;
  return count;
}

RecoveryAction Recovery::climbLadder(Ladder ladder, RecoveryAction entry,
                                     const SolverCapabilities& capabilities) {
  const std::vector<Rung>& rungs = ladderRungs(ladder);

  // enter at entry, unless the ladder already got past it
  Int& rung = current_rung[ladder];
  rung = std::max(rung, rungIndex(rungs, entry));

  // move past the rungs used up or not allowed; the ladder stays on the rung
  // chosen, which can be chosen again while it is not used up
  for (; rung < (Int)rungs.size(); ++rung) {
    const Rung& r = rungs[rung];
    if (allowedBy(r.action, capabilities) &&
        timesProposed(r.action) < r.max_times)
      return r.action;
  }
  return RecoveryAction::kNone;
}

// I, W or G: the linear solver is in trouble
static bool hasNumericalTrouble(const IterationReport& report) {
  return report.has(IssueType::kVeryInaccurateSolve) ||
         report.has(IssueType::kWrongInertia) ||
         report.has(IssueType::kDynamicReg);
}

bool Recovery::waitingForPreviousAction(Int iter) const {
  return !proposed_actions.empty() &&
         iter - proposed_actions.back().iter < kActionPatience;
}

Recovery::ProposedAction Recovery::proposeAction(
    const IterationReport& report, const SolverCapabilities& capabilities) {
  const Int iter = report.iter;

  if (hasNumericalTrouble(report)) {
    numerical_trouble_seen = true;
    last_numerical_trouble_iter = iter;
  }
  const bool recent_numerical_trouble =
      numerical_trouble_seen &&
      iter - last_numerical_trouble_iter < kActionNumericalWindow;

  if (waitingForPreviousAction(iter))
    return {RecoveryAction::kNone, iter, false};

  // Where to enter the stagnation ladder: recentring helps only if the point is
  // poorly centred, and L is too urgent for it.
  const bool try_recentre = report.has(IssueType::kPoorCentrality) &&
                            !report.has(IssueType::kLargeResIncrease);
  const RecoveryAction stagnation_entry =
      try_recentre ? RecoveryAction::kRecentre
                   : RecoveryAction::kSwitchFactorisation;

  // Issues in decreasing order of priority; the first one for which an action
  // can be chosen is acted upon.
  RecoveryAction action = RecoveryAction::kNone;

  // 1, 2: stagnation
  if (report.stagnation()) {
    if (recent_numerical_trouble)
      action = climbLadder(kNumericalLadder, RecoveryAction::kDropBadPivots,
                           capabilities);
    if (action == RecoveryAction::kNone)
      action = climbLadder(kStagnationLadder, stagnation_entry, capabilities);
  }

  // Numerical trouble without stagnation is not acted upon: the ipm is still
  // progressing, and changing the regularisation can stop it.

  // 3: absorbed rows
  if (action == RecoveryAction::kNone && report.has(IssueType::kAbsorbedRows))
    action = climbLadder(kAbsorbedLadder, RecoveryAction::kDecreaseReg,
                         capabilities);

  const ProposedAction proposed{action, iter,
                                actionProperties(action).reset_to_best};
  if (action != RecoveryAction::kNone) proposed_actions.push_back(proposed);
  return proposed;
}

}  // namespace hipo
