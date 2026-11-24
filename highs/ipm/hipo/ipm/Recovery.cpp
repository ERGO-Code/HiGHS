#include "Recovery.h"

namespace hipo {

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

bool Recovery::available(RecoveryAction action,
                         const SolverCapabilities& capabilities) const {
  return allowedBy(action, capabilities) &&
         timesProposed(action) < actionProperties(action).max_times;
}

RecoveryAction Recovery::firstAvailable(
    RecoveryAction entry, const SolverCapabilities& capabilities) const {
  for (RecoveryAction action = entry; action != RecoveryAction::kNone;
       action = actionProperties(action).next) {
    if (available(action, capabilities)) return action;
  }
  return RecoveryAction::kNone;
}

// I or W: the linear solver is in trouble
static bool hasNumericalTrouble(const IterationReport& report) {
  return report.has(IssueType::kVeryInaccurateSolve) ||
         report.has(IssueType::kWrongInertia);
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

  if (report.clean())
    ++streak_of_clean_iter;
  else
    streak_of_clean_iter = 0;

  if (waitingForPreviousAction(iter))
    return {RecoveryAction::kNone, iter, false};

  // Issues in decreasing order of priority; the first one for which an action
  // can be chosen is acted upon.
  RecoveryAction action = RecoveryAction::kNone;

  // 1, 2: stagnation. Recentring helps only if the point is poorly centred,
  // and L is too urgent for it. Numerical trouble with a short step counts as
  // stagnation, without waiting for N.
  const bool short_step_trouble =
      recent_numerical_trouble && report.alpha_min < kShortStepAlpha;
  if (report.stagnation() || short_step_trouble) {
    const bool try_recentre = report.has(IssueType::kPoorCentrality) &&
                              !report.has(IssueType::kLargeResIncrease);
    RecoveryAction entry = RecoveryAction::kSwitchFactorisation;
    if (recent_numerical_trouble)
      entry = RecoveryAction::kDropBadPivots;
    else if (try_recentre)
      entry = RecoveryAction::kRecentre;
    action = firstAvailable(entry, capabilities);

    // try reg+ before stop, e.g. for stagnation without numerical trouble when
    // switch is unavailable or used up
    if (action == RecoveryAction::kTerminate &&
        available(RecoveryAction::kIncreaseReg, capabilities))
      action = RecoveryAction::kIncreaseReg;

    // a short step alone is not enough to stop
    if (action == RecoveryAction::kTerminate && !report.stagnation())
      action = RecoveryAction::kNone;
  } else if (hasNumericalTrouble(report)) {
    // 3: numerical trouble without stagnation. The ipm is still progressing,
    // and changing the regularisation can stop it, but bad pivots can be
    // dropped.
    if (available(RecoveryAction::kDropBadPivots, capabilities))
      action = RecoveryAction::kDropBadPivots;
  }

  // 4: regularisation was increased and iterations have been clean for some
  // time, so decrease regularisation
  if (action == RecoveryAction::kNone &&
      timesProposed(RecoveryAction::kIncreaseReg) > 0 &&
      streak_of_clean_iter >= kCleanIterBeforeRegRestore) {
    if (available(RecoveryAction::kRestoreReg, capabilities))
      action = RecoveryAction::kRestoreReg;
  }

  // 5: absorbed rows
  if (action == RecoveryAction::kNone && report.has(IssueType::kAbsorbedRows))
    action = firstAvailable(RecoveryAction::kDecreaseReg, capabilities);

  const ProposedAction proposed{action, iter,
                                actionProperties(action).reset_to_best};
  if (action != RecoveryAction::kNone) proposed_actions.push_back(proposed);

  if (action == RecoveryAction::kIncreaseReg) streak_of_clean_iter = 0;

  return proposed;
}

}  // namespace hipo
