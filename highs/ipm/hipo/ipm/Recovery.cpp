#include "Recovery.h"

#include <algorithm>

#include "Parameters.h"

namespace hipo {

RecoveryAction Recovery::nextAction(const IterationReport& report) const {
  if (!report.stagnation()) return RecoveryAction::kNone;

  // For now, the same sequence of actions is used for any issue: switch to
  // FactorHighs once, then terminate.
  if (!tried(RecoveryAction::kSwitchFactorisation))
    return RecoveryAction::kSwitchFactorisation;

  return RecoveryAction::kTerminate;
}

void Recovery::applied(RecoveryAction action, const IterationReport& report) {
  records.push_back({action, report.mainStagnationIssue(), report.iter});
}

bool Recovery::tried(RecoveryAction action) const {
  for (const Record& record : records)
    if (record.action == action) return true;
  return false;
}

bool Recovery::tried(RecoveryAction action, IssueType issue) const {
  for (const Record& record : records)
    if (record.action == action && record.issue == issue) return true;
  return false;
}

// A rung of a ladder: its action, and how many times it can be chosen before
// moving to the next rung.
struct Rung {
  RecoveryAction action;
  Int max_uses;
};

// Ladders, ended by kNone.
static const Int kLadderSize = 4;
static const Rung kLadders[][kLadderSize] = {
    // numerical
    {{RecoveryAction::kMoreRefinement, 1},
     {RecoveryAction::kIncreaseReg, kMaxRegChanges},
     {RecoveryAction::kSwitchFactorisation, 1},
     {RecoveryAction::kNone, 0}},
    // stagnation
    {{RecoveryAction::kRecentre, 1},
     {RecoveryAction::kSwitchFactorisation, 1},
     {RecoveryAction::kTerminate, 1},
     {RecoveryAction::kNone, 0}},
    // absorbed rows
    {{RecoveryAction::kDecreaseReg, kMaxRegChanges},
     {RecoveryAction::kNone, 0}},
};

bool Recovery::proposed(RecoveryAction action) const {
  for (const Proposal& proposal : proposals)
    if (proposal.action == action) return true;
  return false;
}

RecoveryAction Recovery::climb(Ladder ladder, Int entry, bool can_switch) {
  Int& rung = next_rung[ladder];
  Int& uses = rung_uses[ladder];
  if (rung < entry) {
    rung = entry;
    uses = 0;
  }
  while (rung < kLadderSize &&
         kLadders[ladder][rung].action != RecoveryAction::kNone) {
    const Rung& current = kLadders[ladder][rung];
    const RecoveryAction action = current.action;
    const bool available =
        uses < current.max_uses &&
        (action != RecoveryAction::kSwitchFactorisation || can_switch) &&
        (actionProperties(action).temporary || !proposed(action));
    if (available) {
      ++uses;
      return action;
    }
    ++rung;
    uses = 0;
  }
  return RecoveryAction::kNone;
}

RecoveryAction Recovery::propose(const IterationReport& report,
                                 bool can_switch) {
  const Int iter = report.iter;

  const bool numerical_now = report.has(IssueType::kVeryInaccurateSolve) ||
                             report.has(IssueType::kWrongInertia) ||
                             report.has(IssueType::kDynamicReg);
  if (numerical_now) {
    numerical_seen = true;
    last_numerical_iter = iter;
  }
  const bool recent_numerical =
      numerical_seen && iter - last_numerical_iter < kActionNumericalWindow;

  // wait for the effect of the previous action
  if (!proposals.empty() && iter - proposals.back().iter < kActionPatience)
    return RecoveryAction::kNone;

  // I and W make the direction unreliable, so refinement is skipped
  const Int numerical_entry = report.has(IssueType::kVeryInaccurateSolve) ||
                                      report.has(IssueType::kWrongInertia)
                                  ? 1
                                  : 0;

  // recentring helps only if the point is poorly centred; L is urgent
  const Int stagnation_entry = report.has(IssueType::kPoorCentrality) &&
                                       !report.has(IssueType::kLargeResIncrease)
                                   ? 0
                                   : 1;

  // issues in decreasing order of priority; the first one for which an action
  // is available is acted upon
  RecoveryAction action = RecoveryAction::kNone;
  if (report.stagnation() && recent_numerical)
    action = climb(kNumericalLadder, numerical_entry, can_switch);
  if (action == RecoveryAction::kNone && report.stagnation())
    action = climb(kStagnationLadder, stagnation_entry, can_switch);
  if (action == RecoveryAction::kNone && numerical_now)
    action = climb(kNumericalLadder, numerical_entry, can_switch);
  if (action == RecoveryAction::kNone && report.has(IssueType::kAbsorbedRows))
    action = climb(kAbsorbedLadder, 0, can_switch);

  if (action != RecoveryAction::kNone) proposals.push_back({action, iter});
  return action;
}

}  // namespace hipo
