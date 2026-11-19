#include "Recovery.h"

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

}  // namespace hipo
