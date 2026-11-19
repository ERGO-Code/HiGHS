#include "IterationHistory.h"

namespace hipo {

void IterationHistory::append() { records.push_back({}); }
IterationRecord& IterationHistory::back() { return records.back(); }

}  // namespace hipo
