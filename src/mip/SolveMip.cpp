#include "mip/SolveMip.h"

NodeIndex NodeStack::chooseBranchingVariable(const Node& node,
                                             const std::vector<int>& integer_variables) {
  if (integer_variables.size() == 0)
    return kNodeIndexError;
  
  return kNoNodeIndex;
}