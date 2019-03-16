#include "mip/SolveMip.h"

#include <cmath>

// For the moment just return first violated.
NodeIndex NodeStack::chooseBranchingVariable(const Node& node) {
  if (node.integer_variables.size() == 0)
    return kNodeIndexError;
  
  assert(node.integer_variables.size() == node.primal_solution.size());

  for (int col=0; col<node.integer_variables.size(); col++) {
    if (!node.integer_variables[col])
      continue;

    double value = node.primal_solution[col];
    if (std::fabs(value - std::floor(value)) > 0.0000001) {
				// This one is violated.
        return NodeIndex(col);
			}
  }

  return kNoNodeIndex;
}

bool NodeStack::branch(Node& node) {
  NodeIndex branch_col = chooseBranchingVariable(node);
  if (branch_col == kNodeIndexError ||
      branch_col == kNoNodeIndex)
    return false;

  // Pop current node from stack before adding children.
  pop();

  // Branch.
  // Create children and add to node.

  // Add to stack.

  return true;
}