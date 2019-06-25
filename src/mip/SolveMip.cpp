#include "mip/SolveMip.h"

#include <cmath>

#include "io/HighsIO.h"

// For the moment just return first violated.
NodeIndex Tree::chooseBranchingVariable(const Node& node) {
  if (node.integer_variables.size() == 0) return kNodeIndexError;

  assert(node.integer_variables.size() == node.primal_solution.size());

  for (int col = 0; col < node.integer_variables.size(); col++) {
    if (!node.integer_variables[col]) continue;

    double value = node.primal_solution[col];
    if (std::fabs(value - std::floor(value)) > 0.0000001) {
      // This one is violated.
      return NodeIndex(col);
    }
  }

  return kNoNodeIndex;
}

struct testn {
  testn* left;
  testn* right;
  testn(int a, int b) : a_(a) {}
  int a_;
};

bool Tree::branch(Node& node) {
  NodeIndex branch_col = chooseBranchingVariable(node);
  if (branch_col == kNodeIndexError) return false;

  if (branch_col == kNoNodeIndex) {
    // All integer variables are feasible. Update best solution.
    // Assuming minimization.
    HighsPrintMessage(ML_ALWAYS, "Updating best solution at node %d.\n",
                      node.id);

    if (node.objective_value < best_objective_) {
      best_objective_ = node.objective_value;
      best_solution_ = node.primal_solution;
    }
    return false;
  }

  int col = static_cast<int>(branch_col);
  double value = node.primal_solution[col];

  // todo: change always to whatever.
  HighsPrintMessage(ML_ALWAYS, "Branching on variable %d\n", col);

  HighsPrintMessage(
      ML_ALWAYS, "%d(%d, %d) left child: %.2f, right child: %.2f\n", node.id,
      num_nodes + 1, num_nodes + 2, std::floor(value), std::ceil(value));

  // Branch.
  // Create children and add to node.
  num_nodes++;
  node.left_child =
      std::unique_ptr<Node>(new Node(node.id, num_nodes, node.level + 1));
  num_nodes++;
  node.right_child =
      std::unique_ptr<Node>(new Node(node.id, num_nodes, node.level + 1));

  // Copy bounds from parent.
  node.left_child->col_lower_bound = node.col_lower_bound;
  node.left_child->col_upper_bound = node.col_upper_bound;
  node.left_child->integer_variables = node.integer_variables;

  node.right_child->col_lower_bound = node.col_lower_bound;
  node.right_child->col_upper_bound = node.col_upper_bound;
  node.right_child->integer_variables = node.integer_variables;

  node.left_child->col_upper_bound[col] = std::floor(value);
  node.right_child->col_lower_bound[col] = std::ceil(value);

  // Add to stack.
  std::reference_wrapper<Node> left(*(node.left_child).get());
  std::reference_wrapper<Node> right(*(node.right_child).get());
  nodes_.push_back(left);
  nodes_.push_back(right);

  return true;
}