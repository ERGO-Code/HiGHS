/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/SolveMip.h"

#include <cmath>

#include "io/HighsIO.h"

// For the moment just return first violated.
NodeIndex Tree::chooseBranchingVariable(const Node& node) {
  assert(node.integer_variables.size() == node.primal_solution.size());

  for (int col = 0; col < (int)node.integer_variables.size(); col++) {
    if (!node.integer_variables[col]) continue;

    double value = node.primal_solution[col];
    if (std::fabs(value - std::floor(value)) > 0.0000001) {
      // This one is violated.
      return NodeIndex(col);
    }
  }

  return kNoNodeIndex;
}

bool Tree::branch(Node& node) {
  NodeIndex branch_col = chooseBranchingVariable(node);
  if (branch_col == kNodeIndexError) return false;

  if (branch_col == kNoNodeIndex) {
    // All integer variables are feasible. Update best solution if node solution
    // is better. Assuming minimization.

    if (node.objective_value < best_objective_) {
      std::cout << "Updating best solution at node " << node.id
                << ". Objective: " << node.objective_value << std::endl;
      best_objective_ = node.objective_value;
      best_solution_ = node.primal_solution;
    }
    return false;
  }

  int col = static_cast<int>(branch_col);
  double value = node.primal_solution[col];

  std::cout << "Branching on variable " << col << std::endl
            << "(" << num_nodes + 1 << "," << num_nodes + 2
            << ") left child ub: " << std::floor(value)
            << " right child lb: " << std::ceil(value) << std::endl;

  // Branch.
  // Create children and add to node.
  num_nodes++;
  node.left_child =
      std::unique_ptr<Node>(new Node(node.id, num_nodes, node.level + 1));
  num_nodes++;
  node.right_child =
      std::unique_ptr<Node>(new Node(node.id, num_nodes, node.level + 1));

  // Copy bounds from parent and set integer variables.
  node.left_child->branch_col = col;
  node.left_child->col_lower_bound = node.col_lower_bound;
  node.left_child->col_upper_bound = node.col_upper_bound;
  node.left_child->col_upper_bound[col] = std::floor(value);
  node.left_child->integer_variables = node.integer_variables;

  node.right_child->branch_col = col;
  node.right_child->col_lower_bound = node.col_lower_bound;
  node.right_child->col_upper_bound = node.col_upper_bound;
  node.right_child->col_lower_bound[col] = std::ceil(value);
  node.right_child->integer_variables = node.integer_variables;

  // Add to stack.
  std::reference_wrapper<Node> left(*(node.left_child).get());
  std::reference_wrapper<Node> right(*(node.right_child).get());
  nodes_.push_back(left);
  nodes_.push_back(right);

  return true;
}
