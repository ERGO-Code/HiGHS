/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_SOLVEMIP_H_
#define MIP_SOLVEMIP_H_

#include <functional>
#include <memory>
#include <stack>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

struct Node {
  int id;
  int parent_id;
  int level;

  Node();
  Node(int parent, int index, int depth)
      : id(index), parent_id(parent), level(depth) {
    left_child = nullptr;
    right_child = nullptr;
  }

  std::vector<int> integer_variables;
  std::vector<double> primal_solution;
  double objective_value;

  // Minimal information about changes. Just bounds for the moment.
  std::vector<double> col_lower_bound;
  std::vector<double> col_upper_bound;

  std::unique_ptr<Node> left_child;
  std::unique_ptr<Node> right_child;
};

using NodeIndex = int;
constexpr NodeIndex kNoNodeIndex = -1;
constexpr NodeIndex kNodeIndexError = -2;

class Tree {
 public:
  Tree(Node& node) {
    std::reference_wrapper<Node> ref(node);
    nodes_.push_back(ref);
  }

  bool branch(FILE* output, const int message_level, Node& node);

  Node& next() { return nodes_[nodes_.size() - 1]; }
  void pop() { nodes_.erase(nodes_.end() - 1); }
  bool empty() { return (nodes_.size() == 0); }

  const std::vector<double>& getBestSolution() const { return best_solution_; }

  double getBestObjective() { return best_objective_; }

 private:
  std::vector<std::reference_wrapper<Node> > nodes_;
  std::vector<double> best_solution_;
  double best_objective_ = HIGHS_CONST_INF;

  NodeIndex chooseBranchingVariable(const Node& node);

  int num_nodes = 0;
};

#endif
