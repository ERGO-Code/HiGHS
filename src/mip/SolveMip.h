#ifndef MIP_SOLVEMIP_H_
#define MIP_SOLVEMIP_H_

#include <memory>
#include <stack>
#include <functional>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

struct Node {
  int id;
  int parent_id;
  int level;

  Node(int parent, int index, int depth) :
    id(index), parent_id(parent), level(depth) {
      left_child = nullptr;
      right_child = nullptr;
    }

  std::vector<int> integer_variables;
  std::vector<double> primal_solution;

  // Minimal information about changes. Just bounds for the moment.
  std::vector<double> col_lower_bound;
  std::vector<double> col_upper_bound;

  Node * left_child;
  Node * right_child;
};

using NodeIndex = int;
constexpr NodeIndex kNoNodeIndex = -1;
constexpr NodeIndex kNodeIndexError = -2;


class Tree {
public:
  bool branch(Node& node);

  Node& next() { return *(nodes_[nodes_.size() - 1]); }
  bool empty() {return (nodes_.size() == 0); }

private:
  std::vector<Node *> nodes_;
  std::vector<double> best_solution_;

  NodeIndex chooseBranchingVariable(const Node& node);
  
  int num_nodes = 0;
};

#endif