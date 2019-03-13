#ifndef MIP_SOLVEMIP_H_
#define MIP_SOLVEMIP_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

class Node {
  // todo: add code from bnb
};

using NodeIndex = int;

class NodeStack {
public:
  NodeIndex chooseBranchingVariable();

private:
  std::vector<int> integer_variables;
  std::vector<double> best_solution;

  // todo: stack of nodes
};

#endif