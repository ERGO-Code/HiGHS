#include "mip/HighsNodeQueue.h"

#include <algorithm>
#include <tuple>

#include "lp_data/HConst.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsSplay.h"

#define ESTIMATE_WEIGHT (1. / 8.)
#define LOWERBOUND_WEIGHT (7. / 8.)

void HighsNodeQueue::link_estim(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lower_bound +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  assert(node != -1);

  highs_splay_link(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_estim(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lower_bound +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  assert(estimroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::link_lower(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].lp_objective, n);
  };

  assert(node != -1);

  highs_splay_link(node, lowerroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_lower(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].lp_objective, n);
  };

  assert(lowerroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, lowerroot, get_left, get_right, get_key);
}

void HighsNodeQueue::link_domchgs(int node) {
  int numchgs = nodes[node].domchgstack.size();
  nodes[node].domchglinks.resize(numchgs);

  for (int i = 0; i != numchgs; ++i) {
    double val = nodes[node].domchgstack[i].boundval;
    int col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::Lower:
        nodes[node].domchglinks[i] = colLowerNodes[col].emplace(val, node);
        break;
      case HighsBoundType::Upper:
        nodes[node].domchglinks[i] = colUpperNodes[col].emplace(val, node);
    }
  }
}

void HighsNodeQueue::unlink_domchgs(int node) {
  int numchgs = nodes[node].domchgstack.size();

  for (int i = 0; i != numchgs; ++i) {
    int col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::Lower:
        colLowerNodes[col].erase(nodes[node].domchglinks[i]);
        break;
      case HighsBoundType::Upper:
        colUpperNodes[col].erase(nodes[node].domchglinks[i]);
    }
  }

  nodes[node].domchglinks.clear();
  nodes[node].domchglinks.shrink_to_fit();
}

void HighsNodeQueue::link(int node) {
  link_estim(node);
  link_lower(node);
  link_domchgs(node);
}

void HighsNodeQueue::unlink(int node) {
  unlink_estim(node);
  unlink_lower(node);
  unlink_domchgs(node);
  freeslots.push(node);
}

void HighsNodeQueue::setNumCol(int numcol) {
  colLowerNodes.resize(numcol);
  colUpperNodes.resize(numcol);
}

void HighsNodeQueue::checkGlobalBounds(int col, double lb, double ub,
                                       double feastol,
                                       HighsCDouble& treeweight) {
  std::set<int> delnodes;
  auto prunestart = colLowerNodes[col].lower_bound(ub + feastol);
  for (auto it = prunestart; it != colLowerNodes[col].end(); ++it)
    delnodes.insert(it->second);

  auto pruneend = colUpperNodes[col].upper_bound(lb - feastol);
  for (auto it = colUpperNodes[col].begin(); it != pruneend; ++it)
    delnodes.insert(it->second);

  for (int delnode : delnodes) {
    treeweight += std::pow(0.5, nodes[delnode].depth - 1);
    unlink(delnode);
  }
}

double HighsNodeQueue::pruneInfeasibleNodes(HighsDomain& globaldomain,
                                            double feastol) {
  size_t numchgs;

  HighsCDouble treeweight = 0.0;

  do {
    if (globaldomain.infeasible()) break;

    numchgs = globaldomain.getDomainChangeStack().size();

    assert(colLowerNodes.size() == globaldomain.colLower_.size());
    int numcol = colLowerNodes.size();
    for (int i = 0; i != numcol; ++i) {
      checkGlobalBounds(i, globaldomain.colLower_[i], globaldomain.colUpper_[i],
                        feastol, treeweight);
    }

    size_t numopennodes = numNodes();
    if (numopennodes == 0) break;

    for (int i = 0; i != numcol; ++i) {
      if (colLowerNodes[i].size() == numopennodes) {
        double globallb = colLowerNodes[i].begin()->first;
        if (globallb > globaldomain.colLower_[i]) {
          globaldomain.changeBound(HighsBoundType::Lower, i, globallb,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }

      if (colUpperNodes[i].size() == numopennodes) {
        double globalub = colUpperNodes[i].rbegin()->first;
        if (globalub < globaldomain.colUpper_[i]) {
          globaldomain.changeBound(HighsBoundType::Upper, i, globalub,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }
    }

    globaldomain.propagate();
  } while (numchgs != globaldomain.getDomainChangeStack().size());

  return double(treeweight);
}

double HighsNodeQueue::performBounding(double upper_limit) {
  if (lowerroot == -1) return 0.0;

  HighsCDouble treeweight = 0.0;

  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) { return std::make_pair(nodes[n].lower_bound, n); };

  // split the lower bound tree along the bounding value
  lowerroot = highs_splay(std::make_pair(upper_limit, 0), lowerroot, get_left,
                          get_right, get_key);
  int delroot;

  if (nodes[lowerroot].lower_bound < upper_limit) {
    delroot = get_right(lowerroot);
    get_right(lowerroot) = -1;
  } else {
    delroot = lowerroot;
    lowerroot = get_left(delroot);
    get_left(delroot) = -1;
  }

  // now unlink all removed nodes from the estimate tree
  if (delroot != -1) {
    std::vector<int> stack;
    stack.reserve(numNodes());
    stack.push_back(delroot);

    while (!stack.empty()) {
      assert(nodes[delroot].lower_bound >= upper_limit);

      delroot = stack.back();
      stack.pop_back();

      // unlink the node from the best estimate tree
      // and add up the tree weight
      unlink_estim(delroot);
      unlink_domchgs(delroot);
      treeweight += std::pow(0.5, nodes[delroot].depth - 1);

      // put the nodes children on the stack for subsequent processing
      if (get_left(delroot) != -1) stack.push_back(get_left(delroot));
      if (get_right(delroot) != -1) stack.push_back(get_right(delroot));

      // release the memory for domain changes
      nodes[delroot].domchgstack.clear();
      nodes[delroot].domchgstack.shrink_to_fit();

      // remember the free position for reuse
      freeslots.push(delroot);
    }
  }

  return double(treeweight);
}

void HighsNodeQueue::emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                                 double lower_bound, double lp_objective,
                                 double estimate, int depth) {
  int pos;

  if (freeslots.empty()) {
    pos = nodes.size();
    nodes.emplace_back(std::move(domchgs), lower_bound, lp_objective, estimate,
                       depth);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    nodes[pos] = OpenNode(std::move(domchgs), lower_bound, lp_objective,
                          estimate, depth);
  }

  assert(nodes[pos].lower_bound == lower_bound);
  assert(nodes[pos].lp_objective == lp_objective);
  assert(nodes[pos].estimate == estimate);
  assert(nodes[pos].depth == depth);

  link(pos);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestNode() {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lp_objective +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  estimroot =
      highs_splay(std::make_tuple(-HIGHS_CONST_INF, -HIGHS_CONST_I_INF, 0),
                  estimroot, get_left, get_right, get_key);
  int bestestimnode = estimroot;

  unlink(bestestimnode);

  return std::move(nodes[bestestimnode]);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestBoundNode() {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].lp_objective, n);
  };

  lowerroot =
      highs_splay(std::make_tuple(-HIGHS_CONST_INF, -HIGHS_CONST_INF, 0),
                  lowerroot, get_left, get_right, get_key);
  int bestboundnode = lowerroot;

  unlink(bestboundnode);

  return std::move(nodes[bestboundnode]);
}

double HighsNodeQueue::getBestLowerBound() {
  if (lowerroot == -1) return HIGHS_CONST_INF;

  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].lp_objective, n);
  };

  lowerroot =
      highs_splay(std::make_tuple(-HIGHS_CONST_INF, -HIGHS_CONST_INF, 0),
                  lowerroot, get_left, get_right, get_key);
  return nodes[lowerroot].lower_bound;
}