/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsNodeQueue.h"

#include <algorithm>
#include <tuple>

#include "lp_data/HConst.h"
#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsSplay.h"

#define ESTIMATE_WEIGHT .5
#define LOWERBOUND_WEIGHT .5

void HighsNodeQueue::link_estim(HighsInt node) {
  auto get_left = [&](HighsInt n) -> HighsInt& {
    return nodes[n].leftestimate;
  };
  auto get_right = [&](HighsInt n) -> HighsInt& {
    return nodes[n].rightestimate;
  };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lower_bound +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  assert(node != -1);

  highs_splay_link(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_estim(HighsInt node) {
  auto get_left = [&](HighsInt n) -> HighsInt& {
    return nodes[n].leftestimate;
  };
  auto get_right = [&](HighsInt n) -> HighsInt& {
    return nodes[n].rightestimate;
  };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lower_bound +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  assert(estimroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::link_lower(HighsInt node) {
  auto get_left = [&](HighsInt n) -> HighsInt& { return nodes[n].leftlower; };
  auto get_right = [&](HighsInt n) -> HighsInt& { return nodes[n].rightlower; };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  assert(node != -1);

  highs_splay_link(node, lowerroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_lower(HighsInt node) {
  auto get_left = [&](HighsInt n) -> HighsInt& { return nodes[n].leftlower; };
  auto get_right = [&](HighsInt n) -> HighsInt& { return nodes[n].rightlower; };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  assert(lowerroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, lowerroot, get_left, get_right, get_key);
}

void HighsNodeQueue::link_domchgs(HighsInt node) {
  HighsInt numchgs = nodes[node].domchgstack.size();
  nodes[node].domchglinks.resize(numchgs);

  for (HighsInt i = 0; i != numchgs; ++i) {
    double val = nodes[node].domchgstack[i].boundval;
    HighsInt col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::kLower:
        nodes[node].domchglinks[i] =
            colLowerNodes[col].emplace(val, node).first;
        break;
      case HighsBoundType::kUpper:
        nodes[node].domchglinks[i] =
            colUpperNodes[col].emplace(val, node).first;
    }
  }
}

void HighsNodeQueue::unlink_domchgs(HighsInt node) {
  HighsInt numchgs = nodes[node].domchgstack.size();

  for (HighsInt i = 0; i != numchgs; ++i) {
    HighsInt col = nodes[node].domchgstack[i].column;
    switch (nodes[node].domchgstack[i].boundtype) {
      case HighsBoundType::kLower:
        colLowerNodes[col].erase(nodes[node].domchglinks[i]);
        break;
      case HighsBoundType::kUpper:
        colUpperNodes[col].erase(nodes[node].domchglinks[i]);
    }
  }

  nodes[node].domchglinks.clear();
  nodes[node].domchglinks.shrink_to_fit();
}

void HighsNodeQueue::link(HighsInt node) {
  link_estim(node);
  link_lower(node);
  link_domchgs(node);
}

void HighsNodeQueue::unlink(HighsInt node) {
  unlink_estim(node);
  unlink_lower(node);
  unlink_domchgs(node);
  freeslots.push(node);
}

void HighsNodeQueue::setNumCol(HighsInt numcol) {
  colLowerNodes.resize(numcol);
  colUpperNodes.resize(numcol);
}

void HighsNodeQueue::checkGlobalBounds(HighsInt col, double lb, double ub,
                                       double feastol,
                                       HighsCDouble& treeweight) {
  std::set<HighsInt> delnodes;
  auto prunestart =
      colLowerNodes[col].lower_bound(std::make_pair(ub + feastol, -1));
  for (auto it = prunestart; it != colLowerNodes[col].end(); ++it)
    delnodes.insert(it->second);

  auto pruneend =
      colUpperNodes[col].upper_bound(std::make_pair(lb - feastol, kHighsIInf));
  for (auto it = colUpperNodes[col].begin(); it != pruneend; ++it)
    delnodes.insert(it->second);

  for (HighsInt delnode : delnodes) {
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

    assert(colLowerNodes.size() == globaldomain.col_lower_.size());
    HighsInt numcol = colLowerNodes.size();
    for (HighsInt i = 0; i != numcol; ++i) {
      checkGlobalBounds(i, globaldomain.col_lower_[i],
                        globaldomain.col_upper_[i], feastol, treeweight);
    }

    size_t numopennodes = numNodes();
    if (numopennodes == 0) break;

    for (HighsInt i = 0; i != numcol; ++i) {
      if (colLowerNodes[i].size() == numopennodes) {
        double globallb = colLowerNodes[i].begin()->first;
        if (globallb > globaldomain.col_lower_[i]) {
          globaldomain.changeBound(HighsBoundType::kLower, i, globallb,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }

      if (colUpperNodes[i].size() == numopennodes) {
        double globalub = colUpperNodes[i].rbegin()->first;
        if (globalub < globaldomain.col_upper_[i]) {
          globaldomain.changeBound(HighsBoundType::kUpper, i, globalub,
                                   HighsDomain::Reason::unspecified());
          if (globaldomain.infeasible()) break;
        }
      }
    }

    globaldomain.propagate();
  } while (numchgs != globaldomain.getDomainChangeStack().size());

  return double(treeweight);
}

double HighsNodeQueue::pruneNode(HighsInt nodeId) {
  double treeweight = std::pow(0.5, nodes[nodeId].depth - 1);
  unlink(nodeId);
  return treeweight;
}

double HighsNodeQueue::performBounding(double upper_limit) {
  if (lowerroot == -1) return 0.0;

  HighsCDouble treeweight = 0.0;

  auto get_left = [&](HighsInt n) -> HighsInt& { return nodes[n].leftlower; };
  auto get_right = [&](HighsInt n) -> HighsInt& { return nodes[n].rightlower; };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  // split the lower bound tree along the bounding value
  lowerroot = highs_splay(std::make_tuple(upper_limit, -kHighsInf, 0),
                          lowerroot, get_left, get_right, get_key);
  HighsInt delroot;

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
    std::vector<HighsInt> stack;
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
      nodes[delroot].branchings.clear();
      nodes[delroot].domchgstack.shrink_to_fit();
      nodes[delroot].branchings.shrink_to_fit();

      // remember the free position for reuse
      freeslots.push(delroot);
    }
  }

  return double(treeweight);
}

void HighsNodeQueue::emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                                 std::vector<HighsInt>&& branchPositions,
                                 double lower_bound, double estimate,
                                 HighsInt depth) {
  HighsInt pos;

  if (freeslots.empty()) {
    pos = nodes.size();
    nodes.emplace_back(std::move(domchgs), std::move(branchPositions),
                       lower_bound, estimate, depth);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    nodes[pos] = OpenNode(std::move(domchgs), std::move(branchPositions),
                          lower_bound, estimate, depth);
  }

  assert(nodes[pos].lower_bound == lower_bound);
  assert(nodes[pos].estimate == estimate);
  assert(nodes[pos].depth == depth);

  link(pos);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestNode() {
  auto get_left = [&](HighsInt n) -> HighsInt& {
    return nodes[n].leftestimate;
  };
  auto get_right = [&](HighsInt n) -> HighsInt& {
    return nodes[n].rightestimate;
  };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(LOWERBOUND_WEIGHT * nodes[n].lower_bound +
                               ESTIMATE_WEIGHT * nodes[n].estimate,
                           -int(nodes[n].domchgstack.size()), n);
  };

  estimroot = highs_splay(std::make_tuple(-kHighsInf, -kHighsIInf, 0),
                          estimroot, get_left, get_right, get_key);
  HighsInt bestestimnode = estimroot;

  unlink(bestestimnode);

  return std::move(nodes[bestestimnode]);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestBoundNode() {
  auto get_left = [&](HighsInt n) -> HighsInt& { return nodes[n].leftlower; };
  auto get_right = [&](HighsInt n) -> HighsInt& { return nodes[n].rightlower; };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  lowerroot = highs_splay(std::make_tuple(-kHighsInf, -kHighsInf, 0), lowerroot,
                          get_left, get_right, get_key);
  HighsInt bestboundnode = lowerroot;

  unlink(bestboundnode);

  return std::move(nodes[bestboundnode]);
}

double HighsNodeQueue::getBestLowerBound() {
  if (lowerroot == -1) return kHighsInf;

  auto get_left = [&](HighsInt n) -> HighsInt& { return nodes[n].leftlower; };
  auto get_right = [&](HighsInt n) -> HighsInt& { return nodes[n].rightlower; };
  auto get_key = [&](HighsInt n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  lowerroot = highs_splay(std::make_tuple(-kHighsInf, -kHighsInf, 0), lowerroot,
                          get_left, get_right, get_key);
  return nodes[lowerroot].lower_bound;
}
