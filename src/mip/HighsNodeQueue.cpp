/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsNodeQueue.h"

#include <algorithm>
#include <tuple>

#include "lp_data/HConst.h"
#include "util/HighsCDouble.h"
#include "util/HighsSplay.h"

void HighsNodeQueue::link_estim(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  assert(node != -1);

  highs_splay_link(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_estim(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  assert(estimroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, estimroot, get_left, get_right, get_key);
}

void HighsNodeQueue::link_lower(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) { return std::make_pair(nodes[n].lower_bound, n); };

  assert(node != -1);

  highs_splay_link(node, lowerroot, get_left, get_right, get_key);
}

void HighsNodeQueue::unlink_lower(int node) {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) { return std::make_pair(nodes[n].lower_bound, n); };

  assert(lowerroot != -1);
  assert(node != -1);

  highs_splay_unlink(node, lowerroot, get_left, get_right, get_key);
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
                                 double lower_bound, double estimate,
                                 int depth) {
  int pos;

  if (freeslots.empty()) {
    pos = nodes.size();
    nodes.emplace_back(std::move(domchgs), lower_bound, estimate, depth);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    nodes[pos] = OpenNode(std::move(domchgs), lower_bound, estimate, depth);
  }

  assert(nodes[pos].lower_bound == lower_bound);
  assert(nodes[pos].estimate == estimate);
  assert(nodes[pos].depth == depth);

  link_estim(pos);
  link_lower(pos);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestNode() {
  auto get_left = [&](int n) -> int& { return nodes[n].leftestimate; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightestimate; };
  auto get_key = [&](int n) {
    return std::make_tuple(nodes[n].lower_bound, nodes[n].estimate, n);
  };

  estimroot =
      highs_splay(std::make_tuple(-HIGHS_CONST_INF, -HIGHS_CONST_INF, 0),
                  estimroot, get_left, get_right, get_key);
  int bestestimnode = estimroot;
  unlink_estim(bestestimnode);
  unlink_lower(bestestimnode);
  freeslots.push(bestestimnode);
  return std::move(nodes[bestestimnode]);
}

HighsNodeQueue::OpenNode HighsNodeQueue::popBestBoundNode() {
  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) { return std::make_pair(nodes[n].lower_bound, n); };

  lowerroot = highs_splay(std::make_pair(-HIGHS_CONST_INF, 0), lowerroot,
                          get_left, get_right, get_key);
  int bestboundnode = lowerroot;
  unlink_estim(bestboundnode);
  unlink_lower(bestboundnode);
  freeslots.push(bestboundnode);
  return std::move(nodes[bestboundnode]);
}

double HighsNodeQueue::getBestLowerBound() {
  if (lowerroot == -1) return HIGHS_CONST_INF;

  auto get_left = [&](int n) -> int& { return nodes[n].leftlower; };
  auto get_right = [&](int n) -> int& { return nodes[n].rightlower; };
  auto get_key = [&](int n) { return std::make_pair(nodes[n].lower_bound, n); };

  lowerroot = highs_splay(std::make_pair(-HIGHS_CONST_INF, 0), lowerroot,
                          get_left, get_right, get_key);
  return nodes[lowerroot].lower_bound;
}