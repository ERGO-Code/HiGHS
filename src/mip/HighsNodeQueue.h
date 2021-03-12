/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_NODE_QUEUE_H_
#define HIGHS_NODE_QUEUE_H_

#include <queue>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomainChange.h"

class HighsNodeQueue {
 public:
  struct OpenNode {
    std::vector<HighsDomainChange> domchgstack;
    double lower_bound;
    double estimate;
    int depth;
    int leftlower;
    int rightlower;
    int leftestimate;
    int rightestimate;

    OpenNode()
        : domchgstack(),
          lower_bound(-HIGHS_CONST_INF),
          estimate(-HIGHS_CONST_INF),
          depth(0),
          leftlower(-1),
          rightlower(-1),
          leftestimate(-1),
          rightestimate(-1) {}

    OpenNode(std::vector<HighsDomainChange>&& domchgstack, double lower_bound,
             double estimate, int depth)
        : domchgstack(domchgstack),
          lower_bound(lower_bound),
          estimate(estimate),
          depth(depth),
          leftlower(-1),
          rightlower(-1),
          leftestimate(-1),
          rightestimate(-1) {}

    OpenNode& operator=(OpenNode&& other) = default;
    OpenNode(OpenNode&&) = default;

    OpenNode& operator=(const OpenNode& other) = delete;
    OpenNode(const OpenNode&) = delete;
  };

 private:
  std::vector<OpenNode> nodes;
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;
  int lowerroot = -1;
  int estimroot = -1;

  void link_estim(int node);

  void unlink_estim(int node);

  void link_lower(int node);

  void unlink_lower(int node);

 public:
  double performBounding(double upper_limit);

  void emplaceNode(std::vector<HighsDomainChange>&& domchgs, double lower_bound,
                   double estimate, int depth);

  OpenNode popBestNode();

  OpenNode popBestBoundNode();

  double getBestLowerBound();

  void clear() {
    nodes.clear();
    decltype(freeslots)().swap(freeslots);
  }

  bool empty() const { return nodes.size() == freeslots.size(); }

  size_t numNodes() const { return nodes.size() - freeslots.size(); }
};

#endif