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

#include <map>
#include <queue>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsCDouble.h"

class HighsDomain;
class HighsLpRelaxation;

class HighsNodeQueue {
 public:
  struct OpenNode {
    std::vector<HighsDomainChange> domchgstack;
    std::vector<std::multimap<double, int>::iterator> domchglinks;
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

  void checkGlobalBounds(int col, double lb, double ub, double feastol,
                         HighsCDouble& treeweight);

 private:
  std::vector<OpenNode> nodes;
  std::vector<std::multimap<double, int>> colLowerNodes;
  std::vector<std::multimap<double, int>> colUpperNodes;
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;
  int lowerroot = -1;
  int estimroot = -1;

  void link_estim(int node);

  void unlink_estim(int node);

  void link_lower(int node);

  void unlink_lower(int node);

  void link_domchgs(int node);

  void unlink_domchgs(int node);

  void link(int node);

  void unlink(int node);

 public:
  double performBounding(double upper_limit);

  void setNumCol(int numcol);

  void emplaceNode(std::vector<HighsDomainChange>&& domchgs, double lower_bound,
                   double estimate, int depth);

  OpenNode popBestNode();

  OpenNode popBestBoundNode();

  OpenNode popRelatedNode(const HighsLpRelaxation& lprelax);

  size_t numNodesUp(int col) const { return colLowerNodes[col].size(); }

  size_t numNodesDown(int col) const { return colUpperNodes[col].size(); }

  size_t numNodesUp(int col, double val) const {
    auto it = colLowerNodes[col].upper_bound(val);
    if (it == colLowerNodes[col].begin()) return colLowerNodes[col].size();
    return std::distance(colLowerNodes[col].upper_bound(val),
                         colLowerNodes[col].end());
  }

  size_t numNodesDown(int col, double val) const {
    auto it = colUpperNodes[col].lower_bound(val);
    if (it == colUpperNodes[col].end()) return colUpperNodes[col].size();
    return std::distance(colUpperNodes[col].begin(), it);
  }

  double pruneInfeasibleNodes(HighsDomain& globaldomain, double feastol);

  double getBestLowerBound();

  void clear() {
    nodes.clear();
    decltype(freeslots)().swap(freeslots);
  }

  bool empty() const { return nodes.size() == freeslots.size(); }

  size_t numNodes() const { return nodes.size() - freeslots.size(); }
};

#endif
