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

#ifndef HIGHS_NODE_QUEUE_H_
#define HIGHS_NODE_QUEUE_H_

#include <cassert>
#include <queue>
#include <set>
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
    std::vector<HighsInt> branchings;
    std::vector<std::set<std::pair<double, HighsInt>>::iterator> domchglinks;
    double lower_bound;
    double estimate;
    HighsInt depth;
    HighsInt leftlower;
    HighsInt rightlower;
    HighsInt leftestimate;
    HighsInt rightestimate;

    OpenNode()
        : domchgstack(),
          branchings(),
          domchglinks(),
          lower_bound(-kHighsInf),
          estimate(-kHighsInf),
          depth(0),
          leftlower(-1),
          rightlower(-1),
          leftestimate(-1),
          rightestimate(-1) {}

    OpenNode(std::vector<HighsDomainChange>&& domchgstack,
             std::vector<HighsInt>&& branchings, double lower_bound,
             double estimate, HighsInt depth)
        : domchgstack(domchgstack),
          branchings(branchings),
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

  void checkGlobalBounds(HighsInt col, double lb, double ub, double feastol,
                         HighsCDouble& treeweight);

 private:
  std::vector<OpenNode> nodes;
  std::vector<std::set<std::pair<double, HighsInt>>> colLowerNodes;
  std::vector<std::set<std::pair<double, HighsInt>>> colUpperNodes;
  std::priority_queue<HighsInt, std::vector<HighsInt>, std::greater<HighsInt>>
      freeslots;
  HighsInt lowerroot = -1;
  HighsInt estimroot = -1;

  void link_estim(HighsInt node);

  void unlink_estim(HighsInt node);

  void link_lower(HighsInt node);

  void unlink_lower(HighsInt node);

  void link_domchgs(HighsInt node);

  void unlink_domchgs(HighsInt node);

  void link(HighsInt node);

  void unlink(HighsInt node);

 public:
  double performBounding(double upper_limit);

  void setNumCol(HighsInt numcol);

  void emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                   std::vector<HighsInt>&& branchings, double lower_bound,
                   double estimate, HighsInt depth);

  OpenNode popBestNode();

  OpenNode popBestBoundNode();

  int64_t numNodesUp(HighsInt col) const { return colLowerNodes[col].size(); }

  int64_t numNodesDown(HighsInt col) const { return colUpperNodes[col].size(); }

  int64_t numNodesUp(HighsInt col, double val) const {
    assert((HighsInt)colLowerNodes.size() > col);
    auto it = colLowerNodes[col].upper_bound(std::make_pair(val, kHighsIInf));
    if (it == colLowerNodes[col].begin()) return colLowerNodes[col].size();
    return std::distance(it, colLowerNodes[col].end());
  }

  int64_t numNodesDown(HighsInt col, double val) const {
    assert((HighsInt)colUpperNodes.size() > col);
    auto it = colUpperNodes[col].lower_bound(std::make_pair(val, -1));
    if (it == colUpperNodes[col].end()) return colUpperNodes[col].size();
    return std::distance(colUpperNodes[col].begin(), it);
  }

  const std::set<std::pair<double, HighsInt>>& getUpNodes(HighsInt col) const {
    return colLowerNodes[col];
  }

  const std::set<std::pair<double, HighsInt>>& getDownNodes(
      HighsInt col) const {
    return colUpperNodes[col];
  }

  double pruneInfeasibleNodes(HighsDomain& globaldomain, double feastol);

  double pruneNode(HighsInt nodeId);

  double getBestLowerBound();

  void clear() {
    HighsNodeQueue nodequeue;
    nodequeue.setNumCol(colUpperNodes.size());
    std::swap(*this, nodequeue);
  }

  bool empty() const { return nodes.size() == freeslots.size(); }

  int64_t numNodes() const { return nodes.size() - freeslots.size(); }
};

#endif
