/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SEARCH_H_
#define HIGHS_SEARCH_H_

#include <cstdint>
#include <queue>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSeparation.h"
#include "util/HighsHash.h"

class HighsMipSolver;
class HighsImplications;
class HighsCliqueTable;

class HighsSearch {
  HighsMipSolver& mipsolver;
  HighsLpRelaxation* lp;
  HighsDomain localdom;
  HighsPseudocost pseudocost;
  HighsRandom random;
  int64_t nnodes;
  int64_t lpiterations;
  int64_t heurlpiterations;
  int64_t sblpiterations;
  double upper_limit;
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  HighsInt depthoffset;
  bool inbranching;
  bool inheuristic;

 public:
  enum class ChildSelectionRule {
    Up,
    Down,
    RootSol,
    Obj,
    Random,
    BestCost,
    WorstCost,
    Disjunction,
  };

  enum class NodeResult {
    BoundExceeding,
    DomainInfeasible,
    LpInfeasible,
    Branched,
    Open,
  };

 private:
  ChildSelectionRule childselrule;

  HighsCDouble treeweight;

  struct NodeData {
    double lower_bound;
    double estimate;
    double branching_point;
    // we store the lp objective separately to the lower bound since the lower
    // bound could be above the LP objective when cuts age out or below when the
    // LP is unscaled dual infeasible and it is not set. We still want to use
    // the objective for pseudocost updates and tiebreaking of best bound node
    // selection
    double lp_objective;
    HighsDomainChange branchingdecision;
    HighsInt domgchgStackPos;
    uint8_t opensubtrees;

    NodeData(double parentlb = -HIGHS_CONST_INF,
             double parentestimate = -HIGHS_CONST_INF)
        : lower_bound(parentlb),
          estimate(parentestimate),
          lp_objective(-HIGHS_CONST_INF),
          domgchgStackPos(-1),
          opensubtrees(2) {}
  };

  std::vector<NodeData> nodestack;
  HighsHashTable<HighsInt, int> reliableatnode;

  bool branchingVarReliableAtNode(HighsInt col) const {
    auto it = reliableatnode.find(col);
    if (it == nullptr || *it != 3) return false;

    return true;
  }

  void markBranchingVarUpReliableAtNode(HighsInt col) {
    reliableatnode[col] |= 1;
  }

  void markBranchingVarDownReliableAtNode(HighsInt col) {
    reliableatnode[col] |= 2;
  }

 public:
  HighsSearch(HighsMipSolver& mipsolver, const HighsPseudocost& pseudocost);

  void setRINSNeighbourhood(const std::vector<double>& basesol,
                            const std::vector<double>& relaxsol);

  void setRENSNeighbourhood(const std::vector<double>& lpsol);

  double getCutoffBound() const;

  void setLpRelaxation(HighsLpRelaxation* lp) { this->lp = lp; }

  double checkSol(const std::vector<double>& sol, bool& integerfeasible) const;

  void createNewNode();

  void cutoffNode();

  void branchDownwards(HighsInt col, double newub, double branchpoint);

  void branchUpwards(HighsInt col, double newlb, double branchpoint);

  void setMinReliable(HighsInt minreliable);

  void setHeuristic(bool inheuristic) { this->inheuristic = inheuristic; }

  void addBoundExceedingConflict();

  void resetLocalDomain();

  int64_t getHeuristicLpIterations() const;

  int64_t getTotalLpIterations() const;

  int64_t getLocalLpIterations() const;

  int64_t getStrongBranchingLpIterations() const;

  bool hasNode() const { return !nodestack.empty(); }

  bool currentNodePruned() const { return nodestack.back().opensubtrees == 0; }

  double getCurrentEstimate() const { return nodestack.back().estimate; }

  double getCurrentLowerBound() const { return nodestack.back().lower_bound; }

  HighsInt getCurrentDepth() const { return nodestack.size() + depthoffset; }

  void openNodesToQueue(HighsNodeQueue& nodequeue);

  void currentNodeToQueue(HighsNodeQueue& nodequeue);

  void flushStatistics();

  void installNode(HighsNodeQueue::OpenNode&& node);

  void addInfeasibleConflict();

  HighsInt selectBranchingCandidate(int64_t maxSbIters);

  void evalUnreliableBranchCands();

  const NodeData* getParentNodeData() const;

  NodeResult evaluateNode();

  NodeResult branch();

  bool backtrack();

  void printDisplayLine(char first, bool header = false);

  NodeResult dive();

  HighsDomain& getLocalDomain() { return localdom; }

  const HighsDomain& getLocalDomain() const { return localdom; }

  HighsPseudocost& getPseudoCost() { return pseudocost; }

  const HighsPseudocost& getPseudoCost() const { return pseudocost; }

  void solveDepthFirst(int64_t maxbacktracks = 1);
};

#endif
