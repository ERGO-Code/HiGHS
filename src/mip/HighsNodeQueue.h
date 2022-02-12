/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
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
#include "util/HighsRbTree.h"

class HighsDomain;
class HighsLpRelaxation;

class HighsNodeQueue {
 public:
  template <typename T>
  class NodesetAllocator {
    template <int S>
    struct ChunkWithSize {
      ChunkWithSize* next;
      typename std::aligned_storage<S, alignof(T)>::type storage;
    };

    using Chunk =
        ChunkWithSize<4096 - offsetof(ChunkWithSize<sizeof(T)>, storage)>;

    union FreelistNode {
      FreelistNode* next;
      typename std::aligned_storage<sizeof(T), alignof(T)>::type storage;
    };

    FreelistNode* freeListHead = nullptr;
    char* currChunkStart = nullptr;
    char* currChunkEnd = nullptr;
    Chunk* chunkListHead = nullptr;

    void clear() noexcept {
      freeListHead = nullptr;
      currChunkStart = nullptr;
      currChunkEnd = nullptr;
      while (chunkListHead) {
        Chunk* delChunk = chunkListHead;
        chunkListHead = delChunk->next;
        delete delChunk;
      }
    }

   public:
    using value_type = T;
    using size_type = std::size_t;
    using propagate_on_container_move_assignment = std::true_type;

    NodesetAllocator() noexcept = default;
    NodesetAllocator(const NodesetAllocator&) noexcept {}
    template <typename U>
    NodesetAllocator(const NodesetAllocator<U>&) noexcept {}
    NodesetAllocator(NodesetAllocator&& other) noexcept
        : freeListHead(other.freeListHead),
          currChunkStart(other.currChunkStart),
          currChunkEnd(other.currChunkEnd),
          chunkListHead(other.chunkListHead) {
      other.currChunkStart = nullptr;
      other.currChunkEnd = nullptr;
      other.chunkListHead = nullptr;
      other.freeListHead = nullptr;
    }

    NodesetAllocator& operator=(const NodesetAllocator&) noexcept {
      return *this;
    }

    NodesetAllocator& operator=(NodesetAllocator&& other) noexcept {
      clear();
      freeListHead = other.freeListHead;
      currChunkStart = other.currChunkStart;
      currChunkEnd = other.currChunkEnd;
      chunkListHead = other.chunkListHead;
      other.chunkListHead = nullptr;
      other.freeListHead = nullptr;
      other.currChunkStart = nullptr;
      other.currChunkEnd = nullptr;
      return *this;
    }

    ~NodesetAllocator() noexcept { clear(); }

    T* allocate(size_type n) {
      if (n == 1) {
        T* ptr = reinterpret_cast<T*>(freeListHead);
        if (ptr) {
          freeListHead = freeListHead->next;
        } else {
          ptr = reinterpret_cast<T*>(currChunkStart);
          currChunkStart += sizeof(FreelistNode);
          if (currChunkStart > currChunkEnd) {
            auto newChunk = new Chunk;
            newChunk->next = chunkListHead;
            chunkListHead = newChunk;
            currChunkStart = reinterpret_cast<char*>(&newChunk->storage);
            currChunkEnd = currChunkStart + sizeof(newChunk->storage);
            ptr = reinterpret_cast<T*>(currChunkStart);
            currChunkStart += sizeof(FreelistNode);
          }
        }
        return ptr;
      }

      return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* ptr, size_type n) noexcept {
      if (n == 1) {
        FreelistNode* node = reinterpret_cast<FreelistNode*>(ptr);
        node->next = freeListHead;
        freeListHead = node;
      } else {
        ::operator delete(ptr);
      }
    }
  };

  using NodeSet = std::set<std::pair<double, int64_t>,
                           std::less<std::pair<double, int64_t>>,
                           NodesetAllocator<std::pair<double, int64_t>>>;

  struct OpenNode {
    std::vector<HighsDomainChange> domchgstack;
    std::vector<HighsInt> branchings;
    std::vector<NodeSet::iterator> domchglinks;
    double lower_bound;
    double estimate;
    HighsInt depth;
    highs::RbTreeLinks<int64_t> lowerLinks;
    highs::RbTreeLinks<int64_t> hybridEstimLinks;

    OpenNode()
        : domchgstack(),
          branchings(),
          domchglinks(),
          lower_bound(-kHighsInf),
          estimate(-kHighsInf),
          depth(0),
          lowerLinks(),
          hybridEstimLinks() {}

    OpenNode(std::vector<HighsDomainChange>&& domchgstack,
             std::vector<HighsInt>&& branchings, double lower_bound,
             double estimate, HighsInt depth)
        : domchgstack(domchgstack),
          branchings(branchings),
          lower_bound(lower_bound),
          estimate(estimate),
          depth(depth),
          lowerLinks(),
          hybridEstimLinks() {}

    OpenNode& operator=(OpenNode&& other) = default;
    OpenNode(OpenNode&&) = default;

    OpenNode& operator=(const OpenNode& other) = delete;
    OpenNode(const OpenNode&) = delete;
  };

  void checkGlobalBounds(HighsInt col, double lb, double ub, double feastol,
                         HighsCDouble& treeweight);

 private:
  class NodeLowerRbTree;
  class NodeHybridEstimRbTree;
  class SuboptimalNodeRbTree;

  NodesetAllocator<std::pair<double, int64_t>> allocator;
  std::vector<OpenNode> nodes;
  std::vector<NodeSet> colLowerNodes;
  std::vector<NodeSet> colUpperNodes;
  std::priority_queue<int64_t, std::vector<int64_t>, std::greater<int64_t>>
      freeslots;
  int64_t lowerRoot = -1;
  int64_t lowerMin = -1;
  int64_t hybridEstimRoot = -1;
  int64_t hybridEstimMin = -1;
  int64_t suboptimalRoot = -1;
  int64_t suboptimalMin = -1;
  int64_t numSuboptimal = 0;
  double optimality_limit = kHighsInf;

  void link_estim(int64_t node);

  void unlink_estim(int64_t node);

  void link_lower(int64_t node);

  void unlink_lower(int64_t node);

  void link_suboptimal(int64_t node);

  void unlink_suboptimal(int64_t node);

  void link_domchgs(int64_t node);

  void unlink_domchgs(int64_t node);

  double link(int64_t node);

  void unlink(int64_t node);

 public:
  void setOptimalityLimit(double optimality_limit) {
    this->optimality_limit = optimality_limit;
  }

  double performBounding(double upper_limit);

  void setNumCol(HighsInt numcol);

  double emplaceNode(std::vector<HighsDomainChange>&& domchgs,
                     std::vector<HighsInt>&& branchings, double lower_bound,
                     double estimate, HighsInt depth);

  OpenNode&& popBestNode();

  OpenNode&& popBestBoundNode();

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

  const NodeSet& getUpNodes(HighsInt col) const { return colLowerNodes[col]; }

  const NodeSet& getDownNodes(HighsInt col) const { return colUpperNodes[col]; }

  double pruneInfeasibleNodes(HighsDomain& globaldomain, double feastol);

  double pruneNode(int64_t nodeId);

  double getBestLowerBound() const;

  HighsInt getBestBoundDomchgStackSize() const;

  void clear() {
    HighsNodeQueue nodequeue;
    nodequeue.setNumCol(colUpperNodes.size());
    *this = std::move(nodequeue);
  }

  int64_t numNodes() const { return nodes.size() - freeslots.size(); }

  int64_t numActiveNodes() const {
    return nodes.size() - freeslots.size() - numSuboptimal;
  }

  bool empty() const { return numActiveNodes() == 0; }
};

template <typename T, typename U>
bool operator==(const HighsNodeQueue::NodesetAllocator<T>&,
                const HighsNodeQueue::NodesetAllocator<U>&) {
  return true;
}

template <typename T, typename U>
bool operator!=(const HighsNodeQueue::NodesetAllocator<T>&,
                const HighsNodeQueue::NodesetAllocator<U>&) {
  return false;
}

#endif
