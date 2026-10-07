#ifndef HIGHS_COLOUR_REFINEMENT_H
#define HIGHS_COLOUR_REFINEMENT_H

#include <cassert>
#include <vector>

#include "DoubleQuantizer.h"
#include "SumClassifier.h"
#include "util/HighsLinkedLists.h"
#include "util/HighsSparseMatrix.h"

namespace highs {
namespace folding {

// Refines separate colourings of the rows and columns of matrix A (col-wise),
// seen as the biadjacency matrix of a weighted bipartite graph. Row colours and
// column colours are numbered independently, each starting from 0 and
// consecutive. Based on "Tight Lower and Upper Bounds for the Complexity of
// Canonical Colour Refinement", Berkholz, Bonsma, Grohe

class ColourRefinement {
  struct IterableStack {
    std::vector<HighsInt> s;
    HighsInt top;
    std::vector<HighsBool> in_stack;

    explicit IterableStack(HighsInt n) : s(n), top{-1}, in_stack(n, 0) {}
    bool belong(HighsInt i) const { return in_stack[i]; }
    bool empty() const { return top < 0; }
    void push(HighsInt i) {
      assert(!in_stack[i]);
      s[++top] = i;
      in_stack[i] = true;
    }
    HighsInt pop() {
      assert(!empty());
      const HighsInt elem = s[top--];
      in_stack[elem] = false;
      return elem;
    }
    HighsInt peek() const {
      assert(!empty());
      return s[top];
    }

    // Iterators in reverse pop order
    std::vector<HighsInt>::const_iterator begin() const { return s.begin(); }
    std::vector<HighsInt>::const_iterator end() const {
      return s.begin() + top + 1;
    }
  };

  class TouchedNodes {
    // Object that stores the nodes of one side that were touched by the
    // refining colour in the current step and the sum of the weights through
    // which each was reached. Nodes that were not reached have sum zero and are
    // not stored. The nodes are stored using slots, numbered in order of first
    // touch. After groupByColour(), the slots of each colour can be iterated:
    //
    //   for (HighsInt slot : touched.slots(colour))
    //     ...
    //
    // Example: the destination side has 6 nodes and two colours:
    //
    //      colour 0: nodes {0, 1, 2}
    //      colour 1: nodes {3, 4, 5}
    //
    // The refining colour reaches the nodes in this order: 4, 1, 3, 0, 5.
    // Node 2 is never reached.
    //
    // After the calls to add(), slots are numbered in order of first touch:
    //
    //      slot:             0   1   2   3   4
    //      node_of_slot_:    4   1   3   0   5
    //      colour_of_slot_:  1   0   1   0   1
    //      sum_of_slot_:     2   1   2   3   2
    //
    //      num_slots_ = 5
    //      slot_of_node_:    3   1  -1   2   0   4
    //      count_of_colour_: 2   3
    //
    //      colours_:         stack [1, 0] <- top (order of first touch)
    //
    // The slots of the two colours are currently interleaved (colours 1, 0, 1,
    // 0, 1). groupByColour() populates slots_in_order_ so that each
    // colour corresponds to a contiguous range of slots in slots_in_order_:
    //
    //      slots_in_order_ = [ 0, 2, 4  |  1, 3 ]
    //                          colour 1 |  colour 0
    //
    // It also populates end_of_colour_ so that it points to the end of the
    // range of a colour:
    //
    //      end_or_colour_ = [3, 5]
    //
    // slots(colour) return an iterable object that lists the slots used by
    // that colour. It can be used to iterate through the nodes touched,
    // grouped by colour:
    //
    //      slots(1) = 0, 2, 4
    //      slots(0) = 1, 3
    //
    // Helper functions node(slot) and sum(slot) can be used to recover the
    // node and sum corresponding to a specific slot.
    // Within a colour the slots keep their order of first touch. Untouched
    // nodes are in no range.
    //
    // The alternative approach would be to use an array of size n to store the
    // sums and a separate collection of linked list to store the nodes touched
    // by colour. Using the slots approach is beneficial when the number of
    // colours touched is much smaller than n, since the active portions of the
    // arrays are small and memory accesses become cheaper. Also, by ordering
    // the slots of each colour consecutively, the pattern of memory access
    // becomes more regular.

    std::vector<HighsInt> slot_of_node_;
    std::vector<HighsInt> node_of_slot_;
    std::vector<double> sum_of_slot_;
    std::vector<HighsInt> colour_of_slot_;
    std::vector<HighsInt> class_id_of_slot_;
    HighsInt num_slots_;

    std::vector<HighsInt> count_of_colour_;
    IterableStack colours_;

    std::vector<HighsInt> end_of_colour_;
    std::vector<HighsInt> slots_in_order_;

   public:
    explicit TouchedNodes(HighsInt num_nodes)
        : slot_of_node_(num_nodes, -1),
          node_of_slot_(num_nodes),
          sum_of_slot_(num_nodes),
          colour_of_slot_(num_nodes),
          class_id_of_slot_(num_nodes),
          num_slots_{0},
          count_of_colour_(num_nodes, 0),
          colours_(num_nodes),
          end_of_colour_(num_nodes, 0),
          slots_in_order_(num_nodes) {}

    void add(HighsInt node, double weight,
             const std::vector<HighsInt>& colour_of_node);
    void groupByColour();
    void clear();
    const IterableStack& colours() const { return colours_; }

    HighsInt colourCount(HighsInt colour) const {
      return count_of_colour_[colour];
    }
    HighsInt node(HighsInt slot) const { return node_of_slot_[slot]; }
    double sum(HighsInt slot) const { return sum_of_slot_[slot]; }

    // Manage class id from classifier
    void setClass(HighsInt slot, HighsInt class_id) {
      class_id_of_slot_[slot] = class_id;
    }
    HighsInt getClass(HighsInt slot) const { return class_id_of_slot_[slot]; }

    struct SlotRange {
      const HighsInt* first;
      const HighsInt* last;
      const HighsInt* begin() const { return first; }
      const HighsInt* end() const { return last; }
      HighsInt size() const { return static_cast<HighsInt>(last - first); }
      HighsInt front() const { return *first; }
    };

    SlotRange slots(HighsInt colour) const {
      const HighsInt* base = slots_in_order_.data();
      const HighsInt end = end_of_colour_[colour];
      return {base + end - count_of_colour_[colour], base + end};
    }
  };

  // Everything that belongs to one side of the bipartite graph
  struct Side {
    HighsInt n;
    std::vector<HighsInt>& colour;
    HighsInt latest_colour;
    HighsLinkedLists node_by_colour;
    IterableStack to_refine;
    TouchedNodes touched;

    explicit Side(std::vector<HighsInt>& c);

    HighsInt numUntouchedNodes(HighsInt c) const {
      return node_by_colour.length(c) - touched.colourCount(c);
    }
    void moveNode(HighsInt node, HighsInt from, HighsInt to) {
      node_by_colour.remove(node, from);
      node_by_colour.append(node, to);
      colour[node] = to;
    }
    HighsInt peekLength() const {
      return to_refine.empty() ? kHighsIInf
                               : node_by_colour.length(to_refine.peek());
    }
    HighsInt coloursUsed() const { return latest_colour + 1; }
  };

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;
  Side rows_;
  Side cols_;
  SumClassifier classifier_;

  void refine(Side& src, Side& dst, const HighsSparseMatrix& M);
  void computeColourSums(HighsInt refining_colour, Side& src, Side& dst,
                         const HighsSparseMatrix& M);
  void splitColours(Side& side);
  bool colourSplits(const Side& side, HighsInt colour) const;
  void splitColour(Side& side, HighsInt colour);

  void countSums(Side& side, HighsInt colour);
  void assignNewColours(Side& side, HighsInt colour);
  void moveNodesToNewColours(Side& side, HighsInt colour);

  void prepareNextIter(Side& side);

  double time_sums_{};
  mutable double time_check_{};
  double time_count_sums_{};
  double time_assign_{};
  double time_new_colour_{};
  double time_prepare_{};
  void printTimes() const;

 public:
  ColourRefinement(const HighsSparseMatrix& A,
                   std::vector<HighsInt>& row_colour,
                   std::vector<HighsInt>& col_colour);
  void run();

  HighsInt rowColoursUsed() const { return rows_.coloursUsed(); }
  HighsInt colColoursUsed() const { return cols_.coloursUsed(); }
};

}  // namespace folding
}  // namespace highs

#endif