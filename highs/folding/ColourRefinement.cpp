#include "ColourRefinement.h"

#include <algorithm>

namespace highs {
namespace folding {

void ColourRefinement::TouchedNodes::add(
    HighsInt node, double weight, const std::vector<HighsInt>& colour_of_node) {
  HighsInt slot = slot_of_node_[node];
  if (slot >= 0) {
    sum_of_slot_[slot] += weight;
    return;
  }
  slot = num_slots_++;
  slot_of_node_[node] = slot;
  node_of_slot_[slot] = node;
  sum_of_slot_[slot] = weight;

  const HighsInt colour = colour_of_node[node];
  colour_of_slot_[slot] = colour;
  if (count_of_colour_[colour]++ == 0) colours_.push(colour);
}

void ColourRefinement::TouchedNodes::groupByColour() {
  HighsInt next_free = 0;
  for (HighsInt colour : colours_) {
    end_of_colour_[colour] = next_free;
    next_free += count_of_colour_[colour];
  }

  for (HighsInt slot = 0; slot < num_slots_; ++slot) {
    const HighsInt colour = colour_of_slot_[slot];
    const HighsInt pos = end_of_colour_[colour];
    end_of_colour_[colour]++;
    slots_in_order_[pos] = slot;
  }
}

void ColourRefinement::TouchedNodes::clear() {
  for (HighsInt slot = 0; slot < num_slots_; ++slot) {
    const HighsInt node = node_of_slot_[slot];
    slot_of_node_[node] = -1;
  }

  num_slots_ = 0;

  while (!colours_.empty()) {
    const HighsInt colour = colours_.pop();
    count_of_colour_[colour] = 0;
  }
}

ColourRefinement::Side::Side(std::vector<HighsInt>& c)
    : n{static_cast<HighsInt>(c.size())},
      colour{c},
      latest_colour{-1},
      to_refine(n),
      touched(n) {
  node_by_colour.init(n, n);

  for (HighsInt i = 0; i < n; ++i) node_by_colour.append(i, colour[i]);
  if (n > 0) latest_colour = *std::max_element(colour.begin(), colour.end());
  for (HighsInt c = 0; c <= latest_colour; ++c) to_refine.push(c);
}

ColourRefinement::ColourRefinement(const HighsSparseMatrix& A,
                                   std::vector<HighsInt>& row_colour,
                                   std::vector<HighsInt>& col_colour)
    : A_{A}, rows_(row_colour), cols_(col_colour) {
  assert(A.isColwise());
  assert(rows_.n == A.num_row_ && cols_.n == A.num_col_);

  A_.buildOppositeFormat(At_);
}

void ColourRefinement::computeColourSums(HighsInt refining_colour, Side& src,
                                         Side& dst,
                                         const HighsSparseMatrix& M) {
  // Sum of the weights between each dst node and the nodes of refining_colour
  for (HighsInt v : src.node_by_colour.list(refining_colour)) {
    for (HighsInt k = M.start_[v]; k < M.start_[v + 1]; ++k) {
      dst.touched.add(M.index_[k], M.value_[k], dst.colour);
    }
  }
  dst.touched.groupByColour();
}

void ColourRefinement::splitColours(Side& side) {
  for (HighsInt colour : side.touched.colours())
    if (colourSplits(side, colour)) splitColour(side, colour);
}

// Check if the nodes all fall in the same quantization bucket. In this way, the
// unordered_map is not cleared and populated unless it is needed.
bool ColourRefinement::colourSplits(const Side& side, HighsInt colour) const {
  const auto slots = side.touched.slots(colour);
  const bool has_untouched_nodes = side.numUntouchedNodes(colour) > 0;
  const double quantized_reference = DoubleQuantizer::bucket(
      has_untouched_nodes ? 0.0 : side.touched.sum(slots.front()));

  bool colour_split = false;
  for (HighsInt slot : slots) {
    if (DoubleQuantizer::bucket(side.touched.sum(slot)) !=
        quantized_reference) {
      colour_split = true;
      break;
    }
  }

  return colour_split;
}

void ColourRefinement::countSums(Side& side, HighsInt colour) {
  // Classify the nodes of the colour by quantized sum, counting how many nodes
  // are in each class. The class of each touched node is saved, so it does not
  // need to be looked up again when the nodes are moved.

  const auto slots = side.touched.slots(colour);

  // At most one class per touched node, plus the class of untouched nodes
  classifier_.reset(slots.size() + 1);

  const HighsInt nodes_with_zero_sum = side.numUntouchedNodes(colour);
  if (nodes_with_zero_sum > 0) classifier_.add(0.0, nodes_with_zero_sum);

  for (HighsInt slot : slots) {
    const HighsInt new_class_id = classifier_.add(side.touched.sum(slot));
    side.touched.setClass(slot, new_class_id);
  }
}

void ColourRefinement::assignNewColours(Side& side, HighsInt colour) {
  // classifier_ contains the new colour classes in which colour is divided.

  // All colour classes are added to the stack for later refinement, apart from
  // the largest one. One class can reuse the current colour: we use the class
  // with sum zero, if it exists, or the smallest class.
  const HighsInt largest = classifier_.largest();
  const HighsInt reusing_colour = classifier_.zeroClass() >= 0
                                      ? classifier_.zeroClass()
                                      : classifier_.smallest();

  // Add colours to the stack
  const bool colour_in_stack = side.to_refine.belong(colour);

  for (HighsInt j = 0; j < classifier_.size(); ++j) {
    SumData& data = classifier_[j];
    if (j == reusing_colour) {
      data.colour = colour;
      if (!colour_in_stack && j != largest) side.to_refine.push(colour);

    } else {
      ++side.latest_colour;
      data.colour = side.latest_colour;
      if (colour_in_stack || j != largest) side.to_refine.push(data.colour);
    }
  }
}

void ColourRefinement::moveNodesToNewColours(Side& side, HighsInt colour) {
  for (HighsInt slot : side.touched.slots(colour)) {
    const HighsInt node = side.touched.node(slot);
    const HighsInt class_id = side.touched.getClass(slot);
    const HighsInt new_colour = classifier_[class_id].colour;
    if (new_colour != colour) {
      side.moveNode(node, colour, new_colour);
    }
  }
}

void ColourRefinement::splitColour(Side& side, HighsInt colour) {
  countSums(side, colour);
  assignNewColours(side, colour);
  moveNodesToNewColours(side, colour);
}

void ColourRefinement::prepareNextIter(Side& side) { side.touched.clear(); }

// Use top colour of src to refine the colours of dst.
// M must be indexed by nodes of src and hold indices of nodes of dst.
// - if src is rows and dst is cols, then M must be At
// - if src is cols and dst is rows, then M must be A
//
void ColourRefinement::refine(Side& src, Side& dst,
                              const HighsSparseMatrix& M) {
  const HighsInt refining_colour = src.to_refine.pop();
  computeColourSums(refining_colour, src, dst, M);
  splitColours(dst);
  prepareNextIter(dst);
}

HighsInt ColourRefinement::run() {
  while (!rows_.to_refine.empty() || !cols_.to_refine.empty()) {
    if (rows_.peekLength() < cols_.peekLength())
      refine(rows_, cols_, At_);
    else
      refine(cols_, rows_, A_);

    if (rows_.coloursUsed() > rows_.n * kFractionStop &&
        cols_.coloursUsed() > cols_.n * kFractionStop) {
      return 1;
    }
  }

  return 0;
}

}  // namespace folding
}  // namespace highs