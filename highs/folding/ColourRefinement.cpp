#include "ColourRefinement.h"

#include <algorithm>

#include "ipm/hipo/auxiliary/Auxiliary.h"

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
  hipo::Clock clock;
  // Sum of the weights between each dst node and the nodes of refining_colour
  for (HighsInt v : src.node_by_colour.list(refining_colour)) {
    for (HighsInt k = M.start_[v]; k < M.start_[v + 1]; ++k) {
      dst.touched.add(M.index_[k], M.value_[k], dst.colour);
    }
  }
  dst.touched.groupByColour();
  time_sums_ += clock.stop();
}

void ColourRefinement::splitColours(Side& side) {
  for (HighsInt colour : side.touched.colours())
    if (colourSplits(side, colour)) splitColour(side, colour);
}

// Check if the nodes all fall in the same quantization bucket. In this way, the
// unordered_map is not cleared and populated unless it is needed.
bool ColourRefinement::colourSplits(const Side& side, HighsInt colour) const {
  hipo::Clock clock;

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

  time_check_ += clock.stop();

  return colour_split;
}

void ColourRefinement::countSums(const Side& side, HighsInt colour) {
  // Populate sum_map_ with the sums of the nodes of the colour, counting how
  // many nodes have each sum. The map uses the DoubleQuantizer, so sums that
  // are in the same bucket are considered the same.
  hipo::Clock clock;

  sum_map_.clear();

  const HighsInt nodes_with_zero_sum = side.numUntouchedNodes(colour);
  if (nodes_with_zero_sum > 0)
    sum_map_.insert({0.0, SumData{nodes_with_zero_sum, -1}});

  for (HighsInt slot : side.touched.slots(colour)) {
    auto result = sum_map_.insert({side.touched.sum(slot), SumData{1, -1}});
    if (!result.second) {
      SumData& data = result.first->second;
      data.count++;
    }
  }

  time_count_sums_ += clock.stop();
}

void ColourRefinement::assignNewColours(Side& side, HighsInt colour) {
  hipo::Clock clock;

  // sum_map_ contains the new colour classes in which split_colour is divided.
  // Find iterators to the largest and smallest of these classes
  auto it_largest = sum_map_.begin();
  auto it_smallest = sum_map_.begin();

  for (auto it = sum_map_.begin(); it != sum_map_.end(); ++it) {
    const SumData& data = it->second;
    if (data.count > it_largest->second.count) it_largest = it;
    if (data.count < it_smallest->second.count) it_smallest = it;
  }

  // All colour classes are added to the stack for later refinement, apart from
  // the largest one. One class can reuse the current colour: we use the class
  // with sum zero, if it exists, or the smallest class.
  const auto it_zero = sum_map_.find(0.0);
  const auto it_reusing_colour =
      it_zero != sum_map_.end() ? it_zero : it_smallest;

  // Add colours to the stack
  const bool colour_in_stack = side.to_refine.belong(colour);

  for (auto it = sum_map_.begin(); it != sum_map_.end(); ++it) {
    SumData& data = it->second;
    if (it == it_reusing_colour) {
      data.colour = colour;
      if (!colour_in_stack && it != it_largest) side.to_refine.push(colour);

    } else {
      ++side.latest_colour;
      data.colour = side.latest_colour;
      if (colour_in_stack || it != it_largest) side.to_refine.push(data.colour);
    }
  }

  time_assign_ += clock.stop();
}

void ColourRefinement::moveNodesToNewColours(Side& side, HighsInt colour) {
  hipo::Clock clock;

  for (HighsInt slot : side.touched.slots(colour)) {
    const HighsInt node = side.touched.node(slot);
    const double sum = side.touched.sum(slot);
    const HighsInt new_colour = sum_map_.find(sum)->second.colour;
    if (new_colour != colour) {
      side.moveNode(node, colour, new_colour);
    }
  }

  time_new_colour_ += clock.stop();
}

void ColourRefinement::splitColour(Side& side, HighsInt colour) {
  countSums(side, colour);
  assignNewColours(side, colour);
  moveNodesToNewColours(side, colour);
}

void ColourRefinement::prepareNextIter(Side& side) {
  hipo::Clock clock;
  side.touched.clear();
  time_prepare_ += clock.stop();
}

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

void ColourRefinement::run() {
  while (!rows_.to_refine.empty() || !cols_.to_refine.empty()) {
    if (rows_.peekLength() < cols_.peekLength())
      refine(rows_, cols_, At_);
    else
      refine(cols_, rows_, A_);

    if (rows_.coloursUsed() > rows_.n * kFractionStop &&
        cols_.coloursUsed() > cols_.n * kFractionStop) {
      printf("Failed\n");
      break;
    }
  }

  printTimes();
}

void ColourRefinement::printTimes() const {
  const double time_total = time_sums_ + time_check_ + time_count_sums_ +
                            time_assign_ + time_new_colour_ + time_prepare_;

  printf("ColourRefinement timers\n");
  printf("Total         %f\n", time_total);
  printf("  sums        %f\n", time_sums_);
  printf("  check       %f\n", time_check_);
  printf("  count sums  %f\n", time_count_sums_);
  printf("  assign      %f\n", time_assign_);
  printf("  new colour  %f\n", time_new_colour_);
  printf("  prepare     %f\n", time_prepare_);
}

}  // namespace folding
}  // namespace highs