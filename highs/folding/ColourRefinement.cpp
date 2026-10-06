#include "ColourRefinement.h"

#include <algorithm>

#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace highs {
namespace folding {

ColourRefinement::Side::Side(std::vector<HighsInt>& c)
    : n{static_cast<HighsInt>(c.size())},
      colour{c},
      latest_colour{-1},
      sums(n, 0.0),
      node_touched(n, false),
      to_refine(n),
      colours_touched(n) {
  node_by_colour.init(n, n);
  node_touched_by_colour.init(n, n);

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
      touchNeighbour(k, M, dst);
    }
  }
  time_sums_ += clock.stop();
}

void ColourRefinement::touchNeighbour(HighsInt k, const HighsSparseMatrix& M,
                                      Side& side) {
  const HighsInt w = M.index_[k];
  side.sums[w] += M.value_[k];
  if (!side.node_touched[w]) {
    side.node_touched[w] = true;
    side.node_touched_by_colour.append(w, side.colour[w]);
    if (!side.colours_touched.belong(side.colour[w]))
      side.colours_touched.push(side.colour[w]);
  }
}

void ColourRefinement::splitColours(Side& side) {
  for (HighsInt split_colour : side.colours_touched)
    if (colourSplits(side, split_colour)) splitColour(side, split_colour);
}

// Check if the nodes all fall in the same quantization bucket. In this way, the
// unordered_map is not cleared and populated unless it is needed.
bool ColourRefinement::colourSplits(const Side& side, HighsInt colour) const {
  hipo::Clock clock;

  const bool has_untouched_nodes = side.node_by_colour.length(colour) >
                                   side.node_touched_by_colour.length(colour);

  const double quantized_reference = DoubleQuantizer::bucket(
      has_untouched_nodes
          ? 0.0
          : side.sums[side.node_touched_by_colour.list(colour).front()]);

  bool colour_split = false;

  for (HighsInt v : side.node_touched_by_colour.list(colour)) {
    if (DoubleQuantizer::bucket(side.sums[v]) != quantized_reference) {
      colour_split = true;
      break;
    }
  }

  time_check_ += clock.stop();

  return colour_split;
}

void ColourRefinement::splitColour(Side& side, HighsInt split_colour) {
  // Given a colour touched by the refining colour, and the colour sums,
  // populate the map with the sums corresponding to each vertex of that colour.
  // The map uses the DoubleQuantizer, so sums that are in the same bucket are
  // considered the same.

  hipo::Clock clock;

  sum_map_.clear();

  // Insert each sum in side.sums into sum_map_ and count how many times each
  // sum appears.

  const HighsInt nodes_with_zero_sum =
      side.node_by_colour.length(split_colour) -
      side.node_touched_by_colour.length(split_colour);

  if (nodes_with_zero_sum > 0)
    sum_map_.insert({0.0, SumData{nodes_with_zero_sum, -1}});

  for (HighsInt v : side.node_touched_by_colour.list(split_colour)) {
    auto result = sum_map_.insert({side.sums[v], SumData{1, -1}});
    if (!result.second) {
      SumData& data = result.first->second;
      data.count++;
    }
  }

  time_insert_ += clock.stop();
  clock.start();

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
  // the largest one. One class can reuse the current split_colour: we use the
  // class with sum zero, if it exists, or the smallest class.
  const auto it_zero = sum_map_.find(0.0);
  const auto it_reusing_colour =
      it_zero != sum_map_.end() ? it_zero : it_smallest;

  // Add colours to the stack
  const bool split_colour_in_stack = side.to_refine.belong(split_colour);

  for (auto it = sum_map_.begin(); it != sum_map_.end(); ++it) {
    SumData& data = it->second;
    if (it == it_reusing_colour) {
      data.colour = split_colour;
      if (!split_colour_in_stack && it != it_largest)
        side.to_refine.push(split_colour);

    } else {
      ++side.latest_colour;
      data.colour = side.latest_colour;
      if (split_colour_in_stack || it != it_largest)
        side.to_refine.push(data.colour);
    }
  }

  time_add_stack_ += clock.stop();
  clock.start();

  // Assign nodes to new colours
  for (HighsInt v : side.node_touched_by_colour.list(split_colour)) {
    const HighsInt new_colour = sum_map_.find(side.sums[v])->second.colour;
    if (new_colour != split_colour) {
      side.node_by_colour.remove(v, split_colour);
      side.node_by_colour.append(v, new_colour);
      side.colour[v] = new_colour;
    }
  }

  time_new_colour_ += clock.stop();
}

void ColourRefinement::prepareNextIter(Side& side) {
  hipo::Clock clock;

  while (!side.colours_touched.empty()) {
    const HighsInt c = side.colours_touched.pop();
    for (HighsInt v : side.node_touched_by_colour.list(c)) {
      side.sums[v] = 0.0;
      side.node_touched[v] = 0;
    }
    side.node_touched_by_colour.clear(c);
  }

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
    if (!rows_.to_refine.empty())
      refine(rows_, cols_, At_);
    else
      refine(cols_, rows_, A_);
  }

  printTimes();
}

void ColourRefinement::printTimes() const {
  const double time_total = time_sums_ + time_check_ + time_insert_ +
                            time_add_stack_ + time_new_colour_ + time_prepare_;

  printf("ColourRefinement timers\n");
  printf("Total         %f\n", time_total);
  printf("  sums        %f\n", time_sums_);
  printf("  check       %f\n", time_check_);
  printf("  insert      %f\n", time_insert_);
  printf("  add stack   %f\n", time_add_stack_);
  printf("  new colour  %f\n", time_new_colour_);
  printf("  prepare     %f\n", time_prepare_);
}

}  // namespace folding
}  // namespace highs