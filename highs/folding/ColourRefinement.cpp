#include "ColourRefinement.h"

#include <algorithm>

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

  At_ = A_;
  At_.ensureRowwise();
}

void ColourRefinement::computeColourSums(HighsInt refining_colour, Side& src,
                                         Side& dst,
                                         const HighsSparseMatrix& M) {
  // Sum of the weights between each dst node and the nodes of refining_colour
  for (HighsInt v : src.node_by_colour.list(refining_colour)) {
    for (HighsInt k = M.start_[v]; k < M.start_[v + 1]; ++k) {
      touchNeighbour(k, M, dst);
    }
  }
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
  const bool has_untouched_nodes = side.node_by_colour.length(colour) >
                                   side.node_touched_by_colour.length(colour);

  const double quantized_reference = DoubleQuantizer::bucket(
      has_untouched_nodes
          ? 0.0
          : side.sums[side.node_touched_by_colour.list(colour).front()]);

  for (HighsInt v : side.node_touched_by_colour.list(colour)) {
    if (DoubleQuantizer::bucket(side.sums[v]) != quantized_reference)
      return true;
  }

  return false;
}

void ColourRefinement::splitColour(Side& side, HighsInt split_colour) {
  // Given a colour touched by the refining colour, and the colour sums,
  // populate the map with the sums corresponding to each vertex of that colour.
  // The map uses the DoubleQuantizer, so sums that are in the same bucket are
  // considered the same.

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

  // Assign nodes to new colours
  for (HighsInt v : side.node_touched_by_colour.list(split_colour)) {
    const HighsInt new_colour = sum_map_.find(side.sums[v])->second.colour;
    if (new_colour != split_colour) {
      side.node_by_colour.remove(v, split_colour);
      side.node_by_colour.append(v, new_colour);
      side.colour[v] = new_colour;
    }
  }
}

void ColourRefinement::prepareNextIter(Side& side) {
  while (!side.colours_touched.empty()) {
    const HighsInt c = side.colours_touched.pop();
    for (HighsInt v : side.node_touched_by_colour.list(c)) {
      side.sums[v] = 0.0;
      side.node_touched[v] = 0;
    }
    side.node_touched_by_colour.clear(c);
  }
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
}

}  // namespace folding
}  // namespace highs