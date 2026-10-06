#ifndef HIGHS_COLOUR_REFINEMENT_H
#define HIGHS_COLOUR_REFINEMENT_H

#include <cassert>
#include <vector>

#include "DoubleQuantizer.h"
#include "util/HighsLinkedLists.h"
#include "util/HighsSparseMatrix.h"

namespace highs {
namespace folding {

// Refines separate colourings of the rows and columns of matrix A, seen as the
// biadjacency matrix of a weighted bipartite graph. Row colours and column
// colours are numbered independently, each starting from 0 and consecutive.
// Based on "Tight Lower and Upper Bounds for the Complexity of Canonical Colour
// Refinement", Berkholz, Bonsma, Grohe

class ColourRefinement {
  struct IterableStack {
    std::vector<HighsInt> s;
    HighsInt top;
    std::vector<HighsBool> in_stack;

    IterableStack(HighsInt n) : s(n), top{-1}, in_stack(n, 0) {}
    bool belong(HighsInt i) const { return in_stack[i]; }
    bool empty() const { return top < 0; }
    void push(HighsInt i) {
      assert(!in_stack[i]);
      s[++top] = i;
      in_stack[i] = 1;
    }
    HighsInt pop() {
      assert(!empty());
      const HighsInt elem = s[top--];
      in_stack[elem] = 0;
      return elem;
    }

    // Iterators in reverse pop order
    std::vector<HighsInt>::const_iterator begin() const { return s.begin(); }
    std::vector<HighsInt>::const_iterator end() const {
      return s.begin() + (top + 1);
    }
  };

  // Everything that belongs to one side of the bipartite graph
  struct Side {
    HighsInt n;
    std::vector<HighsInt>& colour;
    HighsInt latest_colour;
    std::vector<double> sums;
    std::vector<HighsBool> node_touched;
    HighsLinkedLists node_by_colour;
    HighsLinkedLists node_touched_by_colour;
    IterableStack to_refine;
    IterableStack colours_touched;

    explicit Side(std::vector<HighsInt>& c);
  };

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;
  Side rows_;
  Side cols_;

  struct SumData {
    HighsInt count;
    HighsInt colour;
  };
  QuantizedMap<double, SumData> sum_map_;

  void refine(Side& src, Side& dst, const HighsSparseMatrix& M);
  void computeColourSums(HighsInt refining_colour, Side& src, Side& dst,
                         const HighsSparseMatrix& M);
  void touchNeighbour(HighsInt k, const HighsSparseMatrix& M, Side& side);
  void splitColours(Side& side);
  bool colourSplits(const Side& side, HighsInt colour) const;
  void splitColour(Side& side, HighsInt colour);
  void prepareNextIter(Side& side);

 public:
  ColourRefinement(const HighsSparseMatrix& A,
                   std::vector<HighsInt>& row_colour,
                   std::vector<HighsInt>& col_colour);
  void run();

  HighsInt rowColoursUsed() const { return rows_.latest_colour + 1; }
  HighsInt colColoursUsed() const { return cols_.latest_colour + 1; }
};

}  // namespace folding
}  // namespace highs

#endif