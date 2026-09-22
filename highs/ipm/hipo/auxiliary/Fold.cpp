#include "Fold.h"

#include <algorithm>
#include <stack>

#include "util/HighsSparseMatrix.h"

namespace hipo {

/*
Taken from "Tight Lower and Upper Bounds for the Complexity
of Canonical Colour Refinement", Berkholz, Bonsma, Grohe

Pass a graph adjacency structure in CSC format as matrix A, and the initial
colouring. The colouring must be formed of consecutive colours 0...k.

If the graph is bipartite, set the bipartite flag to true, and set A to be the
structure of the bi-adjacency matrix (rectangular).

If the graph is not bipartite, set the bipartite flag to false, and set A to be
the structure of the adjacency matrix (symmetric).

*/

ColourRefinement::ColourRefinement(const HighsSparseMatrix& A,
                                   const std::vector<Int>& colour,
                                   bool bipartite)
    : n_{bipartite ? A.num_row_ + A.num_col_ : A.num_row_},
      bipartite_{bipartite},
      A_{A},
      colour_{colour},
      colour_degree_(n_, 0),
      max_colour_degree_(n_, 0),
      min_colour_degree_(n_, 0),
      colours_touched_(n_, 0),
      in_colours_touched_(n_, 0),
      colours_split_(n_, 0),
      in_stack_(n_, 0) {
  colour_classes_.init(n_, n_);
  colour_classes_touched_.init(n_, n_);

  assert(colour_.size() == n_);

  for (Int i = 0; i < n_; ++i) colour_classes_.append(i, colour_[i]);

  latest_colour_ = *std::max_element(colour_.begin(), colour_.end());
  for (Int i = 0; i <= latest_colour_; ++i) {
    stack_refine_.push(i);
    in_stack_[i] = 1;
  }

  if (bipartite_) {
    At_ = A_;
    At_.ensureRowwise();
  } else {
    assert(A_.num_row_ == A_.num_col_);
  }
}

void ColourRefinement::forEachNeighbour(Int v,
                                        const std::function<void(int)>& f) {
  if (bipartite_)
    forEachNeighbourBipartite(v, f);
  else
    forEachNeighbourNonBipartite(v, f);
}

void ColourRefinement::forEachNeighbourNonBipartite(
    Int v, const std::function<void(int)>& f) {
  for (Int el = A_.start_[v]; el < A_.start_[v + 1]; ++el) {
    const Int w = A_.index_[el];
    f(w);
  }
}

void ColourRefinement::forEachNeighbourBipartite(
    Int v, const std::function<void(int)>& f) {
  if (v < A_.num_row_) {
    for (Int el = At_.start_[v]; el < At_.start_[v + 1]; ++el) {
      const Int w = At_.index_[el] + A_.num_row_;
      f(w);
    }
  } else {
    for (Int el = A_.start_[v - A_.num_row_];
         el < A_.start_[v - A_.num_row_ + 1]; ++el) {
      const Int w = A_.index_[el];
      f(w);
    }
  }
}

void ColourRefinement::chooseRefiningColour() {
  refining_colour_ = stack_refine_.top();
  stack_refine_.pop();
  in_stack_[refining_colour_] = 0;
}

void ColourRefinement::computeColourDegrees() {
  Int v = colour_classes_.head(refining_colour_);
  while (colour_classes_.cont(v)) {
    forEachNeighbour(v, [this](Int w) {
      colour_degree_[w]++;
      if (colour_degree_[w] == 1) colour_classes_touched_.append(w, colour_[w]);

      if (!in_colours_touched_[colour_[w]]) {
        colours_touched_[top_touched_] = colour_[w];
        top_touched_++;
        in_colours_touched_[colour_[w]] = 1;
      }

      if (colour_degree_[w] > max_colour_degree_[colour_[w]])
        max_colour_degree_[colour_[w]] = colour_degree_[w];
    });

    v = colour_classes_.next(v);
  }

  for (Int el = 0; el < top_touched_; ++el) {
    const Int c = colours_touched_[el];
    if (colour_classes_.length(c) != colour_classes_touched_.length(c))
      min_colour_degree_[c] = 0;
    else {
      min_colour_degree_[c] = max_colour_degree_[c];
      Int v = colour_classes_touched_.head(c);
      while (colour_classes_touched_.cont(v)) {
        if (colour_degree_[v] < min_colour_degree_[c])
          min_colour_degree_[c] = colour_degree_[v];
        v = colour_classes_touched_.next(v);
      }
    }
  }
}

void ColourRefinement::findSplitColours() {
  top_split_ = 0;
  for (Int el = 0; el < top_touched_; ++el) {
    const Int c = colours_touched_[el];
    if (min_colour_degree_[c] < max_colour_degree_[c]) {
      colours_split_[top_split_] = c;
      top_split_++;
    }
  }
  std::sort(colours_split_.begin(), colours_split_.begin() + top_split_);
}

void ColourRefinement::splitColours() {
  for (Int el = 0; el < top_split_; ++el) {
    const Int s = colours_split_[el];
    splitColour(s);
  }
}

void ColourRefinement::splitColour(const Int s) {
  const Int max_degree = max_colour_degree_[s];
  std::vector<Int> degree_count(max_degree + 1, 0);
  degree_count[0] =
      colour_classes_.length(s) - colour_classes_touched_.length(s);

  Int v = colour_classes_touched_.head(s);
  while (colour_classes_touched_.cont(v)) {
    degree_count[colour_degree_[v]]++;
    v = colour_classes_touched_.next(v);
  }

  Int max_degree_count_index = 0;
  for (Int i = 1; i <= max_degree; ++i) {
    if (degree_count[i] > degree_count[max_degree_count_index])
      max_degree_count_index = i;
  }

  std::vector<Int> new_colour(max_degree + 1, 0);
  for (Int i = 0; i <= max_degree; ++i) {
    if (degree_count[i] >= 1) {
      if (i == min_colour_degree_[s]) {
        new_colour[i] = s;
        if (!in_stack_[s] && max_degree_count_index != i) {
          stack_refine_.push(new_colour[i]);
          in_stack_[new_colour[i]] = 1;
        }

      } else {
        latest_colour_++;
        new_colour[i] = latest_colour_;
        if (in_stack_[s] || i != max_degree_count_index) {
          stack_refine_.push(new_colour[i]);
          in_stack_[new_colour[i]] = 1;
        }
      }
    }
  }

  v = colour_classes_touched_.head(s);
  while (colour_classes_touched_.cont(v)) {
    if (new_colour[colour_degree_[v]] != s) {
      colour_classes_.remove(v, s);
      colour_classes_.append(v, new_colour[colour_degree_[v]]);
      colour_[v] = new_colour[colour_degree_[v]];
    }
    v = colour_classes_touched_.next(v);
  }
}

void ColourRefinement::prepareNextIter() {
  for (Int el = 0; el < top_touched_; ++el) {
    const Int c = colours_touched_[el];
    Int v = colour_classes_touched_.head(c);
    while (colour_classes_touched_.cont(v)) {
      colour_degree_[v] = 0;
      v = colour_classes_touched_.next(v);
    }
    max_colour_degree_[c] = 0;
    colour_classes_touched_.clear(c);
    in_colours_touched_[c] = 0;
  }
  top_touched_ = 0;
}

void ColourRefinement::run() {
  while (!stack_refine_.empty()) {
    chooseRefiningColour();
    computeColourDegrees();
    findSplitColours();
    splitColours();
    prepareNextIter();
  }
}

const std::vector<Int>& ColourRefinement::getColour() const { return colour_; }

void test_folding() {
  HighsSparseMatrix A1;
  A1.start_ = {0, 3, 7, 11, 15, 18, 21, 24, 27, 30};
  A1.index_ = {4, 5, 6, 4, 5, 7, 8, 4, 6, 7, 8, 5, 6, 7, 8,
               0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3, 1, 2, 3};
  A1.value_.resize(A1.index_.size());
  A1.num_row_ = A1.start_.size() - 1;
  A1.num_col_ = A1.start_.size() - 1;
  std::vector<Int> initial_colour1(A1.start_.size() - 1, 0);

  ColourRefinement CR1(A1, initial_colour1, false);
  CR1.run();
  const std::vector<Int> colour1 = CR1.getColour();

  printf("\n\n");
  for (Int c : colour1) printf("%d", c);
  printf("\n");

  //
  // bipartite

  HighsSparseMatrix A2;
  A2.start_ = {0, 3, 6, 9, 12, 15};
  A2.index_ = {0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3, 1, 2, 3};
  A2.value_.resize(A2.index_.size());
  A2.num_row_ = *std::max_element(A2.index_.begin(), A2.index_.end()) + 1;
  A2.num_col_ = A2.start_.size() - 1;
  std::vector<Int> initial_colour2(A2.num_row_, 0);
  initial_colour2.insert(initial_colour2.end(), A2.num_col_, 1);

  ColourRefinement CR2(A2, initial_colour2, true);
  CR2.run();
  const std::vector<Int> colour2 = CR2.getColour();

  printf("\n\n");
  for (Int c : colour2) printf("%d", c);
  printf("\n");
}

}  // namespace hipo