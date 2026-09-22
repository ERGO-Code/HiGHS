#include "Fold.h"

#include <algorithm>
#include <stack>

namespace hipo {

/*
Taken from "Tight Lower and Upper Bounds for the Complexity
of Canonical Colour Refinement", Berkholz, Bonsma, Grohe

Pass a graph adjacency structure in CSC format (ptr,adj), and the initial
colouring. The colouring must be formed of consecutive colours 0...k.
*/

ColourRefinement::ColourRefinement(const std::vector<Int>& ptr,
                                   const std::vector<Int>& adj,
                                   const std::vector<Int>& colour)
    : ptr_{ptr},
      adj_{adj},
      n_{static_cast<Int>(ptr_.size() - 1)},
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

  for (Int i = 0; i < n_; ++i) colour_classes_.append(i, colour_[i]);

  latest_colour_ = *std::max_element(colour_.begin(), colour_.end());
  for (Int i = 0; i <= latest_colour_; ++i) {
    stack_refine_.push(i);
    in_stack_[i] = 1;
  }
}

void ColourRefinement::chooseRefiningColour() {
  refining_colour_ = stack_refine_.top();
  stack_refine_.pop();
  in_stack_[refining_colour_] = 0;

  printf("\n");
  for (Int c : colour_) printf("%d", c);
  printf("\n");
  printf("Refine with r = %d\n", refining_colour_);
}

void ColourRefinement::computeColourDegrees() {
  Int v = colour_classes_.head(refining_colour_);
  while (colour_classes_.cont(v)) {
    for (Int el = ptr_[v]; el < ptr_[v + 1]; ++el) {
      const Int w = adj_[el];
      colour_degree_[w]++;
      if (colour_degree_[w] == 1) colour_classes_touched_.append(w, colour_[w]);

      if (!in_colours_touched_[colour_[w]]) {
        colours_touched_[top_touched_] = colour_[w];
        top_touched_++;
        in_colours_touched_[colour_[w]] = 1;
      }

      if (colour_degree_[w] > max_colour_degree_[colour_[w]])
        max_colour_degree_[colour_[w]] = colour_degree_[w];
    }
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
  printf("\tSplit %d\n", s);

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
  const std::vector<Int> ptr = {0,  3,  6,  9,  19, 22, 25,
                                28, 31, 34, 37, 40, 41, 42};
  const std::vector<Int> adj = {1, 2, 3, 0, 2,  3, 0, 1, 3, 0, 1,  2,  4,  5,
                                6, 7, 8, 9, 10, 3, 7, 9, 3, 6, 8,  3,  5,  8,
                                3, 4, 9, 3, 5,  6, 3, 4, 7, 3, 11, 12, 10, 10};

  std::vector<Int> initial_colour(ptr.size() - 1, 0);
  initial_colour = {1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0};

  ColourRefinement CR(ptr, adj, initial_colour);
  CR.run();
  const std::vector<Int> colour = CR.getColour();

  printf("\n\n");
  for (Int c : colour) printf("%d", c);
  printf("\n");
}

}  // namespace hipo