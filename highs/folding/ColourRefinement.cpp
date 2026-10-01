#include "ColourRefinement.h"

#include <algorithm>
#include <map>

#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "util/HighsSparseMatrix.h"

namespace highs {

namespace folding {

/*
Taken from "Tight Lower and Upper Bounds for the Complexity
of Canonical Colour Refinement", Berkholz, Bonsma, Grohe

Matrix A is interpreted as the biadjancency matrix of a bipartite graph.
The initial colouring must be formed of consecutive colours 0...k.

*/

ColourRefinement::ColourRefinement(const HighsSparseMatrix& A,
                                   std::vector<HighsInt>& colour)
    : n_{A.num_row_ + A.num_col_},
      A_{A},
      colour_{colour},
      colour_degree_(n_, 0),
      colour_sums_(n_, 0.0),
      node_touched_(n_, false),
      max_colour_degree_(n_, 0),
      min_colour_degree_(n_, 0),
      stack_(n_),
      colours_touched_(n_) {
  hipo::Clock clock;

  node_by_colour_.init(n_, n_);
  node_touched_by_colour_.init(n_, n_);

  assert(colour_.size() == n_);

  for (HighsInt i = 0; i < n_; ++i) node_by_colour_.append(i, colour_[i]);

  latest_colour_ = *std::max_element(colour_.begin(), colour_.end());
  for (HighsInt i = 0; i <= latest_colour_; ++i) stack_.push(i);

  At_ = A_;
  At_.ensureRowwise();

  colours_split_.reserve(n_);

  time_setup_ = clock.stop();
}

void ColourRefinement::chooseRefiningColour() {
  refining_colour_ = stack_.pop();
}

void ColourRefinement::touchNeighbour(HighsInt w) {
  colour_degree_[w]++;
  if (colour_degree_[w] == 1) node_touched_by_colour_.append(w, colour_[w]);

  colours_touched_.pushIfNotPresent(colour_[w]);

  if (colour_degree_[w] > max_colour_degree_[colour_[w]])
    max_colour_degree_[colour_[w]] = colour_degree_[w];
}

void ColourRefinement::computeColourDegrees() {
  hipo::Clock clock;

  for (HighsInt v : node_by_colour_.list(refining_colour_)) {
    for (Neighbour neigh : neighbours(v)) {
      touchNeighbour(neigh.index);
    }
  }

  for (HighsInt c : colours_touched_) {
    if (node_by_colour_.length(c) != node_touched_by_colour_.length(c))
      min_colour_degree_[c] = 0;
    else {
      min_colour_degree_[c] = max_colour_degree_[c];

      for (HighsInt v : node_touched_by_colour_.list(c)) {
        if (colour_degree_[v] < min_colour_degree_[c])
          min_colour_degree_[c] = colour_degree_[v];
      }
    }
  }

  time_degrees_ += clock.stop();
}

void ColourRefinement::touchNeighbour(const Neighbour& neigh) {
  const HighsInt w = neigh.index;
  colour_sums_[w] += neigh.value;
  if (!node_touched_[w]) {
    node_touched_[w] = true;
    node_touched_by_colour_.append(w, colour_[w]);
  }
  colours_touched_.pushIfNotPresent(colour_[w]);
}

void ColourRefinement::computeColourSums() {
  for (HighsInt v : node_by_colour_.list(refining_colour_)) {
    for (const Neighbour& neigh : neighbours(v)) {
      touchNeighbour(neigh);
    }
  }
}

void ColourRefinement::findSplitColours() {
  hipo::Clock clock;

  colours_split_.clear();
  for (HighsInt c : colours_touched_) {
    if (min_colour_degree_[c] < max_colour_degree_[c]) {
      colours_split_.push_back(c);
    }
  }
  std::sort(colours_split_.begin(), colours_split_.end());

  time_find_split_ += clock.stop();
}

void ColourRefinement::splitColours() {
  hipo::Clock clock;

  for (HighsInt s : colours_split_) splitColour(s);

  time_split_ += clock.stop();
}

void ColourRefinement::splitColour(const HighsInt s) {
  const HighsInt max_degree = max_colour_degree_[s];
  std::vector<HighsInt> degree_count(max_degree + 1, 0);
  degree_count[0] =
      node_by_colour_.length(s) - node_touched_by_colour_.length(s);

  for (HighsInt v : node_touched_by_colour_.list(s)) {
    degree_count[colour_degree_[v]]++;
  }

  HighsInt max_degree_count_index = 0;
  for (HighsInt i = 1; i <= max_degree; ++i) {
    if (degree_count[i] > degree_count[max_degree_count_index])
      max_degree_count_index = i;
  }

  std::vector<HighsInt> new_colour(max_degree + 1, 0);
  for (HighsInt i = 0; i <= max_degree; ++i) {
    if (degree_count[i] >= 1) {
      if (i == min_colour_degree_[s]) {
        new_colour[i] = s;
        if (!stack_.belong(s) && max_degree_count_index != i) {
          stack_.push(new_colour[i]);
        }

      } else {
        latest_colour_++;
        new_colour[i] = latest_colour_;
        if (stack_.belong(s) || i != max_degree_count_index) {
          stack_.push(new_colour[i]);
        }
      }
    }
  }

  for (HighsInt v : node_touched_by_colour_.list(s)) {
    if (new_colour[colour_degree_[v]] != s) {
      node_by_colour_.remove(v, s);
      node_by_colour_.append(v, new_colour[colour_degree_[v]]);
      colour_[v] = new_colour[colour_degree_[v]];
    }
  }
}

void ColourRefinement::prepareNextIter() {
  hipo::Clock clock;

  while (!colours_touched_.empty()) {
    HighsInt c = colours_touched_.pop();
    for (HighsInt v : node_touched_by_colour_.list(c)) {
      colour_degree_[v] = 0;
      colour_sums_[v] = 0.0;
      node_touched_[v] = false;
    }
    max_colour_degree_[c] = 0;
    node_touched_by_colour_.clear(c);
  }

  time_prepare_ += clock.stop();
}

void ColourRefinement::run() {
  hipo::Clock clock;

  while (!stack_.empty()) {
    chooseRefiningColour();
    computeColourDegrees();
    findSplitColours();
    splitColours();
    prepareNextIter();
  }

  printf("ColourRefinement timers\n");
  printf("Total     %f\n", clock.stop());
  printf("  setup   %f\n", time_setup_);
  printf("  degrees %f\n", time_degrees_);
  printf("  find    %f\n", time_find_split_);
  printf("  split   %f\n", time_split_);
  printf("  prepare %f\n", time_prepare_);
}
}  // namespace folding
}  // namespace highs