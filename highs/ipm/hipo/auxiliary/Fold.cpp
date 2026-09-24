#include "Fold.h"

#include <algorithm>
#include <map>

#include "ipm/hipo/auxiliary/Auxiliary.h"
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
                                   std::vector<Int>& colour, bool bipartite)
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
  Clock clock;

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

  time_setup_ = clock.stop();
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
  Clock clock;

  refining_colour_ = stack_refine_.top();
  stack_refine_.pop();
  in_stack_[refining_colour_] = 0;

  time_choose_ += clock.stop();
}

void ColourRefinement::computeColourDegrees() {
  Clock clock;

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

  time_degrees_ += clock.stop();
}

void ColourRefinement::findSplitColours() {
  Clock clock;

  top_split_ = 0;
  for (Int el = 0; el < top_touched_; ++el) {
    const Int c = colours_touched_[el];
    if (min_colour_degree_[c] < max_colour_degree_[c]) {
      colours_split_[top_split_] = c;
      top_split_++;
    }
  }
  std::sort(colours_split_.begin(), colours_split_.begin() + top_split_);

  time_find_split_ += clock.stop();
}

void ColourRefinement::splitColours() {
  Clock clock;

  for (Int el = 0; el < top_split_; ++el) {
    const Int s = colours_split_[el];
    splitColour(s);
  }

  time_split_ += clock.stop();
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
  Clock clock;

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

  time_prepare_ += clock.stop();
}

void ColourRefinement::run() {
  Clock clock;

  while (!stack_refine_.empty()) {
    chooseRefiningColour();
    computeColourDegrees();
    findSplitColours();
    splitColours();
    prepareNextIter();
  }

  printf("ColourRefinement timers\n");
  printf("Total     %f\n", clock.stop());
  printf("  setup   %f\n", time_setup_);
  printf("  choose  %f\n", time_choose_);
  printf("  degrees %f\n", time_degrees_);
  printf("  find    %f\n", time_find_split_);
  printf("  split   %f\n", time_split_);
  printf("  prepare %f\n", time_prepare_);
}

void test_folding() {
  std::vector<double> w0{0, 1, 1, 0, 3, 0, 1, 5, 3, 0, 0, 0, 1};
  std::vector<Int> colour(13, 0);
  ColourRefinementVector CRV(colour.size());
  CRV.run(w0, colour);

  std::vector<double> w1{5, 3, 7, 5, 8, 5, 7, 7, 10, 4, 5, 4, 3};
  CRV.run(w1, colour);

  exit(1);
}

void test_folding(const HighsLp& lp) {
  std::vector<Int> colour_rows(lp.row_lower_.size(), 0);
  ColourRefinementVector CRV_rows(lp.row_lower_.size());
  CRV_rows.run(lp.row_lower_, colour_rows);
  Int colours_used_rows = CRV_rows.run(lp.row_upper_, colour_rows);

  std::vector<Int> colour_cols(lp.col_cost_.size(), 0);
  ColourRefinementVector CRV_cols(lp.col_cost_.size());
  CRV_cols.run(lp.col_cost_, colour_cols);
  CRV_cols.run(lp.col_lower_, colour_cols);
  Int colours_used_cols = CRV_cols.run(lp.col_upper_, colour_cols);

  printf("Rows: used %d out of %zu\n", colours_used_rows, lp.row_lower_.size());
  printf("Cols: used %d out of %zu\n", colours_used_cols, lp.col_cost_.size());

  std::vector<Int> colour(lp.a_matrix_.num_row_, 0);
  colour.insert(colour.end(), lp.a_matrix_.num_col_, 1);
  ColourRefinement CR(lp.a_matrix_, colour, true);
  CR.run();

  if (lp.a_matrix_.num_row_ + lp.a_matrix_.num_col_ < 200) {
    printf("\n\n");
    for (Int c : colour) printf("%d-", c);
  }
  printf("\nUsed %d colours for %d vertices\n\n", CR.coloursUsed(),
         lp.a_matrix_.num_row_ + lp.a_matrix_.num_col_);
}

ColourRefinementVector::ColourRefinementVector(Int n) : n_{n} {
  colour_classes_.init(n_, n_);
}

Int ColourRefinementVector::run(const std::vector<double>& w,
                                std::vector<Int>& colour) {
  assert(w.size() == n_ && colour.size() == n_);

  colour_classes_.clear();
  for (Int i = 0; i < n_; ++i) colour_classes_.append(i, colour[i]);
  Int latest_colour = *std::max_element(colour.begin(), colour.end());
  Int max_initial_colour = latest_colour;

  using ValueColourPair = std::map<double, Int>;
  std::vector<ValueColourPair> info_by_colour;

  for (Int r = 0; r <= max_initial_colour; ++r) {
    info_by_colour.push_back({});
    bool r_used = false;

    Int v = colour_classes_.head(r);
    while (colour_classes_.cont(v)) {
      auto it = info_by_colour[r].find(w[v]);
      if (it == info_by_colour[r].end()) {
        if (r_used) {
          ++latest_colour;
          info_by_colour[r].insert({w[v], latest_colour});
        } else {
          r_used = true;
          info_by_colour[r].insert({w[v], r});
        }
      }

      colour[v] = info_by_colour[r][w[v]];
      v = colour_classes_.next(v);
    }
  }

  // for (double d : w) printf("%7.2f ", d);
  // printf("\n");

  // for (Int i : colour) printf("%8d", i);
  // printf("\n");

  return latest_colour + 1;
}

}  // namespace hipo