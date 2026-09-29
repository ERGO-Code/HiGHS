#include "Fold.h"

#include <algorithm>
#include <map>

#include "DoubleQuantizer.h"
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
      max_colour_degree_(n_, 0),
      min_colour_degree_(n_, 0),
      colours_touched_(n_, 0),
      in_colours_touched_(n_, 0),
      colours_split_(n_, 0),
      in_stack_(n_, 0) {
  hipo::Clock clock;

  colour_classes_.init(n_, n_);
  colour_classes_touched_.init(n_, n_);

  assert(colour_.size() == n_);

  for (HighsInt i = 0; i < n_; ++i) colour_classes_.append(i, colour_[i]);

  latest_colour_ = *std::max_element(colour_.begin(), colour_.end());
  for (HighsInt i = 0; i <= latest_colour_; ++i) {
    stack_refine_.push(i);
    in_stack_[i] = 1;
  }

  At_ = A_;
  At_.ensureRowwise();

  time_setup_ = clock.stop();
}

void ColourRefinement::chooseRefiningColour() {
  hipo::Clock clock;

  refining_colour_ = stack_refine_.top();
  stack_refine_.pop();
  in_stack_[refining_colour_] = 0;

  time_choose_ += clock.stop();
}

void ColourRefinement::computeColourDegrees() {
  hipo::Clock clock;

  for (HighsInt v : colour_classes_.list(refining_colour_)) {
    for (HighsInt w : neighbours(v)) {
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
  }

  for (HighsInt el = 0; el < top_touched_; ++el) {
    const HighsInt c = colours_touched_[el];
    if (colour_classes_.length(c) != colour_classes_touched_.length(c))
      min_colour_degree_[c] = 0;
    else {
      min_colour_degree_[c] = max_colour_degree_[c];

      for (HighsInt v : colour_classes_touched_.list(c)) {
        if (colour_degree_[v] < min_colour_degree_[c])
          min_colour_degree_[c] = colour_degree_[v];
      }
    }
  }

  time_degrees_ += clock.stop();
}

void ColourRefinement::findSplitColours() {
  hipo::Clock clock;

  top_split_ = 0;
  for (HighsInt el = 0; el < top_touched_; ++el) {
    const HighsInt c = colours_touched_[el];
    if (min_colour_degree_[c] < max_colour_degree_[c]) {
      colours_split_[top_split_] = c;
      top_split_++;
    }
  }
  std::sort(colours_split_.begin(), colours_split_.begin() + top_split_);

  time_find_split_ += clock.stop();
}

void ColourRefinement::splitColours() {
  hipo::Clock clock;

  for (HighsInt el = 0; el < top_split_; ++el) {
    const HighsInt s = colours_split_[el];
    splitColour(s);
  }

  time_split_ += clock.stop();
}

void ColourRefinement::splitColour(const HighsInt s) {
  const HighsInt max_degree = max_colour_degree_[s];
  std::vector<HighsInt> degree_count(max_degree + 1, 0);
  degree_count[0] =
      colour_classes_.length(s) - colour_classes_touched_.length(s);

  for (HighsInt v : colour_classes_touched_.list(s)) {
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

  for (HighsInt v : colour_classes_touched_.list(s)) {
    if (new_colour[colour_degree_[v]] != s) {
      colour_classes_.remove(v, s);
      colour_classes_.append(v, new_colour[colour_degree_[v]]);
      colour_[v] = new_colour[colour_degree_[v]];
    }
  }
}

void ColourRefinement::prepareNextIter() {
  hipo::Clock clock;

  for (HighsInt el = 0; el < top_touched_; ++el) {
    const HighsInt c = colours_touched_[el];
    for (HighsInt v : colour_classes_touched_.list(c)) {
      colour_degree_[v] = 0;
    }
    max_colour_degree_[c] = 0;
    colour_classes_touched_.clear(c);
    in_colours_touched_[c] = 0;
  }
  top_touched_ = 0;

  time_prepare_ += clock.stop();
}

void ColourRefinement::run() {
  hipo::Clock clock;

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

void test_folding(const HighsLp& lp) {
  Folder folder(lp);
  std::vector<HighsInt> colour_rows, colour_cols;
  HighsInt colours_used_rows = folder.findInitialRowColour(colour_rows);
  HighsInt colours_used_cols = folder.findInitialColColour(colour_cols);

  printf("Rows: used %d out of %zu\n", colours_used_rows, lp.row_lower_.size());
  printf("Cols: used %d out of %zu\n", colours_used_cols, lp.col_cost_.size());

  for (HighsInt& i : colour_cols) i += colours_used_rows;

  std::vector<HighsInt> colour_matrix = std::move(colour_rows);
  colour_matrix.insert(colour_matrix.end(), colour_cols.begin(),
                       colour_cols.end());

  ColourRefinement CR(lp.a_matrix_, colour_matrix);
  CR.run();

  printf("\nMatrix: used %d out of %d\n\n", CR.coloursUsed(),
         lp.a_matrix_.num_row_ + lp.a_matrix_.num_col_);

  exit(1);
}

template <typename Data>
HighsInt Folder::findInitialColour(HighsInt n, std::vector<HighsInt>& colour) {
  colour.assign(n, 0);
  QuantizedMap<Data, HighsInt> map;
  HighsInt next_colour = 0;

  for (HighsInt i = 0; i < n; ++i) {
    Data data(lp_, i, colour[i]);
    auto result = map.insert({data, next_colour});
    if (result.second) next_colour++;
    colour[i] = result.first->second;
  }

  return next_colour;
}

HighsInt Folder::findInitialRowColour(std::vector<HighsInt>& colour) {
  return findInitialColour<RowData>(lp_.num_row_, colour);
}

HighsInt Folder::findInitialColColour(std::vector<HighsInt>& colour) {
  return findInitialColour<ColData>(lp_.num_col_, colour);
}

}  // namespace folding

}  // namespace highs