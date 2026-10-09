#include "Fold.h"

#include "DoubleQuantizer.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace highs {

namespace folding {

void test_paper_example() {
  HighsLp lp;

  lp.col_cost_ = {2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1, 1, .5, .5, .5, .5};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
                   kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
                   kHighsInf, kHighsInf, kHighsInf};
  lp.row_lower_ = {1, 1, 1, 1, 1, 1, 1};
  lp.row_upper_ = {1, 1, 1, 1, 1, 1, 1};

  lp.a_matrix_.start_ = {0, 6, 13, 19, 24, 29, 34, 39, 41, 43, 48, 53, 58, 63};

  lp.a_matrix_.index_ = {0, 1, 2, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2,
                         3, 4, 5, 0, 1, 2, 4, 6, 0, 1, 2, 3, 5, 0, 1, 2,
                         4, 6, 0, 1, 2, 3, 5, 3, 4, 5, 6, 0, 1, 2, 3, 5,
                         0, 1, 2, 4, 6, 0, 1, 2, 3, 5, 0, 1, 2, 4, 6};

  lp.a_matrix_.value_ = {3,   -1,  1,       1.0 / 3, 1.0 / 3, 2.0 / 3, -1,
                         1,   3,   1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1,
                         3,   -1,  2.0 / 3, 1.0 / 3, 1.0 / 3, .25,     .25,
                         .25, 1.5, 1.5,     .25,     .25,     .25,     1.5,
                         1.5, .25, .25,     .25,     1.5,     1.5,     .25,
                         .25, .25, 1.5,     1.5,     2,       2,       2,
                         2,   3,   -2,      .5,      1,       -1,      -2,
                         3,   .5,  1,       -1,      .5,      .5,      .5,
                         -1,  1,   .5,      .5,      .5,      -1,      1};

  lp.num_row_ = 7;
  lp.num_col_ = 13;

  lp.a_matrix_.num_row_ = 7;
  lp.a_matrix_.num_col_ = 13;

  LpFolder folder(lp);
  folder.run();
  folder.print();

  HighsLp folded_lp = folder.getFoldedLp();
  LpFolder folder2(folded_lp);
  folder2.run();
  folder2.print();

  exit(1);
}

void test_folding(const HighsLp& lp) {
  // test_paper_example();

  LpFolder folder(lp);
  folder.run();
  folder.print();

  if (!folder.isPartitionCorrect()) {
    printf("Check failed\n");
  }

  exit(1);
}

LpFolder::LpFolder(const HighsLp& lp)
    : lp_{lp},
      row_colour_(lp.num_row_),
      col_colour_(lp.num_col_),
      A_{lp_.a_matrix_} {
  assert(A_.isColwise());
  A_.buildOppositeFormat(At_);
}

template <typename Data>
HighsInt LpFolder::findInitialColour(HighsInt n, HighsInt* colour,
                                     HighsInt start) {
  QuantizedMap<Data, HighsInt> map;
  HighsInt next_colour = start;

  for (HighsInt i = 0; i < n; ++i) {
    Data data(lp_, i);
    auto result = map.insert({data, next_colour});
    if (result.second) next_colour++;
    colour[i] = result.first->second;
  }

  return next_colour - start;
}

void LpFolder::findInitialColour() {
  hipo::Clock clock;

  initial_row_colours_ =
      findInitialColour<RowQuantizedData>(lp_.num_row_, row_colour_.data(), 0);

  initial_col_colours_ =
      findInitialColour<ColQuantizedData>(lp_.num_col_, col_colour_.data(), 0);

  initial_time_ = clock.stop();
}

void LpFolder::foldMatrix() {
  hipo::Clock clock;
  ColourRefinement CR(A_, At_, row_colour_, col_colour_);
  ctor_time_ = clock.stop();

  clock.start();
  HighsInt status = CR.run();
  run_time_ = clock.stop();

  if (status) {
    fold_error_ = 1;
    printf(" === Folding failed === \n");
  } else {
    fold_error_ = 0;

    num_row_colours_ = CR.rowColoursUsed();
    num_col_colours_ = CR.colColoursUsed();

    row_colours_size_.assign(num_row_colours_, 0);
    col_colours_size_.assign(num_col_colours_, 0);
    for (HighsInt row = 0; row < lp_.num_row_; ++row)
      ++row_colours_size_[row_colour_[row]];
    for (HighsInt col = 0; col < lp_.num_col_; ++col)
      ++col_colours_size_[col_colour_[col]];
  }
}

HighsInt LpFolder::run() {
  hipo::Clock clock;

  findInitialColour();
  foldMatrix();
  assert(fold_error_ || isPartitionCorrect());

  fold_time_ = clock.stop();

  return fold_error_;
}

HighsLp LpFolder::getFoldedLp() const {
  // Given the partition of rows and columns produced by colour refinement,
  // generate a folded LP with num_row equal to the number of row colours and
  // num_col equal to the number of col colours.
  // For each row colour and each col colour, the folded lp has (up to the
  // quantization tolerance):
  // - row bounds equal to the bounds of any row of that colour.
  // - col bounds equal to the bounds of any col of that colour.
  // - col cost equal to the sum of the costs of all cols of that colour.
  // - matrix entry (i, j) equal to the sum, over the columns of col colour j,
  //   of the entries of any row of row colour i. Entries that quantize to zero
  //   are not stored.
  // To achieve this efficiently, the folded matrix is built row-wise, from one
  // representative row per row colour (the one with the fewest nonzeros) and
  // it is then converted to col-wise.

  hipo::Clock clock;

  if (fold_error_) return HighsLp{};

  // Representative row for each colour:
  // For each colour, find the row with the least number of nonzero entries
  auto nz_in_row = [&](HighsInt row) {
    return At_.start_[row + 1] - At_.start_[row];
  };

  std::vector<HighsInt> row_of_colour(num_row_colours_, -1);
  for (HighsInt row = 0; row < A_.num_row_; ++row) {
    const HighsInt colour = row_colour_[row];
    const HighsInt current_best_row = row_of_colour[colour];
    if (current_best_row < 0 || nz_in_row(row) < nz_in_row(current_best_row))
      row_of_colour[colour] = row;
  }

  for (HighsInt colour = 0; colour < num_row_colours_; ++colour)
    assert(row_of_colour[colour] >= 0);

  // Representative col for each colour:
  // For each colour, find the first column
  std::vector<HighsInt> col_of_colour(num_col_colours_, -1);
  for (HighsInt col = 0; col < A_.num_col_; ++col) {
    const HighsInt colour = col_colour_[col];
    if (col_of_colour[colour] < 0) col_of_colour[colour] = col;
  }

  for (HighsInt colour = 0; colour < num_col_colours_; ++colour)
    assert(col_of_colour[colour] >= 0);

  // Build the folded lp
  HighsLp folded_lp;
  folded_lp.sense_ = lp_.sense_;
  folded_lp.offset_ = lp_.offset_;

  folded_lp.num_row_ = num_row_colours_;
  folded_lp.num_col_ = num_col_colours_;

  // Row bounds: all rows of a colour have the same bounds, up to tolerance, so
  // use the representative
  folded_lp.row_lower_.resize(num_row_colours_);
  folded_lp.row_upper_.resize(num_row_colours_);
  for (HighsInt colour = 0; colour < num_row_colours_; ++colour) {
    const HighsInt row = row_of_colour[colour];
    folded_lp.row_lower_[colour] = lp_.row_lower_[row];
    folded_lp.row_upper_[colour] = lp_.row_upper_[row];
  }

  // Col bounds: all cols of a colour have the same bounds, up to tolerance, so
  // use the representative
  folded_lp.col_lower_.resize(num_col_colours_);
  folded_lp.col_upper_.resize(num_col_colours_);
  for (HighsInt colour = 0; colour < num_col_colours_; ++colour) {
    const HighsInt col = col_of_colour[colour];
    folded_lp.col_lower_[colour] = lp_.col_lower_[col];
    folded_lp.col_upper_[colour] = lp_.col_upper_[col];
  }

  // Col costs: costs are summed over all variables with the same colour
  folded_lp.col_cost_.assign(num_col_colours_, 0.0);
  for (HighsInt col = 0; col < lp_.num_col_; ++col) {
    const HighsInt colour = col_colour_[col];
    folded_lp.col_cost_[colour] += lp_.col_cost_[col];
  }

  // Build the folded At
  HighsSparseMatrix Aft;
  Aft.num_row_ = num_row_colours_;
  Aft.num_col_ = num_col_colours_;
  Aft.format_ = MatrixFormat::kRowwise;
  Aft.start_.reserve(num_row_colours_ + 1);
  // Aft.start_ is already initialised with a leading zero
  assert(Aft.start_.size() == 1 && Aft.start_[0] == 0);

  // For each row colour, accumulate the sum of entries with the same col colour
  std::vector<double> dense_row(num_col_colours_, 0.0);
  std::vector<HighsBool> is_col_used(num_col_colours_, false);
  std::vector<HighsInt> col_used;
  for (HighsInt r_colour = 0; r_colour < num_row_colours_; ++r_colour) {
    const HighsInt row = row_of_colour[r_colour];

    for (HighsInt el = At_.start_[row]; el < At_.start_[row + 1]; ++el) {
      const HighsInt col = At_.index_[el];
      const HighsInt c_colour = col_colour_[col];

      if (!is_col_used[c_colour]) {
        is_col_used[c_colour] = true;
        col_used.push_back(c_colour);
      }
      dense_row[c_colour] += At_.value_[el];
    }

    for (HighsInt c_colour : col_used) {
      const double val = dense_row[c_colour];

      // Skip sums that are quantized to zero
      if (DoubleQuantizer::bucket(val) != 0.0) {
        Aft.index_.push_back(c_colour);
        Aft.value_.push_back(val);
      }

      dense_row[c_colour] = 0.0;
      is_col_used[c_colour] = false;
    }
    col_used.clear();

    Aft.start_.push_back(static_cast<HighsInt>(Aft.index_.size()));
  }

  Aft.buildOppositeFormat(folded_lp.a_matrix_);

  folded_lp_time_ = clock.stop();

  return folded_lp;
}

void LpFolder::print() const {
  printf("\nInitial partition:\n");
  printf("\tRows: %9d out of %9d\n", initial_row_colours_, lp_.num_row_);
  printf("\tCols: %9d out of %9d\n", initial_col_colours_, lp_.num_col_);
  printf("Final:\n");
  printf("\tRows: %9d out of %9d (%.1f%%)\n", num_row_colours_, lp_.num_row_,
         (double)num_row_colours_ / lp_.num_row_ * 100);
  printf("\tCols: %9d out of %9d (%.1f%%)\n", num_col_colours_, lp_.num_col_,
         (double)num_col_colours_ / lp_.num_col_ * 100);
  printf("\tMatr: %9d out of %9d\n", num_row_colours_ + num_col_colours_,
         lp_.num_row_ + lp_.num_col_);

  printf("\n");
  printf("Fold time        %f\n", fold_time_);
  printf("Initial          %f\n", initial_time_);
  printf("CR ctor          %f\n", ctor_time_);
  printf("CR run           %f\n", run_time_);
  printf("Folded lp        %f\n", folded_lp_time_);
}

class ColourClasses {
  std::vector<HighsInt> start_;
  std::vector<HighsInt> nodes_;

 public:
  ColourClasses(const std::vector<HighsInt>& colour, HighsInt num_colours)
      : start_(num_colours + 1, 0), nodes_(colour.size()) {
    for (HighsInt c : colour) ++start_[c + 1];
    for (HighsInt c = 0; c < num_colours; ++c) start_[c + 1] += start_[c];

    std::vector<HighsInt> next(start_.begin(), start_.end() - 1);
    for (HighsInt node = 0; node < static_cast<HighsInt>(colour.size()); ++node)
      nodes_[next[colour[node]]++] = node;
  }
  HighsInt size(HighsInt c) const { return start_[c + 1] - start_[c]; }
  HighsInt node(HighsInt c, HighsInt i) const { return nodes_[start_[c] + i]; }
};

bool LpFolder::isPartitionCorrect() const {
  if (fold_error_) return false;

  for (HighsInt c : row_colour_)
    if (c < 0 || c >= num_row_colours_) return false;
  for (HighsInt c : col_colour_)
    if (c < 0 || c >= num_col_colours_) return false;

  // Check that bounds and costs are partitioned correctly
  std::vector<double> reference;
  std::vector<HighsBool> touched;
  auto is_vector_partition_correct = [&](const std::vector<double>& v,
                                         const std::vector<HighsInt>& c,
                                         HighsInt num_colours) -> bool {
    assert(v.size() == c.size());
    const HighsInt n = static_cast<HighsInt>(v.size());
    reference.assign(num_colours, 0.0);
    touched.assign(num_colours, false);

    for (HighsInt i = 0; i < n; ++i) {
      const HighsInt colour = c[i];
      const double val = DoubleQuantizer::bucket(v[i]);
      if (!touched[colour]) {
        reference[colour] = val;
        touched[colour] = true;
      } else if (val != reference[colour]) {
        return false;
      }
    }
    return true;
  };

  // check that partition of A is equitable
  auto is_matrix_partition_correct =
      [](const HighsSparseMatrix& M, HighsInt num_major_colours,
         HighsInt num_minor_colours, const std::vector<HighsInt>& major_colours,
         const std::vector<HighsInt>& minor_colours) -> bool {
    ColourClasses classes(major_colours, num_major_colours);

    std::vector<double> sum(num_minor_colours, 0.0);
    std::vector<HighsBool> is_touched(num_minor_colours, false);
    std::vector<HighsInt> touched;

    std::vector<double> reference(num_minor_colours, 0.0);
    std::vector<HighsInt> reference_touched;

    for (HighsInt major_colour = 0; major_colour < num_major_colours;
         ++major_colour) {
      for (HighsInt i = 0; i < classes.size(major_colour); ++i) {
        const HighsInt major = classes.node(major_colour, i);

        for (HighsInt el = M.start_[major]; el < M.start_[major + 1]; ++el) {
          const HighsInt minor = M.index_[el];
          const HighsInt minor_colour = minor_colours[minor];

          if (!is_touched[minor_colour]) {
            is_touched[minor_colour] = true;
            touched.push_back(minor_colour);
          }

          sum[minor_colour] += M.value_[el];
        }

        if (i == 0) {
          // The first node of each colour defines the reference
          for (HighsInt j : touched) {
            const double q = DoubleQuantizer::bucket(sum[j]);
            if (q != 0.0) {
              reference[j] = q;
              reference_touched.push_back(j);
            }
          }
        } else {
          // The other nodes compare to the reference
          for (HighsInt j : touched)
            if (DoubleQuantizer::bucket(sum[j]) != reference[j]) return false;
          for (HighsInt j : reference_touched)
            if (DoubleQuantizer::bucket(sum[j]) != reference[j]) return false;
        }

        for (HighsInt j : touched) {
          is_touched[j] = false;
          sum[j] = 0.0;
        }
        touched.clear();
      }

      for (HighsInt j : reference_touched) reference[j] = 0.0;
      reference_touched.clear();
    }

    return true;
  };

  if (!is_vector_partition_correct(lp_.row_lower_, row_colour_,
                                   num_row_colours_) ||
      !is_vector_partition_correct(lp_.row_upper_, row_colour_,
                                   num_row_colours_) ||
      !is_vector_partition_correct(lp_.col_cost_, col_colour_,
                                   num_col_colours_) ||
      !is_vector_partition_correct(lp_.col_lower_, col_colour_,
                                   num_col_colours_) ||
      !is_vector_partition_correct(lp_.col_upper_, col_colour_,
                                   num_col_colours_) ||
      !is_matrix_partition_correct(At_, num_row_colours_, num_col_colours_,
                                   row_colour_, col_colour_) ||
      !is_matrix_partition_correct(A_, num_col_colours_, num_row_colours_,
                                   col_colour_, row_colour_))
    return false;

  return true;
}

HighsSolution LpFolder::unfold(const HighsSolution& folded) const {
  // Given partition matrices C and D for row and column partition respectively,
  // and given the operation of scaling the columns and taking the transpose
  // denoted as .^S, the mapping from folded solution (x',y',z') to original
  // solution (x,y,z) is:
  // - col value:   x  = D      x'
  // - col dual :   z  = D^{ST} z'
  // - row value:   Ax = C      A'x'
  // - row dual :   y  = C^{ST} y'

  HighsSolution sol;
  if (fold_error_) return sol;

  if (folded.value_valid) {
    assert(folded.col_value.size() == static_cast<size_t>(num_col_colours_));
    assert(folded.row_value.size() == static_cast<size_t>(num_row_colours_));
    sol.col_value.resize(lp_.num_col_);
    sol.row_value.resize(lp_.num_row_);
    for (HighsInt c = 0; c < lp_.num_col_; ++c) {
      const HighsInt c_colour = col_colour_[c];
      sol.col_value[c] = folded.col_value[c_colour];
    }
    for (HighsInt r = 0; r < lp_.num_row_; ++r) {
      const HighsInt r_colour = row_colour_[r];
      sol.row_value[r] = folded.row_value[r_colour];
    }
  }

  if (folded.dual_valid) {
    assert(folded.col_dual.size() == static_cast<size_t>(num_col_colours_));
    assert(folded.row_dual.size() == static_cast<size_t>(num_row_colours_));
    sol.col_dual.resize(lp_.num_col_);
    sol.row_dual.resize(lp_.num_row_);
    for (HighsInt c = 0; c < lp_.num_col_; ++c) {
      const HighsInt c_colour = col_colour_[c];
      sol.col_dual[c] = folded.col_dual[c_colour] / col_colours_size_[c_colour];
    }
    for (HighsInt r = 0; r < lp_.num_row_; ++r) {
      const HighsInt r_colour = row_colour_[r];
      sol.row_dual[r] = folded.row_dual[r_colour] / row_colours_size_[r_colour];
    }
  }

  sol.value_valid = folded.value_valid;
  sol.dual_valid = folded.dual_valid;

  return sol;
}

}  // namespace folding

}  // namespace highs