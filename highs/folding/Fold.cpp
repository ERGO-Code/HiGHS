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

  Folder folder(lp);
  folder.run();
  folder.print();
  folder.printFoldedLp();

  HighsLp folded_lp = folder.getFoldedLp();
  Folder folder2(folded_lp);
  folder2.run();
  folder2.print();
  folder2.printFoldedLp();

  exit(1);
}

void test_folding(const HighsLp& lp) {
  // test_paper_example();

  Folder folder(lp);
  folder.run();
  folder.print();

  if (!folder.checkCorrect()) {
    printf("Check failed\n");
  }

  exit(1);
}

Folder::Folder(const HighsLp& lp)
    : lp_{lp},
      row_colour_(lp.num_row_),
      col_colour_(lp.num_col_),
      A_{lp_.a_matrix_} {
  A_.buildOppositeFormat(At_);
}

template <typename Data>
HighsInt Folder::findInitialColour(HighsInt n, HighsInt* colour,
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

void Folder::findInitialColour() {
  hipo::Clock clock;

  initial_row_colours_ =
      findInitialColour<RowQuantizedData>(lp_.num_row_, row_colour_.data(), 0);

  initial_col_colours_ =
      findInitialColour<ColQuantizedData>(lp_.num_col_, col_colour_.data(), 0);

  initial_time_ = clock.stop();
}

HighsInt Folder::foldMatrix() {
  hipo::Clock clock;
  CR_.reset(new ColourRefinement(A_, At_, row_colour_, col_colour_));
  ctor_time_ = clock.stop();

  clock.start();
  HighsInt status = CR_->run();
  run_time_ = clock.stop();

  if (status) printf(" === Folding failed === \n");
  return status;
}

HighsInt Folder::run() {
  hipo::Clock clock;

  findInitialColour();
  HighsInt status = foldMatrix();
  if (status == 0) fold_success_ = true;

  fold_time_ = clock.stop();

  return status;
}

HighsLp Folder::getFoldedLp() const {
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

  if (!fold_success_) return HighsLp{};

  // Representative row for each colour:
  // For each colour, find the row with the least number of nonzero entries
  auto nz_in_row = [&](HighsInt row) {
    return At_.start_[row + 1] - At_.start_[row];
  };

  std::vector<HighsInt> row_of_colour(CR_->rowColoursUsed(), -1);
  for (HighsInt row = 0; row < A_.num_row_; ++row) {
    const HighsInt colour = row_colour_[row];
    const HighsInt current_best_row = row_of_colour[colour];
    if (current_best_row < 0 || nz_in_row(row) < nz_in_row(current_best_row))
      row_of_colour[colour] = row;
  }

  for (HighsInt colour = 0; colour < CR_->rowColoursUsed(); ++colour)
    assert(row_of_colour[colour] >= 0);

  // Representative col for each colour:
  // For each colour, find the first column
  std::vector<HighsInt> col_of_colour(CR_->colColoursUsed(), -1);
  for (HighsInt col = 0; col < A_.num_col_; ++col) {
    const HighsInt colour = col_colour_[col];
    if (col_of_colour[colour] < 0) col_of_colour[colour] = col;
  }

  for (HighsInt colour = 0; colour < CR_->colColoursUsed(); ++colour)
    assert(col_of_colour[colour] >= 0);

  // Build the folded lp
  HighsLp folded_lp;
  folded_lp.sense_ = lp_.sense_;
  folded_lp.offset_ = lp_.offset_;

  const HighsInt num_row_colours = CR_->rowColoursUsed();
  const HighsInt num_col_colours = CR_->colColoursUsed();
  folded_lp.num_row_ = num_row_colours;
  folded_lp.num_col_ = num_col_colours;

  // Row bounds: they are the average of row bounds of a given colour, so just
  // use the representative
  folded_lp.row_lower_.resize(num_row_colours);
  folded_lp.row_upper_.resize(num_row_colours);
  for (HighsInt colour = 0; colour < num_row_colours; ++colour) {
    const HighsInt row = row_of_colour[colour];
    folded_lp.row_lower_[colour] = lp_.row_lower_[row];
    folded_lp.row_upper_[colour] = lp_.row_upper_[row];
  }

  // Col bounds: they are the average of col bounds of a given colour, so just
  // use the representative
  folded_lp.col_lower_.resize(num_col_colours);
  folded_lp.col_upper_.resize(num_col_colours);
  for (HighsInt colour = 0; colour < num_col_colours; ++colour) {
    const HighsInt col = col_of_colour[colour];
    folded_lp.col_lower_[colour] = lp_.col_lower_[col];
    folded_lp.col_upper_[colour] = lp_.col_upper_[col];
  }

  // Col costs: costs are summed over all variables with the same colour
  folded_lp.col_cost_.assign(num_col_colours, 0.0);
  for (HighsInt col = 0; col < lp_.num_col_; ++col) {
    const HighsInt colour = col_colour_[col];
    folded_lp.col_cost_[colour] += lp_.col_cost_[col];
  }

  // Build the folded At
  HighsSparseMatrix Aft;
  Aft.num_row_ = num_row_colours;
  Aft.num_col_ = num_col_colours;
  Aft.format_ = MatrixFormat::kRowwise;
  Aft.start_.reserve(num_row_colours + 1);
  // Aft.start_ is already initialised with a leading zero

  // For each row colour, accumulate the sum of entries with the same col colour
  std::vector<double> dense_row(num_col_colours, 0.0);
  std::vector<HighsBool> is_col_used(num_col_colours, false);
  std::vector<HighsInt> col_used;
  for (HighsInt r_colour = 0; r_colour < num_row_colours; ++r_colour) {
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

void Folder::print() const {
  printf("\nInitial partition:\n");
  printf("\tRows: %9d out of %9d\n", initial_row_colours_, lp_.num_row_);
  printf("\tCols: %9d out of %9d\n", initial_col_colours_, lp_.num_col_);
  printf("Final:\n");
  printf("\tRows: %9d out of %9d (%.1f%%)\n", CR_->rowColoursUsed(),
         lp_.num_row_, (double)CR_->rowColoursUsed() / lp_.num_row_ * 100);
  printf("\tCols: %9d out of %9d (%.1f%%)\n", CR_->colColoursUsed(),
         lp_.num_col_, (double)CR_->colColoursUsed() / lp_.num_col_ * 100);
  printf("\tMatr: %9d out of %9d\n",
         CR_->rowColoursUsed() + CR_->colColoursUsed(),
         lp_.num_row_ + lp_.num_col_);

  printf("\n");
  printf("Fold time        %f\n", fold_time_);
  printf("Initial          %f\n", initial_time_);
  printf("CR ctor          %f\n", ctor_time_);
  printf("CR run           %f\n", run_time_);
  printf("Folded lp        %f\n", folded_lp_time_);

  /*
  for (HighsInt i = 0; i < lp_.num_row_; ++i) printf("%d ", colour_[i]);
  printf("\n");
  for (HighsInt i = 0; i < lp_.num_col_; ++i)
    printf("%d ", colour_[lp_.num_row_ + i]);
  printf("\n");
  */
}

void Folder::printFoldedLp() const {
  HighsLp flp = getFoldedLp();

  printf("row l: ");
  for (double d : flp.row_lower_) printf("%.1f ", d);
  printf("\n");

  printf("row u: ");
  for (double d : flp.row_upper_) printf("%.1f ", d);
  printf("\n");

  printf("col c: ");
  for (double d : flp.col_cost_) printf("%.1f ", d);
  printf("\n");

  printf("col l: ");
  for (double d : flp.col_lower_) printf("%.1f ", d);
  printf("\n");

  printf("col u: ");
  for (double d : flp.col_upper_) printf("%.1f ", d);
  printf("\n");

  printf("start: ");
  for (HighsInt i : flp.a_matrix_.start_) printf("%d ", i);
  printf("\n");

  printf("index: ");
  for (HighsInt i : flp.a_matrix_.index_) printf("%d ", i);
  printf("\n");

  printf("value: ");
  for (double d : flp.a_matrix_.value_) printf("%.1f ", d);
  printf("\n");
}

bool Folder::checkCorrect() const {
  if (!fold_success_) return false;

  const HighsInt num_row_colours = CR_->rowColoursUsed();
  const HighsInt num_col_colours = CR_->colColoursUsed();
  bool success = true;

  {
    std::vector<double> reference;
    std::vector<HighsBool> seen;

    auto check_vector = [&](const std::vector<double>& v,
                            const std::vector<HighsInt>& c,
                            HighsInt num_colours) -> bool {
      assert(v.size() == c.size());
      const HighsInt n = static_cast<HighsInt>(v.size());
      reference.assign(num_colours, 0.0);
      seen.assign(num_colours, false);

      for (HighsInt i = 0; i < n; ++i) {
        const HighsInt colour = c[i];
        const double val = DoubleQuantizer::bucket(v[i]);
        if (!seen[colour]) {
          reference[colour] = val;
          seen[colour] = true;
        } else if (val != reference[colour]) {
          return false;
        }
      }
      return true;
    };

    success = check_vector(lp_.row_lower_, row_colour_, num_row_colours) &&
              check_vector(lp_.row_upper_, row_colour_, num_row_colours) &&
              check_vector(lp_.col_cost_, col_colour_, num_col_colours) &&
              check_vector(lp_.col_lower_, col_colour_, num_col_colours) &&
              check_vector(lp_.col_upper_, col_colour_, num_col_colours);
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
      for (HighsInt node = 0; node < static_cast<HighsInt>(colour.size());
           ++node)
        nodes_[next[colour[node]]++] = node;
    }
    HighsInt size(HighsInt c) const { return start_[c + 1] - start_[c]; }
    HighsInt node(HighsInt c, HighsInt i) const {
      return nodes_[start_[c] + i];
    }
  };

  std::vector<double> reference;
  std::vector<double> sum;

  auto checkMatrix = [&](const HighsSparseMatrix& M, HighsInt num_major_colours,
                         HighsInt num_minor_colours,
                         const std::vector<HighsInt>& major_colour,
                         const std::vector<HighsInt>& minor_colour) -> bool {
    ColourClasses classes(major_colour, num_major_colours);
    for (HighsInt r_colour = 0; r_colour < num_major_colours; ++r_colour) {
      reference.assign(num_minor_colours, 0.0);

      for (HighsInt node = 0; node < classes.size(r_colour); ++node) {
        sum.assign(num_minor_colours, 0.0);

        const HighsInt row = classes.node(r_colour, node);
        for (HighsInt el = M.start_[row]; el < M.start_[row + 1]; ++el) {
          const HighsInt col = M.index_[el];
          const HighsInt c_colour = minor_colour[col];
          if (node == 0)
            reference[c_colour] += M.value_[el];
          else
            sum[c_colour] += M.value_[el];
        }

        if (node > 0) {
          for (HighsInt i = 0; i < num_minor_colours; ++i)
            if (DoubleQuantizer::bucket(reference[i]) !=
                DoubleQuantizer::bucket(sum[i]))
              return false;
        }
      }
    }

    return true;
  };

  success = success &&
            checkMatrix(At_, num_row_colours, num_col_colours, row_colour_,
                        col_colour_) &&
            checkMatrix(A_, num_col_colours, num_row_colours, col_colour_,
                        row_colour_);

  return success;
}

}  // namespace folding

}  // namespace highs