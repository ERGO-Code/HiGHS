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

  exit(1);
}

void test_folding(const HighsLp& lp) {
  // test_paper_example();

  Folder folder(lp);
  folder.run();
  folder.print();

  exit(1);
}

Folder::Folder(const HighsLp& lp)
    : lp_{lp}, row_colour_(lp.num_row_), col_colour_(lp.num_col_) {}

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
  initial_row_colours_ =
      findInitialColour<RowQuantizedData>(lp_.num_row_, row_colour_.data(), 0);

  initial_col_colours_ =
      findInitialColour<ColQuantizedData>(lp_.num_col_, col_colour_.data(), 0);
}

void Folder::foldMatrix() {
  CR_.reset(new ColourRefinement(lp_.a_matrix_, row_colour_, col_colour_));
  CR_->run();
}

void Folder::run() {
  hipo::Clock clock;

  findInitialColour();
  foldMatrix();

  printf("Fold:         %f\n", clock.stop());
}

void Folder::print() const {
  printf("Initial:\n");
  printf("Rows: %d out of %d\n", initial_row_colours_, lp_.num_row_);
  printf("Cols: %d out of %d\n", initial_col_colours_, lp_.num_col_);
  printf("Final:\n");
  printf("Rows: %d out of %d\n", CR_->rowColoursUsed(), lp_.num_row_);
  printf("Cols: %d out of %d\n", CR_->colColoursUsed(), lp_.num_col_);
  printf("Matr: %d out of %d\n", CR_->rowColoursUsed() + CR_->colColoursUsed(),
         lp_.num_row_ + lp_.num_col_);

  /*
  for (HighsInt i = 0; i < lp_.num_row_; ++i) printf("%d ", colour_[i]);
  printf("\n");
  for (HighsInt i = 0; i < lp_.num_col_; ++i)
    printf("%d ", colour_[lp_.num_row_ + i]);
  printf("\n");
  */
}

}  // namespace folding

}  // namespace highs