#include "Fold.h"

#include "DoubleQuantizer.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp) {
  Folder folder(lp);
  folder.run();
  folder.print();

  exit(1);
}

Folder::Folder(const HighsLp& lp) : lp_{lp} {
  colour_.resize(lp_.num_row_ + lp_.num_col_);
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
  initial_row_colours_ =
      findInitialColour<RowQuantizedData>(lp_.num_row_, colour_.data(), 0);

  initial_col_colours_ = findInitialColour<ColQuantizedData>(
      lp_.num_col_, colour_.data() + lp_.num_row_, initial_row_colours_);
}

void Folder::foldMatrix() {
  ColourRefinement CR(lp_.a_matrix_, colour_);
  CR.runNew();
  matrix_colours_ = CR.coloursUsed();
}

void Folder::run() {
  findInitialColour();
  foldMatrix();
}

void Folder::print() const {
  printf("Rows: %d out of %d\n", initial_row_colours_, lp_.num_row_);
  printf("Cols: %d out of %d\n", initial_col_colours_, lp_.num_col_);
  printf("Matr: %d out of %d\n", matrix_colours_, lp_.num_row_ + lp_.num_col_);
}

}  // namespace folding

}  // namespace highs