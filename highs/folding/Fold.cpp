#include "Fold.h"

#include "DoubleQuantizer.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp) {
  Folder folder(lp);
  std::vector<HighsInt> colour(lp.num_row_ + lp.num_col_);
  HighsInt* colour_rows = colour.data();
  HighsInt* colour_cols = colour.data() + lp.num_row_;

  HighsInt colours_used_rows = folder.findInitialRowColour(colour_rows, 0);
  HighsInt colours_used_cols =
      folder.findInitialColColour(colour_cols, colours_used_rows);

  printf("Rows: used %d out of %zu\n", colours_used_rows, lp.row_lower_.size());
  printf("Cols: used %d out of %zu\n", colours_used_cols, lp.col_cost_.size());

  ColourRefinement CR(lp.a_matrix_, colour);
  CR.run();

  printf("\nMatrix: used %d out of %d\n\n", CR.coloursUsed(),
         lp.a_matrix_.num_row_ + lp.a_matrix_.num_col_);

  exit(1);
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

HighsInt Folder::findInitialRowColour(HighsInt* colour, HighsInt start) {
  return findInitialColour<RowQuantizedData>(lp_.num_row_, colour, start);
}

HighsInt Folder::findInitialColColour(HighsInt* colour, HighsInt start) {
  return findInitialColour<ColQuantizedData>(lp_.num_col_, colour, start);
}

}  // namespace folding

}  // namespace highs