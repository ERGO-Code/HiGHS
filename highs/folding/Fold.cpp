#include "Fold.h"

#include "DoubleQuantizer.h"

namespace highs {

namespace folding {

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
  return findInitialColour<RowQuantizedData>(lp_.num_row_, colour);
}

HighsInt Folder::findInitialColColour(std::vector<HighsInt>& colour) {
  return findInitialColour<ColQuantizedData>(lp_.num_col_, colour);
}

}  // namespace folding

}  // namespace highs