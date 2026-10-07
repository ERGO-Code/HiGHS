#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include "ColourRefinement.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp);

class Folder {
  const HighsLp& lp_;
  std::vector<HighsInt> row_colour_;
  std::vector<HighsInt> col_colour_;
  std::unique_ptr<ColourRefinement> CR_;

  HighsInt initial_row_colours_;
  HighsInt initial_col_colours_;

  template <typename Data>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);

  void findInitialColour();
  void foldMatrix();

 public:
  Folder(const HighsLp& lp);
  void run();
  void print() const;
};

}  // namespace folding

}  // namespace highs

#endif