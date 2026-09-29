#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include "ColourRefinement.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp);

class Folder {
  const HighsLp& lp_;

  template <typename Data>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);

 public:
  Folder(const HighsLp& lp) : lp_{lp} {}

  HighsInt findInitialRowColour(HighsInt* colour, HighsInt start);
  HighsInt findInitialColColour(HighsInt* colour, HighsInt start);
};

}  // namespace folding

}  // namespace highs

#endif