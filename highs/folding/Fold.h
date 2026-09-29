#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include "ColourRefinement.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp);

class Folder {
  const HighsLp& lp_;

  template <typename Data>
  HighsInt findInitialColour(HighsInt num, std::vector<HighsInt>& colour);

 public:
  Folder(const HighsLp& lp) : lp_{lp} {}

  HighsInt findInitialRowColour(std::vector<HighsInt>& colour);
  HighsInt findInitialColColour(std::vector<HighsInt>& colour);
};

}  // namespace folding

}  // namespace highs

#endif