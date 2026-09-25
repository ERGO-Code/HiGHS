#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include <functional>
#include <stack>

#include "CollectionLinkedLists.h"
#include "lp_data/HighsLp.h"

namespace highs {

namespace folding {

class ColourRefinement {
  const HighsInt n_;

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;

  std::vector<HighsInt>& colour_;

  std::vector<HighsInt> colour_degree_;
  std::vector<HighsInt> max_colour_degree_;
  std::vector<HighsInt> min_colour_degree_;

  std::vector<HighsInt> colours_touched_;
  std::vector<HighsBool> in_colours_touched_;
  HighsInt top_touched_{};

  std::vector<HighsInt> colours_split_;
  HighsInt top_split_{};

  CollectionLinkedLists colour_classes_;
  CollectionLinkedLists colour_classes_touched_;

  HighsInt latest_colour_;
  HighsInt refining_colour_;

  std::stack<HighsInt> stack_refine_;
  std::vector<HighsBool> in_stack_;

  double time_setup_{};
  double time_choose_{};
  double time_degrees_{};
  double time_find_split_{};
  double time_split_{};
  double time_prepare_{};

  void chooseRefiningColour();
  void computeColourDegrees();
  void findSplitColours();
  void splitColours();
  void splitColour(HighsInt split_colour);
  void prepareNextIter();

  void forEachNeighbour(HighsInt v, const std::function<void(int)>& f);

 public:
  ColourRefinement(const HighsSparseMatrix& A, std::vector<HighsInt>& colour);
  void run();
  HighsInt coloursUsed() const { return latest_colour_ + 1; }
};

void test_folding();
void test_folding(const HighsLp& lp);

class ColourRefinementVector {
  const HighsInt n_;
  CollectionLinkedLists colour_classes_;

 public:
  ColourRefinementVector(HighsInt n);
  HighsInt run(const std::vector<double>& w, std::vector<HighsInt>& colour);
};

}  // namespace folding

}  // namespace highs

#endif