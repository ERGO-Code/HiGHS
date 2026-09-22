#ifndef HIPO_FOLDING_H
#define HIPO_FOLDING_H

#include "CollectionLinkedLists.h"

namespace hipo {

class ColourRefinement {
  const std::vector<Int>& ptr_;
  const std::vector<Int>& adj_;
  const Int n_;

  std::vector<Int> colour_;
  std::vector<Int> colour_degree_;
  std::vector<Int> max_colour_degree_;
  std::vector<Int> min_colour_degree_;

  std::vector<Int> colours_touched_;
  std::vector<HighsBool> in_colours_touched_;
  Int top_touched_{};

  std::vector<Int> colours_split_;
  Int top_split_{};

  CollectionLinkedLists colour_classes_;
  CollectionLinkedLists colour_classes_touched_;

  Int latest_colour_{};

  std::stack<Int> stack_refine_;
  std::vector<HighsBool> in_stack_;

  Int chooseRefiningColour();
  void computeColourDegrees(Int refining_colour);
  void findSplitColours();
  void splitColours();
  void splitColour(Int split_colour);
  void prepareNextIter();

 public:
  ColourRefinement(const std::vector<Int>& ptr, const std::vector<Int>& adj);
  void run();
  const std::vector<Int>& getColour() const;
};

void test_folding();

}  // namespace hipo

#endif