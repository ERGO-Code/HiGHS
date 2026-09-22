#ifndef HIPO_FOLDING_H
#define HIPO_FOLDING_H

#include <functional>
#include <stack>

#include "CollectionLinkedLists.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

class ColourRefinement {
  const Int n_;
  const bool bipartite_ = false;

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;

  std::vector<Int>& colour_;

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

  Int latest_colour_;
  Int refining_colour_;

  std::stack<Int> stack_refine_;
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
  void splitColour(Int split_colour);
  void prepareNextIter();

  void forEachNeighbour(Int v, const std::function<void(int)>& f);
  void forEachNeighbourNonBipartite(Int v, const std::function<void(int)>& f);
  void forEachNeighbourBipartite(Int v, const std::function<void(int)>& f);

 public:
  ColourRefinement(const HighsSparseMatrix& A, std::vector<Int>& colour,
                   bool bipartite);
  void run();
  Int coloursUsed() const { return latest_colour_; }
};

void test_folding();
void test_folding(const HighsSparseMatrix& A);

}  // namespace hipo

#endif