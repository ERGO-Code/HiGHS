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

  LinkedLists colour_classes_;
  LinkedLists colour_classes_touched_;

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

  // Define iterator for range-based loop:
  //  for (HighsInt i : neighbours(j))
  //    ...
  // returns the neighbours of j in the bipartite graph, selecting correctly
  // from rows or columns, and applying the correct offset.
  struct Neighbours {
    const HighsInt* first;
    const HighsInt* last;
    const HighsInt offset;

    struct Iterator {
      const HighsInt* idx;
      const HighsInt offset;

      HighsInt operator*() const { return *idx + offset; }
      Iterator& operator++() {
        ++idx;
        return *this;
      }
      bool operator!=(const Iterator& o) const { return idx != o.idx; }
      bool operator==(const Iterator& o) const { return idx == o.idx; }
    };

    Iterator begin() const { return {first, offset}; }
    Iterator end() const { return {last, offset}; }
  };

  Neighbours neighbours(HighsInt i) const {
    if (i < A_.num_row_) {
      return {At_.index_.data() + At_.start_[i],
              At_.index_.data() + At_.start_[i + 1], A_.num_row_};
    } else {
      return {A_.index_.data() + A_.start_[i - A_.num_row_],
              A_.index_.data() + A_.start_[i + 1 - A_.num_row_], 0};
    }
  }

  void chooseRefiningColour();
  void computeColourDegrees();
  void findSplitColours();
  void splitColours();
  void splitColour(HighsInt split_colour);
  void prepareNextIter();

 public:
  ColourRefinement(const HighsSparseMatrix& A, std::vector<HighsInt>& colour);
  void run();
  HighsInt coloursUsed() const { return latest_colour_ + 1; }
};

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