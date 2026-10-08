#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include <memory>

#include "ColourRefinement.h"

namespace highs {

namespace folding {

void test_folding(const HighsLp& lp);

class Folder {
  const HighsLp& lp_;
  std::vector<HighsInt> row_colour_;
  std::vector<HighsInt> col_colour_;
  std::unique_ptr<ColourRefinement> CR_;

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;

  HighsInt initial_row_colours_;
  HighsInt initial_col_colours_;

  bool fold_success_ = false;

  double initial_time_;
  double ctor_time_;
  double run_time_;
  double fold_time_;
  mutable double folded_lp_time_;

  template <typename Data>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);

  void findInitialColour();
  HighsInt foldMatrix();

 public:
  explicit Folder(const HighsLp& lp);

  // Folder holds references to the lp, so it cannot bind to a temporary HighsLp
  Folder(const HighsLp&&) = delete;

  // Avoid copies because colour refinement holds reference to the Folder
  Folder(const Folder&) = delete;
  Folder& operator=(const Folder&) = delete;

  HighsInt run();
  HighsLp getFoldedLp() const;
  void print() const;
  void printFoldedLp() const;

  bool isPartitionCorrect() const;
};

}  // namespace folding

}  // namespace highs

#endif