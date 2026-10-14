#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include "ColourRefinement.h"
#include "model/HighsModel.h"

namespace highs {

namespace folding {

class LpFolder {
  const HighsModel& model_;
  const HighsLp& lp_;
  const HighsHessian& Q_;

  std::vector<HighsInt> row_colour_;
  std::vector<HighsInt> col_colour_;

  const HighsSparseMatrix& A_;
  HighsSparseMatrix At_;

  HighsInt initial_row_colours_{};
  HighsInt initial_col_colours_{};
  HighsInt num_row_colours_{};
  HighsInt num_col_colours_{};

  std::vector<HighsInt> row_colours_size_;
  std::vector<HighsInt> col_colours_size_;

  HighsInt fold_error_ = -1;

  double initial_time_{};
  double ctor_time_{};
  double run_time_{};
  double fold_time_{};
  mutable double folded_lp_time_{};

  template <typename Data>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);

  void findInitialColour();
  void foldMatrix();

 public:
  explicit LpFolder(const HighsModel& model);

  // LpFolder holds references to the lp, so it cannot bind to a temporary
  // HighsLp
  LpFolder(const HighsModel&&) = delete;

  // Avoid copies because colour refinement holds reference to the LpFolder
  LpFolder(const LpFolder&) = delete;
  LpFolder& operator=(const LpFolder&) = delete;

  HighsInt run();
  HighsModel getFoldedModel() const;
  HighsSolution unfold(const HighsSolution& folded) const;

  void print() const;

  bool isPartitionCorrect() const;
};

}  // namespace folding

}  // namespace highs

#endif