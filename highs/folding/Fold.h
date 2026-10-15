#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include "ColourRefinement.h"
#include "model/HighsModel.h"

namespace highs {

namespace folding {

// LP or convex QP folding
//
// Finds a partition of the rows and columns of the model with colour
// refinement, and builds a smaller model by merging the rows and columns that
// belong to the same class.
//
// Usage:
//    Folder F(model);
//    if (F.run() == 0) {
//      HighsModel folded_model = F.getFoldedModel();
//      ...   // solve the folded model
//      HighsSolution unfolded_sol = F.unfold(folded_sol);
//    }
//
// Requirements:
//  - the model is a minimization problem with continuous variables (integrality
//    is not supported) and, if it has a Hessian, a symmetric positive
//    semidefinite one (not checked);
//  - the matrix of the model is col-wise;
//  - the Hessian can be either square or triangular;
//  - the model must outlive the Folder and must not be modified meanwhile. The
//    Folder cannot be copied, moved, or constructed from a temporary.
//
// Result:
//  - run() returns 0 on success, and nonzero if folding was abandoned because
//    the partition was too fine to be worthwhile. In that case, the folded
//    model and the unfolded solution are empty.
//  - the folded model has one row per row class and one column per column
//    class. A folded variable stands for the common value of its class. Bounds
//    come from a representative of the class, costs and Hessian entries are
//    summed over the class, and the offset is kept. All of this holds up to the
//    quantization tolerance.
//  - the matrix of the folded model is col-wise, its Hessian is lower
//    triangular.
//  - unfold() copies primal values to all the members of a class, and divides
//    the dual values by the size of the class. The unfolded solution is
//    optimal for the original model (up to the tolerance), but no basis is
//    produced.
//
// NB: the unfolded solution is not guaranteed to be a vertex of the original
// polyhedron, even if the folded solution is a vertex of the folded one. It is
// guaranteed to lie on the optimal face and it is a vertex if the original
// problem has a unique optimal solution (for LPs).

class Folder {
  const HighsModel& model_;
  const HighsLp& lp_;

  std::vector<HighsInt> row_colour_;
  std::vector<HighsInt> col_colour_;

  const HighsSparseMatrix& A_;
  const HighsSparseMatrix At_;
  const HighsHessian Q_;

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
  void foldModel();

 public:
  explicit Folder(const HighsModel& model);

  // Folder holds references to the lp, so it cannot bind to a temporary
  // HighsModel
  Folder(const HighsModel&&) = delete;

  // Avoid copies because colour refinement holds reference to the Folder
  Folder(const Folder&) = delete;
  Folder& operator=(const Folder&) = delete;

  HighsInt run();
  HighsModel getFoldedModel() const;
  HighsSolution unfold(const HighsSolution& folded) const;

  void print() const;

  bool isPartitionCorrect() const;
};

}  // namespace folding

}  // namespace highs

#endif