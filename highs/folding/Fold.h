#ifndef HIGHS_FOLDING_H
#define HIGHS_FOLDING_H

#include <memory>

#include "ColourRefinement.h"
#include "model/HighsModel.h"
#include "util/HighsTimer.h"

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
//    F.run();
//    if (F.getInfo().status == Status::kOk) {
//      HighsModel folded_model = F.getFoldedModel();
//      ...   // solve the folded model
//      HighsSolution unfolded_sol = F.unfold(folded_sol);
//    }
//
// Requirements:
//  - the model has only continuous variables (integrality is not supported)
//    and, if it has a Hessian, a symmetric positive semidefinite one (not
//    checked);
//  - the matrix of the model is col-wise;
//  - the Hessian can be either square or triangular;
//  - the model must outlive the Folder and must not be modified meanwhile. The
//    Folder cannot be copied, moved, or constructed from a temporary.
//
// Result:
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

struct FoldInfo {
  HighsInt initial_row_colours{};
  HighsInt initial_col_colours{};
  HighsInt num_row_colours{};
  HighsInt num_col_colours{};

  Status status = Status::kNotRun;

  double initial_time{};
  double ctor_time{};
  double run_time{};
  double fold_time{};
  mutable double folded_lp_time{};
};

class Folder {
  const HighsLp& lp_;
  const HighsSparseMatrix& A_;
  const HighsHessian& hessian_;

  struct Data {
    const HighsSparseMatrix At;
    const HighsHessian Q;
    std::vector<HighsInt> row_colour;
    std::vector<HighsInt> col_colour;

    Data(const HighsLp& lp, const HighsHessian& hessian)
        : At(lp.a_matrix_.buildOppositeFormat()),
          Q(hessian.toSquare()),
          row_colour(lp.num_row_),
          col_colour(lp.num_col_) {}
  };
  std::unique_ptr<Data> data_;

  std::vector<HighsInt> row_colours_size_;
  std::vector<HighsInt> col_colours_size_;

  const HighsTimer* timer_ = nullptr;
  double time_limit_ = kHighsInf;

  FoldInfo info_;

  template <typename T>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);
  void findInitialColour();
  void foldModel();
  bool isPartitionCorrect() const;
  bool success() const { return info_.status == Status::kOk; }
  bool checkTimeLimit();

 public:
  explicit Folder(const HighsModel& model)
      : lp_{model.lp_}, A_{model.lp_.a_matrix_}, hessian_{model.hessian_} {}

  // Folder holds references to the lp, so it cannot bind to a temporary
  // HighsModel
  Folder(const HighsModel&&) = delete;

  // Avoid copies because colour refinement holds reference to the Folder
  Folder(const Folder&) = delete;
  Folder& operator=(const Folder&) = delete;

  void setTimeLimit(const HighsTimer& timer, double time_limit);
  void run();
  HighsModel getFoldedModel() const;
  HighsSolution unfold(const HighsSolution& folded) const;
  const FoldInfo& getInfo() const { return info_; }

  void print() const;
};

}  // namespace folding

}  // namespace highs

#endif