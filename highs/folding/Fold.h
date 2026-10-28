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
//    Folder F;
//    F.run(model);
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
//  - the model is not modified by other threads while run() is executing.
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
  HighsInt initial_row_colours;
  HighsInt initial_col_colours;
  HighsInt num_row_colours;
  HighsInt num_col_colours;
  Status status;

  FoldInfo() { clear(); }
  void clear() {
    initial_row_colours = 0;
    initial_col_colours = 0;
    num_row_colours = 0;
    num_col_colours = 0;
    status = Status::kNotRun;
  }
};

class Folder {
  struct Data {
    // This data exists only during run()
    const HighsLp& lp;
    const HighsSparseMatrix At;
    const HighsHessian Q;

    Data(const HighsModel& model)
        : lp{model.lp_},
          At(lp.a_matrix_.getRowwise()),
          Q(model.hessian_.toSquare()) {}
  };
  std::unique_ptr<Data> data_;

  HighsInt num_row_{};
  HighsInt num_col_{};

  std::vector<HighsInt> row_colour_;
  std::vector<HighsInt> col_colour_;

  std::vector<HighsInt> row_colours_size_;
  std::vector<HighsInt> col_colours_size_;

  const HighsTimer* timer_ = nullptr;
  double time_limit_ = kHighsInf;

  FoldInfo info_;

  HighsModel folded_model_;

  template <typename T>
  HighsInt findInitialColour(HighsInt num, HighsInt* colour, HighsInt start);
  void findInitialColour();
  void foldModel();
  bool isPartitionCorrect() const;
  bool success() const { return info_.status == Status::kOk; }
  bool checkTimeLimit();
  void buildFoldedModel();
  void discardResult();

 public:
  Folder() = default;

  Folder(const Folder&) = delete;
  Folder& operator=(const Folder&) = delete;
  Folder(Folder&&) = default;
  Folder& operator=(Folder&&) = default;

  void setTimeLimit(const HighsTimer& timer, double time_limit);
  void run(const HighsModel& model);
  const HighsModel& getFoldedModel() const { return folded_model_; }
  HighsModel takeFoldedModel() { return std::move(folded_model_); }
  HighsSolution unfold(const HighsSolution& folded) const;
  const FoldInfo& getInfo() const { return info_; }
};

}  // namespace folding

}  // namespace highs

#endif