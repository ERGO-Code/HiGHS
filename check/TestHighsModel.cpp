#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

// No commas in test case name.
TEST_CASE("HighsModel", "[highs_model]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  HighsStatus status;
  Highs highs;
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  HighsModel model;
  model.lp_ = lp;
  // Add an identity Hessian
  HighsHessian& hessian = model.hessian_;
  const HighsInt dim = lp.numCol_;
  hessian.dim_ = dim;
  hessian.q_start_.resize(dim + 1);
  hessian.q_index_.resize(dim);
  hessian.q_value_.resize(dim);
  for (int iCol = 0; iCol < dim; iCol++) {
    hessian.q_start_[iCol] = iCol;
    hessian.q_index_[iCol] = iCol;
    hessian.q_value_[iCol] = 1.0;
  }
  hessian.q_start_[dim] = dim;
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
}
