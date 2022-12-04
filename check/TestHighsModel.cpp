#include <cassert>

#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("HighsModel", "[highs_model]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  HighsStatus status;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  HighsModel model;
  model.lp_ = lp;
  // Add an identity Hessian
  HighsHessian& hessian = model.hessian_;
  HighsInt dim = lp.num_col_;
  hessian.dim_ = dim;
  hessian.start_.resize(dim + 1);
  hessian.index_.resize(dim);
  hessian.value_.resize(dim);
  for (int iCol = 0; iCol < dim; iCol++) {
    hessian.start_[iCol] = iCol;
    hessian.index_[iCol] = iCol;
    hessian.value_[iCol] = 1.0;
  }
  hessian.start_[dim] = dim;
  status = highs.passModel(model);
  REQUIRE(status == HighsStatus::kOk);
  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  highs.clear();
  model.clear();

  SpecialLps special_lps;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(model.lp_, require_model_status,
                             optimal_objective);
  dim = model.lp_.num_col_;
  // A Hessian with dimesion but no nonzeros (ie identically zero)
  // should be ignored OK
  hessian.dim_ = dim;
  hessian.start_.resize(dim + 1);
  hessian.start_[0] = 0;
  hessian.start_[dim] = 0;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.passModel(model);
  REQUIRE(status == HighsStatus::kOk);

  // A Hessian with small diagonal entries should cause an error. The
  // first is eliminated when normalising the Hessian - due to being
  // too small - and the second when checking the diagonal entries -
  // for being less than small_matrix_value.
  const double illegal_small_hessian_diagonal_entry = 1e-12;
  const double illegal_negative_hessian_diagonal_entry = -1;
  assert(dim == 2);
  hessian.start_[1] = 1;
  hessian.start_[2] = 2;
  hessian.index_.resize(dim);
  hessian.index_[0] = 0;
  hessian.index_[1] = 1;
  hessian.value_.resize(dim);
  hessian.value_[0] = illegal_small_hessian_diagonal_entry;
  hessian.value_[1] = illegal_negative_hessian_diagonal_entry;
  status = highs.passModel(model);
  REQUIRE(status == HighsStatus::kOk);
  status = highs.run();
  REQUIRE(status == HighsStatus::kError);
}
