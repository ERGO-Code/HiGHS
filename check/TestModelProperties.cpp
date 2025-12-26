#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;

TEST_CASE("simplest-ill-conditioning", "[highs_model_properties]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  const double epsilon = 1e-4;
  const double ill_conditioning_bound = 1;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {2, 2 + epsilon};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {2, 2 + epsilon};
  lp.row_upper_ = {inf, inf};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, 1 + epsilon};
  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsIllConditioning ill_conditioning;

  const bool constraint = true;
  highs.getIllConditioning(ill_conditioning, constraint);
  REQUIRE(ill_conditioning.record.size() == 2);
  // Both multipliers should be large
  for (HighsInt iX = 0; iX < 2; iX++) {
    REQUIRE(std::fabs(ill_conditioning.record[iX].multiplier) > 0.45);
    REQUIRE(std::fabs(ill_conditioning.record[iX].multiplier) < 0.55);
  }
  highs.getIllConditioning(ill_conditioning, !constraint);

  REQUIRE(highs.getIllConditioning(ill_conditioning, constraint, 1, 0.1) ==
          HighsStatus::kOk);
  REQUIRE(highs.getIllConditioning(ill_conditioning, constraint, 1,
                                   ill_conditioning_bound) == HighsStatus::kOk);
  REQUIRE(highs.getIllConditioning(ill_conditioning, constraint, 1, 10) ==
          HighsStatus::kOk);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("simple-ill-conditioning", "[highs_model_properties]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  const double epsilon = 1e-4;
  lp.num_col_ = 3;
  lp.num_row_ = 3;
  lp.col_cost_ = {3, 2, 3 + epsilon};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {inf, inf, inf};
  lp.row_lower_ = {3, 2, 3 + epsilon};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 5, 8};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 1, 1, 1, 1 + epsilon};

  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsIllConditioning ill_conditioning;
  const bool constraint = true;
  highs.getIllConditioning(ill_conditioning, constraint);
  REQUIRE(ill_conditioning.record.size() == 3);
  // First two multipliers should be the large ones
  for (HighsInt iX = 0; iX < 2; iX++) {
    REQUIRE(std::fabs(ill_conditioning.record[iX].multiplier) > 0.45);
    REQUIRE(std::fabs(ill_conditioning.record[iX].multiplier) < 0.55);
  }
  highs.getIllConditioning(ill_conditioning, !constraint);
}

TEST_CASE("afiro-ill-conditioning", "[highs_model_properties]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/afiro.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  HighsInt num_nz;
  std::vector<HighsInt> index(lp.num_col_);
  std::vector<double> value(lp.num_col_);
  highs.run();
  const HighsBasis& highs_basis = highs.getBasis();
  lp.a_matrix_.getRow(0, num_nz, index.data(), value.data());
  if (dev_run) {
    for (HighsInt iEl = 0; iEl < num_nz; iEl++) {
      HighsInt iCol = index[iEl];
      printf("%s: %d %19.12g (%s)\n", lp.col_names_[iCol].c_str(), int(iCol),
             value[iEl],
             highs.basisStatusToString(highs_basis.col_status[iCol]).c_str());
    }
  }
  value[0] += 1e-4;

  const bool negate_bad_row = false;
  if (negate_bad_row)
    for (HighsInt iEl = 0; iEl < num_nz; iEl++) value[iEl] *= -1;

  highs.addRow(0, 0, num_nz, index.data(), value.data());
  HighsInt bad_row = lp.num_row_ - 1;
  highs.passRowName(bad_row, "R09bad");
  // Find a nonbasic row to replace bad row in the basis
  HighsInt nonbasic_row = -1;
  for (HighsInt iRow = 1; iRow < lp.num_row_; iRow++) {
    if (dev_run) {
      printf("Row %d (%s) has status %s\n", int(iRow),
             lp.row_names_[iRow].c_str(),
             highs.basisStatusToString(highs_basis.row_status[iRow]).c_str());
    }
    if (highs_basis.row_status[iRow] != HighsBasisStatus::kBasic) {
      nonbasic_row = iRow;
      break;
    }
  }
  // Bad row should be basic - since it's been added
  REQUIRE(highs_basis.row_status[bad_row] == HighsBasisStatus::kBasic);
  REQUIRE(nonbasic_row >= 0);
  HighsBasis basis = highs_basis;
  // Make the bad row nonbasic at lower bound - it's FX - and make the
  // nonbasic_row basic
  basis.row_status[bad_row] = HighsBasisStatus::kLower;
  basis.row_status[nonbasic_row] = HighsBasisStatus::kBasic;
  highs.setBasis(basis);
  if (dev_run) {
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      printf("Row %d (%s) has status %s\n", int(iRow),
             lp.row_names_[iRow].c_str(),
             highs.basisStatusToString(highs_basis.row_status[iRow]).c_str());
    }
  }
  HighsIllConditioning ill_conditioning;
  const bool constraint = true;
  highs.getIllConditioning(ill_conditioning, constraint);
  highs.getIllConditioning(ill_conditioning, !constraint);
}
