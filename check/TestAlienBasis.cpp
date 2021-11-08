#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

TEST_CASE("AlienBasis-avgas", "[highs_test_alien_basis]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = lp.num_col_;
  const HighsInt num_row = lp.num_row_;
  HighsBasis basis;
  basis.col_status.resize(num_col);
  basis.row_status.resize(num_row);
  basis.debug_origin_name = "AlienBasis-avgas";
  // Create a full-dimension basis using struturals and then enough logicals
  HighsBasisStatus status = HighsBasisStatus::kBasic;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (iCol >= num_row) status = HighsBasisStatus::kNonbasic;
    basis.col_status[iCol] = status;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (num_col + iRow >= num_row) status = HighsBasisStatus::kNonbasic;
    basis.row_status[iRow] = status;
  }
  REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
  highs.run();
  // Create a rectangular basis using just struturals
  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    basis.col_status[iCol] = HighsBasisStatus::kBasic;
  assert(num_col < num_row);
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
  REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
  highs.run();
}
