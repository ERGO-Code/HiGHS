#include <sstream>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
TEST_CASE("highs-names", "[highs_names]") {
  const std::string model = "avgas";
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  const HighsLp& lp = highs.getLp();

  // Change all names to distinct new names
  REQUIRE(highs.passColName(-1, "FRED") == HighsStatus::kError);
  REQUIRE(highs.passColName(lp.num_col_, "FRED") == HighsStatus::kError);
  REQUIRE(highs.passColName(0, "") == HighsStatus::kError);
  std::string col0_name;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    std::stringstream ss;
    ss.str(std::string());
    ss << model << "_col_" << iCol << "\0";
    const std::string name = ss.str();
    if (iCol == 0) col0_name = name;
    if (dev_run) printf("Col %d name is to be %s\n", int(iCol), name.c_str());
    REQUIRE(highs.passColName(iCol, name) == HighsStatus::kOk);
  }
  REQUIRE(highs.passRowName(-1, "FRED") == HighsStatus::kError);
  REQUIRE(highs.passRowName(lp.num_row_, "FRED") == HighsStatus::kError);
  REQUIRE(highs.passRowName(0, "") == HighsStatus::kError);
  std::string row0_name;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    std::stringstream ss;
    ss.str(std::string());
    ss << model << "_row_" << iRow << "\0";
    const std::string name = ss.str();
    if (iRow == 0) row0_name = name;
    if (dev_run) printf("Row %d name is to be %s\n", int(iRow), name.c_str());
    REQUIRE(highs.passRowName(iRow, name) == HighsStatus::kOk);
  }
  highs.run();
  REQUIRE(highs.writeModel("") == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", 1);

  // Change name of column num_col/2 to be the same as column 0
  std::string name;
  REQUIRE(highs.getColName(0, name) == HighsStatus::kOk);
  REQUIRE(name == col0_name);
  HighsInt iCol = lp.num_col_ / 2;
  std::string iCol_name;
  REQUIRE(highs.getColName(iCol, iCol_name) == HighsStatus::kOk);
  REQUIRE(highs.passColName(iCol, col0_name) == HighsStatus::kOk);
  // Model can't be written
  REQUIRE(highs.writeModel("") == HighsStatus::kError);
  if (dev_run) highs.writeSolution("", 1);

  // Reinstate name and model writes OK
  REQUIRE(highs.passColName(iCol, iCol_name) == HighsStatus::kOk);
  REQUIRE(highs.writeModel("") == HighsStatus::kOk);

  // Change name of row num_row/2 to be the same as row 0
  REQUIRE(highs.getRowName(0, name) == HighsStatus::kOk);
  REQUIRE(name == row0_name);
  HighsInt iRow = lp.num_row_ / 2;
  REQUIRE(highs.passRowName(iRow, row0_name) == HighsStatus::kOk);
  // Model can't be written
  REQUIRE(highs.writeModel("") == HighsStatus::kError);
  if (dev_run) highs.writeSolution("", 1);

  // Now work with a name-less model
  HighsLp local_lp = lp;
  local_lp.col_names_.clear();
  local_lp.row_names_.clear();
  highs.passModel(local_lp);
  const std::string solution_file = "temp.sol";
  REQUIRE(highs.writeSolution(solution_file, 1) == HighsStatus::kOk);

  // Cannot get name of column or row 0
  REQUIRE(highs.getColName(0, name) == HighsStatus::kError);
  REQUIRE(highs.getRowName(0, name) == HighsStatus::kError);

  std::remove(solution_file.c_str());
}
