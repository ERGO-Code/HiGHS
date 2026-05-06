#include <sstream>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
TEST_CASE("highs-names", "[highs_names]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string solution_file = test_name + ".sol";
  std::string name;
  const std::string model = "avgas";
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  const HighsLp& lp = highs.getLp();

  HighsInt iCol, iRow;
  HighsStatus status;

  HighsInt ck_iCol;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    status = highs.getColName(iCol, name);
    REQUIRE(status == HighsStatus::kOk);
    status = highs.getColByName(name, ck_iCol);
    REQUIRE(ck_iCol == iCol);
  }
  HighsInt ck_iRow;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    status = highs.getRowName(iRow, name);
    REQUIRE(status == HighsStatus::kOk);
    status = highs.getRowByName(name, ck_iRow);
    REQUIRE(ck_iRow == iRow);
  }

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

  status = highs.getColByName(col0_name, iCol);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(iCol == 0);
  status = highs.getRowByName(row0_name, iRow);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(iRow == 0);

  // Change name of column num_col/2 to be the same as column 0
  REQUIRE(highs.getColName(0, name) == HighsStatus::kOk);
  REQUIRE(name == col0_name);
  iCol = lp.num_col_ / 2;
  std::string iCol_name;
  REQUIRE(highs.getColName(iCol, iCol_name) == HighsStatus::kOk);
  REQUIRE(highs.passColName(iCol, col0_name) == HighsStatus::kOk);

  // column num_col/2 is no longer called iCol_name
  status = highs.getColByName(iCol_name, iCol);
  REQUIRE(status == HighsStatus::kError);

  status = highs.getColByName(col0_name, iCol);
  REQUIRE(status == HighsStatus::kError);

  // Model can (since duplicates lead to generic names in fix-2887) be
  // written
  REQUIRE(highs.writeModel("") == HighsStatus::kWarning);
  if (dev_run) highs.writeSolution("", 1);

  // Reinstate name and model writes OK
  REQUIRE(highs.passColName(iCol, iCol_name) == HighsStatus::kOk);
  REQUIRE(highs.writeModel("") == HighsStatus::kOk);

  // Change name of row num_row/2 to be the same as row 0
  REQUIRE(highs.getRowName(0, name) == HighsStatus::kOk);
  REQUIRE(name == row0_name);
  iRow = lp.num_row_ / 2;
  REQUIRE(highs.passRowName(iRow, row0_name) == HighsStatus::kOk);
  // Model can (since duplicates lead to generic names in fix-2887) be
  // written
  REQUIRE(highs.writeModel("") == HighsStatus::kWarning);
  if (dev_run) highs.writeSolution("", 1);

  // Now work with a name-less model
  HighsLp local_lp = lp;
  local_lp.col_names_.clear();
  local_lp.row_names_.clear();
  highs.passModel(local_lp);
  REQUIRE(highs.writeSolution(solution_file, 1) == HighsStatus::kWarning);

  std::remove(solution_file.c_str());
}

TEST_CASE("highs-model-name", "[model_names]") {
  Highs highs;
  const HighsLp& lp = highs.getLp();

  std::string name = lp.model_name_;
  REQUIRE(name == "");

  highs.passModelName("new_name");
  name = lp.model_name_;
  REQUIRE(name == "new_name");
}

TEST_CASE("highs-illegal-col-row-name", "[model_names]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string lp_file = test_name + ".lp";
  const std::string mps_file = test_name + ".mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {1, 2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {1, 1};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {5};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};

  HighsStatus mps_write_return = HighsStatus::kWarning;
  HighsStatus lp_write_return = HighsStatus::kWarning;
  HighsStatus status = HighsStatus::kOk;
  for (HighsInt k = 0; k < 3; k++) {
    if (k == 0) {
      // Repacing space with "_" is OK, but "^" is illegal for LP
      lp.col_names_ = {"Col 0", "Col^1"};
      lp.row_names_ = {"Row 0"};
    } else if (k == 1) {
      // Replacing blank with "c_ekk0" is OK, but "^" is illegal for
      // LP
      lp.col_names_ = {"", "Col^1"};
      lp.row_names_ = {"Row 0"};
    } else {
      // Repacing space yields duplicate, and "^" is illegal for
      // LP
      lp.col_names_ = {"Col 0", "Col_0"};
      lp.row_names_ = {"Row^0"};
    }
    status = h.passModel(lp);
    REQUIRE(status != HighsStatus::kError);

    status = h.writeModel(mps_file);
    REQUIRE(status == mps_write_return);

    status = h.writeModel(lp_file);
    REQUIRE(status == lp_write_return);
  }
  status = h.readModel(mps_file);
  REQUIRE(status == HighsStatus::kOk);

  status = h.readModel(lp_file);
  REQUIRE(status == HighsStatus::kOk);

  std::remove(mps_file.c_str());
  std::remove(lp_file.c_str());
}
