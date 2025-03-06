#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderEms.h"
#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"

const bool dev_run = false;
const double inf = kHighsInf;

TEST_CASE("filereader-edge-cases", "[highs_filereader]") {
  std::string model = "";
  std::string model_file;
  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  // Switching off for debugging
  const bool run_first_tests = true;

  const bool test_garbage_mps = true;
  const bool test_garbage_ems = true;
  const bool test_garbage_lp = true;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  if (run_first_tests) {
    // Try to run HiGHS with default options. No model loaded so OK
    run_status = highs.run();
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kModelEmpty);
    REQUIRE(run_status == HighsStatus::kOk);

    // Load a non-existent MPS file and try to run HiGHS
    model = "";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    return_status = highs.readModel(model_file);
    REQUIRE(return_status == HighsStatus::kError);
    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);

    // Load a non-existent LP file and try to run HiGHS
    model = "";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
    return_status = highs.readModel(model_file);
    REQUIRE(return_status == HighsStatus::kError);
    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);

    // Load a non-supported file type and try to run HiGHS
    model = "model";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".xyz";
    return_status = highs.readModel(model_file);
    REQUIRE(return_status == HighsStatus::kError);
    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);

    // Load an LP file that does not end with a newline
    model = "no-newline-eof";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
    return_status = highs.readModel(model_file);
    REQUIRE(return_status == HighsStatus::kOk);

    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);

    // Load an existing MPS file and run HiGHS
    model = "adlittle";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    return_status = highs.readModel(model_file);
    REQUIRE(return_status == HighsStatus::kOk);

    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(info.simplex_iteration_count == 87);

    model = "garbage";
    if (test_garbage_mps) {
      if (dev_run) printf("\ngarbage.mps\n");
      model_file =
          std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
      read_status = highs.readModel(model_file);
      REQUIRE(read_status == HighsStatus::kError);
    }

    if (test_garbage_ems) {
      if (dev_run) printf("\ngarbage.ems\n");
      model_file =
          std::string(HIGHS_DIR) + "/check/instances/" + model + ".ems";
      read_status = highs.readModel(model_file);
      REQUIRE(read_status == HighsStatus::kError);
    }

    if (test_garbage_lp) {
      if (dev_run) printf("\ngarbage.lp\n");
      model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
      read_status = highs.readModel(model_file);
      REQUIRE(read_status == HighsStatus::kError);
    }
  }

  // blah blah not a good file
  model = "1448";
  if (dev_run) printf("\n%s.mps\n", model.c_str());
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kError);

  // Gurobi cannot read
  //
  // Minimize a subject to a >= 1 bounds a <= 0
  model = "1449a";
  if (dev_run) printf("\n%s.mps\n", model.c_str());
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kError);

  // Gurobi can read
  //
  // Minimize a subject to a >= 1
  //
  // However, requiring "end" checks for file corruption
  model = "1449b";
  if (dev_run) printf("\n%s.mps\n", model.c_str());
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kError);

  model = "1451";
  if (dev_run) printf("\n%s.lp\n", model.c_str());
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kOk);
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getInfo().objective_function_value == 2);
}

void freeFixedModelTest(const std::string model_name) {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/" + model_name + ".mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);

  HighsModel model_free = highs.getModel();

  REQUIRE(highs.setOptionValue("mps_parser_type_free", false) ==
          HighsStatus::kOk);

  REQUIRE(highs.readModel(filename) == HighsStatus::kWarning);

  HighsModel model_fixed = highs.getModel();

  bool are_the_same = model_free == model_fixed;
  REQUIRE(are_the_same);
}

TEST_CASE("filereader-free-format-parser-qp", "[highs_filereader]") {
  freeFixedModelTest("qjh");
  freeFixedModelTest("qjh_quadobj");
  // This test can't be used since fixed format reader can't handle
  // QMATRIX section
  //
  //  freeFixedModelTest("qjh_qmatrix");
}

TEST_CASE("filereader-free-format-parser-lp", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  HighsStatus status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_free = highs.getLp();

  status = highs.setOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_fixed = highs.getLp();

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("filereader-read-mps-ems-lp", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);
  HighsLp lp_mps = highs.getLp();

  // Write lp
  std::string filename_lp = "adlittle.lp";
  status = highs.writeModel(filename_lp);
  REQUIRE(status == HighsStatus::kOk);

  /*
  bool are_the_same;
  // Write ems
  std::string filename_ems = "adlittle.ems";
  status = highs.writeModel(filename_ems);
  REQUIRE(status == HighsStatus::kOk);

  // Read ems and compare with mps
  std::cout << "Reading " << filename_ems << std::endl;
  status = highs.readModel(filename_ems);
  REQUIRE(status == HighsStatus::kOk);

  std::cout << "Compare LP from .ems and .mps" << std::endl;
  are_the_same = lp_mps == highs.getLp();
  REQUIRE(are_the_same);

  std::remove(filename_ems.c_str());
  */

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  const HighsInfo& info = highs.getInfo();
  double mps_objective_function_value = info.objective_function_value;

  // Read lp and compare objective with mps
  if (dev_run) std::cout << "Reading " << filename_lp << std::endl;
  status = highs.readModel(filename_lp);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  double delta =
      std::fabs(mps_objective_function_value - info.objective_function_value);
  if (dev_run && delta)
    printf("TestFilereader: LP-MPS |objective difference| of %g\n", delta);
  REQUIRE(delta < 1e-8);

  std::remove(filename_lp.c_str());
}

TEST_CASE("filereader-integrality-constraints", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";

  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<HighsVarType> kIntegers{
      HighsVarType::kContinuous, HighsVarType::kContinuous,
      HighsVarType::kInteger,    HighsVarType::kInteger,
      HighsVarType::kContinuous, HighsVarType::kContinuous,
      HighsVarType::kContinuous, HighsVarType::kContinuous};

  HighsStatus status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_free = highs.getLp();

  REQUIRE((int)lp_free.integrality_.size() == lp_free.num_col_);
  REQUIRE(lp_free.integrality_ == kIntegers);

  // Read mps with fixed format parser.
  status = highs.setOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_fixed = highs.getLp();

  REQUIRE((int)lp_fixed.integrality_.size() == lp_fixed.num_col_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

/*
TEST_CASE("filereader-nan", "[highs_filereader]") {
  // Check that if
  std::string model_file;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  model_file = std::string(HIGHS_DIR) + "/check/instances/nan0.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kError);
  model_file = std::string(HIGHS_DIR) + "/check/instances/nan1.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kError);
  model_file = std::string(HIGHS_DIR) + "/check/instances/nan2.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kError);
}
*/

TEST_CASE("filereader-fixed-integer", "[highs_filereader]") {
  double objective_value;
  const double optimal_objective_value = 0;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/fixed-binary.lp";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  objective_value = highs.getInfo().objective_function_value;
  REQUIRE(objective_value == optimal_objective_value);
}

TEST_CASE("filereader-dD2e", "[highs_filereader]") {
  // dD2e.mps is min -x1 - 2x2 with upper bounds 1.0D3 and 1.0d3
  //
  // If read correctly, the optimal objective value is -3000
  double objective_value;
  const double optimal_objective_value = -3000;
  std::string model_file = std::string(HIGHS_DIR) + "/check/instances/dD2e.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  objective_value = highs.getInfo().objective_function_value;
  REQUIRE(objective_value == optimal_objective_value);
}

// TEST_CASE("filereader-comment", "[highs_filereader]") {
//   // Check that comments - either whole line with * in first column,
//   // or rest of line following */$ are handled correctly
//   const double optimal_objective_value = -4;
//   std::string model_file =
//       std::string(HIGHS_DIR) + "/check/instances/comment.mps";
//   Highs highs;
//   highs.setOptionValue("output_flag", dev_run);
//   REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
//   REQUIRE(highs.run() == HighsStatus::kOk);
//   double objective_value = highs.getInfo().objective_function_value;
//   REQUIRE(objective_value == optimal_objective_value);
// }

TEST_CASE("writeLocalModel", "[highs_filereader]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  std::string write_model_file = "foo.mps";
  HighsModel model;
  HighsLp& lp = model.lp_;
  ;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {8, 10};
  lp.row_lower_ = {7, 12, 6};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {2, 3, 2, 2, 4, 1};

  if (dev_run) printf("\nModel with no column lower or upper bounds\n");
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);
  lp.col_lower_ = {0, 0};

  if (dev_run) printf("\nModel with no column upper bounds\n");
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);
  lp.col_upper_ = {inf, inf};

  // Model has no dimensions for a_matrix_, but these are set in
  // writeLocalModel.
  if (dev_run) printf("\nModel with no column or row names\n");
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kWarning);
  lp.col_names_ = {"C0", "C1"};
  lp.row_names_ = {"R0", "R1", "R2"};

  if (dev_run) printf("\nModel with column and row names\n");
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kOk);

  // Introduce illegal start
  if (dev_run) printf("\nModel with start entry > num_nz\n");
  lp.a_matrix_.start_ = {0, 7, 6};
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);

  // Introduce illegal start
  if (dev_run) printf("\nModel with start entry -1\n");
  lp.a_matrix_.start_ = {0, -1, 6};
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);
  lp.a_matrix_.start_ = {0, 3, 6};

  // Introduce illegal index
  if (dev_run) printf("\nModel with index entry -1\n");
  lp.a_matrix_.index_ = {0, -1, 2, 0, 1, 2};
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);

  // Introduce illegal index
  if (dev_run) printf("\nModel with index entry 3 >= num_row\n");
  lp.a_matrix_.index_ = {0, 1, 3, 0, 1, 2};
  REQUIRE(h.writeLocalModel(model, write_model_file) == HighsStatus::kError);

  std::remove(write_model_file.c_str());
}

TEST_CASE("write-MI-bound-model", "[highs_filereader]") {
  std::string write_model_file = "temp.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.addCol(1, -kHighsInf, 1, 0, nullptr, nullptr);
  h.changeColIntegrality(0, HighsVarType::kInteger);
  h.passColName(0, "x");
  std::vector<HighsInt> index = {0};
  std::vector<double> value = {1};
  h.addRow(-10, kHighsInf, 1, index.data(), value.data());
  h.passRowName(0, "r");
  h.run();
  double required_objective_value = h.getInfo().objective_function_value;
  // writeModel must ensure that there is a line
  //
  // MI BOUND x
  h.writeModel(write_model_file);
  h.readModel(write_model_file);
  h.run();
  REQUIRE(required_objective_value == h.getInfo().objective_function_value);
  std::remove(write_model_file.c_str());
}

TEST_CASE("mps-warnings", "[highs_filereader]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/warnings.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  HighsStatus return_status = h.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kWarning);
}

TEST_CASE("mps-silly-names", "[highs_filereader]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/silly-names.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  HighsStatus return_status = h.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);
}
