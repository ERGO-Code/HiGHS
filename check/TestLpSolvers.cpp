#include <fstream>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

void testDualObjective(const std::string model) {
  HighsStatus return_status;

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  highs.readModel(model_file);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  double dual_objective;
  return_status = highs.getDualObjectiveValue(dual_objective);
  REQUIRE(return_status == HighsStatus::kOk);
  double primal_objective = highs.getInfo().objective_function_value;
  double relative_primal_dual_gap =
      std::fabs(primal_objective - dual_objective) /
      std::max(1.0, std::fabs(primal_objective));
  REQUIRE(relative_primal_dual_gap < 1e-12);
}

TEST_CASE("mip-with-lp-solver", "[highs_lp_solver]") {
  // When solving the relaxation of a MIP. Exposed #1406
  HighsStatus status;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);
  highs.setOptionValue("solver", kIpmString);
  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
}

TEST_CASE("dual-objective-upper-bound", "[highs_lp_solver]") {
  std::string filename;
  HighsStatus status;
  HighsModelStatus model_status;
  bool bool_status;
  const double min_objective_function_value = -11.6389290663705;
  const double max_objective_function_value = 111.650960689315;
  const double smaller_min_objective_bound = -110.0;
  const double larger_min_objective_bound = -45.876;
  const double use_max_objective_bound = 150.0;
  double save_objective_bound;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  //  status = highs.setOptionValue("log_dev_level",
  //  kHighsLogDevLevelVerbose);

  double error;
  filename = std::string(HIGHS_DIR) + "/check/instances/e226.mps";
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  // Solve vanilla
  if (dev_run) printf("\nSolving vanilla LP\n");
  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  error = fabs((info.objective_function_value - min_objective_function_value) /
               min_objective_function_value);
  if (dev_run) printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-14);

  // Set dual objective value upper bound after saving the default value
  status = highs.getOptionValue("objective_bound", save_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.setOptionValue("objective_bound", larger_min_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  if (dev_run)
    printf(
        "\nSolving LP with presolve and dual objective value upper bound of "
        "%g\n",
        larger_min_objective_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  // Switch off presolve
  status = highs.setOptionValue("presolve", "off");
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  // This larger dual objective value upper bound is satisfied during phase 2
  if (dev_run)
    printf(
        "\nSolving LP without presolve and larger dual objective value upper "
        "bound of %g\n",
        larger_min_objective_bound);
  status = highs.clearSolver();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kObjectiveBound);

  // Solve again
  // This smaller dual objective value upper bound is satisfied at the start of
  // phase 2
  if (dev_run)
    printf(
        "\nSolving LP without presolve and smaller dual objective value upper "
        "bound of %g\n",
        smaller_min_objective_bound);
  status = highs.setOptionValue("objective_bound", smaller_min_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kObjectiveBound);

  // Solve as maximization and ensure that the dual objective value upper bound
  // isn't used
  bool_status =
      highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk;
  REQUIRE(bool_status);

  status = highs.setOptionValue("objective_bound", use_max_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  if (dev_run)
    printf(
        "\nSolving LP as maximization without presolve and dual objective "
        "value "
        "upper bound of %g\n",
        use_max_objective_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  error = fabs((info.objective_function_value - max_objective_function_value) /
               max_objective_function_value);
  if (dev_run) printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-10);
}

TEST_CASE("blending-lp-ipm", "[highs_lp_solver]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {-8, -10};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {-kHighsInf, -kHighsInf};
  lp.row_upper_ = {80, 120};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 2, 4};
  highs.passModel(lp);
  highs.setOptionValue("solver", kIpmString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.run();
  HighsInfo info = highs.getInfo();
  if (dev_run) {
    printf("Num primal infeasibilities = %d\n",
           int(info.num_primal_infeasibilities));
    printf("Max primal infeasibility   = %g\n", info.max_primal_infeasibility);
    printf("Sum primal infeasibilities = %g\n",
           info.sum_primal_infeasibilities);
    printf("Num   dual infeasibilities = %d\n",
           int(info.num_dual_infeasibilities));
    printf("Max   dual infeasibility   = %g\n", info.max_dual_infeasibility);
    printf("Sum   dual infeasibilities = %g\n", info.sum_dual_infeasibilities);
  }
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
}

TEST_CASE("dual-objective-max", "[highs_lp_solver]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 10;
  lp.col_cost_ = {8, 10};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {-kHighsInf, -kHighsInf};
  lp.row_upper_ = {80, 120};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 2, 4};
  highs.passModel(lp);
  highs.run();
  double dual_objective;
  HighsStatus return_status = highs.getDualObjectiveValue(dual_objective);
  REQUIRE(return_status == HighsStatus::kOk);
  double primal_objective = highs.getInfo().objective_function_value;
  double relative_primal_dual_gap =
      std::fabs(primal_objective - dual_objective) /
      std::max(1.0, std::fabs(primal_objective));
  REQUIRE(relative_primal_dual_gap < 1e-12);
}

TEST_CASE("dual-objective", "[highs_lp_solver]") {
  testDualObjective("avgas");
  testDualObjective("adlittle");
  testDualObjective("etamacro");
  testDualObjective("stair");
}

void testStandardForm(const HighsLp& lp) {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsInt sense = HighsInt(lp.sense_);
  highs.passModel(lp);
  highs.run();
  //  highs.writeSolution("", kSolutionStylePretty);
  double required_objective_function_value =
      highs.getInfo().objective_function_value;

  HighsInt num_col;
  HighsInt num_row;
  HighsInt num_nz;
  double offset;
  REQUIRE(highs.getStandardFormLp(num_col, num_row, num_nz, offset) ==
          HighsStatus::kOk);

  std::vector<double> cost(num_col);
  std::vector<double> rhs(num_row);
  std::vector<HighsInt> start(num_col + 1);
  std::vector<HighsInt> index(num_nz);
  std::vector<double> value(num_nz);
  REQUIRE(highs.getStandardFormLp(num_col, num_row, num_nz, offset, cost.data(),
                                  rhs.data(), start.data(), index.data(),
                                  value.data()) == HighsStatus::kOk);

  HighsLp standard_form_lp;
  standard_form_lp.num_col_ = num_col;
  standard_form_lp.num_row_ = num_row;
  standard_form_lp.offset_ = offset;
  standard_form_lp.col_cost_ = cost;
  standard_form_lp.col_lower_.assign(num_col, 0);
  standard_form_lp.col_upper_.assign(num_col, kHighsInf);
  standard_form_lp.row_lower_ = rhs;
  standard_form_lp.row_upper_ = rhs;
  standard_form_lp.a_matrix_.start_ = start;
  standard_form_lp.a_matrix_.index_ = index;
  standard_form_lp.a_matrix_.value_ = value;
  REQUIRE(highs.passModel(standard_form_lp) == HighsStatus::kOk);
  // highs.writeModel("");
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  highs.writeSolution("", kSolutionStylePretty);
  double objective_function_value =
      sense * highs.getInfo().objective_function_value;
  double objective_difference =
      std::fabs(objective_function_value - required_objective_function_value) /
      std::max(1.0, std::fabs(required_objective_function_value));
  REQUIRE(objective_difference < 1e-10);
  const bool look_at_presolved_lp = false;
  if (look_at_presolved_lp) {
    // Strange that presolve doesn't convert the constraints
    //
    // Ax-s = b; s >= 0 into Ax >= b
    REQUIRE(highs.passModel(standard_form_lp) == HighsStatus::kOk);
    REQUIRE(highs.presolve() == HighsStatus::kOk);
    HighsLp presolved_lp = highs.getPresolvedLp();
    REQUIRE(highs.passModel(presolved_lp) == HighsStatus::kOk);
    highs.writeModel("");
  }
}

void testStandardFormModel(const std::string model) {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  ;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  HighsLp lp = highs.getLp();
  testStandardForm(lp);
}

TEST_CASE("standard-form-mps", "[highs_lp_solver]") {
  testStandardFormModel("avgas");
  testStandardFormModel("afiro");
}

TEST_CASE("standard-form-lp", "[highs_lp_solver]") {
  HighsLp lp;
  lp.offset_ = -0.5;
  lp.num_col_ = 4;
  lp.num_row_ = 3;
  lp.col_cost_ = {1, 1, 1, -1};
  lp.col_lower_ = {1, -kHighsInf, -kHighsInf, -1};
  lp.col_upper_ = {kHighsInf, kHighsInf, 2, 3};
  lp.row_lower_ = {0, 1, -kHighsInf};
  lp.row_upper_ = {4, kHighsInf, 4};
  lp.a_matrix_.start_ = {0, 2, 4, 6, 8};
  lp.a_matrix_.index_ = {0, 2, 0, 1, 1, 2, 0, 2};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 1, 1, 1, 1};

  testStandardForm(lp);
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  std::vector<HighsInt> index;
  std::vector<double> value;
  // Add a fixed column and a fixed row, and maximize
  highs.passModel(lp);
  index = {0, 1, 2};
  value = {-1, 1, -1};
  REQUIRE(highs.addCol(-2.0, 1.0, 1.0, 3, index.data(), value.data()) ==
          HighsStatus::kOk);
  index = {0, 1, 2, 3};
  value = {-2, -1, 1, 3};
  REQUIRE(highs.addRow(1.0, 1.0, 4, index.data(), value.data()) ==
          HighsStatus::kOk);
  REQUIRE(highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk);
  if (dev_run)
    printf(
        "\nNow test by adding a fixed column and a fixed row, and "
        "maximizing\n");
  testStandardForm(highs.getLp());
}

TEST_CASE("simplex-stats", "[highs_lp_solver]") {
  HighsStatus return_status;

  Highs h;
  const HighsSimplexStats& simplex_stats = h.getSimplexStats();
  h.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  REQUIRE(simplex_stats.valid);
  REQUIRE(simplex_stats.iteration_count == 0);
  REQUIRE(simplex_stats.num_invert == 1);
  REQUIRE(simplex_stats.last_invert_num_el > 0);
  REQUIRE(simplex_stats.last_factored_basis_num_el > 0);
  REQUIRE(simplex_stats.col_aq_density == 0);
  REQUIRE(simplex_stats.row_ep_density == 0);
  REQUIRE(simplex_stats.row_ap_density == 0);
  REQUIRE(simplex_stats.row_DSE_density == 0);
  if (dev_run) h.reportSimplexStats(stdout);

  h.clearSolver();
  h.setOptionValue("presolve", kHighsOffString);
  REQUIRE(h.run() == HighsStatus::kOk);
  REQUIRE(simplex_stats.valid);
  REQUIRE(simplex_stats.iteration_count > 0);
  REQUIRE(simplex_stats.num_invert > 0);
  REQUIRE(simplex_stats.last_invert_num_el > 0);
  REQUIRE(simplex_stats.last_factored_basis_num_el > 0);
  REQUIRE(simplex_stats.col_aq_density > 0);
  REQUIRE(simplex_stats.row_ep_density > 0);
  REQUIRE(simplex_stats.row_ap_density > 0);
  REQUIRE(simplex_stats.row_DSE_density > 0);
  if (dev_run) h.reportSimplexStats(stdout);
}

TEST_CASE("use_warm_start", "[highs_lp_solver]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);

  h.run();
  HighsInt required_iteration_count = h.getInfo().simplex_iteration_count;
  h.setOptionValue("use_warm_start", false);
  h.run();
  HighsInt iteration_count = h.getInfo().simplex_iteration_count;
  REQUIRE(iteration_count == required_iteration_count);
}

bool fileExists(const std::string& file_name) {
  std::ifstream infile(file_name);
  return static_cast<bool>(infile.good());
}

TEST_CASE("highs-files-lp", "[highs_lp_solver]") {
  Highs h;
  std::string write_solution_file = "temp.sol";
  std::string write_basis_file = "temp.bas";
  std::string write_model_file = "temp.mps";
  h.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);

  h.setOptionValue("solution_file", write_solution_file);
  h.setOptionValue("write_basis_file", write_basis_file);
  h.setOptionValue("write_model_file", write_model_file);

  h.run();

  REQUIRE(fileExists(write_model_file));
  REQUIRE(fileExists(write_solution_file));
  REQUIRE(fileExists(write_basis_file));

  h.setOptionValue("solution_file", "");
  h.setOptionValue("write_basis_file", "");
  h.setOptionValue("write_model_file", "");

  REQUIRE(h.readModel(write_model_file) == HighsStatus::kOk);

  h.setOptionValue("read_basis_file", write_basis_file);
  h.run();
  REQUIRE(h.getInfo().simplex_iteration_count == 0);

  std::remove(write_model_file.c_str());
  std::remove(write_solution_file.c_str());
  std::remove(write_basis_file.c_str());
}

TEST_CASE("highs-files-mip", "[highs_lp_solver]") {
  Highs h;
  std::string write_solution_file = "temp.sol";
  std::string write_basis_file = "temp.bas";
  std::string write_model_file = "temp.mps";
  h.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);

  h.setOptionValue("solution_file", write_solution_file);
  h.setOptionValue("write_model_file", write_model_file);

  h.run();

  // Ideally we'd check that the files have been created, but this
  // breaks the meson build
  //
  // REQUIRE(fileExists(write_model_file));
  // REQUIRE(fileExists(write_solution_file));

  // However, std::remove returning zero is a test for existence
  REQUIRE(std::remove(write_model_file.c_str()) == 0);
  REQUIRE(std::remove(write_solution_file.c_str()) == 0);
  REQUIRE(std::remove(write_basis_file.c_str()) != 0);
}
