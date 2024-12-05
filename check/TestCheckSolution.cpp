// #include <cstdio>
#include <iostream>

#include "HCheckConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

void runWriteReadCheckSolution(Highs& highs, const std::string model,
                               const HighsModelStatus require_model_status,
                               const HighsInt write_solution_style);

void runSetLpSolution(const std::string model);

TEST_CASE("check-solution", "[highs_check_solution]") {
  std::string model = "";
  std::string model_file;
  HighsStatus read_status;
  HighsStatus require_read_status;
  HighsModelStatus require_model_status;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  //  const HighsInfo& info = highs.getInfo();

  HighsInt write_solution_style = kSolutionStyleRaw;
  for (HighsInt pass = 0; pass < 2; pass++) {
    const bool test_st_test23 = false;
    if (test_st_test23) {
      model = "st-test23";
      model_file = "st-test23.lp";
      require_read_status = HighsStatus::kWarning;
    } else {
      model = "avgas";  // 25fv47";
      model_file =
          std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
      require_read_status = HighsStatus::kOk;
    }

    read_status = highs.readModel(model_file);
    REQUIRE(read_status == require_read_status);

    require_model_status = HighsModelStatus::kOptimal;
    runWriteReadCheckSolution(highs, model, require_model_status,
                              write_solution_style);
    SpecialLps special_lps;
    HighsLp lp;
    double optimal_objective;

    model = "distillation";
    special_lps.distillationMip(lp, require_model_status, optimal_objective);
    highs.passModel(lp);
    runWriteReadCheckSolution(highs, model, require_model_status,
                              write_solution_style);

    lp.clear();
    model = "primalDualInfeasible1Lp";
    special_lps.primalDualInfeasible1Lp(lp, require_model_status);
    highs.passModel(lp);
    runWriteReadCheckSolution(highs, model, require_model_status,
                              write_solution_style);
    // Second pass uses sparse format
    write_solution_style = kSolutionStyleSparse;
  }
}

TEST_CASE("check-set-mip-solution", "[highs_check_solution]") {
  HighsStatus return_status;
  const std::string model = "flugpl";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  if (dev_run) printf("\n********************\nSolving from scratch\n");
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  HighsLp lp = highs.getLp();

  highs.run();
  HighsSolution optimal_solution = highs.getSolution();

  HighsInt scratch_num_nodes = info.mip_node_count;
  if (dev_run) printf("Num nodes = %d\n", int(scratch_num_nodes));

  std::string solution_file = model + ".sol";
  if (dev_run) return_status = highs.writeSolution("");
  return_status = highs.writeSolution(solution_file);
  REQUIRE(return_status == HighsStatus::kOk);

  highs.clear();

  const bool other_tests = true;
  const bool test0 = other_tests;
  bool valid, integral, feasible;
  if (test0) {
    if (dev_run)
      printf("\n***************************\nSolving from saved solution\n");
    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    return_status = highs.setSolution(optimal_solution);
    REQUIRE(return_status == HighsStatus::kOk);

    return_status = highs.assessPrimalSolution(valid, integral, feasible);
    REQUIRE(return_status == HighsStatus::kOk);

    highs.run();
    if (dev_run) printf("Num nodes = %d\n", int(info.mip_node_count));
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }

  const bool test1 = other_tests;
  if (test1) {
    if (dev_run)
      printf("\n***************************\nSolving from solution file\n");
    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    return_status = highs.readSolution(solution_file);
    REQUIRE(return_status == HighsStatus::kOk);

    return_status = highs.assessPrimalSolution(valid, integral, feasible);
    REQUIRE(return_status == HighsStatus::kOk);

    highs.run();
    if (dev_run) printf("Num nodes = %d\n", int(info.mip_node_count));
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }

  const bool test2 = other_tests;
  if (test2) {
    if (dev_run)
      printf(
          "\n***************************\nSolving from saved integer "
          "solution\n");
    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    HighsSolution solution = optimal_solution;
    double solution_dl = 0;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kInteger) continue;
      solution_dl += std::fabs(solution.col_value[iCol]);
      solution.col_value[iCol] = 0;
    }
    REQUIRE(solution_dl);

    return_status = highs.setSolution(solution);
    REQUIRE(return_status == HighsStatus::kOk);

    return_status = highs.assessPrimalSolution(valid, integral, feasible);
    REQUIRE(return_status == HighsStatus::kWarning);

    highs.run();
    if (dev_run) printf("Num nodes = %d\n", int(info.mip_node_count));
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }

  const bool test3 = other_tests;
  if (test3) {
    if (dev_run)
      printf(
          "\n***************************\nSolving from column solution file\n");
    std::string column_solution_file =
        std::string(HIGHS_DIR) + "/check/instances/flugpl_integer.sol";

    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    return_status = highs.readSolution(column_solution_file);
    REQUIRE(return_status == HighsStatus::kOk);

    return_status = highs.assessPrimalSolution(valid, integral, feasible);
    REQUIRE(return_status == HighsStatus::kWarning);

    highs.run();
    if (dev_run) printf("Num nodes = %d\n", int(info.mip_node_count));
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }

  const bool test4 = other_tests;
  if (test4) {
    if (dev_run)
      printf(
          "\n***************************\nSolving from illegal column solution "
          "file\n");
    std::string column_solution_file =
        std::string(HIGHS_DIR) + "/check/instances/flugpl_illegal_integer.sol";

    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    return_status = highs.readSolution(column_solution_file);
    REQUIRE(return_status == HighsStatus::kError);

    highs.clear();
  }

  const bool test5 = other_tests;
  if (test5) {
    HighsSolution starting_solution = optimal_solution;
    if (dev_run)
      printf(
          "\n***************************\nSolving from partial integer "
          "solution\n");
    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);

    HighsInt k = 0;
    const HighsInt max_k = 1;
    // Set a proportion of the integer variables to a fractional value
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] != HighsVarType::kInteger) continue;
      if (k <= max_k) {
        starting_solution.col_value[iCol] = 0.5;
        k++;
      } else {
        k = 0;
      }
    }
    return_status = highs.setSolution(starting_solution);
    REQUIRE(return_status == HighsStatus::kOk);
    highs.run();
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }

  const bool test6 = other_tests;
  if (test6) {
    if (dev_run)
      printf(
          "\n***************************\nSolving from sparse integer "
          "solution\n");
    HighsInt num_integer_variable = 0;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      if (lp.integrality_[iCol] == HighsVarType::kInteger)
        num_integer_variable++;

    highs.setOptionValue("output_flag", dev_run);
    highs.readModel(model_file);
    std::vector<HighsInt> index;
    std::vector<double> value;
    // Check that duplicate values are spotted
    index.push_back(0);
    value.push_back(0);
    index.push_back(1);
    value.push_back(1);
    index.push_back(0);
    value.push_back(2);
    HighsInt num_entries = index.size();
    return_status = highs.setSolution(num_entries, index.data(), value.data());
    REQUIRE(return_status == HighsStatus::kWarning);

    index.clear();
    value.clear();
    std::vector<bool> is_set;
    is_set.assign(lp.num_col_, false);
    HighsInt num_to_set = 2;
    assert(num_to_set > 0);
    HighsRandom random;
    for (HighsInt iSet = 0; iSet < num_to_set;) {
      HighsInt iCol = random.integer(lp.num_col_);
      if (lp.integrality_[iCol] != HighsVarType::kInteger) continue;
      if (is_set[iCol]) continue;
      is_set[iCol] = true;
      index.push_back(iCol);
      value.push_back(optimal_solution.col_value[iCol]);
      iSet++;
    }
    num_entries = index.size();
    assert(num_entries == num_to_set);
    return_status = highs.setSolution(num_entries, index.data(), value.data());
    REQUIRE(return_status == HighsStatus::kOk);
    highs.run();
    REQUIRE(info.mip_node_count != scratch_num_nodes);
    highs.clear();
  }
  assert(other_tests);
  std::remove(solution_file.c_str());
}

TEST_CASE("set-pathological-solution", "[highs_check_solution]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsSolution solution;

  solution.clear();
  highs.clearSolver();
  highs.addCol(1.0, 0, 1, 0, nullptr, nullptr);
  HighsInt index = 0;
  double value = 1.0;
  highs.addRow(0, 1, 1, &index, &value);
  solution.col_value.push_back(0);
  solution.row_value.push_back(0);
  highs.setSolution(solution);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  solution.clear();
  highs.clearModel();
  highs.addCol(1.0, -kHighsInf, kHighsInf, 0, nullptr, nullptr);
  solution.col_value.push_back(0);
  highs.setSolution(solution);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
}

TEST_CASE("check-set-lp-solution", "[highs_check_solution]") {
  //  runSetLpSolution("avgas");
  runSetLpSolution("adlittle");
  runSetLpSolution("shell");
  runSetLpSolution("stair");
}

TEST_CASE("check-set-rowwise-lp-solution", "[highs_check_solution]") {
  const HighsInt num_col = 100;
  std::vector<HighsInt> indices;
  std::vector<double> values;
  indices.resize(num_col);
  values.resize(num_col);
  for (HighsInt i = 0; i < num_col; i++) {
    indices[i] = i;
    values[i] = sin((double)i + 1.0);
  }
  // Round 1
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  for (HighsInt i = 0; i < num_col; i++) {
    highs.addCol(-1.0, 0.0, 1.0, 0, nullptr, nullptr);
    highs.changeColIntegrality(i, HighsVarType::kInteger);
  }
  highs.addRow(0.0, 1.0, num_col, indices.data(), values.data());
  highs.run();
  double objective1 = highs.getInfo().objective_function_value;
  HighsSolution solution = highs.getSolution();
  solution.row_value.clear();
  highs.clear();
  // Round 2
  highs.setOptionValue("output_flag", dev_run);
  for (HighsInt i = 0; i < num_col; i++) {
    highs.addCol(-1.0, 0.0, 1.0, 0, nullptr, nullptr);
    highs.changeColIntegrality(i, HighsVarType::kInteger);
  }
  highs.addRow(0.0, 1.0, num_col, indices.data(), values.data());
  highs.setSolution(solution);
  highs.run();
  double objective2 = highs.getInfo().objective_function_value;
  REQUIRE(fabs(objective1 - objective2) / max(1.0, objective1) < 1e-5);
}

TEST_CASE("check-set-mip-solution-extra-row", "[highs_check_solution]") {
  Highs highs;
  const std::string solution_file_name = "temp.sol";
  highs.setOptionValue("output_flag", dev_run);
  highs.addVar(0, 2);
  highs.addVar(0, 2);
  highs.changeColCost(0, 1);
  highs.changeColCost(1, 10);
  highs.changeColIntegrality(0, HighsVarType::kInteger);
  std::vector<HighsInt> index = {0, 1};
  std::vector<double> value = {1, 1};
  highs.addRow(1, kHighsInf, 2, index.data(), value.data());
  highs.run();
  highs.writeSolution(solution_file_name);
  if (dev_run) highs.writeSolution("", 1);
  highs.clearSolver();
  // Add a constraint that cuts off the optimal solution, but leaves
  // the integer assignment feasible
  value[0] = 1;
  value[1] = 4;
  highs.addRow(4, kHighsInf, 2, index.data(), value.data());
  // Read the original solution - testing that the row section is not
  // used
  REQUIRE(highs.readSolution(solution_file_name) == HighsStatus::kOk);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  std::remove(solution_file_name.c_str());
}

TEST_CASE("check-set-illegal-solution", "[highs_check_solution]") {
  HighsStatus return_status;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  const HighsLp& lp = highs.getLp();
  HighsSolution solution;
  REQUIRE(highs.setSolution(solution) == HighsStatus::kError);
  solution.col_value.assign(lp.num_col_, 0);
  REQUIRE(highs.setSolution(solution) == HighsStatus::kOk);
}

TEST_CASE("read-miplib-solution", "[highs_check_solution]") {
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMaximize;
  lp.col_cost_ = {8, 5, 3, 11, 7};
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {11};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 5};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 4};
  lp.a_matrix_.value_ = {4, 3, 1, 5, 4};
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("presolve", kHighsOffString);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  REQUIRE(h.run() == HighsStatus::kOk);
  //  REQUIRE(h.writeSolution("", kSolutionStylePretty) == HighsStatus::kOk);
  const std::vector<double>& col_value = h.getSolution().col_value;
  std::string miplib_sol_file = "miplib.sol";
  FILE* file = fopen(miplib_sol_file.c_str(), "w");
  REQUIRE(file != 0);
  fprintf(file, "=obj= 22\n");
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    std::string col_name = "c" + std::to_string(int(iCol));
    lp.col_names_.push_back(col_name);
    if (std::fabs(col_value[iCol]) < 1e-2) continue;
    std::string line = col_name + " 1\n";
    fprintf(file, "%s", line.c_str());
  }
  fclose(file);
  // Can't read file yet, as model has no column names
  REQUIRE(h.readSolution(miplib_sol_file) == HighsStatus::kError);

  // Pass model again now that column names have been defined

  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  //  REQUIRE(h.writeModel("miplib.mps") == HighsStatus::kOk);
  REQUIRE(h.readSolution(miplib_sol_file) == HighsStatus::kOk);
  REQUIRE(h.run() == HighsStatus::kOk);
  std::remove(miplib_sol_file.c_str());
}

void runWriteReadCheckSolution(Highs& highs, const std::string model,
                               const HighsModelStatus require_model_status,
                               const HighsInt write_solution_style) {
  HighsStatus run_status;
  HighsStatus return_status;
  std::string solution_file;
  HighsModelStatus status = HighsModelStatus::kNotset;
  if (dev_run) printf("\nSolving model %s from scratch\n", model.c_str());
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  status = highs.getModelStatus();
  REQUIRE(status == require_model_status);

  solution_file = model + ".sol";
  if (dev_run)
    printf("Writing solution in style %d to %s\n", int(write_solution_style),
           solution_file.c_str());
  if (dev_run) return_status = highs.writeSolution("", write_solution_style);
  return_status = highs.writeSolution(solution_file, write_solution_style);
  REQUIRE(return_status == HighsStatus::kOk);

  const bool& value_valid = highs.getSolution().value_valid;
  bool valid, integral, feasible;

  // primalDualInfeasible1Lp has no values in the solution file so,
  // after it's read, HiGHS::solution.value_valid is false
  if (dev_run) printf("Reading solution from %s\n", solution_file.c_str());
  return_status = highs.readSolution(solution_file);
  if (value_valid) {
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    REQUIRE(return_status == HighsStatus::kWarning);
  }

  return_status = highs.assessPrimalSolution(valid, integral, feasible);
  if (value_valid) {
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    REQUIRE(return_status == HighsStatus::kError);
  }
  if (dev_run) printf("Solving model from solution read from file\n");
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  status = highs.getModelStatus();
  REQUIRE(status == require_model_status);

  std::remove(solution_file.c_str());
}

void runSetLpSolution(const std::string model) {
  HighsStatus return_status;
  Highs highs;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  const HighsInfo& info = highs.getInfo();
  if (dev_run) printf("\nSolving %s from scratch\n", model.c_str());
  highs.setOptionValue("output_flag", dev_run);
  if (dev_run) highs.setOptionValue("log_dev_level", 1);

  highs.readModel(model_file);
  highs.run();
  HighsInt simplex_iteration_count0 = info.simplex_iteration_count;
  HighsSolution solution = highs.getSolution();
  highs.clear();
  if (dev_run) printf("\nSolving from saved solution\n");
  highs.setOptionValue("output_flag", dev_run);
  if (dev_run) highs.setOptionValue("log_dev_level", 1);
  highs.readModel(model_file);

  // solution.col_value.assign(highs.getNumCol(), 0);

  return_status = highs.setSolution(solution);
  REQUIRE(return_status == HighsStatus::kOk);

  highs.run();
  // Use a reduction in iteration count as a sanity check that
  // starting from the optimal solution has worked
  HighsInt simplex_iteration_count1 = info.simplex_iteration_count;
  if (dev_run)
    printf(
        "For model %s: iteration counts are %d from logical basis and %d from "
        "optimal solution\n",
        model.c_str(), (int)simplex_iteration_count0,
        (int)simplex_iteration_count1);
  REQUIRE(simplex_iteration_count1 < simplex_iteration_count0);

  // Now write a sparse solution, and read it in to hot start
  HighsInt write_solution_style = kSolutionStyleSparse;
  std::string solution_file = model + ".sol";
  if (dev_run) printf("Writing sparse solution to %s\n", solution_file.c_str());
  if (dev_run) return_status = highs.writeSolution("");
  return_status = highs.writeSolution(solution_file, write_solution_style);
  REQUIRE(return_status == HighsStatus::kOk);

  highs.clear();
  highs.setOptionValue("output_flag", dev_run);
  if (dev_run) highs.setOptionValue("log_dev_level", 1);

  highs.readModel(model_file);
  if (dev_run)
    printf("Reading sparse solution from %s\n", solution_file.c_str());
  return_status = highs.readSolution(solution_file);
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) printf("Solving model from sparse solution read from file\n");
  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  HighsModelStatus status = highs.getModelStatus();
  REQUIRE(status == HighsModelStatus::kOptimal);

  highs.clear();

  std::remove(solution_file.c_str());
}
