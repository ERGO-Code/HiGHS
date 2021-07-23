#include <cstdio>

#include "FilereaderLp.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

TEST_CASE("qpsolver", "[qpsolver]") {
  double required_objective_function_value;
  double required_x0;
  double required_x1;
  double required_x2;
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qptestnw.lp";

  required_objective_function_value = -6.45;
  required_x0 = 1.4;
  required_x1 = 1.7;

  Highs highs;
  const HighsModel& model = highs.getModel();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  const double& objective_function_value = info.objective_function_value;

  if (!dev_run) highs.setOptionValue("output_flag", false);

  HighsStatus return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  double alt_objective_function_value =
      model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);

  // Check with qjh.mps
  filename = std::string(HIGHS_DIR) + "/check/instances/qjh.mps";
  required_objective_function_value = -5.25;
  required_x0 = 0.5;
  required_x1 = 5.0;
  required_x2 = 1.5;

  return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) printf("Objective = %g\n", objective_function_value);

  alt_objective_function_value = model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - required_x2) < double_equal_tolerance);
  REQUIRE(return_status == HighsStatus::kOk);

  // Test writeModel by writing out qjh.mps...
  filename = "qjh.mps";
  highs.writeModel(filename);

  // ... and reading it in again
  return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) printf("Objective = %g\n", objective_function_value);

  alt_objective_function_value = model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - required_x2) < double_equal_tolerance);
}

TEST_CASE("test-qod", "[qpsolver]") {
  HighsStatus return_status;
  HighsModelStatus model_status;
  double required_objective_function_value;
  double required_x0;
  double required_x1;

  HighsModel local_model;
  HighsLp& lp = local_model.lp_;
  HighsHessian& hessian = local_model.hessian_;
  const double inf = kHighsInf;

  // Oscar's edge case
  //
  // min x^2 + x = x(x + 1)

  lp.model_name_ = "qod";
  lp.num_col_ = 1;
  lp.num_row_ = 0;
  lp.col_cost_ = {1.0};
  lp.col_lower_ = {-inf};
  lp.col_upper_ = {inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.25;
  hessian.dim_ = lp.num_col_;
  hessian.q_start_ = {0, 1};
  hessian.q_index_ = {0};
  hessian.q_value_ = {2.0};

  Highs highs;
  const HighsModel& model = highs.getModel();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  const double& objective_function_value = info.objective_function_value;

  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  const bool pretty = true;
  if (dev_run) {
    printf("One variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", pretty);
  }

  required_objective_function_value = 0;
  required_x0 = -0.5;

  double alt_objective_function_value =
      model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);

  // Add a variable x1 with objective x1^2 - x1
  //
  // Add the variable
  highs.addCol(-1, -inf, inf, 0, NULL, NULL);
  if (dev_run) highs.writeModel("");

  // Cannot solve the model until the Hessian has been replaced
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kError);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kModelError);

  // Pass the new Hessian
  hessian.dim_ = 2;
  hessian.q_start_ = {0, 1, 2};
  hessian.q_index_ = {0, 1};
  hessian.q_value_ = {2.0, 2.0};
  return_status = highs.passHessian(hessian);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("Two variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", pretty);
  }

  required_objective_function_value = -0.25;
  required_x1 = 0.5;

  alt_objective_function_value = model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);

  // Illustrate methods for getting and changing the offset by getting
  // the current offset, shifting it by the current objective and
  // checking that the objective value is changed accordingly

  double offset;
  return_status = highs.getObjectiveOffset(offset);
  REQUIRE(return_status == HighsStatus::kOk);
  double dl_offset = -objective_function_value;
  offset += dl_offset;
  return_status = highs.changeObjectiveOffset(offset);
  required_objective_function_value += dl_offset;
  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);

  // Add the constraint 0.5 <= x0 + x1
  lp.a_index_ = {0, 1};
  lp.a_value_ = {1, 1};
  highs.addRow(0.5, inf, 2, &lp.a_index_[0], &lp.a_value_[0]);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("Two variable constrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", pretty);
  }

  required_objective_function_value = 0.125;
  required_x0 = -0.25;
  required_x1 = 0.75;
  alt_objective_function_value = model.objectiveValue(solution.col_value);
  REQUIRE(fabs(objective_function_value - alt_objective_function_value) <
          double_equal_tolerance);

  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);
}

TEST_CASE("test-qjh", "[qpsolver]") {
  // Test passing/reading and solving the problem qjh
  //
  // minimize -x_2 - 3x_3 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
  //
  // subject to x_1 + x_3 <= 2; x>=0
  HighsStatus return_status;
  HighsModelStatus model_status;
  double required_objective_function_value;

  HighsModel local_model;
  HighsLp& lp = local_model.lp_;
  HighsHessian& hessian = local_model.hessian_;
  const double inf = kHighsInf;
  // Start with an unconstrained QP
  lp.model_name_ = "qjh";
  lp.num_col_ = 3;
  lp.num_row_ = 0;
  lp.col_cost_ = {0.0, -1.0, -3.0};
  lp.col_lower_ = {-inf, -inf, -inf};
  lp.col_upper_ = {inf, inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  hessian.dim_ = lp.num_col_;

  //  hessian.format_ = HessianFormat::kSquare;
  //  hessian.q_start_ = {0, 2, 3, 5};
  //  hessian.q_index_ = {0, 2, 1, 0, 2};
  //  hessian.q_value_ = {2.0, -1.0, 0.2, -1.0, 2.0};

  hessian.format_ = HessianFormat::kTriangular;
  hessian.q_start_ = {0, 2, 3, 4};
  hessian.q_index_ = {0, 2, 1, 2};
  hessian.q_value_ = {2.0, -1.0, 0.2, 2.0};

  Highs highs;
  const HighsInfo& info = highs.getInfo();
  const double& objective_function_value = info.objective_function_value;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  required_objective_function_value = -5.50;
  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);

  if (dev_run) printf("Objective = %g\n", objective_function_value);
  if (dev_run) highs.writeSolution("", true);

  // Now with a constraint
  lp.num_row_ = 1;
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {2};
  lp.a_start_ = {0, 1, 1, 2};
  lp.a_index_ = {0, 0};
  lp.a_value_ = {1.0, 1.0};
  lp.format_ = MatrixFormat::kColwise;
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  required_objective_function_value = -5.25;
  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);

  if (dev_run) printf("Objective = %g\n", objective_function_value);
  if (dev_run) highs.writeSolution("", true);

  // Make the problem infeasible
  return_status = highs.changeColBounds(0, 3, inf);
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.changeColBounds(2, 3, inf);
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) highs.writeSolution("", true);
  model_status = highs.getModelStatus();
  if (dev_run)
    printf("Infeasible QP status: %s\n",
           highs.modelStatusToString(model_status).c_str());
  REQUIRE(model_status == HighsModelStatus::kInfeasible);

  return_status = highs.clearModel();

  std::string filename;
  for (HighsInt test_k = 0; test_k < 2; test_k++) {
    if (test_k == 0) {
      filename = std::string(HIGHS_DIR) + "/check/instances/qjh.mps";
    } else if (test_k == 1) {
      filename = std::string(HIGHS_DIR) + "/check/instances/qjh_quadobj.mps";
    } else {
      filename = std::string(HIGHS_DIR) + "/check/instances/qjh_qmatrix.mps";
    }

    return_status = highs.readModel(filename);
    REQUIRE(return_status == HighsStatus::kOk);
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(fabs(objective_function_value - required_objective_function_value) <
            double_equal_tolerance);
    return_status = highs.clearModel();
  }
}
