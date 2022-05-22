#include <cstdio>

#include "FilereaderLp.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

TEST_CASE("qp-unbounded", "[qpsolver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qpunbounded.lp";

  Highs highs;
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
}

TEST_CASE("qp-infeasible", "[qpsolver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qpinfeasible.lp";

  Highs highs;
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}

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
  std::remove(filename.c_str());
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
  hessian.start_ = {0, 1};
  hessian.index_ = {0};
  hessian.value_ = {2.0};

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

  if (dev_run) {
    printf("One variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", kSolutionStylePretty);
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
  hessian.start_ = {0, 1, 2};
  hessian.index_ = {0, 1};
  hessian.value_ = {2.0, 2.0};
  return_status = highs.passHessian(hessian);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("Two variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", kSolutionStylePretty);
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
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1, 1};
  highs.addRow(0.5, inf, 2, &lp.a_matrix_.index_[0], &lp.a_matrix_.value_[0]);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("Two variable constrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", kSolutionStylePretty);
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
  //  hessian.start_ = {0, 2, 3, 5};
  //  hessian.index_ = {0, 2, 1, 0, 2};
  //  hessian.value_ = {2.0, -1.0, 0.2, -1.0, 2.0};

  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 2, 3, 4};
  hessian.index_ = {0, 2, 1, 2};
  hessian.value_ = {2.0, -1.0, 0.2, 2.0};

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
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  // Now with a constraint
  lp.num_row_ = 1;
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {2};
  lp.a_matrix_.start_ = {0, 1, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  required_objective_function_value = -5.25;
  REQUIRE(fabs(objective_function_value - required_objective_function_value) <
          double_equal_tolerance);

  if (dev_run) printf("Objective = %g\n", objective_function_value);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  // Make the problem infeasible
  return_status = highs.changeColBounds(0, 3, inf);
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.changeColBounds(2, 3, inf);
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) highs.writeSolution("", kSolutionStylePretty);
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

TEST_CASE("test-min-negative-definite", "[qpsolver]") {
  HighsModel model;
  model.lp_.model_name_ = "min-negative-definite";
  model.lp_.num_col_ = 1;
  model.lp_.num_row_ = 0;
  model.lp_.col_cost_ = {0};
  model.lp_.col_lower_ = {0};
  model.lp_.col_upper_ = {inf};
  model.hessian_.dim_ = model.lp_.num_col_;
  model.hessian_.start_ = {0, 1};
  model.hessian_.index_ = {0};
  model.hessian_.value_ = {-0.5};
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  // Should load OK
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  // Run should fail since objective is non-convex
  REQUIRE(highs.run() == HighsStatus::kError);
}

TEST_CASE("test-max-negative-definite", "[qpsolver]") {
  HighsLp lp;
  HighsHessian hessian;
  // Construct a QP with negative definite Hessian
  const double col_lower = 0;
  const double col_upper = inf;  // 10.0;

  lp.model_name_ = "max-negative-definite";
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.offset_ = -1.0;
  lp.col_cost_ = {1.0, 1.0, 2.0};
  lp.col_lower_ = {col_lower, col_lower, col_lower};
  lp.col_upper_ = {col_upper, col_upper, col_upper};
  lp.row_lower_ = {4};
  lp.row_upper_ = {inf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3};
  lp.a_matrix_.index_ = {0, 1, 2};
  lp.a_matrix_.value_ = {1.0, 1.0, 1.0};
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  hessian.dim_ = lp.num_col_;
  hessian.start_ = {0, 2, 3, 4};
  hessian.index_ = {0, 2, 1, 2};
  hessian.value_ = {-2.0, -1, -1.0, -1};
  if (dev_run) {
    printf("\nNegative definite Hessian\n");
    hessian.print();
  }
  // Should load OK
  REQUIRE(highs.passHessian(hessian) == HighsStatus::kOk);
  // Make the problem a maximization
  REQUIRE(highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  const double required_objective_function_value = 1.25;
  REQUIRE(fabs(highs.getObjectiveValue() - required_objective_function_value) <
          double_equal_tolerance);
  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 0.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 1.5) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - 2.5) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_dual[0] + 1.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.row_dual[0] + 0.5) < double_equal_tolerance);
}

TEST_CASE("test-semi-definite0", "[qpsolver]") {
  HighsStatus return_status;
  HighsModelStatus model_status;
  double required_objective_function_value;

  HighsModel local_model;
  HighsLp& lp = local_model.lp_;
  HighsHessian& hessian = local_model.hessian_;
  // Construct a QP with positive semi-definite Hessian
  const double col_lower = 0;
  const double col_upper = inf;  // 10.0;

  lp.model_name_ = "semi-definite";
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.col_cost_ = {1.0, 1.0, 2.0};
  lp.col_lower_ = {col_lower, col_lower, col_lower};
  lp.col_upper_ = {col_upper, col_upper, col_upper};
  lp.row_lower_ = {2};
  lp.row_upper_ = {inf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0, 1.0};
  hessian.dim_ = lp.num_col_;
  hessian.start_ = {0, 2, 2, 3};
  hessian.index_ = {0, 2, 2};
  hessian.value_ = {2.0, -1.0, 1.0};

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);

  //  highs.writeModel("semi-definite.mps");

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);
}

TEST_CASE("test-semi-definite1", "[qpsolver]") {
  HighsStatus return_status;
  HighsModelStatus model_status;
  double required_objective_function_value;

  HighsLp lp;
  HighsHessian hessian;

  lp.model_name_ = "semi-definite";
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-2.0, 0.0};
  lp.col_lower_ = {-inf, -inf};
  lp.col_upper_ = {inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {1};
  lp.row_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1.0, 1.0};

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  hessian.dim_ = lp.num_col_;
  hessian.start_ = {0, 2, 3};
  hessian.index_ = {0, 1, 1};
  hessian.value_ = {1.0, -1.0, 1.0};
  REQUIRE(highs.passHessian(hessian) == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);
  const HighsSolution& solution = highs.getSolution();
  const double objective_function_value = highs.getObjectiveValue();
  REQUIRE(fabs(objective_function_value + 1.5) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[0] - 1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1]) < double_equal_tolerance);
}
