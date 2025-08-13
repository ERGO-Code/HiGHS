#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderLp.h"

const bool dev_run = false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

bool okValueDifference(const double& v_test, const double& v_true) {
  double difference = fabs(v_test - v_true) / std::max(1.0, fabs(v_true));
  return difference < double_equal_tolerance;
}

void testPrimalDualObjective(Highs& h,
                             const double required_objective_function_value) {
  const HighsInfo& info = h.getInfo();
  double objective_function_value = info.objective_function_value;
  const HighsSolution& solution = h.getSolution();
  double alt_objective_function_value =
      h.getModel().objectiveValue(solution.col_value);
  double dual_objective_function_value;
  REQUIRE(h.getDualObjectiveValue(dual_objective_function_value) ==
          HighsStatus::kOk);
  double alt_objective_function_value_error =
      fabs(objective_function_value - alt_objective_function_value);
  if (dev_run)
    printf(
        "(Primal, Alt, Dual) objective = (%17.10g, %17.10g, %17.10g) alt error "
        "= %17.10g; P-D error = %17.10g\n",
        objective_function_value, alt_objective_function_value,
        dual_objective_function_value, alt_objective_function_value_error,
        info.primal_dual_objective_error);
  REQUIRE(okValueDifference(objective_function_value,
                            required_objective_function_value));
  REQUIRE(okValueDifference(dual_objective_function_value,
                            required_objective_function_value));
  REQUIRE(okValueDifference(alt_objective_function_value,
                            required_objective_function_value));
  double optimality_tolerance;
  h.getOptionValue("optimality_tolerance", optimality_tolerance);
  REQUIRE(fabs(info.primal_dual_objective_error) < optimality_tolerance);
}

TEST_CASE("qp-unbounded", "[qpsolver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qpunbounded.lp";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("qp-infeasible", "[qpsolver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qpinfeasible.lp";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("qpsolver", "[qpsolver]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  double required_objective_function_value;
  double required_x0;
  double required_x1;
  double required_x2;
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qptestnw.lp";

  required_objective_function_value = -6.45;
  required_x0 = 1.4;
  required_x1 = 1.7;

  const double required_col_dual0 = 0;
  const double required_col_dual1 = 0;
  const double required_row_dual0 = 0.8;
  const double required_row_dual1 = 0;

  // At the optimal solution g-Qx = [0.8, -1.6] with only constraint 0
  // active. It has normal [1, -2], so dual of 0.8 is correct

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsModel& model = highs.getModel();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  const double& objective_function_value = info.objective_function_value;

  HighsStatus return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  testPrimalDualObjective(highs, required_objective_function_value);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);

  REQUIRE(fabs(solution.col_dual[0] - required_col_dual0) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.col_dual[1] - required_col_dual1) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.row_dual[0] - required_row_dual0) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.row_dual[1] - required_row_dual1) <
          double_equal_tolerance);
  REQUIRE(fabs(solution.row_dual[1] - required_row_dual1) <
          double_equal_tolerance);

  if (dev_run) highs.writeSolution("", 1);

  // Check with qjh.mps
  filename = std::string(HIGHS_DIR) + "/check/instances/qjh.mps";
  required_objective_function_value = -5.25;
  required_x0 = 0.5;
  required_x1 = 5.0;
  required_x2 = 1.5;

  return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) printf("Objective = %g\n", objective_function_value);

  testPrimalDualObjective(highs, required_objective_function_value);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - required_x2) < double_equal_tolerance);
  REQUIRE(return_status == HighsStatus::kOk);

  // Test writeModel by writing out qjh.mps...
  filename = test_name + ".mps";
  highs.writeModel(filename);

  // ... and reading it in again
  return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) printf("Objective = %g\n", objective_function_value);

  testPrimalDualObjective(highs, required_objective_function_value);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - required_x2) < double_equal_tolerance);
  std::remove(filename.c_str());

  // Test that attempting to solve MIQP yields error
  HighsInt num_col = highs.getNumCol();
  std::vector<HighsVarType> integrality;
  integrality.assign(num_col, HighsVarType::kInteger);
  REQUIRE(highs.changeColsIntegrality(0, num_col - 1, integrality.data()) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kError);

  // Test that attempting to solve MIQP relaxation is OK
  highs.setOptionValue("solve_relaxation", true);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  highs.resetGlobalScheduler(true);
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
  // min 1/4 + x^2 + x = 1/4 + x(x + 1)
  //
  // x* = -1/2; f* = 0

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
  highs.setOptionValue("output_flag", dev_run);
  const HighsModel& model = highs.getModel();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  const double& objective_function_value = info.objective_function_value;

  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("\nOne variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", kSolutionStylePretty);
  }

  required_objective_function_value = 0;
  required_x0 = -0.5;

  testPrimalDualObjective(highs, required_objective_function_value);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);

  // Add a variable x1 with objective x1^2 - x1
  //
  // Add the variable
  highs.addCol(-1, -inf, inf, 0, NULL, NULL);
  if (dev_run) highs.writeModel("");

  // Can solve the model before the Hessian has been replaced
  if (dev_run)
    printf(
        "\nTwo variable unconstrained QP with semi-definite Hessian - is "
        "unbounded\n");
  // Reinstate the QP regularization since the Hessian is only semi-definite
  REQUIRE(highs.setOptionValue("qp_regularization_value",
                               kHessianRegularizationValue) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kUnbounded);

  // Pass the new Hessian
  hessian.dim_ = 2;
  hessian.start_ = {0, 1, 2};
  hessian.index_ = {0, 1};
  hessian.value_ = {2.0, 2.0};
  return_status = highs.passHessian(hessian);
  REQUIRE(return_status == HighsStatus::kOk);

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("Two variable unconstrained QP: objective = %g; solution:\n",
           objective_function_value);
    highs.writeSolution("", kSolutionStylePretty);
  }

  required_objective_function_value = -0.25;
  required_x1 = 0.5;

  testPrimalDualObjective(highs, required_objective_function_value);

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
  highs.addRow(0.5, inf, 2, lp.a_matrix_.index_.data(),
               lp.a_matrix_.value_.data());
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

  testPrimalDualObjective(highs, required_objective_function_value);

  REQUIRE(fabs(solution.col_value[0] - required_x0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - required_x1) < double_equal_tolerance);

  highs.resetGlobalScheduler(true);
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
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();
  const double& objective_function_value = info.objective_function_value;
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeModel("");
  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  required_objective_function_value = -5.50;
  testPrimalDualObjective(highs, required_objective_function_value);

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
  testPrimalDualObjective(highs, required_objective_function_value);

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

    testPrimalDualObjective(highs, required_objective_function_value);

    return_status = highs.clearModel();
  }

  highs.resetGlobalScheduler(true);
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
  highs.setOptionValue("output_flag", dev_run);

  // Should load OK
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  // Run should fail since objective is non-convex
  REQUIRE(highs.run() == HighsStatus::kError);

  highs.resetGlobalScheduler(true);
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
  highs.setOptionValue("output_flag", dev_run);

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
  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  const double required_objective_function_value = 1.25;

  testPrimalDualObjective(highs, required_objective_function_value);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 0.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 1.5) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - 2.5) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_dual[0] + 1.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.row_dual[0] + 0.5) < double_equal_tolerance);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-semi-definite0", "[qpsolver]") {
  HighsStatus return_status;

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
  highs.setOptionValue("output_flag", dev_run);
  return_status = highs.passModel(local_model);
  REQUIRE(return_status == HighsStatus::kOk);

  //  highs.writeModel("semi-definite.mps");

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-semi-definite1", "[qpsolver]") {
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
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  hessian.dim_ = lp.num_col_;
  hessian.start_ = {0, 2, 3};
  hessian.index_ = {0, 1, 1};
  hessian.value_ = {1.0, -1.0, 1.0};
  REQUIRE(highs.passHessian(hessian) == HighsStatus::kOk);

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  const double required_objective_function_value = -1.5;

  testPrimalDualObjective(highs, required_objective_function_value);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1]) < double_equal_tolerance);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-semi-definite2", "[qpsolver]") {
  HighsLp lp;
  HighsHessian hessian;

  lp.model_name_ = "semi-definite";
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {0.0, -1.0};
  lp.col_lower_ = {-inf, -inf};
  lp.col_upper_ = {inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1.0, 1.0};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Illustrates how final entries of 0 on the diagonal are handled
  hessian.dim_ = lp.num_col_;
  hessian.start_ = {0, 1, 1};
  hessian.index_ = {0};
  hessian.value_ = {1.0};
  REQUIRE(highs.passHessian(hessian) == HighsStatus::kOk);

  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  const double required_objective_function_value = -1.5;

  testPrimalDualObjective(highs, required_objective_function_value);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] + 1) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 2) < double_equal_tolerance);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-qp-modification", "[qpsolver]") {
  //  HighsStatus return_status;
  //  HighsModelStatus model_status;
  //  double required_objective_function_value;

  HighsModel model;

  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {1.0, -1.0};
  lp.col_lower_ = {-inf, -inf};
  lp.col_upper_ = {inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1.0, 1.0};
  hessian.dim_ = 1;
  hessian.start_ = {0, 1};
  hessian.value_ = {1.0};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsModel& incumbent_model = highs.getModel();
  // Cannot have Hessian with index exceeding hessian.dim_-1
  hessian.index_ = {1};
  REQUIRE(highs.passModel(model) == HighsStatus::kError);

  // Correct the Hessian index
  hessian.index_[0] = 0;
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nNow solve the QP\n\n");
    incumbent_model.hessian_.print();
  }
  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);

  highs.run();
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);
  // Add a new variables and ensure that the Hessian dimension is correct
  std::vector<HighsInt> index = {0};
  std::vector<double> value = {1};
  REQUIRE(highs.addCol(-1, 0, 1, 1, index.data(), value.data()) ==
          HighsStatus::kOk);
  REQUIRE((incumbent_model.hessian_.dim_ == 0 ||
           incumbent_model.hessian_.dim_ == incumbent_model.lp_.num_col_));
  if (dev_run) {
    printf("\nNow solve the QP after adding new variable\n\n");
    incumbent_model.hessian_.print();
  }
  highs.run();
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  HighsInt dim = incumbent_model.lp_.num_col_;
  std::vector<double> arg0;
  std::vector<double> arg1;
  std::vector<double> result0;
  std::vector<double> result1;
  arg0.resize(dim);
  HighsRandom random;
  for (HighsInt iCol = 0; iCol < dim; iCol++) arg0[iCol] = random.fraction();
  HighsHessian hessian0 = incumbent_model.hessian_;
  arg1 = arg0;

  // Deleting column 1 removes no nonzeros from the Hessian
  HighsInt delete_col = 1;
  REQUIRE(highs.deleteCols(delete_col, delete_col) == HighsStatus::kOk);
  REQUIRE((incumbent_model.hessian_.dim_ == 0 ||
           incumbent_model.hessian_.dim_ == incumbent_model.lp_.num_col_));
  if (dev_run) {
    printf("\nNow solve the QP after deleting column 1\n\n");
    incumbent_model.hessian_.print();
  }
  highs.run();
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  dim--;
  for (HighsInt iCol = delete_col; iCol < dim; iCol++)
    arg1[iCol] = arg1[iCol + 1];
  arg0[delete_col] = 0;
  hessian0.product(arg0, result0);
  for (HighsInt iCol = delete_col; iCol < dim; iCol++)
    result0[iCol] = result0[iCol + 1];

  arg1.resize(dim);
  incumbent_model.hessian_.product(arg1, result1);
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    REQUIRE(result0[iCol] == result1[iCol]);

  // Deleting column 0 removes only nonzero from the Hessian, so problem is an
  // LP
  delete_col = 0;
  REQUIRE(highs.deleteCols(delete_col, delete_col) == HighsStatus::kOk);
  REQUIRE((incumbent_model.hessian_.dim_ == 0 ||
           incumbent_model.hessian_.dim_ == incumbent_model.lp_.num_col_));
  if (dev_run) {
    printf("\nNow solve the LP after deleting column 0\n\n");
    incumbent_model.hessian_.print();
  }
  highs.run();
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-qp-delete-col", "[qpsolver]") {
  HighsModel model;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.col_cost_ = {-2, -1, 0, 1, 2};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {inf, inf, inf, inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 5};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 4};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 1};
  hessian.dim_ = 5;
  hessian.start_ = {0, 4, 7, 10, 11, 12};
  hessian.index_ = {0, 1, 3, 4, 1, 2, 4, 2, 3, 4, 3, 4};
  hessian.value_ = {11, 21, 41, 51, 22, 32, 52, 33, 43, 53, 44, 55};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsModel& incumbent_model = highs.getModel();
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  if (dev_run) incumbent_model.hessian_.print();

  HighsInt dim = incumbent_model.lp_.num_col_;
  std::vector<double> arg0;
  std::vector<double> arg1;
  std::vector<double> result0;
  std::vector<double> result1;
  arg0.resize(dim);
  HighsRandom random;
  for (HighsInt iCol = 0; iCol < dim; iCol++) arg0[iCol] = random.fraction();
  HighsHessian hessian0 = incumbent_model.hessian_;
  arg1 = arg0;

  std::vector<HighsInt> set = {1, 3};
  REQUIRE(highs.deleteCols(2, set.data()) == HighsStatus::kOk);
  if (dev_run) incumbent_model.hessian_.print();

  arg1[1] = arg1[2];
  arg1[2] = arg1[4];
  arg0[1] = 0;
  arg0[3] = 0;
  hessian0.product(arg0, result0);
  result0[1] = result0[2];
  result0[2] = result0[4];

  dim = 3;
  arg1.resize(dim);
  incumbent_model.hessian_.product(arg1, result1);
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    REQUIRE(result0[iCol] == result1[iCol]);

  dim = 100;
  lp.clear();
  lp.num_col_ = dim;
  lp.num_row_ = 1;

  lp.col_cost_.resize(dim);
  lp.col_lower_.resize(dim);
  lp.col_upper_.resize(dim);
  lp.a_matrix_.index_.resize(dim);
  lp.a_matrix_.value_.resize(dim);

  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    lp.col_cost_[iCol] = double(iCol);
    lp.col_lower_[iCol] = 0;
    lp.col_upper_[iCol] = inf;
    lp.a_matrix_.index_[iCol] = iCol;
    lp.a_matrix_.value_[iCol] = 1;
  }
  lp.a_matrix_.start_.push_back(dim);

  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;

  hessian.clear();
  hessian.dim_ = dim;
  std::vector<double> hessian_col(dim);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    hessian_col.assign(dim, 0);
    HighsInt kmax = std::max(HighsInt(1), (dim - iCol) / 3);
    hessian_col[iCol] = dim * iCol + iCol;
    if (iCol < dim - 1) {
      for (HighsInt k = 0; k < kmax; k++) {
        HighsInt iRow = iCol + 1 + random.integer(dim - iCol - 1);
        assert(iRow >= 0);
        assert(iRow < dim);
        hessian_col[iRow] = dim * iCol + iRow;
      }
    }
    for (HighsInt iRow = iCol; iRow < dim; iRow++) {
      if (hessian_col[iRow]) {
        hessian.index_.push_back(iRow);
        hessian.value_.push_back(hessian_col[iRow]);
      }
    }
    hessian.start_.push_back(HighsInt(hessian.index_.size()));
  }

  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  if (dev_run && dim < 20) incumbent_model.hessian_.print();

  arg0.resize(dim);
  for (HighsInt iCol = 0; iCol < dim; iCol++) arg0[iCol] = random.fraction();
  hessian0 = incumbent_model.hessian_;
  arg1 = arg0;

  std::vector<HighsInt> mask;
  mask.assign(dim, 0);

  HighsInt kmax = std::max(HighsInt(1), dim / 3);
  for (HighsInt k = 0; k < kmax; k++) {
    HighsInt iRow = random.integer(dim);
    assert(iRow >= 0);
    assert(iRow < dim);
    mask[iRow] = 1;
  }
  highs.deleteCols(mask.data());
  if (dev_run && dim < 20) incumbent_model.hessian_.print();

  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iRow = mask[iCol];
    if (iRow < 0) {
      arg0[iCol] = 0;
    } else {
      arg1[iRow] = arg1[iCol];
    }
  }
  hessian0.product(arg0, result0);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iRow = mask[iCol];
    if (iRow >= 0) result0[iRow] = result0[iCol];
  }
  dim = incumbent_model.hessian_.dim_;
  arg1.resize(dim);
  incumbent_model.hessian_.product(arg1, result1);

  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    REQUIRE(fabs(result0[iCol] - result1[iCol]) < 1e-8);
  }
}

TEST_CASE("test-qp-hot-start", "[qpsolver]") {
  // Test hot start
  HighsStatus return_status;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();

  double required_objective_function_value = 0;
  for (HighsInt k = 0; k < 2; k++) {
    if (dev_run)
      printf(
          "\n"
          "===================\n"
          "Hot start test %d\n"
          "===================\n",
          int(k));
    if (k == 1) {
      const std::string filename =
          std::string(HIGHS_DIR) + "/check/instances/primal1.mps";
      REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
      required_objective_function_value = -0.035012965733477348;
    } else if (k == 2) {
      // Not currently tested
      const std::string filename =
          std::string(HIGHS_DIR) + "/check/instances/qptestnw.lp";
      REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
    } else {
      HighsModel model;
      model.lp_.num_col_ = 2;
      model.lp_.num_row_ = 1;
      model.lp_.col_cost_ = {-2, -2};
      model.lp_.col_lower_ = {-inf, -inf};
      model.lp_.col_upper_ = {inf, inf};
      model.lp_.row_lower_ = {1};
      model.lp_.row_upper_ = {inf};
      model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
      model.lp_.a_matrix_.start_ = {0, 2};
      model.lp_.a_matrix_.index_ = {0, 1};
      model.lp_.a_matrix_.value_ = {1, 1};
      model.hessian_.dim_ = 2;
      model.hessian_.start_ = {0, 1, 2};
      model.hessian_.index_ = {0, 1};
      model.hessian_.value_ = {2, 2};
      REQUIRE(highs.passModel(model) == HighsStatus::kOk);
      required_objective_function_value = -2;
    }
    // Zero the QP regularization so "true" solution is obtained
    REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
            HighsStatus::kOk);
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);

    testPrimalDualObjective(highs, required_objective_function_value);

    if (dev_run) highs.writeSolution("", 1);

    HighsBasis basis = highs.getBasis();
    HighsSolution solution = highs.getSolution();
    if (dev_run) printf("Saved basis has validity = %d\n", basis.valid);

    if (dev_run)
      printf(
          "================\n"
          "Hot start re-run\n"
          "================\n");
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(info.qp_iteration_count == 0);

    if (dev_run)
      printf(
          "===========================\n"
          "Hot start using saved basis\n"
          "===========================\n");
    highs.setBasis(basis);
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(info.qp_iteration_count == 0);

    // QP Hot start needs a saved solution as well as a basis after
    // clearSolver()
    if (dev_run)
      printf(
          "==============================================================\n"
          "Hot start using saved basis and solution after clearing solver\n"
          "==============================================================\n");
    highs.clearSolver();
    highs.setSolution(solution);
    highs.setBasis(basis);
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(info.qp_iteration_count == 0);
    /*
        if (dev_run)
    printf("=================================================\n"
           "Hot start using saved basis after clearing solver\n"
           "=================================================\n");
    highs.clearSolver();
    highs.setBasis(basis);
    return_status = highs.run();
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(info.qp_iteration_count == 0);
    */
    // QP Hot start needs a saved solution as well as a basis after
    // clearSolver()
    if (dev_run)
      printf(
          "==============================================================\n"
          "Hot start using alien basis and solution after clearing solver\n"
          "==============================================================\n");
    highs.clearSolver();
    highs.setSolution(solution);
    basis.alien = true;
    highs.setBasis(basis);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(info.qp_iteration_count == 0);
    if (k == 0) {
      // Modify the constraint so that the solution and basis are not
      // feasible and one iteration is needed
      REQUIRE(highs.changeCoeff(0, 1, 2.0) == HighsStatus::kOk);
      REQUIRE(highs.changeRowBounds(0, 4.0, kHighsInf) == HighsStatus::kOk);
      highs.clearSolver();
      basis.alien = false;
      highs.setBasis(basis);
      highs.setSolution(solution);
      return_status = highs.run();
      REQUIRE(info.qp_iteration_count == 1);
    }
  }

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-qp-terminations", "[qpsolver]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/qptestnw.lp";
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);

  highs.setOptionValue("qp_iteration_limit", 1);
  REQUIRE(highs.run() == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kIterationLimit);
  highs.clearSolver();
  highs.setOptionValue("qp_iteration_limit", kHighsIInf);

  highs.setOptionValue("time_limit", 0);
  REQUIRE(highs.run() == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kTimeLimit);
  highs.setOptionValue("time_limit", kHighsInf);

  filename = std::string(HIGHS_DIR) + "/check/instances/primal1.mps";
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);

  highs.setOptionValue("qp_nullspace_limit", 1);
  REQUIRE(highs.run() == HighsStatus::kError);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kSolveError);
  highs.setOptionValue("qp_nullspace_limit", 4000);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("rowless-qp", "[qpsolver]") {
  HighsModel model;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  lp.num_col_ = 2;
  lp.num_row_ = 0;
  lp.col_cost_ = {0, -3};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  hessian.dim_ = 2;
  hessian.start_ = {0, 2, 3};
  hessian.index_ = {0, 1, 1};
  hessian.value_ = {2, 1, 2};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  // Zero the QP regularization so "true" solution is obtained
  REQUIRE(highs.setOptionValue("qp_regularization_value", 0) ==
          HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.writeSolution("", kSolutionStylePretty) ==
          HighsStatus::kWarning);

  const double required_objective_function_value = -2.25;

  testPrimalDualObjective(highs, required_objective_function_value);

  const std::vector<double>& col_value = highs.getSolution().col_value;
  if (dev_run)
    printf("Solution (%24.18g, %24.18g)\n", col_value[0], col_value[1]);
  double dl_solution = std::fabs(col_value[0]);
  REQUIRE(dl_solution < 1e-6);
  dl_solution = std::fabs(col_value[1] - 1.5);
  REQUIRE(dl_solution < 1e-6);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("2489", "[qpsolver]") {
  // This QP is
  //
  // Min, x^2/2 + x
  //
  // 0*x + 0*y <= 1
  //
  // -10 <= (x, y) <= 10
  //
  // Hence it has a constraint, but its coefficients are zero
  Highs h;
  //  h.setOptionValue("output_flag", dev_run);
  assert(h.setOptionValue("log_dev_level", 3) == HighsStatus::kOk);
  assert(h.addCol(1.0, -10.0, 10.0, 0, NULL, NULL) == HighsStatus::kOk);
  assert(h.addCol(0.0, -10.0, 10.0, 0, NULL, NULL) == HighsStatus::kOk);
  assert(h.addRow(0.0, 0.0, 0, NULL, NULL) == HighsStatus::kOk);
  HighsHessian hessian;
  hessian.dim_ = 1;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 1};
  hessian.index_ = {0};
  hessian.value_ = {1.0};
  assert(h.passHessian(hessian) == HighsStatus::kOk);
  assert(h.run() == HighsStatus::kOk);

  h.resetGlobalScheduler(true);
}
