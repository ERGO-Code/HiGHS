#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

bool smallDoubleDifference(double v0, double v1) {
  double difference = std::fabs(v0 - v1);
  //  printf("smallDoubleDifference = %g\n", difference);
  return difference < 1e-4;
}

TEST_CASE("multi-objective", "[util]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {-kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {18, 8, 14};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {3, 1, 1, 1, 1, 2};
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.passModel(lp);

  HighsLinearObjective linear_objective;
  std::vector<HighsLinearObjective> linear_objectives;

  // Begin with an illegal linear objective
  if (dev_run) printf("\nPass illegal linear objective\n");
  linear_objective.weight = -1;
  linear_objective.offset = -1;
  linear_objective.coefficients = {2, 1, 0};
  linear_objective.abs_tolerance = 0.0;
  linear_objective.rel_tolerance = 1.0;
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kError);

  // Now legalise the linear objective so LP has nonunique optimal
  // solutions on the line joining (2, 6) and (5, 3)
  if (dev_run) printf("\nPass legal linear objective\n");
  linear_objective.coefficients = {1, 1};
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getInfo().objective_function_value, -7));
  // Save the linear objective for the next
  linear_objectives.push_back(linear_objective);

  // Add a second linear objective with a very small minimization
  // weight that should push the optimal solution to (2, 6)
  if (dev_run) printf("\nPass second linear objective\n");
  linear_objective.weight = 1e-4;
  linear_objective.offset = 0;
  linear_objective.coefficients = {1, 0};
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));
  linear_objectives.push_back(linear_objective);

  if (dev_run) printf("\nClear and pass two linear objectives\n");
  REQUIRE(h.clearLinearObjectives() == HighsStatus::kOk);
  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);
  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));

  // Set illegal priorities - that can be passed OK since
  // blend_multi_objectives = true
  if (dev_run)
    printf(
        "\nSetting priorities that will be illegal when using lexicographic "
        "optimization\n");
  linear_objectives[0].priority = 0;
  linear_objectives[1].priority = 0;
  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);

  // Now test lexicographic optimization
  h.setOptionValue("blend_multi_objectives", false);

  if (dev_run) printf("\nLexicographic using illegal priorities\n");
  REQUIRE(h.run() == HighsStatus::kError);

  if (dev_run)
    printf(
        "\nSetting priorities that are illegal now blend_multi_objectives = "
        "false\n");
  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kError);

  if (dev_run)
    printf("\nSetting legal priorities for blend_multi_objectives = false\n");
  linear_objectives[0].priority = 10;
  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);

  if (dev_run) printf("\nLexicographic using existing multi objective data\n");
  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));

  // Back to blending
  h.setOptionValue("blend_multi_objectives", true);
  //  h.setOptionValue("output_flag", true);
  REQUIRE(h.clearLinearObjectives() == HighsStatus::kOk);
  linear_objectives[0].coefficients = {1.0001, 1};
  linear_objectives[0].abs_tolerance = -1e-5;
  linear_objectives[0].rel_tolerance = 0.95;

  if (dev_run) printf("\nBlending: Illegal abs_tolerance \n");
  REQUIRE(h.passLinearObjectives(1, linear_objectives.data()) ==
          HighsStatus::kError);
  linear_objectives[0].abs_tolerance = 1e-5;

  if (dev_run) printf("\nBlending: Illegal rel_tolerance \n");
  REQUIRE(h.passLinearObjectives(1, linear_objectives.data()) ==
          HighsStatus::kError);
  linear_objectives[0].rel_tolerance = 1.05;

  linear_objectives[1].weight = 1e-3;
  if (dev_run)
    printf(
        "\nBlending: first solve objective just giving unique optimal "
        "solution\n");
  REQUIRE(h.passLinearObjectives(1, linear_objectives.data()) ==
          HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);

  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);

  // Back to lexicographic optimization
  h.setOptionValue("blend_multi_objectives", false);

  if (dev_run) printf("\nLexicographic using non-trivial tolerances\n");
  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);

  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 4.9));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 3.1));

  linear_objectives[0].abs_tolerance = kHighsInf;

  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);

  //  printf("Solution = [%23.18g, %23.18g]\n", h.getSolution().col_value[0],
  //  h.getSolution().col_value[1]);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 1.30069));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6.34966));
}
