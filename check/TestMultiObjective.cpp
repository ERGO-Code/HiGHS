#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

bool smallDoubleDifference(double v0, double v1) {
  return std::fabs(v0 - v1) < 1e-12;
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
  printf("\nPass illegal linear objective\n");
  linear_objective.weight = -1;
  linear_objective.offset = -1;
  linear_objective.coefficients = {2, 1, 0};
  linear_objective.abs_tolerance = 0.0;
  linear_objective.rel_tolerance = 1.0;
  linear_objective.priority = 10;
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kError);

  // Now legalise the linear objective so LP has nonunique optimal
  // solutions on the line joining (2, 6) and (5, 3)
  printf("\nPass legal linear objective\n");
  linear_objective.coefficients = {1, 1};
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getInfo().objective_function_value, -7));
  // Save the linear objective for the next
  linear_objectives.push_back(linear_objective);

  // Add a second linear objective with a very small minimization
  // weight that should push the optimal solution to (2, 6)
  printf("\nPass second linear objective\n");
  linear_objective.weight = 1e-4;
  linear_objective.offset = 0;
  linear_objective.coefficients = {1, 0};
  linear_objective.priority = 0;
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));
  linear_objectives.push_back(linear_objective);

  printf("\nClear and pass two linear objectives\n");
  REQUIRE(h.clearLinearObjectives() == HighsStatus::kOk);
  REQUIRE(h.passLinearObjectives(2, linear_objectives.data()) ==
          HighsStatus::kOk);
  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));

  // Now test lexicographic optimization
  h.setOptionValue("blend_multi_objectives", false);
  printf("\nLexicographic using existing multi objective data\n");
  REQUIRE(h.run() == HighsStatus::kOk);
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[0], 2));
  REQUIRE(smallDoubleDifference(h.getSolution().col_value[1], 6));
}
