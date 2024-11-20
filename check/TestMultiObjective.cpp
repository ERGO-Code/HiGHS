#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

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
  linear_objective.weight = -1;
  linear_objective.offset = -1;
  linear_objective.coefficients = {2, 1, 0};
  linear_objective.abs_tolerance = 0.0;
  linear_objective.rel_tolerance = 1.0;
  linear_objective.priority = 0;
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kError);

  // Now legalise the linear objective so LP has nonunique optimal
  // solutions on the line joining (2, 6) and (5, 3)
  linear_objective.coefficients = {1, 1};
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  h.run();
  h.writeSolution("", kSolutionStylePretty);
  REQUIRE(h.getInfo().objective_function_value == -7);
  // Save the linear objective for the next
  linear_objectives.push_back(linear_objective);

  // Add a second linear objective with a very small minimization
  // weight that should push the optimal solution to (2, 6)
  linear_objective.weight = 1e-4;
  linear_objective.offset = 0;
  linear_objective.coefficients = {-1, 0};
  REQUIRE(h.addLinearObjective(linear_objective) == HighsStatus::kOk);

  h.run();
  h.writeSolution("", kSolutionStylePretty);
  //  REQUIRE(h.getSolution().col_value[0] == 2);
  //  REQUIRE(h.getSolution().col_value[1] == 6);
}
