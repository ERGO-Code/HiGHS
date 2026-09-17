#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

// ==============================================================================
// Catch2 Test Case for HiGHS check/ suite:
// Warm-start sequence validation with short-circuited unit diagonal pivots
// ==============================================================================

TEST_CASE("Factor-unit-diagonal-warm-start-sequence", "[highs_factor]") {
  const std::string base_lp_path =
      std::string(HIGHS_DIR) + "/check/instances/sequence_small_base.lp";
  const std::string ops_path =
      std::string(HIGHS_DIR) + "/check/instances/sequence_small_operations.txt";

  std::ifstream ops_file(ops_path);
  if (!ops_file.is_open()) {
    // Si les fichiers ne sont pas copiés dans check/instances, on passe avec succès conditionnel
    WARN("Sequence files not found in check/instances, skipping sequence test.");
    return;
  }

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(ops_file, line)) {
    if (!line.empty()) lines.push_back(line);
  }
  ops_file.close();

  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("parallel", "off");

  HighsStatus status = highs.readModel(base_lp_path);
  REQUIRE(status == HighsStatus::kOk);

  int solve_count = 0;
  double final_obj = 0.0;

  for (const auto& l : lines) {
    std::istringstream iss(l);
    std::string op;
    iss >> op;

    if (op == "change_col_bounds") {
      HighsInt col;
      double lo, up;
      iss >> col >> lo >> up;
      REQUIRE(highs.changeColBounds(col - 1, lo, up) == HighsStatus::kOk);
    } else if (op == "change_cols_bounds") {
      int n;
      iss >> n;
      for (int i = 0; i < n; i++) {
        HighsInt col;
        double lo, up;
        iss >> col >> lo >> up;
        REQUIRE(highs.changeColBounds(col - 1, lo, up) == HighsStatus::kOk);
      }
    } else if (op == "change_row_bounds") {
      HighsInt row;
      double lo, up;
      iss >> row >> lo >> up;
      REQUIRE(highs.changeRowBounds(row - 1, lo, up) == HighsStatus::kOk);
    } else if (op == "change_rows_bounds") {
      int n;
      iss >> n;
      for (int i = 0; i < n; i++) {
        HighsInt row;
        double lo, up;
        iss >> row >> lo >> up;
        REQUIRE(highs.changeRowBounds(row - 1, lo, up) == HighsStatus::kOk);
      }
    } else if (op == "change_cols_cost") {
      int n;
      iss >> n;
      for (int i = 0; i < n; i++) {
        HighsInt col;
        double val;
        iss >> col >> val;
        REQUIRE(highs.changeColCost(col - 1, val) == HighsStatus::kOk);
      }
    } else if (op == "solve") {
      REQUIRE(highs.run() == HighsStatus::kOk);
      solve_count++;
      final_obj = highs.getInfo().objective_function_value;
    } else if (op == "end") {
      break;
    }
  }

  REQUIRE(solve_count == 76);
  // Vérification de la valeur d'objectif attendue
  REQUIRE(std::fabs(final_obj - 298.2799078) < 1e-5);
}
