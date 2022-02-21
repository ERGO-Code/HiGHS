#include "Avgas.h"
#include "Highs.h"
#include "HighsLpUtils.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;

// No commas in test case name.
TEST_CASE("LP-dimension-validation", "[highs_data]") {
  // Create an LP with lots of illegal values
  const HighsInt true_num_col = 2;
  const HighsInt true_num_row = 3;
  HighsLp lp;
  lp.num_col_ = -1;
  lp.num_row_ = -1;
  lp.col_cost_.resize(1);
  lp.col_lower_.resize(1);
  lp.col_upper_.resize(1);
  lp.a_matrix_.format_ = MatrixFormat::kRowwisePartitioned;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_.resize(2);
  lp.a_matrix_.start_[0] = -1;
  lp.a_matrix_.start_[1] = 2;
  lp.a_matrix_.index_.resize(2);
  lp.a_matrix_.index_[0] = -1;
  lp.a_matrix_.index_[1] = 0;
  lp.a_matrix_.value_.resize(2);
  lp.row_lower_.resize(1);
  lp.row_upper_.resize(1);

  // Set up invalid scale data once scaling can be imported
  //  lp.scale_.strategy = -1;
  //  lp.scale_.num_col = 1;
  //  lp.scale_.num_row = 1;
  //  lp.scale_.has_scaling = false;
  //  lp.scale_.row.resize(1);
  //  lp.scale_.col.resize(1);
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  HighsStatus return_status;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid number of columns\n");
  lp.num_col_ = true_num_col;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid number of rows\n");
  lp.num_row_ = true_num_row;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid col_cost.size()\n");
  lp.col_cost_.resize(true_num_col);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid col_lower.size()\n");
  lp.col_lower_.resize(true_num_col);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid col_upper.size()\n");
  lp.col_upper_.resize(true_num_col);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid a_matrix_.format_\n");
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid a_matrix_.start.size()\n");
  lp.a_matrix_.start_.resize(true_num_row + 1);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid row_lower.size()\n");
  lp.row_lower_.resize(true_num_row);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid row_upper.size()\n");
  lp.row_upper_.resize(true_num_row);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  /*
  if (dev_run) printf("Give valid scale_.strategy\n");
  lp.scale_.strategy = kSimplexScaleStrategyOff;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.num_col\n");
  lp.scale_.num_col = 0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.num_row\n");
  lp.scale_.num_row = 0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.col.size()\n");
  lp.scale_.col.resize(0);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.row.size()\n");
  lp.scale_.row.resize(0);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  if (dev_run) printf("Set scale_.strategy =
  kSimplexScaleStrategyMaxValue015\n"); lp.scale_.strategy =
  kSimplexScaleStrategyMaxValue015; REQUIRE(highs.passModel(lp) ==
  HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.num_col\n");
  lp.scale_.num_col = true_num_col;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.num_row\n");
  lp.scale_.num_row = true_num_row;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.col.size()\n");
  lp.scale_.col.resize(true_num_col);
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid scale_.row.size()\n");
  lp.scale_.row.resize(true_num_row);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  */
}

TEST_CASE("LP-validation", "[highs_data]") {
  // Create an empty LP
  HighsLp lp;
  HighsOptions options;
  HighsStatus return_status;
  options.log_dev_level = kHighsLogDevLevelVerbose;
  if (!dev_run) options.output_flag = false;

  Avgas avgas;
  const HighsInt avgas_num_col = 8;
  const HighsInt avgas_num_row = 10;
  HighsInt num_row = 0;
  HighsInt num_row_nz = 0;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<HighsInt> ARstart;
  vector<HighsInt> ARindex;
  vector<double> ARvalue;

  for (HighsInt row = 0; row < avgas_num_row; row++) {
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);
  }

  HighsInt num_col = 0;
  HighsInt num_col_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<HighsInt> Astart;
  vector<HighsInt> Aindex;
  vector<double> Avalue;
  for (HighsInt col = 0; col < avgas_num_col; col++) {
    avgas.col(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
              Aindex, Avalue);
  }

  return_status = assessLp(lp, options);
  REQUIRE(return_status == HighsStatus::kOk);
  //  reportLp(lp, HighsLogType::kVerbose);

  const double my_infinity = 1e30;
  Highs highs;
  highs.passOptions(options);

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  return_status =
      highs.addRows(num_row, &rowLower[0], &rowUpper[0], 0, NULL, NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status =
      highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                    num_col_nz, &Astart[0], &Aindex[0], &Avalue[0]);
  REQUIRE(return_status == HighsStatus::kOk);

  // Create an empty column
  HighsInt XnumNewCol = 1;
  HighsInt XnumNewNZ = 0;
  vector<double> XcolCost;
  XcolCost.resize(XnumNewCol);
  XcolCost[0] = 1;
  vector<double> XcolLower;
  XcolLower.resize(XnumNewCol);
  XcolLower[0] = 0;
  vector<double> XcolUpper;
  XcolUpper.resize(XnumNewCol);
  XcolUpper[0] = 1e25;
  vector<HighsInt> XAstart;
  XAstart.resize(XnumNewCol);
  vector<HighsInt> XAindex;
  vector<double> XAvalue;
  // Add an empty column
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);
  XcolUpper[0] = my_infinity;
  //  reportLp(lp, HighsLogType::kVerbose);

  // Try to add a column with illegal cost
  HighsStatus require_return_status;
  if (kHighsAllowInfiniteCosts) {
    require_return_status = HighsStatus::kOk;
  } else {
    require_return_status = HighsStatus::kError;
  }
  XcolCost[0] = my_infinity;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == require_return_status);

  XcolCost[0] = -my_infinity;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == require_return_status);

  // Reset to a legitimate cost
  XcolCost[0] = 1;

  // Add a column with bound inconsistency due to upper
  XcolLower[0] = 0;
  XcolUpper[0] = -1;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kWarning);

  // Add a column with bound inconsistency due to lower
  XcolLower[0] = 1;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kWarning);

  // Add a column with illegal bound due to lower
  XcolLower[0] = my_infinity;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kError);

  // Add a column with illegal bound due to upper
  XcolLower[0] = 0;
  XcolUpper[0] = -my_infinity;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kError);

  // Add a legitimate column
  XcolLower[0] = 0;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);

  //  reportLp(lp, HighsLogType::kVerbose);

  // Add a couple of non-empty columns with some small and large values
  XnumNewCol = 2;
  XnumNewNZ = 7;
  XcolCost.resize(XnumNewCol);
  XcolCost[0] = 1;
  XcolCost[1] = 2;
  XcolLower.resize(XnumNewCol);
  XcolLower[0] = 0;
  XcolLower[1] = 0;
  XcolUpper.resize(XnumNewCol);
  XcolUpper[0] = 1;
  XcolUpper[1] = 1;
  XAstart.resize(XnumNewCol + 1);
  XAindex.resize(XnumNewNZ);
  XAstart[1] = 4;
  XAstart[2] = XnumNewNZ;
  XAindex[0] = 0;
  XAindex[1] = 2;
  XAindex[2] = 3;
  XAindex[3] = 9;
  XAindex[4] = 1;
  XAindex[5] = 3;
  XAindex[6] = 8;
  XAvalue.resize(XnumNewNZ);
  XAvalue[0] = 1;
  XAvalue[1] = 1e-12;
  XAvalue[2] = -1e-20;
  XAvalue[3] = -1;
  XAvalue[4] = -1e60;
  XAvalue[5] = 1e100;
  XAvalue[6] = -1;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], &XAindex[0], &XAvalue[0]);
  REQUIRE(return_status == HighsStatus::kError);

  // Legitimise large matrix entries. Small entries now cause warning
  XAvalue[4] = -1;
  XAvalue[5] = 1;
  return_status =
      highs.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                    XnumNewNZ, &XAstart[0], &XAindex[0], &XAvalue[0]);
  REQUIRE(return_status == HighsStatus::kWarning);

  if (!dev_run) highs.setOptionValue("output_flag", false);

  const HighsLp& internal_lp = highs.getLp();
  double check_value;
  REQUIRE(highs.getCoeff(-1, 0, check_value) == HighsStatus::kError);
  REQUIRE(highs.getCoeff(0, -1, check_value) == HighsStatus::kError);
  REQUIRE(highs.getCoeff(internal_lp.num_row_, 0, check_value) ==
          HighsStatus::kError);
  REQUIRE(highs.getCoeff(0, internal_lp.num_col_, check_value) ==
          HighsStatus::kError);

  const HighsInt check_col = 4;
  const HighsInt check_row = 7;
  REQUIRE(highs.getCoeff(check_col, check_row, check_value) ==
          HighsStatus::kOk);
  REQUIRE(check_value == 0);

  const double value = -3;
  REQUIRE(highs.getCoeff(check_row, check_col, check_value) ==
          HighsStatus::kOk);
  REQUIRE(check_value == value);

  // This is a highly anomalous LP. It has two pairs of inconsistent
  // bounds (cols 11 and 12) but also has costs of 1e+30 and -1e+30
  // for columns 9 and 10.

  // LP is found to be unbounded by presolve, but is primal
  // infeasible. With isBoundInfeasible check in solveLp,
  // infeasiblility is identified before reaching a solver, so
  // presolve isn't called
  HighsStatus run_status;
  HighsModelStatus model_status;
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kInfeasible);

  REQUIRE(highs.changeCoeff(-1, 0, check_value) == HighsStatus::kError);
  REQUIRE(highs.changeCoeff(0, -1, check_value) == HighsStatus::kError);
  REQUIRE(highs.changeCoeff(internal_lp.num_row_, 0, check_value) ==
          HighsStatus::kError);
  REQUIRE(highs.changeCoeff(0, internal_lp.num_col_, check_value) ==
          HighsStatus::kError);

  const double to_value = 99;
  REQUIRE(highs.changeCoeff(check_row, check_col, to_value) ==
          HighsStatus::kOk);
  REQUIRE(highs.getCoeff(check_row, check_col, check_value) ==
          HighsStatus::kOk);
  REQUIRE(check_value == to_value);
}

TEST_CASE("LP-extreme-coefficient", "[highs_data]") {
  HighsStatus return_status;
  std::string filename;
  Highs highs;
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {1};
  lp.col_lower_ = {0};
  lp.col_upper_ = {1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1e300};
  lp.row_lower_ = {1};
  lp.row_upper_ = {inf};
  if (!dev_run) highs.setOptionValue("output_flag", false);
  // Reading a model with a large value results in HighsStatus::kError
  return_status = highs.passModel(lp);
  if (dev_run)
    printf("highs.passModel(filename); returns %d\n", (int)return_status);
  REQUIRE(return_status == HighsStatus::kError);
  // Solving a model with a large value results in HighsStatus::kError
  return_status = highs.run();
  if (dev_run) printf("highs.run(); returns %d\n", (int)return_status);
  REQUIRE(return_status == HighsStatus::kError);
  lp.a_matrix_.value_[0] = 1e-300;
  return_status = highs.passModel(lp);
  if (dev_run)
    printf("highs.passModel(filename); returns %d\n", (int)return_status);
  REQUIRE(return_status == HighsStatus::kWarning);
  // Solving a model with a small value results in HighsStatus::kOk
  return_status = highs.run();
  if (dev_run) printf("highs.run(); returns %d\n", (int)return_status);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}

TEST_CASE("LP-change-coefficient", "[highs_data]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsInt change_coefficient_col = 4;
  const HighsInt zero_coefficient_row = 4;
  const HighsInt add_coefficient_row = 5;
  const double zero_coefficient_value = 1e-10;
  const double add_coefficient_value = -3;
  double required_objective_value;
  double delta_objective_value;
  double original_value;
  REQUIRE(highs.getCoeff(zero_coefficient_row, change_coefficient_col,
                         original_value) == HighsStatus::kOk);
  REQUIRE(original_value == -1);
  // Change coefficient (4, 4) from -1 to a small value - zeroing it
  REQUIRE(highs.changeCoeff(zero_coefficient_row, change_coefficient_col,
                            zero_coefficient_value) == HighsStatus::kOk);
  required_objective_value = -8.6666666667;
  highs.run();
  delta_objective_value = std::fabs(required_objective_value -
                                    highs.getInfo().objective_function_value);
  REQUIRE(delta_objective_value < 1e-8);

  // Restore coefficient (4, 4)
  REQUIRE(highs.changeCoeff(zero_coefficient_row, change_coefficient_col,
                            original_value) == HighsStatus::kOk);
  required_objective_value = -7.75;
  highs.run();
  delta_objective_value = std::fabs(required_objective_value -
                                    highs.getInfo().objective_function_value);
  REQUIRE(delta_objective_value < 1e-8);

  // Change coefficient (5, 4) to -3
  REQUIRE(highs.changeCoeff(add_coefficient_row, change_coefficient_col,
                            add_coefficient_value) == HighsStatus::kOk);
  required_objective_value = -7.5;
  highs.run();
  delta_objective_value = std::fabs(required_objective_value -
                                    highs.getInfo().objective_function_value);
  REQUIRE(delta_objective_value < 1e-8);
}
