#include "Avgas.h"
#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"

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
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid a_matrix_.start_[0]\n");
  lp.a_matrix_.start_[0] = 0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run)
    printf("Give valid a_matrix_.start_[2] and a_matrix_.start_[3]\n");
  lp.a_matrix_.start_[2] = 2;
  lp.a_matrix_.start_[3] = 2;
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid a_matrix_.index_[0]\n");
  // Yields duplicate index, but values are still zero, so both are
  // discarded and a warning is returned
  lp.a_matrix_.index_[0] = 0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  if (dev_run)
    printf("Give nonzero a_matrix_.value_[0] and a_matrix_.value_[1]\n");
  // Yields duplicate index, but values are still zero, so both are
  // discarded and a warning is returned
  lp.a_matrix_.value_[0] = 1;
  lp.a_matrix_.value_[1] = 1;
  // Now the duplicate indices yield an erorr
  REQUIRE(highs.passModel(lp) == HighsStatus::kError);

  if (dev_run) printf("Give valid a_matrix_.index_[1]\n");
  lp.a_matrix_.index_[1] = 1;
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
    avgas.addRow(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
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
    avgas.addCol(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
                 Aindex, Avalue);
  }

  return_status = assessLp(lp, options);
  REQUIRE(return_status == HighsStatus::kOk);
  //  reportLp(lp, HighsLogType::kVerbose);

  const double my_infinity = 1e30;
  Highs highs;
  highs.passOptions(options);

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  return_status = highs.addRows(num_row, rowLower.data(), rowUpper.data(), 0,
                                NULL, NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status =
      highs.addCols(num_col, colCost.data(), colLower.data(), colUpper.data(),
                    num_col_nz, Astart.data(), Aindex.data(), Avalue.data());
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(!highs.getLp().has_infinite_cost_);

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
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(!highs.getLp().has_infinite_cost_);

  XcolUpper[0] = my_infinity;
  //  reportLp(lp, HighsLogType::kVerbose);

  // Try to add a column with infinite cost
  XcolCost[0] = my_infinity;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  XcolCost[0] = -my_infinity;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  // Reset to a legitimate cost
  XcolCost[0] = 1;

  // Add a column with bound inconsistency due to upper
  XcolLower[0] = 0;
  XcolUpper[0] = -1;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kWarning);

  // Add a column with bound inconsistency due to lower
  XcolLower[0] = 1;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kWarning);

  // Add a column with illegal bound due to lower
  XcolLower[0] = my_infinity;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kError);

  // Add a column with illegal bound due to upper
  XcolLower[0] = 0;
  XcolUpper[0] = -my_infinity;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
  REQUIRE(return_status == HighsStatus::kError);

  // Add a legitimate column
  XcolLower[0] = 0;
  XcolUpper[0] = 0;
  return_status =
      highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                    XcolUpper.data(), XnumNewNZ, XAstart.data(), NULL, NULL);
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
  return_status = highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                                XcolUpper.data(), XnumNewNZ, XAstart.data(),
                                XAindex.data(), XAvalue.data());
  REQUIRE(return_status == HighsStatus::kError);

  // Legitimise large matrix entries. Small entries now cause warning
  XAvalue[4] = -1;
  XAvalue[5] = 1;
  return_status = highs.addCols(XnumNewCol, XcolCost.data(), XcolLower.data(),
                                XcolUpper.data(), XnumNewNZ, XAstart.data(),
                                XAindex.data(), XAvalue.data());
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
  // infeasible. With infeasibleBoundsOk check in solveLp,
  // infeasiblility is identified before reaching a solver, so
  // presolve isn't called
  HighsStatus run_status;
  HighsModelStatus model_status;
  run_status = highs.run();
  // LP has free variable with infinite costs so cannot be solved
  REQUIRE(run_status == HighsStatus::kError);
  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kUnknown);

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

TEST_CASE("LP-row-index-duplication", "[highs_data]") {
  HighsStatus return_status;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsInt num_col = 10;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) highs.addVar(0, 1);
  std::vector<HighsInt> start = {0, 6, 8};
  std::vector<HighsInt> index = {0, 1, 4, 5, 6, 7, 0, 1, 0, 1, 4, 5, 6, 7, 4};
  std::vector<double> value = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> lower = {0, 0, 0};
  std::vector<double> upper = {inf, inf, inf};
  HighsInt num_nz = index.size();
  REQUIRE(highs.addRows(3, lower.data(), upper.data(), num_nz, start.data(),
                        index.data(), value.data()) == HighsStatus::kError);
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

TEST_CASE("LP-inf-cost", "[highs_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  const double my_infinite_bound = 1e30;
  const double my_infinite_cost = kHighsInf;
  const double my_large_cost = 1e20;
  highs.setOptionValue("infinite_cost", my_infinite_cost);
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {-3, -2, -my_infinite_cost};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {my_infinite_bound, my_infinite_bound, 1};
  lp.row_lower_ = {-my_infinite_bound, 12};
  lp.row_upper_ = {7.0, 12.0};
  lp.a_matrix_.start_ = {0, 2, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {1.0, 4.0, 1.0, 2.0, 1, 1.0};

  lp.integrality_.resize(lp.num_col_);
  lp.integrality_[0] = HighsVarType::kContinuous;
  lp.integrality_[0] = HighsVarType::kContinuous;
  lp.integrality_[0] = HighsVarType::kInteger;

  HighsStatus status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  // Now just make the cost large
  status = highs.changeColCost(2, -my_large_cost);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(!highs.getLp().has_infinite_cost_);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  // Now transform to LP and make cost infinite
  status = highs.changeColCost(2, -my_infinite_cost);
  status = highs.changeColIntegrality(2, HighsVarType::kContinuous);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  highs.clearModel();

  // Formulate min -inf x s.t. x <= 2; 1 <= x <= 3
  //
  // Fixing x at 3 is infeasible, so kUnknown status should be
  // returned with kWarning
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {-my_infinite_cost};
  lp.col_lower_ = {1};
  lp.col_upper_ = {3};
  lp.row_lower_ = {-my_infinite_bound};
  lp.row_upper_ = {2};
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1};

  lp.integrality_.clear();

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);

  highs.changeObjectiveSense(ObjSense::kMaximize);
  // Switching to maximization, fixing x at 1 is feasible, so kOk
  // status should be returned with kOptimal and objective = -inf

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getInfo().objective_function_value == -my_infinite_cost);

  // Formulate min -inf x s.t. x <= 2; 0.5 <= x <= 3.5; x integer
  //
  // Fixing x at 3 is infeasible, so kUnknown status should be
  // returned with kWarning

  highs.clearModel();

  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {-my_infinite_cost};
  lp.col_lower_ = {0.5};
  lp.col_upper_ = {3.5};
  lp.row_lower_ = {-my_infinite_bound};
  lp.row_upper_ = {2};
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1};
  lp.integrality_ = {HighsVarType::kInteger};
  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getLp().has_infinite_cost_);

  status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);

  highs.changeObjectiveSense(ObjSense::kMaximize);
  // Switching to maximization, fixing x at 1 is feasible, so kOk
  // status should be returned with kOptimal and objective = -inf

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getInfo().objective_function_value == -my_infinite_cost);
  // Check that x was fixed at 1, not 0.5
  REQUIRE(highs.getSolution().col_value[0] == 1);
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

TEST_CASE("LP-illegal-empty-start-ok", "[highs_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 0;
  lp.num_row_ = 1;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {1};
  lp.a_matrix_.start_ = {1};
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  REQUIRE(highs.getLp().a_matrix_.start_[0] == 0);
}

TEST_CASE("LP-row-wise", "[highs_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.sense_ = ObjSense::kMaximize;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {10, 25};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 2, 1, 4};
  lp.row_lower_ = {-inf, -inf};
  lp.row_upper_ = {80, 120};
  highs.passModel(lp);
  highs.run();
}

TEST_CASE("LP-infeasible-bounds", "[highs_data]") {
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  double epsilon = 1e-10;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.sense_ = ObjSense::kMaximize;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {10, 25};
  lp.col_lower_ = {1, 2.5 + epsilon};
  lp.col_upper_ = {1 - epsilon, 2.5 - 2 * epsilon};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 2, 1, 4};
  lp.row_lower_ = {6, -inf};
  lp.row_upper_ = {6 - epsilon, 11 - epsilon};
  highs.passModel(lp);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  if (dev_run) highs.writeSolution("", 1);

  highs.changeColBounds(0, 0, -1);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}
