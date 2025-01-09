#include "Avgas.h"
#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("LP-orientation", "[lp_orientation]") {
  Avgas avgas;
  HighsInt num_row = 0;
  HighsInt num_row_nz = 0;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<HighsInt> ARstart;
  vector<HighsInt> ARindex;
  vector<double> ARvalue;

  for (HighsInt row = 0; row < avgas_num_row; row++)
    avgas.addRow(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);

  ARstart.push_back(num_row_nz);

  HighsInt num_col = 0;
  HighsInt num_col_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<HighsInt> Astart;
  vector<HighsInt> Aindex;
  vector<double> Avalue;
  for (HighsInt col = 0; col < avgas_num_col; col++)
    avgas.addCol(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
              Aindex, Avalue);
  Astart.push_back(num_col_nz);

  assert(num_col_nz == num_row_nz);

  double optimal_objective_function_value = -7.75;
  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  } else {
    highs.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  }
  const HighsLp& highs_lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();

  REQUIRE(highs_lp.a_matrix_.format_ == MatrixFormat::kColwise);

  // Set up the LP externally
  HighsLp lp;
  lp.num_col_ = num_col;
  lp.num_row_ = num_row;
  lp.col_cost_ = colCost;
  lp.col_lower_ = colLower;
  lp.col_upper_ = colUpper;
  lp.row_lower_ = rowLower;
  lp.row_upper_ = rowUpper;
  lp.a_matrix_.start_ = Astart;
  lp.a_matrix_.index_ = Aindex;
  lp.a_matrix_.value_ = Avalue;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  // Have to set the matrix dimension so that setFormat can be used
  lp.a_matrix_.num_col_ = num_col;
  lp.a_matrix_.num_row_ = num_row;
  // Pass the LP
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Make the external LP row-wise then pass and solve it
  lp.ensureRowwise();
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Make the external LP col-wise then pass and solve it
  lp.ensureColwise();
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Clear the internal LP
  highs.clearModel();
  REQUIRE(highs.addCols(num_col, colCost.data(), colLower.data(),
                        colUpper.data(), 0, NULL, NULL,
                        NULL) == HighsStatus::kOk);
  REQUIRE(highs.addRows(num_row, rowLower.data(), rowUpper.data(), num_row_nz,
                        ARstart.data(), ARindex.data(),
                        ARvalue.data()) == HighsStatus::kOk);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Clear the internal LP
  highs.clearModel();
  highs.addCols(num_col, colCost.data(), colLower.data(), colUpper.data(), 0,
                NULL, NULL, NULL);
  double one_row_lower;
  double one_row_upper;
  vector<HighsInt> one_row_index;
  vector<double> one_row_value;
  HighsInt one_row_num_nz;
  for (HighsInt row = 0; row < avgas_num_row; row++) {
    avgas.getRow(row, one_row_lower, one_row_upper, one_row_index, one_row_value);
    one_row_num_nz = one_row_index.size();
    REQUIRE(highs.addRow(one_row_lower, one_row_upper, one_row_num_nz, 
			 one_row_index.data(),
			 one_row_value.data()) == HighsStatus::kOk);
  }
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  //  return_status = highs.writeModel("");
}
