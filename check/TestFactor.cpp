#include "Highs.h"
#include "catch.hpp"
#include "util/HFactor.h"

const bool dev_run = true;

HVector rhs;
HVector col_aq;
HVector row_ep;
std::vector<HighsInt> basic_set;
std::vector<double> solution;
HighsLp lp;
HighsInt num_col;
HighsInt num_row;
HFactor factor;

bool iterate(const HighsInt row_out, const HighsInt variable_out, const HighsInt variable_in);
bool testSolve();

TEST_CASE("Factor-get-set-invert", "[highs_test_factor]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
   Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  lp = highs.getLp();
  num_col = lp.num_col_;
  num_row = lp.num_row_;
  std::vector<HighsInt> row_out = {8, 1, 7, 4, 0, 6};
  std::vector<HighsInt> variable_out = {16, 9, 15, 12, 8, 14};
  std::vector<HighsInt> variable_in = {5, 2, 0, 4, 3, 6};
  basic_set = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  HighsRandom random;
  solution.resize(num_row);
  for (HighsInt iRow=0; iRow<num_row; iRow++) 
    solution[iRow] = random.fraction();
  rhs.setup(num_row);
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  factor.setup(lp.a_matrix_, basic_set);
  factor.build();
  HighsInt from_basis_change = 0;
  HighsInt to_basis_change = 3;
  for (HighsInt basis_change = from_basis_change; basis_change < to_basis_change; basis_change++)
    REQUIRE(iterate(row_out[basis_change], variable_out[basis_change], variable_in[basis_change]));
  std::vector<HighsInt> get_basic_set = basic_set;
  InvertibleRepresentation invert = factor.getInvert();
  std::vector<InvertibleRepresentation> invert_set;
  from_basis_change = to_basis_change;
  to_basis_change = row_out.size();
  for (HighsInt basis_change = from_basis_change; basis_change < to_basis_change; basis_change++) {
    REQUIRE(iterate(row_out[basis_change], variable_out[basis_change], variable_in[basis_change]));
    invert_set.push_back(factor.getInvert());
  }
  basic_set = get_basic_set;
  factor.setInvert(invert);
  REQUIRE(testSolve());

  for (HighsInt basis_change = from_basis_change; basis_change < to_basis_change; basis_change++)
    REQUIRE(iterate(row_out[basis_change], variable_out[basis_change], variable_in[basis_change]));
 
}

bool iterate(const HighsInt row_out, const HighsInt variable_out, const HighsInt variable_in) {
  assert(basic_set[row_out] == variable_out);
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = row_out;
  row_ep.array[row_out] = 1;
  row_ep.packFlag = true;
  factor.btranCall(row_ep, 1);
    
    
  col_aq.clear();
  col_aq.packFlag = true;
  lp.a_matrix_.collectAj(col_aq, variable_in, 1);
  factor.ftranCall(col_aq, 1);
    
  basic_set[row_out] = variable_in;
  HighsInt rebuild_reason = 0;
  HighsInt lc_row_out = row_out;
  factor.update(&col_aq, &row_ep, &lc_row_out, &rebuild_reason);
  if (rebuild_reason) return false;
    
  return testSolve();
}

bool testSolve() {
  // FTRAN
  rhs.clear();
  for (HighsInt iCol=0; iCol<num_row; iCol++)
    lp.a_matrix_.collectAj(rhs, basic_set[iCol], solution[iCol]);
  factor.ftranCall(rhs, 1);
  double error_norm = 0;
  for (HighsInt iRow=0; iRow<num_row; iRow++) 
    error_norm = std::max(std::fabs(solution[iRow] - rhs.array[iRow]), error_norm);
  if (error_norm > 1e-4) return false;
  // BTRAN
  rhs.clear();
  for (HighsInt iCol=0; iCol<num_row; iCol++) {
    rhs.array[iCol] = lp.a_matrix_.computeDot(solution, basic_set[iCol]);
    if (rhs.array[iCol]) rhs.index[rhs.count++] = iCol;
  }
  factor.btranCall(rhs, 1);
  error_norm = 0;
  for (HighsInt iRow=0; iRow<num_row; iRow++) 
    error_norm = std::max(std::fabs(solution[iRow] - rhs.array[iRow]), error_norm);
  return error_norm < 1e-4;
}
