// HiGHS is designed to solve linear optimization problems of the form
//
// Min (1/2)x^TQx + c^Tx + d subject to L <= Ax <= U; l <= x <= u
//
// where A is a matrix with m rows and n columns, and Q is either zero
// or positive definite. If Q is zero, HiGHS can determine the optimal
// integer-valued solution.
//
// The scalar n is num_col_
// The scalar m is num_row_
//
// The vector c is col_cost_
// The scalar d is offset_
// The vector l is col_lower_
// The vector u is col_upper_
// The vector L is row_lower_
// The vector U is row_upper_
//
// The matrix A is represented in packed vector form, either
// row-wise or column-wise: only its nonzeros are stored
//
// * The number of nonzeros in A is num_nz
//
// * The indices of the nonnzeros in the vectors of A are stored in a_index
//
// * The values of the nonnzeros in the vectors of A are stored in a_value
//
// * The position in a_index/a_value of the index/value of the first
// nonzero in each vector is stored in a_start
//
// Note that a_start[0] must be zero
//
// The matrix Q is represented in packed column form
//
// * The dimension of Q is dim_
//
// * The number of nonzeros in Q is hessian_num_nz
//
// * The indices of the nonnzeros in the vectors of A are stored in q_index
//
// * The values of the nonnzeros in the vectors of A are stored in q_value
//
// * The position in q_index/q_value of the index/value of the first
// nonzero in each column is stored in q_start
//
// Note
//
// * By default, Q is zero. This is indicated by dim_ being initialised to zero.
//
// * q_start[0] must be zero
//
#include "Highs.h"
#include <cassert>

using std::cout;
using std::endl;

int main() {
  // Create and populate a HighsModel instance for the LP
  //
  // Min    f  =  x_0 +  x_1 + 3
  // s.t.                x_1 <= 7
  //        5 <=  x_0 + 2x_1 <= 15
  //        6 <= 3x_0 + 2x_1
  // 0 <= x_0 <= 4; 1 <= x_1
  //
  // Although the first constraint could be expressed as an upper
  // bound on x_1, it serves to illustrate a non-trivial packed
  // column-wise matrix.
  //
  HighsModel model;
  model.lp_.num_col_ = 2;
  model.lp_.num_row_ = 3;
  model.lp_.sense_ = ObjSense::kMinimize;
  model.lp_.offset_ = 3;
  model.lp_.col_cost_ = {1.0, 1.0};
  model.lp_.col_lower_ = {0.0, 1.0};
  model.lp_.col_upper_ = {4.0, 1.0e30};
  model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
  model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};
  //
  // Here the orientation of the matrix is column-wise
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  // a_start_ has num_col_1 entries, and the last entry is the number
  // of nonzeros in A, allowing the number of nonzeros in the last
  // column to be defined
  model.lp_.a_matrix_.start_ = {0, 2, 5};
  model.lp_.a_matrix_.index_ = {1, 2, 0, 1, 2};
  model.lp_.a_matrix_.value_ = {1.0, 3.0, 1.0, 2.0, 2.0};
  //
  // Create a Highs instance
  Highs highs;
  HighsStatus return_status;
  //
  // Pass the model to HiGHS
  return_status = highs.passModel(model);
  assert(return_status==HighsStatus::kOk);
  //
  // Get a const reference to the LP data in HiGHS
  const HighsLp& lp = highs.getLp();
  //
  // Solve the model
  return_status = highs.run();
  assert(return_status==HighsStatus::kOk);
  //
  // Get the model status
  const HighsModelStatus& model_status = highs.getModelStatus();
  assert(model_status==HighsModelStatus::kOptimal);
  cout << "Model status: " << highs.modelStatusToString(model_status) << endl;
  //
  // Get the solution information
  const HighsInfo& info = highs.getInfo();
  cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
  cout << "Objective function value: " << info.objective_function_value << endl;
  cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << endl;
  cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << endl;
  cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;
  const bool has_values = info.primal_solution_status;
  const bool has_duals = info.dual_solution_status;
  const bool has_basis = info.basis_validity;
  //
  // Get the solution values and basis
  const HighsSolution& solution = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  //
  // Report the primal and solution values and basis
  for (int col=0; col < lp.num_col_; col++) {
    cout << "Column " << col;
    if (has_values) cout << "; value = " << solution.col_value[col];
    if (has_duals) cout << "; dual = " << solution.col_dual[col];
    if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.col_status[col]);
    cout << endl;
  }
  for (int row=0; row < lp.num_row_; row++) {
    cout << "Row    " << row;
    if (has_values) cout << "; value = " << solution.row_value[row];
    if (has_duals) cout << "; dual = " << solution.row_dual[row];
    if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.row_status[row]);
    cout << endl;
  }

  // Now indicate that all the variables must take integer values
  model.lp_.integrality_.resize(lp.num_col_);
  for (int col=0; col < lp.num_col_; col++)
    model.lp_.integrality_[col] = HighsVarType::kInteger;

  highs.passModel(model);
  // Solve the model
  return_status = highs.run();
  assert(return_status==HighsStatus::kOk);
  // Report the primal solution values
  for (int col=0; col < lp.num_col_; col++) {
    cout << "Column " << col;
    if (info.primal_solution_status) cout << "; value = " << solution.col_value[col];
    cout << endl;
  }
  for (int row=0; row < lp.num_row_; row++) {
    cout << "Row    " << row;
    if (info.primal_solution_status) cout << "; value = " << solution.row_value[row];
    cout << endl;
  }

  return 0;
}
