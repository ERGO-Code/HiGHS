#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// gcc call_highs_from_c.c -o highstest -I install_folder/include/ -L install_folder/lib/ -lhighs

void minimal_api() {
  // This illustrates the use of Highs_lpCall, the simple C interface to
  // HiGHS. It's designed to solve the general LP problem
  //
  // Min c^Tx + d subject to L <= Ax <= U; l <= x <= u
  //
  // where A is a matrix with m rows and n columns
  //
  // The scalar n is num_col
  // The scalar m is num_row
  //
  // The vector c is col_cost
  // The scalar d is offset
  // The vector l is col_lower
  // The vector u is col_upper
  // The vector L is row_lower
  // The vector U is row_upper
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
  // After a successful call to Highs_lpCall, the primal and dual
  // solution, and the simplex basis are returned as follows
  //
  // The vector x is col_value
  // The vector Ax is row_value
  // The vector of dual values for the variables x is col_dual
  // The vector of dual values for the variables Ax is row_dual
  // The basic/nonbasic status of the variables x is col_basis_status
  // The basic/nonbasic status of the variables Ax is row_basis_status
  //
  // The status of the solution obtained is model_status
  //
  // The use of Highs_lpCall is illustrated for the LP
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
  const int num_col = 2;
  const int num_row = 3;
  const int num_nz = 5;
  // Define the optimization sense and objective offset
  int sense = kHighsObjSenseMinimize;
  const double offset = 3;

  // Define the column costs, lower bounds and upper bounds
  const double col_cost[2] = {1.0, 1.0};
  const double col_lower[2] = {0.0, 1.0};
  const double col_upper[2] = {4.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  const double row_lower[3] = {-1.0e30, 5.0, 6.0};
  const double row_upper[3] = {7.0, 15.0, 1.0e30};
  // Define the constraint matrix column-wise
  const int a_format = 1;
  const int a_start[2] = {0, 2};
  const int a_index[5] = {1, 2, 0, 1, 2};
  const double a_value[5] = {1.0, 3.0, 1.0, 2.0, 2.0};

  double objective_value;
  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* col_dual = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);
  double* row_dual = (double*)malloc(sizeof(double) * num_row);

  int* col_basis_status = (int*)malloc(sizeof(int) * num_col);
  int* row_basis_status = (int*)malloc(sizeof(int) * num_row);

  int model_status;
  int run_status;

  run_status = Highs_lpCall(num_col, num_row, num_nz, a_format,
			   sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
			   a_start, a_index, a_value,
			   col_value, col_dual, row_value, row_dual,
			   col_basis_status, row_basis_status,
			   &model_status);
  // The run must be successful, and the model status optimal
  assert(run_status == kHighsStatusOk);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  objective_value = offset;
  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf; dual = %lf; status = %d\n", i, col_value[i], col_dual[i], col_basis_status[i]);
    objective_value += col_value[i]*col_cost[i];
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf; dual = %lf; status = %d\n", i, row_value[i], row_dual[i], row_basis_status[i]);
  }
  printf("Optimal objective value = %g\n", objective_value);

  // Switch the sense to maximization and solve the LP again
  sense = kHighsObjSenseMaximize;
  run_status = Highs_lpCall(num_col, num_row, num_nz, a_format,
			   sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
			   a_start, a_index, a_value,
			   col_value, col_dual, row_value, row_dual,
			   col_basis_status, row_basis_status,
			   &model_status);
  // The run must be successful, and the model status optimal
  assert(run_status == kHighsStatusOk);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  // Compute the objective value
  objective_value = offset;
  for (int i = 0; i < num_col; i++) objective_value += col_value[i]*col_cost[i];
  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf; dual = %lf; status = %d\n", i, col_value[i], col_dual[i], col_basis_status[i]);
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf; dual = %lf; status = %d\n", i, row_value[i], row_dual[i], row_basis_status[i]);
  }
  printf("Optimal objective value = %g\n", objective_value);
  // 
  // Indicate that the optimal solution for both columns must be
  // integer valued and solve the model as a MIP
  int integrality[2] = {1, 1};
  run_status = Highs_mipCall(num_col, num_row, num_nz, a_format,
			     sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
			     a_start, a_index, a_value,
			     integrality,
			     col_value, row_value, 
			     &model_status);
  // The run must be successful, and the model status optimal
  assert(run_status == kHighsStatusOk);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  // Compute the objective value
  objective_value = offset;
  for (int i = 0; i < num_col; i++) objective_value += col_value[i]*col_cost[i];
  // Report the column primal values
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf\n", i, col_value[i]);
  }
  // Report the row primal values
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf\n", i, row_value[i]);
  }
  printf("Optimal objective value = %g\n", objective_value);

  free(col_value);
  free(col_dual);
  free(row_value);
  free(row_dual);
  free(col_basis_status);
  free(row_basis_status);
}

void minimal_api_qp() {
  // Illustrate the solution of a QP
  //
  // minimize -x_2 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
  //
  // subject to x_1 + x_2 + x_3 >= 1; x>=0
  
  const int num_col = 3;
  const int num_row = 1;
  const int num_nz = 3;
  // Define the optimization sense and objective offset
  int sense = kHighsObjSenseMinimize;
  const double offset = 0;

  // Define the column costs, lower bounds and upper bounds
  const double col_cost[3] = {0.0, -1.0, 0.0};
  const double col_lower[3] = {0.0, 0.0, 0.0};
  const double col_upper[3] = {1.0e30, 1.0e30, 1.0e30};
  // Define the row lower bounds and upper bounds
  const double row_lower[1] = {1};
  const double row_upper[1] = {1.0e30};
  // Define the constraint matrix row-wise
  const int a_format = kHighsMatrixFormatRowwise;
  const int a_start[2] = {0, 3};
  const int a_index[3] = {0, 1, 2};
  const double a_value[3] = {1.0, 1.0, 1.0};

  const int q_format = kHighsHessianFormatTriangular;
  const int q_num_nz = 5;
  const int q_start[3] = {0, 2, 3};
  const int q_index[5] = {0, 2, 1, 0, 2};
  const double q_value[5] = {2.0, -1.0, 0.2, -1.0, 2.0};  

  double objective_value;
  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* col_dual = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);
  double* row_dual = (double*)malloc(sizeof(double) * num_row);

  int* col_basis_status = (int*)malloc(sizeof(int) * num_col);
  int* row_basis_status = (int*)malloc(sizeof(int) * num_row);

  int model_status;
  int run_status;

  run_status = Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, 
			   sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
			   a_start, a_index, a_value,
			   q_start, q_index, q_value,
			   col_value, col_dual, row_value, row_dual,
			   col_basis_status, row_basis_status,
			   &model_status);
  // The run must be successful, and the model status optimal
  assert(run_status == kHighsStatusOk);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  // Compute the objective value
  objective_value = offset;
  for (int i = 0; i < num_col; i++) objective_value += col_value[i]*col_cost[i];
  for (int i = 0; i < num_col; i++) {
    int from_el = q_start[i];
    int to_el;
    if (i+1<num_col) {
      to_el = q_start[i+1];
    } else {
      to_el = q_num_nz;
    }
    for (int el = from_el; el < to_el; el++) {
      int j = q_index[el];
      objective_value += 0.5*col_value[i]*col_value[j]*q_value[el];
    }
  }

  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    //    printf("Col%d = %lf; dual = %lf; status = %d\n", i, col_value[i], col_dual[i], col_basis_status[i]);
    printf("Col%d = %lf; dual = %lf\n", i, col_value[i], col_dual[i]);
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    //    printf("Row%d = %lf; dual = %lf; status = %d\n", i, row_value[i], row_dual[i], row_basis_status[i]);
    printf("Row%d = %lf; dual = %lf\n", i, row_value[i], row_dual[i]);
  }
  printf("Optimal objective value = %g\n", objective_value);

  free(col_value);
  free(col_dual);
  free(row_value);
  free(row_dual);
  free(col_basis_status);
  free(row_basis_status);
}

void minimal_api_mps() {
  // Illustrate the minimal interface for reading an mps file. Assumes
  // that the model file is check/instances/avgas.mps
  
  const char* filename = "../HiGHS/check/instances/avgas.mps";
  // Create a Highs instance
  void* highs = Highs_create();
  int run_status;
  run_status = Highs_readModel(highs, filename);
  assert(run_status == kHighsStatusOk);
  run_status = Highs_run(highs);
  int model_status = Highs_getModelStatus(highs);
  // The run must be successful, and the model status optimal
  assert(run_status == kHighsStatusOk);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  double objective_function_value;
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  printf("Optimal objective value = %g\n", objective_function_value);
  assert(abs(objective_function_value+7.75)<1e-5);
}

void full_api() {
  // This example does exactly the same as the minimal example above,
  // but illustrates the full C API.  It first forms and solves the LP
  //
  // Min    f  =  x_0 +  x_1 + 3
  // s.t.                x_1 <= 7
  //        5 <=  x_0 + 2x_1 <= 15
  //        6 <= 3x_0 + 2x_1
  // 0 <= x_0 <= 4; 1 <= x_1
  //
  // It then solves it as a maximization, then as a MIP.
  //
  const int num_col = 2;
  const int num_row = 3;
  const int num_nz = 5;
  // Define the optimization sense and objective offset
  int sense = kHighsObjSenseMinimize;
  const double offset = 3;

  // Define the column costs, lower bounds and upper bounds
  const double col_cost[2] = {1.0, 1.0};
  const double col_lower[2] = {0.0, 1.0};
  const double col_upper[2] = {4.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  const double row_lower[3] = {-1.0e30, 5.0, 6.0};
  const double row_upper[3] = {7.0, 15.0, 1.0e30};
  // Define the constraint matrix column-wise
  const int a_format = kHighsHessianFormatTriangular;
  const int a_start[2] = {0, 2};
  const int a_index[5] = {1, 2, 0, 1, 2};
  const double a_value[5] = {1.0, 3.0, 1.0, 2.0, 2.0};

  int run_status;
  int model_status;
  double objective_function_value;
  int simplex_iteration_count;
  int64_t mip_node_count;
  int primal_solution_status;
  int dual_solution_status;
  int basis_validity;

  // Create a Highs instance
  void* highs = Highs_create();

  // Pass the LP to HiGHS
  run_status = Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset,
			    col_cost, col_lower, col_upper,
			    row_lower, row_upper,
			    a_start, a_index, a_value);
  assert(run_status == kHighsStatusOk);

  // Solve the incumbent model
  run_status = Highs_run(highs);
  // The run must be successful
  assert(run_status == kHighsStatusOk);
  // Get the model status - which must be optimal
  model_status = Highs_getModelStatus(highs);
  assert(model_status == kHighsModelStatusOptimal);

  printf("Run status = %d; Model status = %d\n", run_status, model_status);

  // Get scalar information about the solution
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
  Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
  Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
  Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

  // The primal and dual solution status values should indicate feasibility
  assert(primal_solution_status == kHighsSolutionStatusFeasible);
  assert(dual_solution_status == kHighsSolutionStatusFeasible);
  // The basis status should indicate validity
  assert(basis_validity == kHighsBasisValidityValid);

  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* col_dual = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);
  double* row_dual = (double*)malloc(sizeof(double) * num_row);

  int* col_basis_status = (int*)malloc(sizeof(int) * num_col);
  int* row_basis_status = (int*)malloc(sizeof(int) * num_row);

  // Get the primal and dual solution
  Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
  // Get the basis
  Highs_getBasis(highs, col_basis_status, row_basis_status);

  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf; dual = %lf; status = %d; \n", i, col_value[i], col_dual[i], col_basis_status[i]);
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf; dual = %lf; status = %d; \n", i, row_value[i], row_dual[i], row_basis_status[i]);
  }
  printf("Objective value = %g; Iteration count = %d\n", objective_function_value, simplex_iteration_count);

  // Illustrate extraction of model data
  int check_sense;
  run_status = Highs_getObjectiveSense(highs, &check_sense);
  assert(run_status==0);
  printf("LP problem has objective sense = %d\n", check_sense);
  assert(check_sense == sense);
 
  // Illustrate change of model data
  Highs_changeObjectiveSense(highs, kHighsObjSenseMaximize);

  // Solve the incumbent model
  run_status = Highs_run(highs);
  // The run must be successful
  assert(run_status == kHighsStatusOk);
  // Get the model status - which must be optimal
  model_status = Highs_getModelStatus(highs);
  assert(model_status == kHighsModelStatusOptimal);

  printf("Run status = %d; Model status = %d\n", run_status, model_status);

  // Get scalar information about the solution
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
  Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
  Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
  Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);


  // The primal and dual solution status values should indicate feasibility
  assert(primal_solution_status == kHighsSolutionStatusFeasible);
  assert(dual_solution_status == kHighsSolutionStatusFeasible);
  // The basis status should indicate validity
  assert(basis_validity == kHighsBasisValidityValid);

  // Get the primal and dual solution
  Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
  // Get the basis
  Highs_getBasis(highs, col_basis_status, row_basis_status);

  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf; dual = %lf; status = %d; \n", i, col_value[i], col_dual[i], col_basis_status[i]);
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf; dual = %lf; status = %d; \n", i, row_value[i], row_dual[i], row_basis_status[i]);
  }
  printf("Objective value = %g; Iteration count = %d\n", objective_function_value, simplex_iteration_count);
  
  // Now illustrate how LPs can be built within HiGHS by constructing
  // the same maximization problem. First clear the incumbent model
  Highs_clearModel(highs);

  // Demonstrate that the incumbent model is empty
  const int check_num_col = Highs_getNumCol(highs);
  const int check_num_row = Highs_getNumRow(highs);
  const int check_num_nz = Highs_getNumNz(highs);
  assert(check_num_col == 0);
  assert(check_num_row == 0);
  assert(check_num_nz == 0);
  printf("\nCleared model has %d columns, %d rows and %d nonzeros\n",
	 check_num_col, check_num_row, check_num_nz);

  // Define the constraint matrix row-wise, as it is added to the LP
  // with the rows
  const int ar_start[3] = {0, 1, 3};
  const int ar_index[5] = {1, 0, 1, 0, 1};
  const double ar_value[5] = {1.0, 1.0, 2.0, 3.0, 2.0};

  // Add two columns to the empty LP
  run_status = Highs_addCols(highs, num_col, col_cost, col_lower, col_upper, 0, NULL, NULL, NULL);
  assert(run_status==0);
  // Add three rows to the 2-column LP
  run_status = Highs_addRows(highs, num_row, row_lower, row_upper, num_nz, ar_start, ar_index, ar_value);
  assert(run_status==0);

  // By default, the optimization sense is minimization, and the
  // objective offset is zero, so these need to be changed
  Highs_changeObjectiveSense(highs, kHighsObjSenseMaximize);
  Highs_changeObjectiveOffset(highs, offset);

  // Illustrate extraction of model data
  run_status = Highs_getObjectiveSense(highs, &check_sense);
  assert(run_status==0);
  printf("LP problem has objective sense = %d\n", check_sense);
  assert(check_sense == kHighsObjSenseMaximize);

  // Illustrate a change of option

  // Check what type of option value you should provide
  int option_type;
  const char* option_string = "primal_feasibility_tolerance";
  run_status = Highs_getOptionType(highs, option_string, &option_type);
  printf("Option %s is of type %d\n", option_string, option_type);
  assert(run_status == kHighsStatusOk);
  assert(option_type == 2);

  double primal_feasibility_tolerance;
  Highs_getDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  printf("primal_feasibility_tolerance = %g: setting it to 1e-6\n", primal_feasibility_tolerance);
  primal_feasibility_tolerance = 1e-6;
  Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", primal_feasibility_tolerance);

  // Illustrate how HiGHS can run quietly
  Highs_setBoolOptionValue(highs, "output_flag", 0);
  printf("Running quietly...\n");
  // Solve the incumbent model
  run_status = Highs_run(highs);
  printf("Running loudly...\n");
  Highs_setBoolOptionValue(highs, "output_flag", 1);

  assert(run_status==0);
  // Get the model status - which must be optimal
  model_status = Highs_getModelStatus(highs);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
  Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
  Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
  Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

  // The primal and dual solution status values should indicate feasibility
  assert(primal_solution_status == kHighsSolutionStatusFeasible);
  assert(dual_solution_status == kHighsSolutionStatusFeasible);
  // The basis status should indicate validity
  assert(basis_validity == kHighsBasisValidityValid);

  // Get the primal
  Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
  Highs_getBasis(highs, col_basis_status, row_basis_status);

  // Report the column primal and dual values, and basis status
  for (int i = 0; i < num_col; i++) {
    printf("Col%d = %lf; dual = %lf; status = %d; \n", i, col_value[i], col_dual[i], col_basis_status[i]);
  }
  // Report the row primal and dual values, and basis status
  for (int i = 0; i < num_row; i++) {
    printf("Row%d = %lf; dual = %lf; status = %d; \n", i, row_value[i], row_dual[i], row_basis_status[i]);
  }
  printf("Objective value = %g; Iteration count = %d\n", objective_function_value, simplex_iteration_count);
  
  // 
  // Indicate that the optimal solution for both columns must be
  // integer valued and solve the model as a MIP
  int integrality[2] = {1, 1};
  Highs_changeColsIntegralityByRange(highs, 0, 1, integrality);

  // Solve the incumbent model quietly
  Highs_setBoolOptionValue(highs, "output_flag", 0);
  run_status = Highs_run(highs);
  Highs_setBoolOptionValue(highs, "output_flag", 1);

  assert(run_status == kHighsStatusOk);
  // Get the model status - which must be optimal
  model_status = Highs_getModelStatus(highs);
  assert(model_status == kHighsModelStatusOptimal);

  printf("\nRun status = %d; Model status = %d\n", run_status, model_status);

  // Get scalar information about the solution
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
  Highs_getInt64InfoValue(highs, "mip_node_count", &mip_node_count);
  Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
  Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
  Highs_getIntInfoValue(highs, "basis_validity", &basis_validity);

  // The primal solution status value should indicate feasibility
  assert(primal_solution_status == kHighsSolutionStatusFeasible);
  // The dual solution status value should indicate no solution
  assert(dual_solution_status == kHighsSolutionStatusNone);
  // The basis status value should indicate invalidity
  assert(basis_validity == kHighsBasisValidityInvalid);
  assert(mip_node_count == 1);

  // Get the primal and dual solution
  Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);

  // Report the column primal values
  for (int i = 0; i < num_col; i++)
    printf("Col%d = %lf\n", i, col_value[i]);
  // Report the row primal values
  for (int i = 0; i < num_row; i++)
    printf("Row%d = %lf\n", i, row_value[i]);
  printf("Objective value = %g; Iteration count = %d\n", objective_function_value, simplex_iteration_count);

  free(col_value);
  free(col_dual);
  free(row_value);
  free(row_dual);
  free(col_basis_status);
  free(row_basis_status);

  Highs_destroy(highs);
}

int main() {
  minimal_api();
  minimal_api_qp();
  minimal_api_mps();
  full_api();
  return 0;
}
