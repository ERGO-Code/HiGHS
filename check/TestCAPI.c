#include "interfaces/highs_c_api.h"

#include "HCheckConfig.h"
#include <stdio.h>
#include <stdlib.h>
// Force asserts to be checked always.
#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>

const HighsInt dev_run = 0;
const double double_equal_tolerance = 1e-5;

void checkGetCallbackDataOutPointer(const HighsCallbackDataOut* data_out, const char* name, HighsInt valid) {
  const void* name_p = Highs_getCallbackDataOutItem(data_out, name);
  if (valid) {
    if (!name_p) printf("checkGetCallbackDataOutItem fail for %s (valid = %d)\n", name, (int)valid);
    assert(name_p);
  } else {
    if (name_p) printf("checkGetCallbackDataOutItem fail for %s (valid = %d)\n",
			name, (int)valid);
    assert(!name_p);
  }
}
    
void checkGetCallbackDataOutHighsInt(const HighsCallbackDataOut* data_out, const char* name, HighsInt value) {
  const void* name_p = Highs_getCallbackDataOutItem(data_out, name);
  if (!name_p) {
    printf("checkGetCallbackDataOutItem fail for %s\n", name);
    assert(name_p);
  } else {
    HighsInt check_value = *(HighsInt*)(name_p);
    HighsInt value_ok = check_value == value;
    if (!value_ok) printf("checkGetCallbackDataOutItem fail for %s (%d = check_value != value = %d)\n",
			  name, (int)check_value, (int)value);
    assert(value_ok);
  }
}
    
void checkGetCallbackDataOutInt(const HighsCallbackDataOut* data_out, const char* name, int value) {
  const void* name_p = Highs_getCallbackDataOutItem(data_out, name);
  if (!name_p) {
    printf("checkGetCallbackDataOutInt fail for %s\n", name);
    assert(name_p);
  } else {
    int check_value = *(int*)(name_p);
    int value_ok = check_value == value;
    if (!value_ok) printf("checkGetCallbackDataOutInt fail for %s (%d = check_value != value = %d)\n",
			  name, check_value, value);
    assert(value_ok);
  }
}
    
void checkGetCallbackDataOutInt64(const HighsCallbackDataOut* data_out, const char* name, int64_t value) {
  const void* name_p = Highs_getCallbackDataOutItem(data_out, name);
  if (!name_p) {
    printf("checkGetCallbackDataOutInt64 fail for %s\n", name);
    assert(name_p);
  } else {
    int64_t check_value = *(int*)(name_p);
    int value_ok = check_value == value;
    if (!value_ok) printf("checkGetCallbackDataOutInt64 fail for %s (%d = check_value != value = %d)\n",
			  name, (int)check_value, (int)value);
    assert(value_ok);
  }
}
    
void checkGetCallbackDataOutDouble(const HighsCallbackDataOut* data_out, const char* name, double value) {
  const void* name_p = Highs_getCallbackDataOutItem(data_out, name);
  if (!name_p) {
    printf("checkGetCallbackDataOutDouble fail for %s\n", name);
    assert(name_p);
  } else {
    double check_value = *(double*)(name_p);
    double value_ok = check_value == value;
    if (!value_ok) printf("checkGetCallbackDataOutDouble fail for %s (%g = check_value != value = %g)\n",
			  name, check_value, value);
    assert(value_ok);
  }
}
    
static void userCallback(const int callback_type, const char* message,
			 const HighsCallbackDataOut* data_out,
			 HighsCallbackDataIn* data_in,
			 void* user_callback_data) {
  // Extract the double value pointed to from void* user_callback_data
  const double local_callback_data = user_callback_data == NULL ? -1 : *(double*)user_callback_data;

  if (callback_type == kHighsCallbackLogging) {
    if (dev_run) printf("userCallback(%11.4g): %s\n", local_callback_data, message);
  } else if (callback_type == kHighsCallbackMipImprovingSolution) {
    // Test the accessor function for data_out
    //
    // Check that passing an valid name returns a non-null pointer,
    // and that the corresponding value is the same as obtained using
    // the struct
    const void* objective_function_value_p =
	Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutObjectiveFunctionValueName);
    assert(objective_function_value_p);
    double objective_function_value = *(double*)(objective_function_value_p);
    assert(objective_function_value == data_out->objective_function_value);
    if (dev_run) printf("userCallback(%11.4g): improving solution with objective = %g\n",
			local_callback_data, objective_function_value);
    // Now test all more simply
    checkGetCallbackDataOutInt(data_out,
			       kHighsCallbackDataOutLogTypeName, -1);
    checkGetCallbackDataOutDouble(data_out,
				  kHighsCallbackDataOutRunningTimeName,
				  data_out->running_time);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutSimplexIterationCountName,
				    data_out->simplex_iteration_count);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutIpmIterationCountName,
				    data_out->ipm_iteration_count);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutPdlpIterationCountName,
				    data_out->pdlp_iteration_count);
    checkGetCallbackDataOutDouble(data_out,
				  kHighsCallbackDataOutObjectiveFunctionValueName,
				  data_out->objective_function_value);
    checkGetCallbackDataOutInt64(data_out,
				 kHighsCallbackDataOutMipNodeCountName,
				 data_out->mip_node_count);
    checkGetCallbackDataOutDouble(data_out,
				  kHighsCallbackDataOutMipPrimalBoundName,
				  data_out->mip_primal_bound);
    checkGetCallbackDataOutDouble(data_out,
				  kHighsCallbackDataOutMipDualBoundName,
				  data_out->mip_dual_bound);
    checkGetCallbackDataOutDouble(data_out,
				  kHighsCallbackDataOutMipGapName,
				  data_out->mip_gap);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutCutpoolNumColName, 0);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutCutpoolNumCutName, 0);
    checkGetCallbackDataOutHighsInt(data_out,
				    kHighsCallbackDataOutCutpoolNumNzName, 0);

    // Check that passing an unrecognised name returns NULL
    const void* foo_p = Highs_getCallbackDataOutItem(data_out, "foo");
    assert(!foo_p);
    // Check that passing the name of an assigned vector returns
    // non-NULL, and that the corresponding value is the same as
    // obtained using the struct
    const void* mip_solution_void_p =
      Highs_getCallbackDataOutItem(data_out,
				   kHighsCallbackDataOutMipSolutionName);
    assert(mip_solution_void_p);
    double mip_solution0 = *(double*)(mip_solution_void_p);
    assert(mip_solution0 == *(data_out->mip_solution));
    if (dev_run) printf("userCallback(%11.4g): improving solution with value[0] = %g\n",
			local_callback_data, mip_solution0);
    // Check that passing names of the unassigned vectors returns NULL
    assert(!Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutCutpoolStartName));
    assert(!Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutCutpoolIndexName));
    assert(!Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutCutpoolValueName));
    assert(!Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutCutpoolLowerName));
    assert(!Highs_getCallbackDataOutItem(data_out, kHighsCallbackDataOutCutpoolUpperName));
  } else if (callback_type == kHighsCallbackMipLogging) {
    if (dev_run) printf("userCallback(%11.4g): MIP logging\n", local_callback_data);
    data_in->user_interrupt = 1;
  } else if (callback_type == kHighsCallbackMipInterrupt) {
    if (dev_run) printf("userCallback(%11.4g): MIP interrupt\n", local_callback_data);
    data_in->user_interrupt = 1;
  }
}

HighsInt highsIntArraysEqual(const HighsInt dim, const HighsInt* array0, const HighsInt* array1) {
  for (HighsInt ix = 0; ix < dim; ix++) if (array0[ix] != array1[ix]) return 0;
  return 1;
}

HighsInt doubleArraysEqual(const double dim, const double* array0, const double* array1) {
  for (HighsInt ix = 0; ix < dim; ix++) if (array0[ix] != array1[ix]) return 0;
  return 1;
}

void assertDoubleValuesEqual(const char* name, const double is, const double should_be) {
  const double dl = fabs(is-should_be);
  if (dl > double_equal_tolerance) {
    printf("Value %s = %g differs from %g by %g but should be equal\n", name, is, should_be, dl);
    assert(1==0);
  }
}

void assertIntValuesEqual(const char* name, const HighsInt is, const HighsInt should_be) {
  if (is != should_be) {
    printf("Value %s = %"HIGHSINT_FORMAT" should be %"HIGHSINT_FORMAT"\n", name, is, should_be);
    assert(1==0);
  }
}

void assertLogical(const char* name, const HighsInt is) {
  if (is == 0) {
    printf("Value %s = %"HIGHSINT_FORMAT" should not be 0\n", name, is);
    assert(1==0);
  }
}

void version_api() {
  if (dev_run) {
    printf("HiGHS version %s\n", Highs_version());
    printf("HiGHS version major %"HIGHSINT_FORMAT"\n", Highs_versionMajor());
    printf("HiGHS version minor %"HIGHSINT_FORMAT"\n", Highs_versionMinor());
    printf("HiGHS version patch %"HIGHSINT_FORMAT"\n", Highs_versionPatch());
    printf("HiGHS githash: %s\n", Highs_githash());
    // Compilation date is deprecated.
    // printf("HiGHS compilation date %s\n", Highs_compilationDate());
  }
}

void minimal_api_lp() {
  // This illustrates the use of Highs_call, the simple C interface to
  // HiGHS. It's designed to solve the general LP problem
  //
  // Min c^Tx subject to L <= Ax <= U; l <= x <= u
  //
  // where A is a matrix with m rows and n columns
  //
  // The scalar n is num_col
  // The scalar m is num_row
  //
  // The vector c is col_cost
  // The vector l is col_lower
  // The vector u is col_upper
  // The vector L is row_lower
  // The vector U is row_upper
  //
  // The matrix A is represented in packed column-wise form: only its
  // nonzeros are stored
  //
  // * The number of nonzeros in A is num_nz
  //
  // * The row indices of the nonnzeros in A are stored column-by-column
  // in a_index
  //
  // * The values of the nonnzeros in A are stored column-by-column in
  // a_value
  //
  // * The position in a_index/a_value of the index/value of the first
  // nonzero in each column is stored in a_start
  //
  // Note that a_start[0] must be zero
  //
  // After a successful call to Highs_call, the primal and dual
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
  // To solve maximization problems, the values in c must be negated
  //
  // The use of Highs_lpCall is illustrated for the LP
  //
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  const HighsInt num_col = 2;
  const HighsInt num_row = 3;
  const HighsInt num_nz = 5;
  HighsInt a_format = kHighsMatrixFormatColwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;

  // Define the column costs, lower bounds and upper bounds
  double col_cost[2] = {2.0, 3.0};
  double col_lower[2] = {0.0, 1.0};
  double col_upper[2] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double row_lower[3] = {-1.0e30, 10.0, 8.0};
  double row_upper[3] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix column-wise
  HighsInt a_start[2] = {0, 2};
  HighsInt a_index[5] = {1, 2, 0, 1, 2};
  double a_value[5] = {1.0, 2.0, 1.0, 2.0, 1.0};

  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* col_dual = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);
  double* row_dual = (double*)malloc(sizeof(double) * num_row);

  HighsInt* col_basis_status = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* row_basis_status = (HighsInt*)malloc(sizeof(HighsInt) * num_row);

  HighsInt model_status;

  HighsInt return_status = Highs_lpCall(num_col, num_row, num_nz, a_format,
					sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
					a_start, a_index, a_value,
					col_value, col_dual, row_value, row_dual,
					col_basis_status, row_basis_status,
					&model_status);

  assert( return_status == kHighsStatusOk );

  if (dev_run) {
    printf("Run status = %"HIGHSINT_FORMAT"; Model status = %"HIGHSINT_FORMAT"\n", return_status, model_status);
    
    HighsInt i;
    if (model_status == kHighsModelStatusOptimal) {
      double objective_value = 0;
      // Report the column primal and dual values, and basis status
      for (i = 0; i < num_col; i++) {
	printf("Col%"HIGHSINT_FORMAT" = %lf; dual = %lf; status = %"HIGHSINT_FORMAT"; \n",
	       i, col_value[i], col_dual[i], col_basis_status[i]);
	objective_value += col_value[i]*col_cost[i];
      }
      // Report the row primal and dual values, and basis status
      for (i = 0; i < num_row; i++) {
	printf("Row%"HIGHSINT_FORMAT" = %lf; dual = %lf; status = %"HIGHSINT_FORMAT"; \n",
	       i, row_value[i], row_dual[i], row_basis_status[i]);
      }
      printf("Optimal objective value = %g\n", objective_value);
    }
  }

  free(col_value);
  free(col_dual);
  free(row_value);
  free(row_dual);
  free(col_basis_status);
  free(row_basis_status);
}

void minimal_api_mip() {
  // The use of Highs_mipCall is illustrated for the MIP
  //
  // Min    f  = -3x_0 - 2x_1 - x_2
  // s.t.          x_0 +  x_1 + x_2 <=  7
  //              4x_0 + 2x_1 + x_2  = 12
  //              x_0 >=0; x_1 >= 0; x_2 binary

  const HighsInt num_col = 3;
  const HighsInt num_row = 2;
  const HighsInt num_nz = 6;
  HighsInt a_format = kHighsMatrixFormatColwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;

  // Define the column costs, lower bounds and upper bounds
  double col_cost[3] = {-3.0, -2.0, -1.0};
  double col_lower[3] = {0.0, 0.0, 0.0};
  double col_upper[3] = {1.0e30, 1.0e30, 1.0};
  // Define the row lower bounds and upper bounds
  double row_lower[2] = {-1.0e30, 12.0};
  double row_upper[2] = { 7.0,    12.0};
  // Define the constraint matrix column-wise
  HighsInt a_start[3] = {0, 2, 4};
  HighsInt a_index[6] = {0, 1, 0, 1, 0, 1};
  double a_value[6] = {1.0, 4.0, 1.0, 2.0, 1.0, 1.0};

  // Give an illegal value to an entry in integrality
  HighsInt integrality[3] = {kHighsVarTypeContinuous, kHighsVarTypeContinuous, -1};

  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);

  HighsInt model_status;
  HighsInt return_status;

  return_status = Highs_mipCall(num_col, num_row, num_nz, a_format,
				sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
				a_start, a_index, a_value,
				integrality,
				col_value, row_value,
				&model_status);
  // Should return error, with model status not set
  assert( return_status == kHighsStatusError );
  assert( model_status == kHighsModelStatusNotset );

  // Correct integrality
  integrality[num_col-1] = kHighsVarTypeInteger;

  return_status = Highs_mipCall(num_col, num_row, num_nz, a_format,
				sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
				a_start, a_index, a_value,
				integrality,
				col_value, row_value,
				&model_status);
  // Should return OK
  assert( return_status == kHighsStatusOk );

  if (dev_run) {
    printf("Run status = %"HIGHSINT_FORMAT"; Model status = %"HIGHSINT_FORMAT"\n", return_status, model_status);
    
    HighsInt i;
    if (model_status == kHighsModelStatusOptimal) {
      double objective_value = 0;
      // Report the column primal values
      for (i = 0; i < num_col; i++) {
	printf("Col%"HIGHSINT_FORMAT" = %lf; \n", i, col_value[i]);
	objective_value += col_value[i]*col_cost[i];
      }
      // Report the row primal values
      for (i = 0; i < num_row; i++) {
	printf("Row%"HIGHSINT_FORMAT" = %lf; \n", i, row_value[i]);
      }
      printf("Optimal objective value = %g\n", objective_value);
    }
  }

  free(col_value);
  free(row_value);

}

void minimal_api_qp() {
  // Test solving the problem qjh
  //
  // minimize -x_2 - 3x_3 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
  //
  // subject to x_1 + x_3 <= 2; x>=0
  const double inf = 1e30;
  HighsInt num_col = 3;
  HighsInt num_row = 1;
  HighsInt num_nz = 2;
  HighsInt q_num_nz = 4;
  HighsInt a_format = kHighsMatrixFormatColwise;
  HighsInt q_format = kHighsHessianFormatTriangular;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;
  double col_cost[3] = {0.0, -1.0, -3.0};
  double col_lower[3] = {-inf, -inf, -inf};
  double col_upper[3] = {inf, inf, inf};
  double row_lower[1] = {-inf};
  double row_upper[1] = {2};
  HighsInt a_start[3] = {0, 1, 1};
  HighsInt a_index[2] = {0, 0};
  double a_value[2] = {1.0, 1.0};
  HighsInt q_start[3] = {0, 2, 3};
  HighsInt q_index[4] = {0, 2, 1, 2};
  double q_value[4] = {2.0, -1.0, 0.2, 2.0};
  
  double* col_value = (double*)malloc(sizeof(double) * num_col);
  HighsInt model_status;
  HighsInt return_status = Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset,
					col_cost, col_lower, col_upper, row_lower, row_upper,
					a_start, a_index, a_value, q_start, q_index, q_value,
					col_value, NULL, NULL, NULL, NULL, NULL, &model_status);
  assert( return_status == kHighsStatusOk );
  assertIntValuesEqual("Model status for QP qph", model_status, kHighsModelStatusOptimal);
  double required_x[3] = {0.5, 5.0, 1.5};
  if (dev_run) {
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      printf("x%d1 = %g\n", (int)iCol, col_value[iCol]);
      assertDoubleValuesEqual("Solution value for QP qph", col_value[iCol], required_x[iCol]);
    }
  }
  free(col_value);
}

void minimal_api_illegal_lp() {
  const double inf = 1e30;
  HighsInt num_col = 2;
  HighsInt num_row = 1;
  HighsInt num_nz = 2;
  HighsInt a_format = kHighsMatrixFormatRowwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;
  double col_cost[2] = {0.0, -1.0};
  double col_lower[2] = {-inf, -inf};
  double col_upper[2] = {inf, inf};
  double row_lower[1] = {-inf};
  double row_upper[1] = {2};
  HighsInt a_start[1] = {0};
  HighsInt a_index[2] = {0, -1}; // Illegal index
  double a_value[2] = {1.0, 1.0};

  HighsInt model_status;
  HighsInt return_status = Highs_lpCall(num_col, num_row, num_nz, a_format,
					sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
					a_start, a_index, a_value,
					NULL, NULL, NULL, NULL, 
					NULL, NULL, 
					&model_status);
  // Should return error, with model status not set
  assert( return_status == kHighsStatusError );
  assert( model_status == kHighsModelStatusNotset );
}

void full_api() {
  void* highs = Highs_create();

  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);

  HighsInt num_col = 2;
  HighsInt num_row = 2;
  HighsInt num_nz = 4;
  HighsInt a_format = kHighsMatrixFormatRowwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;
  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  HighsInt a_start[3] = {0, 2, 4};
  HighsInt a_index[4] = {0, 1, 0, 1};
  double a_value[4] = {1.0, 2.0, 1.0, 3.0};

  assert( Highs_addCols(highs, 2, cc, cl, cu, 0, NULL, NULL, NULL) == 0);
  assert( Highs_addRows(highs, 2, rl, ru,  4, a_start, a_index, a_value) == 0);

  assert( Highs_getNumCols(highs) == num_col);
  assert( Highs_getNumRows(highs) == num_row);
  assert( Highs_getNumNz(highs) == num_nz);
  assert( Highs_getHessianNumNz(highs) == 0);

  HighsInt ck_num_col;
  HighsInt ck_num_row;
  HighsInt ck_num_nz;
  HighsInt ck_hessian_num_nz;
  HighsInt ck_rowwise;
  HighsInt ck_sense;
  double ck_offset;
  double ck_cc[2];
  double ck_cl[2];
  double ck_cu[2];
  double ck_rl[2];
  double ck_ru[2];
  HighsInt ck_a_start[3];
  HighsInt ck_a_index[4];
  double ck_a_value[4];
  HighsInt return_status;
  return_status = Highs_getModel(highs, a_format, 0,
				 &ck_num_col, &ck_num_row, &ck_num_nz, NULL,
				 &ck_sense, &ck_offset,
				 ck_cc, ck_cl, ck_cu, ck_rl, ck_ru,
				 ck_a_start, ck_a_index, ck_a_value,
				 NULL, NULL, NULL, NULL);
  assert( return_status == kHighsStatusOk );

  assert( ck_num_col == num_col );
  assert( ck_num_row == num_row );
  assert( ck_num_nz == num_nz );
  assert( ck_sense == sense );
  assert( ck_offset == offset );
  assert( doubleArraysEqual(num_col, ck_cc, cc) );
  assert( doubleArraysEqual(num_col, ck_cl, cl) );
  assert( doubleArraysEqual(num_col, ck_cu, cu) );
  assert( doubleArraysEqual(num_row, ck_rl, rl) );
  assert( doubleArraysEqual(num_row, ck_ru, ru) );
  assert( highsIntArraysEqual(num_col, ck_a_start, a_start) );
  assert( highsIntArraysEqual(num_nz, ck_a_index, a_index) );
  assert( doubleArraysEqual(num_nz, ck_a_value, a_value) );

  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );

  char* col_prefix = "Col";
  char* row_prefix = "Row";
  // Check index out of bounds
  return_status = Highs_passColName(highs, -1, col_prefix);
  assert( return_status == kHighsStatusError );
  return_status = Highs_passColName(highs, num_col, col_prefix);
  assert( return_status == kHighsStatusError );
  return_status = Highs_passRowName(highs, -1, row_prefix);
  assert( return_status == kHighsStatusError );
  return_status = Highs_passRowName(highs, num_row, row_prefix);
  assert( return_status == kHighsStatusError );

  // Define all column names to be the same
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    return_status = Highs_passColName(highs, iCol, col_prefix);
    assert( return_status == kHighsStatusOk );
  }
  return_status = Highs_writeModel(highs, "");
  assert( return_status == kHighsStatusError );

  // Define all column names to be different
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    const char suffix = iCol + '0';
    const char* suffix_p = &suffix;
    char name[5];  // 3 chars prefix, 1 char iCol, 1 char 0-terminator
    sprintf(name, "%s%" HIGHSINT_FORMAT "", col_prefix, iCol);
    const char* name_p = name;
    return_status = Highs_passColName(highs, iCol, name_p);
    assert( return_status == kHighsStatusOk );
  }
  return_status = Highs_writeModel(highs, "");
  assert( return_status == kHighsStatusOk );

  // Check that the columns can be found by name
  HighsInt ck_iCol;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    char name[5];
    return_status = Highs_getColName(highs, iCol, name);
    assert( return_status == kHighsStatusOk );
    return_status = Highs_getColByName(highs, name, &ck_iCol);
    assert( return_status == kHighsStatusOk );
    assert( ck_iCol == iCol );
  }
  return_status = Highs_getColByName(highs, "FRED", &ck_iCol);
  assert( return_status == kHighsStatusError );
  
  // Define all row names to be the same
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    return_status = Highs_passRowName(highs, iRow, row_prefix);
    assert( return_status == kHighsStatusOk );
  }
  return_status = Highs_writeModel(highs, "");
  assert( return_status == kHighsStatusError );
  
  // Define all row names to be different
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const char suffix = iRow + '0';
    const char* suffix_p = &suffix;
    char name[5];  // 3 chars prefix, 1 char iCol, 1 char 0-terminator
    sprintf(name, "%s%" HIGHSINT_FORMAT "", row_prefix, iRow);
    const char* name_p = name;
    return_status = Highs_passRowName(highs, iRow, name_p);
    assert( return_status == kHighsStatusOk );
  }
  return_status = Highs_writeModel(highs, "");
  assert( return_status == kHighsStatusOk );
  
  // Check that the rows can be found by name
  HighsInt ck_iRow;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    char name[5];
    return_status = Highs_getRowName(highs, iRow, name);
    assert( return_status == kHighsStatusOk );
    return_status = Highs_getRowByName(highs, name, &ck_iRow);
    assert( return_status == kHighsStatusOk );
    assert( ck_iRow == iRow );
  }
  return_status = Highs_getRowByName(highs, "FRED", &ck_iRow);
  assert( return_status == kHighsStatusError );
  
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    char name[5];
    char* name_p = name;
    return_status = Highs_getColName(highs, iCol, name_p);
    assert( return_status == kHighsStatusOk );
    if (dev_run) printf("Column %" HIGHSINT_FORMAT " has name %s\n", iCol, name_p);
  }
  
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    char name[5];
    char* name_p = name;
    return_status = Highs_getRowName(highs, iRow, name_p);
    assert( return_status == kHighsStatusOk );
    if (dev_run) printf("Row    %" HIGHSINT_FORMAT " has name %s\n", iRow, name_p);
  }

  Highs_destroy(highs);
}

void full_api_options() {
  void* highs;

  highs = Highs_create();
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);

  const double kHighsInf = Highs_getInfinity(highs);
  HighsInt simplex_scale_strategy;
  HighsInt return_status;
  return_status = Highs_getIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("simplex_scale_strategy = %"HIGHSINT_FORMAT": setting it to 3\n", simplex_scale_strategy);
  simplex_scale_strategy = 3;
  return_status = Highs_setIntOptionValue(highs, "simplex_scale_strategy", simplex_scale_strategy);

  const HighsInt presolve_index = 0;
  char* name = NULL;
  return_status = Highs_getOptionName(highs, presolve_index, &name);
  if (dev_run) printf("option %"HIGHSINT_FORMAT" has name %s\n", presolve_index, name);
  const char* presolve = "presolve";
  assert( *name == *presolve );
  free(name);

  HighsInt check_simplex_scale_strategy;
  HighsInt min_simplex_scale_strategy;
  HighsInt max_simplex_scale_strategy;
  HighsInt default_simplex_scale_strategy;
  return_status = Highs_getIntOptionValues(highs, "scale_strategy", NULL, NULL, NULL, NULL);
  assert( return_status == kHighsStatusError );
  return_status = Highs_getDoubleOptionValues(highs, "simplex_scale_strategy", NULL, NULL, NULL, NULL);
  assert( return_status == kHighsStatusError );
  return_status = Highs_getIntOptionValues(highs, "simplex_scale_strategy",
					   &check_simplex_scale_strategy,
					   &min_simplex_scale_strategy,
					   &max_simplex_scale_strategy,
					   &default_simplex_scale_strategy);
  assert( return_status == kHighsStatusOk );
  assert( check_simplex_scale_strategy == simplex_scale_strategy );
  assert( min_simplex_scale_strategy ==  0 );
  assert( max_simplex_scale_strategy ==  5 );
  assert( default_simplex_scale_strategy == 1 );


  // There are some functions to check what type of option value you should
  // provide.
  HighsInt option_type;
  return_status = Highs_getOptionType(highs, "simplex_scale_strategy", &option_type);
  assert( return_status == kHighsStatusOk );
  assert( option_type == kHighsOptionTypeInt );
  return_status = Highs_getOptionType(highs, "bad_option", &option_type);
  assert( return_status == kHighsStatusError );

  double primal_feasibility_tolerance;
  return_status = Highs_getDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("primal_feasibility_tolerance = %g: setting it to 1e-6\n", primal_feasibility_tolerance);
  primal_feasibility_tolerance = 1e-6;
  return_status = Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", primal_feasibility_tolerance);
  assert( return_status == kHighsStatusOk );

  double check_primal_feasibility_tolerance;
  return_status = Highs_getDoubleOptionValues(highs, "primal_feasibility_tolerance",
  					      &check_primal_feasibility_tolerance, NULL, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  assert( check_primal_feasibility_tolerance == primal_feasibility_tolerance );
  double default_primal_feasibility_tolerance;
  double min_primal_feasibility_tolerance;
  double max_primal_feasibility_tolerance;
  return_status = Highs_getDoubleOptionValues(highs, "primal_feasibility_tolerance",
  					      &check_primal_feasibility_tolerance,
  					      &min_primal_feasibility_tolerance,
  					      &max_primal_feasibility_tolerance,
  					      &default_primal_feasibility_tolerance);
  assert( min_primal_feasibility_tolerance == 1e-10 );
  assert( max_primal_feasibility_tolerance == kHighsInf );
  assert( default_primal_feasibility_tolerance ==  1e-7 );

  Highs_setStringOptionValue(highs, "presolve", "off");

  return_status = Highs_getStringOptionValues(highs, "pre-solve", NULL, NULL);
  assert( return_status == kHighsStatusError );
  //  char check_presolve_value[kHighsMaximumStringLength];
  char check_presolve_value[512];
  return_status = Highs_getStringOptionValues(highs, "presolve", check_presolve_value, NULL);
  assert( return_status == kHighsStatusOk );

  // const HighsInt output_flag = 1;
  // return_status = Highs_setBoolOptionValue(highs, "output_flag", output_flag);
  return_status = Highs_setBoolOptionValue(highs, "output_flag", 1);

  assert( return_status == kHighsStatusOk );

  HighsInt check_output_flag, default_output_flag;
  return_status = Highs_getBoolOptionValues(highs, "output_flag", NULL, NULL);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_getBoolOptionValues(highs, "output_flag",
					    &check_output_flag, NULL);
  assert( return_status == kHighsStatusOk );
  //    assert( check_output_flag == output_flag );
  assert( check_output_flag == 1 );
  return_status = Highs_getBoolOptionValues(highs, "output_flag",
					    &check_output_flag,
					    &default_output_flag);
  assert( return_status == kHighsStatusOk );
  //    assert( default_output_flag == output_flag );
  assert( default_output_flag == 1 );
  
  HighsInt num_string_option = 0;
  char* option = NULL;
  HighsInt type;
  HighsInt num_options = Highs_getNumOptions(highs);
  char current_string_value[512];
 
  if (dev_run)
    printf("\nString options are:\n");
  for (HighsInt index = 0; index < num_options; index++) {
    Highs_getOptionName(highs, index, &option);
    Highs_getOptionType(highs, option, &type);
    if (type != kHighsOptionTypeString) {
      free(option);
      continue;
    }
    Highs_getStringOptionValues(highs, option, current_string_value, NULL);
    num_string_option++;
    if (dev_run)
      printf("%"HIGHSINT_FORMAT": %-24s \"%s\"\n",
	     num_string_option, option, current_string_value);    
    free(option);
  }

  Highs_destroy(highs);

}

void full_api_lp() {
  // Form and solve the LP
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  void* highs;

  highs = Highs_create();
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);

  const HighsInt num_col = 2;
  const HighsInt num_row = 3;
  const HighsInt num_nz = 5;

  // Define the column costs, lower bounds and upper bounds
  double col_cost[2] = {2.0, 3.0};
  double col_lower[2] = {0.0, 1.0};
  double col_upper[2] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double row_lower[3] = {-1.0e30, 10.0, 8.0};
  double row_upper[3] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix row-wise, as it is added to the LP
  // with the rows
  HighsInt arstart[3] = {0, 1, 3};
  HighsInt arindex[5] = {1, 0, 1, 0, 1};
  double arvalue[5] = {1.0, 1.0, 2.0, 2.0, 1.0};

  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* col_dual = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);
  double* row_dual = (double*)malloc(sizeof(double) * num_row);

  HighsInt* col_basis_status = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* row_basis_status = (HighsInt*)malloc(sizeof(HighsInt) * num_row);

  // Add two columns to the empty LP
  assert( Highs_addCols(highs, num_col, col_cost, col_lower, col_upper, 0, NULL, NULL, NULL) == 0);
  // Add three rows to the 2-column LP
  assert( Highs_addRows(highs, num_row, row_lower, row_upper, num_nz, arstart, arindex, arvalue) == 0);

  HighsInt sense;
  HighsInt return_status;
  return_status = Highs_getObjectiveSense(highs, &sense);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("LP problem has objective sense = %"HIGHSINT_FORMAT"\n", sense);
  assert( sense == kHighsObjSenseMinimize );

  sense *= -1;
  return_status = Highs_changeObjectiveSense(highs, sense);
  assert( return_status == kHighsStatusOk );
  assert( sense == kHighsObjSenseMaximize );

  sense *= -1;
  return_status = Highs_changeObjectiveSense(highs, sense);
  assert( return_status == kHighsStatusOk );

  return_status = Highs_getObjectiveSense(highs, &sense);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("LP problem has old objective sense = %"HIGHSINT_FORMAT"\n", sense);
  assert( sense == kHighsObjSenseMinimize );



  // fetch column data (just first column)
  {
    const HighsInt get_col = 0;
    const HighsInt num_get_col = 1;
    HighsInt get_num_col = 0;
    double* get_costs = (double*)malloc(sizeof(double) * num_get_col);
    double* get_lower = (double*)malloc(sizeof(double) * num_get_col);
    double* get_upper = (double*)malloc(sizeof(double) * num_get_col);
    HighsInt get_num_nz = 0;

    return_status = Highs_getColsByRange(highs, get_col, get_col,
					 &get_num_col, get_costs, get_lower, get_upper, &get_num_nz,
					 NULL, NULL, NULL);
    assert( return_status == kHighsStatusOk );

    assertIntValuesEqual("getCols get_num_col", get_num_col, num_get_col);
    assertDoubleValuesEqual("getCols get_costs", get_costs[0], col_cost[get_col]);
    assertDoubleValuesEqual("getCols get_lower", get_lower[0], col_lower[get_col]);
    assertDoubleValuesEqual("getCols get_upper", get_upper[0], col_upper[get_col]);
    assertIntValuesEqual("getCols get_num_nz", get_num_nz, 2);

      // could also check coefficients by calling again...

    free(get_upper);
    free(get_lower);
    free(get_costs);
  }

  // fetch row data (just 2nd row: 10 <=  x_0 + 2x_1 <= 14)
  {
    const HighsInt get_row = 1;
    const HighsInt num_get_row = 1;
    HighsInt get_num_row = 0;
    double* get_lower = (double*)malloc(sizeof(double) * num_get_row);
    double* get_upper = (double*)malloc(sizeof(double) * num_get_row);
    HighsInt get_num_nz = 0;

    assertIntValuesEqual("getNumRows", Highs_getNumRows(highs), num_row);

    return_status = Highs_getRowsByRange(highs, get_row, get_row,
					 &get_num_row, get_lower, get_upper, &get_num_nz,
					 NULL, NULL, NULL);
    assert( return_status == kHighsStatusOk );

    assertIntValuesEqual("getRows get_num_row", get_num_row, num_get_row);
    assertDoubleValuesEqual("getRows get_lower", get_lower[0], row_lower[get_row]);
    assertDoubleValuesEqual("getRows get_upper", get_upper[0], row_upper[get_row]);
    assertIntValuesEqual("getRows get_num_nz", get_num_nz, 2);

      // could also check coefficients by calling again...

    free(get_upper);
    free(get_lower);
  }

  



  return_status = Highs_setBoolOptionValue(highs, "output_flag", 0);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("Running quietly...\n");
  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );
  if (dev_run)
    printf("Running loudly...\n");

  // Get the model status
  HighsInt model_status = Highs_getModelStatus(highs);

  if (dev_run)
    printf("Run status = %"HIGHSINT_FORMAT"; Model status = %"HIGHSINT_FORMAT"\n", return_status, model_status);

  double objective_function_value;
  return_status = Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  assert( return_status == kHighsStatusOk );
  HighsInt simplex_iteration_count;
  return_status = Highs_getIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);
  assert( return_status == kHighsStatusOk );
  HighsInt primal_solution_status;
  return_status = Highs_getIntInfoValue(highs, "primal_solution_status", &primal_solution_status);
  assert( return_status == kHighsStatusOk );
  HighsInt dual_solution_status;
  return_status = Highs_getIntInfoValue(highs, "dual_solution_status", &dual_solution_status);
  assert( return_status == kHighsStatusOk );

  if (dev_run) {
    printf("Objective value = %g; Iteration count = %"HIGHSINT_FORMAT"\n",
	   objective_function_value, simplex_iteration_count);
    if (model_status == kHighsModelStatusOptimal) {
      // Get the primal and dual solution
      return_status = Highs_getSolution(highs, col_value, col_dual, row_value, row_dual);
      assert( return_status == kHighsStatusOk );
      // Get the basis
      return_status = Highs_getBasis(highs, col_basis_status, row_basis_status);
      assert( return_status == kHighsStatusOk );
      // Report the column primal and dual values, and basis status
      for (HighsInt iCol = 0; iCol < num_col; iCol++)
	printf("Col%"HIGHSINT_FORMAT" = %lf; dual = %lf; status = %"HIGHSINT_FORMAT"; \n",
	       iCol, col_value[iCol], col_dual[iCol], col_basis_status[iCol]);
      // Report the row primal and dual values, and basis status
      for (HighsInt iRow = 0; iRow < num_row; iRow++)
	printf("Row%"HIGHSINT_FORMAT" = %lf; dual = %lf; status = %"HIGHSINT_FORMAT"; \n",
	       iRow, row_value[iRow], row_dual[iRow], row_basis_status[iRow]);
    }
  }
  free(col_value);
  free(col_dual);
  free(row_value);
  free(row_dual);
  free(col_basis_status);
  free(row_basis_status);

  Highs_destroy(highs);

  // Define the constraint matrix to pass to the LP
  HighsInt a_format = kHighsMatrixFormatColwise;
  sense = kHighsObjSenseMinimize;
  double offset = 0;
  HighsInt a_start[2] = {0, 2};
  HighsInt a_index[5] = {1, 2, 0, 1, 2};
  double a_value[5] = {1.0, 2.0, 1.0, 2.0, 1.0};
  highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);
  return_status = Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset,
			   col_cost, col_lower, col_upper,
			   row_lower, row_upper,
			   a_start, a_index, a_value);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );
  model_status = Highs_getModelStatus(highs);
  assert( model_status == kHighsModelStatusOptimal );
  if (dev_run)
    printf("Run status = %"HIGHSINT_FORMAT"; Model status = %"HIGHSINT_FORMAT"\n", return_status, model_status);

  Highs_destroy(highs);
}

void full_api_mip() {
  // The use of the full HiGHS API is illustrated for the MIP
  //
  // Min    f  = -3x_0 - 2x_1 - x_2
  // s.t.          x_0 +  x_1 + x_2 <=  7
  //              4x_0 + 2x_1 + x_2  = 12
  //              x_0 >=0; x_1 >= 0; x_2 binary

  const HighsInt num_col = 3;
  const HighsInt num_row = 2;
  const HighsInt num_nz = 6;
  HighsInt a_format = kHighsMatrixFormatColwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;

  // Define the column costs, lower bounds and upper bounds
  double col_cost[3] = {-3.0, -2.0, -1.0};
  double col_lower[3] = {0.0, 0.0, 0.0};
  double col_upper[3] = {1.0e30, 1.0e30, 1.0};
  // Define the row lower bounds and upper bounds
  double row_lower[2] = {-1.0e30, 12.0};
  double row_upper[2] = { 7.0,    12.0};
  // Define the constraint matrix column-wise
  HighsInt a_start[3] = {0, 2, 4};
  HighsInt a_index[6] = {0, 1, 0, 1, 0, 1};
  double a_value[6] = {1.0, 4.0, 1.0, 2.0, 1.0, 1.0};

  HighsInt integrality[3] = {kHighsVarTypeInteger, kHighsVarTypeInteger, kHighsVarTypeInteger};

  double* col_value = (double*)malloc(sizeof(double) * num_col);
  double* row_value = (double*)malloc(sizeof(double) * num_row);

  HighsInt model_status;
  HighsInt return_status;

  void* highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);
  return_status = Highs_passMip(highs, num_col, num_row, num_nz, a_format,
				sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper,
				a_start, a_index, a_value,
				integrality);
  assert(return_status == kHighsStatusOk);
  Highs_setStringOptionValue(highs, "presolve", "off");
  return_status = Highs_run(highs);
  // mip_node_count is always int64_t, so the following should be an
  // error depending on whether HIGHSINT64 is set
  HighsInt mip_node_count_int;
  HighsInt required_return_status = kHighsStatusError;
#ifdef HIGHSINT64
  required_return_status = kHighsStatusOk;
#endif
  return_status = Highs_getIntInfoValue(highs, "mip_node_count", &mip_node_count_int);
  assert(return_status == required_return_status);
  int64_t mip_node_count;
  return_status = Highs_getInt64InfoValue(highs, "mip_node_count", &mip_node_count);
  assert( return_status == kHighsStatusOk );
  assert( mip_node_count == 1 );

  // Test Highs_getColIntegrality
  HighsInt col_integrality;
  return_status = Highs_getColIntegrality(highs, -1, &col_integrality);
  assert( return_status == kHighsStatusError );
  return_status = Highs_getColIntegrality(highs, num_col, &col_integrality);
  assert( return_status == kHighsStatusError );
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    return_status = Highs_getColIntegrality(highs, iCol, &col_integrality);
    assert( return_status == kHighsStatusOk );
    assert( col_integrality == 1 );
  }

  Highs_destroy(highs);

  free(col_value);
  free(row_value);

}

void full_api_qp() {
  double required_objective_function_value;
  double required_x0;
  double required_x1;
  double objective_function_value;
  HighsInt model_status;
  HighsInt return_status;
  void* highs = Highs_create();
  const double inf = Highs_getInfinity(highs);
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);

  // Oscar's edge case
  //
  // min x^2 + x = x(x + 1)

  HighsInt num_col = 0;
  return_status = Highs_addCol(highs, 1.0, -inf, inf, 0, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  num_col++;

  double offset = 0.25;
  return_status = Highs_changeObjectiveOffset(highs, offset);
  assert( return_status == kHighsStatusOk );

  HighsInt q_dim = 1;
  HighsInt q_num_nz = 1;
  HighsInt q_format = kHighsHessianFormatTriangular;
  HighsInt* q_start = (HighsInt*)malloc(sizeof(HighsInt) * q_dim);
  HighsInt* q_index = (HighsInt*)malloc(sizeof(HighsInt) * q_num_nz);
  double* q_value = (double*)malloc(sizeof(double) * q_num_nz);
  q_start[0] = 0;
  q_index[0] = 0;
  q_value[0] = 2.0;
  return_status = Highs_passHessian(highs, q_dim, q_num_nz, q_format, q_start, q_index, q_value);
  assert( return_status == kHighsStatusOk );
  if (dev_run) Highs_writeModel(highs, "");
  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );
  model_status = Highs_getModelStatus(highs);
  assertIntValuesEqual("Model status for 1-d QP", model_status, kHighsModelStatusOptimal);

  required_objective_function_value = 0;
  required_x0 = -0.5;
  objective_function_value = Highs_getObjectiveValue(highs);
  assertDoubleValuesEqual("Objective", objective_function_value, required_objective_function_value);

  double* col_solution = (double*)malloc(sizeof(double) * num_col);

  return_status = Highs_getSolution(highs, col_solution, NULL, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  assertDoubleValuesEqual("x0", col_solution[0], required_x0);

  if (dev_run) Highs_writeSolutionPretty(highs, "");
  // Add a variable x1 with objective x1^2 - x1
  //
  // Add the variable
  return_status = Highs_addCol(highs, -1.0, -inf, inf, 0, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  num_col++;
  // Cannot solve the model until the Hessian has been replaced
  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusError );
  assertIntValuesEqual("Run status for 2-d QP with illegal Hessian", return_status, -1);

  model_status = Highs_getModelStatus(highs);
  assertIntValuesEqual("Model status for 2-d QP with illegal Hessian", model_status, 2);

  free(q_start);
  free(q_index);
  free(q_value);

  // Pass the new Hessian
  q_dim = 2;
  q_num_nz = 2;
  q_start = (HighsInt*)malloc(sizeof(HighsInt) * q_dim);
  q_index = (HighsInt*)malloc(sizeof(HighsInt) * q_num_nz);
  q_value = (double*)malloc(sizeof(double) * q_num_nz);
  q_start[0] = 0;
  q_index[0] = 0;
  q_value[0] = 2.0;
  q_start[1] = 1;
  q_index[1] = 1;
  q_value[1] = 2.0;
  return_status = Highs_passHessian(highs, q_dim, q_num_nz, q_format, q_start, q_index, q_value);
  assert( return_status == kHighsStatusOk );
  if (dev_run) Highs_writeModel(highs, "");

  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );

  model_status = Highs_getModelStatus(highs);
  assertIntValuesEqual("Model status for 2-d QP", model_status, kHighsModelStatusOptimal);
  
  required_objective_function_value = -0.25;
  required_x1 = 0.5;
  objective_function_value = Highs_getObjectiveValue(highs);
  assertDoubleValuesEqual("Objective", objective_function_value, required_objective_function_value);

  free(col_solution);
  col_solution = (double*)malloc(sizeof(double) * num_col);

  return_status = Highs_getSolution(highs, col_solution, NULL, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  assertDoubleValuesEqual("x0", col_solution[0], required_x0);
  assertDoubleValuesEqual("x1", col_solution[1], required_x1);

  // Illustrate methods for getting and changing the offset by getting
  // the current offset, shifting it by the current objective and
  // checking that the objective value is changed accordingly

  double check_offset;
  return_status = Highs_getObjectiveOffset(highs, &check_offset);
  assert( return_status == kHighsStatusOk );
  assertDoubleValuesEqual("Offset", check_offset, offset);

  double dl_offset = -objective_function_value;
  offset += dl_offset;

  return_status = Highs_changeObjectiveOffset(highs, offset);
  assert( return_status == kHighsStatusOk );
  required_objective_function_value += dl_offset;
  objective_function_value = Highs_getObjectiveValue(highs);
  assertDoubleValuesEqual("Objective with new offset", objective_function_value, required_objective_function_value);

  // Add the constraint 0.5 <= x0 + x1
  HighsInt a_index[2] = {0, 1};
  double a_value[2] = {1, 1};
  return_status = Highs_addRow(highs, 0.5, inf, 2, a_index, a_value);
  assert( return_status == kHighsStatusOk );
  if (dev_run) Highs_writeModel(highs, "");

  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );
  assertIntValuesEqual("Run status for 2-d QP with constraint", return_status, kHighsStatusOk);
  
  model_status = Highs_getModelStatus(highs);
  assertIntValuesEqual("Model status for 2-d QP with constraint", model_status, kHighsModelStatusOptimal);

  required_objective_function_value = 0.125;
  required_x0 = -0.25;
  required_x1 = 0.75;

  objective_function_value = Highs_getObjectiveValue(highs);
  assertDoubleValuesEqual("Objective", objective_function_value, required_objective_function_value);

  return_status = Highs_getSolution(highs, col_solution, NULL, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  assertDoubleValuesEqual("x0", col_solution[0], required_x0);
  assertDoubleValuesEqual("x1", col_solution[1], required_x1);

  // Add bounds to make the QP infeasible
  return_status = Highs_changeColBounds(highs, 0, -inf, 0);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_changeColBounds(highs, 1, -inf, 0);
  assert( return_status == kHighsStatusOk );
 
  if (dev_run) Highs_writeModel(highs, "");

  return_status = Highs_run(highs);
  assert( return_status == kHighsStatusOk );
  assertIntValuesEqual("Run status for infeasible 2-d QP", return_status, 0);
  
  model_status = Highs_getModelStatus(highs);
  assertIntValuesEqual("Model status for infeasible 2-d QP", model_status, 8);
  assert( model_status == kHighsModelStatusInfeasible );

  Highs_destroy(highs);

  free(q_start);
  free(q_index);
  free(q_value);
  free(col_solution);

}

void pass_presolve_get_lp() {
  // Form and solve the LP
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  void* highs;

  highs = Highs_create();
  const double kHighsInf = Highs_getInfinity(highs);
  HighsInt model_status;
  HighsInt return_status;
  
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);
  HighsInt a_format = kHighsMatrixFormatColwise;
  HighsInt sense = kHighsObjSenseMinimize;
  double offset = 0;
  // Define the column costs, lower bounds and upper bounds

  const HighsInt num_col = 2;
  const HighsInt num_row = 3;
  const HighsInt num_nz = 5;

  double col_cost[2] = {2.0, 3.0};
  double col_lower[2] = {0.0, 1.0};
  double col_upper[2] = {3.0, kHighsInf};
  // Define the row lower bounds and upper bounds
  double row_lower[3] = {-kHighsInf, 10.0, 8.0};
  double row_upper[3] = {6.0, 14.0, kHighsInf};
  HighsInt a_start[2] = {0, 2};
  HighsInt a_index[5] = {1, 2, 0, 1, 2};
  double a_value[5] = {1.0, 2.0, 1.0, 2.0, 1.0};

  return_status = Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset,
			       col_cost, col_lower, col_upper,
			       row_lower, row_upper,
			       a_start, a_index, a_value);
  assert( return_status == kHighsStatusOk );

  return_status = Highs_presolve(highs);
  assert( return_status == kHighsStatusOk );
  for (HighsInt k = 0; k < 2; k++) {
    // Loop twice: once for col-wise; once for row-wise
    HighsInt presolved_num_col = Highs_getPresolvedNumCol(highs);
    HighsInt presolved_num_row = Highs_getPresolvedNumRow(highs);
    HighsInt presolved_num_nz = Highs_getPresolvedNumNz(highs);
    HighsInt presolved_a_format = k == 0 ? kHighsMatrixFormatColwise : kHighsMatrixFormatRowwise;
    HighsInt presolved_sense;
    double presolved_offset;
    double* presolved_col_cost = (double*)malloc(sizeof(double) * presolved_num_col);
    double* presolved_col_lower = (double*)malloc(sizeof(double) * presolved_num_col);
    double* presolved_col_upper = (double*)malloc(sizeof(double) * presolved_num_col);
    double* presolved_row_lower = (double*)malloc(sizeof(double) * presolved_num_row);
    double* presolved_row_upper = (double*)malloc(sizeof(double) * presolved_num_row);
    HighsInt* presolved_a_start = (HighsInt*)malloc(sizeof(HighsInt) * (presolved_num_col+1));
    HighsInt* presolved_a_index = (HighsInt*)malloc(sizeof(HighsInt) * presolved_num_nz);
    double* presolved_a_value = (double*)malloc(sizeof(double) * presolved_num_nz);
  
    return_status = Highs_getPresolvedLp(highs, presolved_a_format,
					 &presolved_num_col, &presolved_num_row, &presolved_num_nz,
					 &presolved_sense, &presolved_offset,
					 presolved_col_cost, presolved_col_lower, presolved_col_upper,
					 presolved_row_lower, presolved_row_upper,
					 presolved_a_start, presolved_a_index, presolved_a_value, NULL);
    assert( return_status == kHighsStatusOk );
    // Solve the presolved LP within a local version of HiGHS
    void* local_highs;
    local_highs = Highs_create();
    Highs_setBoolOptionValue(local_highs, "output_flag", dev_run);
    Highs_setStringOptionValue(local_highs, "presolve", "off");
    return_status = Highs_passLp(local_highs,
				 presolved_num_col, presolved_num_row, presolved_num_nz,
				 presolved_a_format, presolved_sense, presolved_offset,
				 presolved_col_cost, presolved_col_lower, presolved_col_upper,
				 presolved_row_lower, presolved_row_upper,
				 presolved_a_start, presolved_a_index, presolved_a_value);
    assert( return_status == kHighsStatusOk );
    return_status = Highs_run(local_highs);
    
    double* col_value = (double*)malloc(sizeof(double) * num_col);
    double* col_dual = (double*)malloc(sizeof(double) * num_col);
    double* row_dual = (double*)malloc(sizeof(double) * num_row);

    return_status = Highs_getSolution(local_highs, col_value, col_dual, NULL, row_dual);
    assert( return_status == kHighsStatusOk );

    return_status = Highs_postsolve(highs, col_value, col_dual, row_dual);
    assert( return_status == kHighsStatusOk );

    model_status = Highs_getModelStatus(highs);
    assert( model_status == kHighsModelStatusOptimal );
    
    // With just the primal solution, optimality cannot be determined

    return_status = Highs_postsolve(highs, col_value, NULL, NULL);
    assert( return_status == kHighsStatusWarning );

    model_status = Highs_getModelStatus(highs);
    assert( model_status == kHighsModelStatusUnknown );

    free(presolved_col_cost);
    free(presolved_col_lower);
    free(presolved_col_upper);
    free(presolved_row_lower);
    free(presolved_row_upper);
    free(presolved_a_start);
    free(presolved_a_index);
    free(presolved_a_value);
    free(col_value);
    free(col_dual);
    free(row_dual);

    
  }
}

void options() {
  void* highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);

  HighsInt simplex_scale_strategy;
  HighsInt return_status;
  return_status = Highs_setIntOptionValue(highs, "simplex_scale_strategy", 0);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_getIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
  assert( return_status == kHighsStatusOk );
  assert( simplex_scale_strategy == 0 );

  double primal_feasibility_tolerance;
  return_status = Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", 2.0);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_getDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  assert( return_status == kHighsStatusOk );
  assert( primal_feasibility_tolerance == 2.0 );

  Highs_destroy(highs);

}

void test_getColsByRange() {
  void* highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);
  HighsInt return_status;
  return_status = Highs_addCol(highs, -1.0, 0.0, 1.0, 0, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  return_status = Highs_addCol(highs, -1.0, 0.0, 1.0, 0, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  HighsInt a_index[2] = {0, 1};
  double a_value[2] = {1.0, -1.0};
  return_status = Highs_addRow(highs, 0.0, 0.0, 2, a_index, a_value);
  assert( return_status == kHighsStatusOk );
  HighsInt num_cols;
  HighsInt num_nz;
  HighsInt matrix_start[2] = {-1, -1};
  return_status = Highs_getColsByRange(highs, 0, 1, &num_cols, NULL, NULL, NULL, &num_nz,
		       matrix_start, NULL, NULL);
  assert( return_status == kHighsStatusOk );
  assert( num_cols == 2 );
  assert( num_nz == 2 );
  assert( matrix_start[0] == 0 );
  assert( matrix_start[1] == 1 );
  HighsInt matrix_index[2] = {-1, -1};
  double matrix_value[2] = {0.0, 0.0};
  return_status = Highs_getColsByRange(highs, 0, 1, &num_cols, NULL, NULL, NULL, &num_nz,
		       matrix_start, matrix_index, matrix_value);
  assert( return_status == kHighsStatusOk );
  assert( matrix_index[0] == 0 );
  assert( matrix_index[1] == 0 );
  assert( matrix_value[0] == 1.0 );
  assert( matrix_value[1] == -1.0 );

  Highs_destroy(highs);
}

void test_passHessian() {
  void* highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);
  Highs_addCol(highs, 2.0, 0.0, 2.0, 0, NULL, NULL);
  Highs_changeObjectiveSense(highs, kHighsObjSenseMaximize);
  HighsInt start[1] = {0};
  HighsInt index[1] = {0};
  double value[1] = {-2.0};
  HighsInt return_status;
  return_status = Highs_passHessian(highs, 1, 1, 1, start, index, value);
  assertIntValuesEqual("Return of passHessian", return_status, kHighsStatusOk);
  Highs_run(highs);
  // Solving max -x^2 + 2x
  const double optimal_objective_value = 1;
  const double primal = 1;
  const double dual = 0;
  assertIntValuesEqual("Status", Highs_getModelStatus(highs), kHighsModelStatusOptimal);  // kOptimal
  double col_value[1] = {-123.0};
  double col_dual[1] = {0.0};
  Highs_getSolution(highs, col_value, col_dual, NULL, NULL);
  double objective_value = Highs_getObjectiveValue(highs);
  assertDoubleValuesEqual("Objective", objective_value, optimal_objective_value);
  assertDoubleValuesEqual("Primal", col_value[0], primal);
  assertDoubleValuesEqual("Dual", col_dual[0], dual);

  Highs_destroy(highs);
}

void test_ranging() {

  void* highs = Highs_create();
  if (!dev_run) Highs_setBoolOptionValue(highs, "output_flag", 0);
  //
  // Set up
  //        min y
  //        s.t.
  //        -x + y >= 2
  //        x + y >= 0
  //
  double inf = Highs_getInfinity(highs);
  Highs_addVar(highs, -inf, inf);
  Highs_addVar(highs, -inf, inf);
  Highs_changeColCost(highs, 0, 0);
  Highs_changeColCost(highs, 1, 1);
  HighsInt index[2] = {0.0, 1.0};
  double value[2] = {-1, 1};
  Highs_addRow(highs, 2, inf, 2, index, value);
  value[0] = 1.0;
  Highs_addRow(highs, 0, inf, 2, index, value);
  // Cost ranging
  // c0 2 -1 1 0
  // c1 0 0 inf inf
  // 
  // Bound ranging
  // Columns
  // c0 1 -inf inf 1
  // c1 1 1 inf 1
  // Rows
  // r0 -inf -inf inf inf
  // r1 -inf -inf inf inf
  Highs_run(highs);
  HighsInt num_col = Highs_getNumCol(highs);
  HighsInt num_row = Highs_getNumRow(highs);
  double* col_cost_up_value = (double*)malloc(sizeof(double) * num_col);
  double* col_cost_up_objective = (double*)malloc(sizeof(double) * num_col);
  HighsInt* col_cost_up_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* col_cost_up_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  double* col_cost_dn_value = (double*)malloc(sizeof(double) * num_col);
  double* col_cost_dn_objective = (double*)malloc(sizeof(double) * num_col);
  HighsInt* col_cost_dn_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* col_cost_dn_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  double* col_bound_up_value = (double*)malloc(sizeof(double) * num_col);
  double* col_bound_up_objective = (double*)malloc(sizeof(double) * num_col);
  HighsInt* col_bound_up_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* col_bound_up_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  double* col_bound_dn_value = (double*)malloc(sizeof(double) * num_col);
  double* col_bound_dn_objective = (double*)malloc(sizeof(double) * num_col);
  HighsInt* col_bound_dn_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  HighsInt* col_bound_dn_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_col);
  double* row_bound_up_value = (double*)malloc(sizeof(double) * num_row);
  double* row_bound_up_objective = (double*)malloc(sizeof(double) * num_row);
  HighsInt* row_bound_up_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_row);
  HighsInt* row_bound_up_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_row);
  double* row_bound_dn_value = (double*)malloc(sizeof(double) * num_row);
  double* row_bound_dn_objective = (double*)malloc(sizeof(double) * num_row);
  HighsInt* row_bound_dn_in_var = (HighsInt*)malloc(sizeof(HighsInt) * num_row);
  HighsInt* row_bound_dn_ou_var = (HighsInt*)malloc(sizeof(HighsInt) * num_row);
  HighsInt status = 
  Highs_getRanging(highs,
		   //
		   col_cost_up_value, col_cost_up_objective, col_cost_up_in_var, col_cost_up_ou_var, 
		   col_cost_dn_value, col_cost_dn_objective, col_cost_dn_in_var, col_cost_dn_ou_var, 
		   col_bound_up_value, col_bound_up_objective, col_bound_up_in_var, col_bound_up_ou_var, 
		   col_bound_dn_value, col_bound_dn_objective, col_bound_dn_in_var, col_bound_dn_ou_var, 
		   row_bound_up_value, row_bound_up_objective, row_bound_up_in_var, row_bound_up_ou_var, 
		   row_bound_dn_value, row_bound_dn_objective, row_bound_dn_in_var, row_bound_dn_ou_var);
  assert(status == kHighsStatusOk);

  assertDoubleValuesEqual("col_cost_dn_objective[0]", col_cost_dn_objective[0], 2);
  assertDoubleValuesEqual("col_cost_dn_value[0]", col_cost_dn_value[0], -1);
  assertDoubleValuesEqual("col_cost_up_value[0]", col_cost_up_value[0], 1);
  assertDoubleValuesEqual("col_cost_up_objective[0]", col_cost_up_objective[0], 0);
  assertDoubleValuesEqual("col_cost_dn_objective[1]", col_cost_dn_objective[1], 0);
  assertDoubleValuesEqual("col_cost_dn_value[1]", col_cost_dn_value[1], 0);
  assertDoubleValuesEqual("col_cost_up_value[1]", col_cost_up_value[1], inf);
  assertDoubleValuesEqual("col_cost_up_objective[1]", col_cost_up_objective[1], inf);

  assertDoubleValuesEqual("col_bound_dn_objective[0]", col_bound_dn_objective[0], 1);
  assertDoubleValuesEqual("col_bound_dn_value[0]", col_bound_dn_value[0], -inf);
  assertDoubleValuesEqual("col_bound_up_value[0]", col_bound_up_value[0], inf);
  assertDoubleValuesEqual("col_bound_up_objective[0]", col_bound_up_objective[0], 1);
  assertDoubleValuesEqual("col_bound_dn_objective[1]", col_bound_dn_objective[1], 1);
  assertDoubleValuesEqual("col_bound_dn_value[1]", col_bound_dn_value[1], 1);
  assertDoubleValuesEqual("col_bound_up_value[1]", col_bound_up_value[1], inf);
  assertDoubleValuesEqual("col_bound_up_objective[1]", col_bound_up_objective[1], 1);

  assertDoubleValuesEqual("row_bound_dn_objective[0]", row_bound_dn_objective[0], -inf);
  assertDoubleValuesEqual("row_bound_dn_value[0]", row_bound_dn_value[0], -inf);
  assertDoubleValuesEqual("row_bound_up_value[0]", row_bound_up_value[0], inf);
  assertDoubleValuesEqual("row_bound_up_objective[0]", row_bound_up_objective[0], inf);
  assertDoubleValuesEqual("row_bound_dn_objective[1]", row_bound_dn_objective[1], -inf);
  assertDoubleValuesEqual("row_bound_dn_value[1]", row_bound_dn_value[1], -inf);
  assertDoubleValuesEqual("row_bound_up_value[1]", row_bound_up_value[1], inf);
  assertDoubleValuesEqual("row_bound_up_objective[1]", row_bound_up_objective[1], inf);

  free(col_cost_up_value);
  free(col_cost_up_objective);
  free(col_cost_up_in_var);
  free(col_cost_up_ou_var);
  free(col_cost_dn_value);
  free(col_cost_dn_objective);
  free(col_cost_dn_in_var);
  free(col_cost_dn_ou_var);
  free(col_bound_up_value);
  free(col_bound_up_objective);
  free(col_bound_up_in_var);
  free(col_bound_up_ou_var);
  free(col_bound_dn_value);
  free(col_bound_dn_objective);
  free(col_bound_dn_in_var);
  free(col_bound_dn_ou_var);
  free(row_bound_up_value);
  free(row_bound_up_objective);
  free(row_bound_up_in_var);
  free(row_bound_up_ou_var);
  free(row_bound_dn_value);
  free(row_bound_dn_objective);
  free(row_bound_dn_in_var);
  free(row_bound_dn_ou_var);

  Highs_destroy(highs);
}

void test_callback() {
  HighsInt num_col = 7;
  HighsInt num_row = 1;
  HighsInt num_nz = num_col;
  HighsInt a_format = kHighsMatrixFormatRowwise;
  HighsInt sense = kHighsObjSenseMaximize;
  double offset = 0;
  double col_cost[7] = {8, 1, 7, 2, 1, 2, 1};
  double col_lower[7] = {0, 0, 0, 0, 0, 0, 0};
  double col_upper[7] = {1, 1, 1, 1, 1, 1, 1};
  double row_lower[1] = {0};
  double row_upper[1] = {28};
  HighsInt a_start[2] = {0, 7};
  HighsInt a_index[7] = {0, 1, 2, 3, 4, 5, 6};
  double a_value[7] = {9, 6, 7, 9, 7, 9, 9};
  HighsInt integrality[7] = {kHighsVarTypeInteger, kHighsVarTypeInteger,
			     kHighsVarTypeInteger, kHighsVarTypeInteger,
			     kHighsVarTypeInteger, kHighsVarTypeInteger,
			     kHighsVarTypeInteger};

  void* highs;
  highs = Highs_create();
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);
  Highs_passMip(highs, num_col, num_row, num_nz, a_format, sense, offset,
		col_cost, col_lower, col_upper,
		row_lower, row_upper,
		a_start, a_index, a_value,
		integrality);
  
  Highs_setCallback(highs, userCallback, NULL);
  Highs_startCallback(highs, kHighsCallbackLogging);
  Highs_startCallback(highs, kHighsCallbackMipInterrupt);
  Highs_run(highs);
  double objective_function_value;
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  double inf = Highs_getInfinity(highs);
  assertDoubleValuesEqual("objective_function_value", objective_function_value, inf);
  Highs_stopCallback(highs, kHighsCallbackMipInterrupt);
  Highs_run(highs);
  Highs_getDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  assertDoubleValuesEqual("objective_function_value", objective_function_value, 17);

  double user_callback_data = inf;
  void* p_user_callback_data = (void*)(&user_callback_data);
  
  Highs_setCallback(highs, userCallback, p_user_callback_data);
  Highs_clearSolver(highs);
  Highs_startCallback(highs, kHighsCallbackMipImprovingSolution);
  Highs_run(highs);
  
  Highs_destroy(highs);
}

void test_getModel() {
  void* highs;
  highs = Highs_create();
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);
  const double inf = Highs_getInfinity(highs);
  
  HighsInt num_col = 2;
  HighsInt num_row = 2;
  HighsInt num_nz = 4;
  HighsInt sense = -1;
  double offset;
  double col_cost[2] = {8, 10};
  double col_lower[2] = {0, 0};
  double col_upper[2] = {inf, inf};
  double row_lower[2] = {-inf, -inf};
  double row_upper[2] = {120, 210};
  HighsInt a_index[4] = {0, 1, 0, 1};
  double a_value[4] = {0.3, 0.5, 0.7, 0.5};
  HighsInt a_start[2] = {0, 2};
  Highs_addVars(highs, num_col, col_lower, col_upper);
  Highs_changeColsCostByRange(highs, 0, num_col-1, col_cost);
  Highs_addRows(highs, num_row, row_lower, row_upper, num_nz, a_start, a_index, a_value);
  Highs_changeObjectiveSense(highs, sense);
  Highs_run(highs);

  HighsInt ck_num_col;
  HighsInt ck_num_row;
  HighsInt ck_num_nz;
  HighsInt ck_sense;
  double ck_offset;

  // Get the model dimensions by passing array pointers as NULL
  Highs_getLp(highs, kHighsMatrixFormatRowwise,
	      &ck_num_col, &ck_num_row, &ck_num_nz,		  
	      &ck_sense, &ck_offset, NULL,
	      NULL, NULL, NULL,
	      NULL, NULL, NULL,
	      NULL, NULL);

  assert( ck_num_col == num_col );
  assert( ck_num_row == num_row );
  assert( ck_num_nz == num_nz );
  // Motivated by #1712, ensure that the correct sense is returned when maximizing
  assert( ck_sense == sense );

  double* ck_col_cost = (double*)malloc(sizeof(double) * ck_num_col);;
  double* ck_col_lower = (double*)malloc(sizeof(double) * ck_num_col);
  double* ck_col_upper = (double*)malloc(sizeof(double) * ck_num_col);
  double* ck_row_lower = (double*)malloc(sizeof(double) * ck_num_row);
  double* ck_row_upper = (double*)malloc(sizeof(double) * ck_num_row);
  HighsInt* ck_a_start = (HighsInt*)malloc(sizeof(HighsInt) * ck_num_col);
  HighsInt* ck_a_index = (HighsInt*)malloc(sizeof(HighsInt) * ck_num_nz);
  double* ck_a_value = (double*)malloc(sizeof(double) * num_nz);
  
  // Get the arrays
  Highs_getLp(highs, kHighsMatrixFormatRowwise,
	      &ck_num_col, &ck_num_row, &ck_num_nz,
	      &ck_sense, &ck_offset, ck_col_cost,
	      ck_col_lower, ck_col_upper, ck_row_lower,
	      ck_row_upper, ck_a_start, ck_a_index,
	      ck_a_value, NULL);

  assert( doubleArraysEqual(num_col, ck_col_cost, col_cost) );
  assert( doubleArraysEqual(num_col, ck_col_lower, col_lower) );
  assert( doubleArraysEqual(num_col, ck_col_upper, col_upper) );
  assert( doubleArraysEqual(num_row, ck_row_lower, row_lower) );
  assert( doubleArraysEqual(num_row, ck_row_upper, row_upper) );
  assert( highsIntArraysEqual(num_col, ck_a_start, a_start) );
  assert( highsIntArraysEqual(num_nz, ck_a_index, a_index) );
  assert( doubleArraysEqual(num_nz, ck_a_value, a_value) );

}

/*
The horrible C in this causes problems in some of the CI tests,
so suppress thius test until the C has been improved

void test_setSolution() {
  void* highs = Highs_create();
  // Perform in C the equivalent of std::string model_file =
  // std::string(HIGHS_DIR) + "/check/instances/shell.mps";

  char* dir = HIGHS_DIR;
  char model_file0[100];
  strcat(model_file0, dir);
  strcat(model_file0, "/check/instances/shell.mps");
  strcat(model_file0, "\0");
  char* substr = model_file0 + 1;
  memmove(model_file0, substr, strlen(substr) + 1);
  HighsInt length = strlen(model_file0) + 1;
  char model_file[length];
  strcpy(model_file, model_file0);
  
  if (dev_run) printf("\nSolving from scratch\n");
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);

  Highs_readModel(highs, model_file);
  Highs_run(highs);
 HighsInt iteration_count0;
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &iteration_count0);
 HighsInt num_col = Highs_getNumCol(highs);
  double* col_value = (double*)malloc(sizeof(double) * num_col);
  Highs_getSolution(highs, col_value, NULL, NULL, NULL);
  Highs_clear(highs);
  if (dev_run) printf("\nSolving from saved solution\n");
  Highs_setBoolOptionValue(highs, "output_flag", dev_run);
  Highs_readModel(highs, model_file);
  Highs_setSolution(highs, col_value, NULL, NULL, NULL);
  Highs_run(highs);
  HighsInt iteration_count1;
  Highs_getIntInfoValue(highs, "simplex_iteration_count", &iteration_count1);
  HighsInt logic = iteration_count0 > iteration_count1;
  printf("Iteration counts are %d and %d\n", iteration_count0, iteration_count1);
  assertLogical("Dual", logic);
  
  Highs_destroy(highs);
}
*/
int main() {
  minimal_api_illegal_lp();
  test_callback();
  version_api();
  full_api();
  minimal_api_lp();
  minimal_api_mip();
  minimal_api_qp();
  full_api_options();
  full_api_lp();
  full_api_mip();
  full_api_qp();
  pass_presolve_get_lp();
  options();
  test_getColsByRange();
  test_passHessian();
  test_ranging();
  test_getModel();
  //  test_setSolution();
  return 0;
}
