/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/CupdlpWrapper.cpp
 * @brief
 * @author Julian Hall
 */
#include "pdlp/CupdlpWrapper.h"

typedef enum CONSTRAINT_TYPE { EQ = 0, LEQ, GEQ, BOUND } constraint_type;

void reportParams(CUPDLPwork *w,
		  cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
		  cupdlp_bool *ifChangeFloatParam,
		  cupdlp_float *floatParam);

void getUserParamsFromOptions(const HighsOptions& options,
			      cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
			      cupdlp_bool *ifChangeFloatParam,
			      cupdlp_float *floatParam);

HighsStatus solveLpCupdlp(HighsLpSolverObject& solver_object) {
  return solveLpCupdlp(solver_object.options_, solver_object.timer_, solver_object.lp_, 
		       solver_object.basis_, solver_object.solution_, 
		       solver_object.model_status_, solver_object.highs_info_,
		       solver_object.callback_);
}


HighsStatus solveLpCupdlp(const HighsOptions& options,
			  HighsTimer& timer,
			  const HighsLp& lp, 
			  HighsBasis& highs_basis,
			  HighsSolution& highs_solution,
			  HighsModelStatus& model_status,
			  HighsInfo& highs_info,
			  HighsCallback& callback) {
  // Indicate that there is no valid primal solution, dual solution or basis
  highs_basis.valid = false;
  highs_solution.value_valid = false;
  highs_solution.dual_valid = false;
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  char *fout = nullptr;

  int nCols;
  int nRows;
  int nEqs;
  int nCols_origin;
  cupdlp_bool ifSaveSol = false;
  cupdlp_bool ifPresolve = false;

  int nnz = 0;
  double *rhs = NULL;
  double *cost = NULL;

  cupdlp_float *lower = NULL;
  cupdlp_float *upper = NULL;

  // -------------------------
  int *csc_beg = NULL, *csc_idx = NULL;
  double *csc_val = NULL;
  double offset =
      0.0;  // true objVal = sig * c'x - offset, sig = 1 (min) or -1 (max)
  double sign_origin = 1;  // 1 (min) or -1 (max)
  int *constraint_new_idx = NULL;
  cupdlp_float *x_origin = cupdlp_NULL;
  cupdlp_float *y_origin = cupdlp_NULL;

  void *model = NULL;
  void *presolvedmodel = NULL;
  void *model2solve = NULL;

  CUPDLPscaling *scaling =
      (CUPDLPscaling *)cupdlp_malloc(sizeof(CUPDLPscaling));

  // claim solvers variables
  // prepare pointers
  CUPDLP_MATRIX_FORMAT src_matrix_format = CSC;
  CUPDLP_MATRIX_FORMAT dst_matrix_format = CSR_CSC;
  CUPDLPcsc *csc_cpu = cupdlp_NULL;
  CUPDLPproblem *prob = cupdlp_NULL;

  // load parameters

  // set solver parameters
  cupdlp_bool ifChangeIntParam[N_INT_USER_PARAM] = {false};
  cupdlp_int intParam[N_INT_USER_PARAM] = {0};
  cupdlp_bool ifChangeFloatParam[N_FLOAT_USER_PARAM] = {false};
  cupdlp_float floatParam[N_FLOAT_USER_PARAM] = {0.0};

  // Transfer from options
  getUserParamsFromOptions(options,
			   ifChangeIntParam, intParam,
			   ifChangeFloatParam, floatParam);

  formulateLP_highs(lp, &cost, &nCols, &nRows, &nnz, &nEqs,
		    &csc_beg, &csc_idx, &csc_val, &rhs, &lower,
		    &upper, &offset, &sign_origin, &nCols_origin,
		    &constraint_new_idx);
 

  Init_Scaling(scaling, nCols, nRows, cost, rhs);
  cupdlp_int ifScaling = 1;

  CUPDLPwork *w = cupdlp_NULL;
  cupdlp_init_work(w, 1);

  problem_create(&prob);

  // currently, only supprot that input matrix is CSC, and store both CSC and
  // CSR
  csc_create(&csc_cpu);
  csc_cpu->nRows = nRows;
  csc_cpu->nCols = nCols;
  csc_cpu->nMatElem = nnz;
  csc_cpu->colMatBeg = (int *)malloc((1 + nCols) * sizeof(int));
  csc_cpu->colMatIdx = (int *)malloc(nnz * sizeof(int));
  csc_cpu->colMatElem = (double *)malloc(nnz * sizeof(double));
  memcpy(csc_cpu->colMatBeg, csc_beg, (nCols + 1) * sizeof(int));
  memcpy(csc_cpu->colMatIdx, csc_idx, nnz * sizeof(int));
  memcpy(csc_cpu->colMatElem, csc_val, nnz * sizeof(double));

  cupdlp_float scaling_time = getTimeStamp();
  PDHG_Scale_Data_cuda(csc_cpu, ifScaling, scaling, cost, lower,
                                   upper, rhs);
  scaling_time = getTimeStamp() - scaling_time;

  cupdlp_float alloc_matrix_time = 0.0;
  cupdlp_float copy_vec_time = 0.0;

  problem_alloc(prob, nRows, nCols, nEqs, cost, offset, sign_origin,
		csc_cpu, src_matrix_format, dst_matrix_format, rhs,
		lower, upper, &alloc_matrix_time, &copy_vec_time);

  w->problem = prob;
  w->scaling = scaling;
  PDHG_Alloc(w);
  w->timers->dScalingTime = scaling_time;
  w->timers->dPresolveTime = 0;//presolve_time;
  cupdlp_copy_vec(w->rowScale, scaling->rowScale, cupdlp_float, nRows);
  cupdlp_copy_vec(w->colScale, scaling->colScale, cupdlp_float, nCols);

  cupdlp_printf("--------------------------------------------------\n");
  cupdlp_printf("enter main solve loop\n");
  cupdlp_printf("--------------------------------------------------\n");
  // CUPDLP_CALL(LP_SolvePDHG(prob, cupdlp_NULL, cupdlp_NULL, cupdlp_NULL,
  // cupdlp_NULL));
  //   CUPDLP_CALL(LP_SolvePDHG(prob, ifChangeIntParam, intParam,
  //                               ifChangeFloatParam, floatParam, fout));

  cupdlp_init_double(x_origin, nCols_origin);
  cupdlp_init_double(y_origin, nRows);
  LP_SolvePDHG(w, ifChangeIntParam, intParam, ifChangeFloatParam,
	       floatParam, fout, x_origin, nCols_origin, y_origin,
	       ifSaveSol, constraint_new_idx);
 
  HighsStatus return_status =
    pdlpSolutionToHighsSolution(x_origin, nCols_origin,
				y_origin, nRows,
				options, lp, highs_solution);
  // Set the status to optimal until other statuses can be identified
  model_status = HighsModelStatus::kOptimal;
  return return_status;
}

int formulateLP_highs(const HighsLp& lp,
		      double **cost, int *nCols,
		      int *nRows, int *nnz, int *nEqs, int **csc_beg,
		      int **csc_idx, double **csc_val, double **rhs,
		      double **lower, double **upper, double *offset,
		      double *sign_origin, int *nCols_origin,
		      int **constraint_new_idx) {
  int retcode = 0;

  // problem size for malloc
  int nCols_clp = lp.num_col_;
  int nRows_clp = lp.num_row_;
  int nnz_clp = lp.a_matrix_.start_[lp.num_col_];
  *nCols_origin = nCols_clp;
  *nRows = nRows_clp;    // need not recalculate
  *nCols = nCols_clp;    // need recalculate
  *nEqs = 0;             // need recalculate
  *nnz = nnz_clp;        // need recalculate
  *offset = lp.offset_;  // need not recalculate
  if (lp.sense_ == ObjSense::kMinimize) {
    *sign_origin = 1.0;
    printf("Minimize\n");
  } else if (lp.sense_ == ObjSense::kMaximize) {
    *sign_origin = -1.0;
    printf("Maximize\n");
  }
  if (*offset != 0.0) {
    printf("Has obj offset %f\n", *offset);
  } else {
    printf("No obj offset\n");
  }
  // allocate buffer memory
  //  constraint_type *constraint_type_clp = NULL;  // the ONLY one need to free
  // int *constraint_original_idx = NULL;  // pass by user is better, for
  // postsolve recovering dual

  const double *lhs_clp = lp.row_lower_.data();
  const double *rhs_clp = lp.row_upper_.data();
  const HighsInt *A_csc_beg = lp.a_matrix_.start_.data();
  const HighsInt *A_csc_idx = lp.a_matrix_.index_.data();
  const double *A_csc_val = lp.a_matrix_.value_.data();
  int has_lower, has_upper;

  std::vector<constraint_type> constraint_type_clp(nRows_clp);

  cupdlp_init_int(*constraint_new_idx, *nRows);

  // recalculate nRows and nnz for Ax - z = 0
  for (int i = 0; i < nRows_clp; i++) {
    has_lower = lhs_clp[i] > -1e20;
    has_upper = rhs_clp[i] < 1e20;

    // count number of equations and rows
    if (has_lower && has_upper && lhs_clp[i] == rhs_clp[i]) {
      constraint_type_clp[i] = EQ;
      (*nEqs)++;
    } else if (has_lower && !has_upper) {
      constraint_type_clp[i] = GEQ;
    } else if (!has_lower && has_upper) {
      constraint_type_clp[i] = LEQ;
    } else if (has_lower && has_upper) {
      constraint_type_clp[i] = BOUND;
      (*nCols)++;
      (*nnz)++;
      (*nEqs)++;
    } else {
      // printf("Error: constraint %d has no lower and upper bound\n", i);
      // retcode = 1;
      // goto exit_cleanup;

      // what if regard free as bounded
      printf("Warning: constraint %d has no lower and upper bound\n", i);
      constraint_type_clp[i] = BOUND;
      (*nCols)++;
      (*nnz)++;
      (*nEqs)++;
    }
  }

  // allocate memory
  cupdlp_init_double(*cost, *nCols);
  cupdlp_init_double(*lower, *nCols);
  cupdlp_init_double(*upper, *nCols);
  cupdlp_init_int(*csc_beg, *nCols + 1);
  cupdlp_init_int(*csc_idx, *nnz);
  cupdlp_init_double(*csc_val, *nnz);
  cupdlp_init_double(*rhs, *nRows);

  // cost, lower, upper
  for (int i = 0; i < nCols_clp; i++) {
    (*cost)[i] = lp.col_cost_[i] * (*sign_origin);
    (*lower)[i] = lp.col_lower_[i];

    (*upper)[i] = lp.col_upper_[i];
  }
  // slack costs
  for (int i = nCols_clp; i < *nCols; i++) {
    (*cost)[i] = 0.0;
  }
  // slack bounds
  for (int i = 0, j = nCols_clp; i < *nRows; i++) {
    if (constraint_type_clp[i] == BOUND) {
      (*lower)[j] = lhs_clp[i];
      (*upper)[j] = rhs_clp[i];
      j++;
    }
  }

  for (int i = 0; i < *nCols; i++) {
    if ((*lower)[i] < -1e20) (*lower)[i] = -INFINITY;
    if ((*upper)[i] > 1e20) (*upper)[i] = INFINITY;
  }

  // permute LP rhs
  // EQ or BOUND first
  for (int i = 0, j = 0; i < *nRows; i++) {
    if (constraint_type_clp[i] == EQ) {
      (*rhs)[j] = lhs_clp[i];
      (*constraint_new_idx)[i] = j;
      j++;
    } else if (constraint_type_clp[i] == BOUND) {
      (*rhs)[j] = 0.0;
      (*constraint_new_idx)[i] = j;
      j++;
    }
  }
  // then LEQ or GEQ
  for (int i = 0, j = *nEqs; i < *nRows; i++) {
    if (constraint_type_clp[i] == LEQ) {
      (*rhs)[j] = -rhs_clp[i];  // multiply -1
      (*constraint_new_idx)[i] = j;
      j++;
    } else if (constraint_type_clp[i] == GEQ) {
      (*rhs)[j] = lhs_clp[i];
      (*constraint_new_idx)[i] = j;
      j++;
    }
  }

  // formulate and permute LP matrix
  // beg remains the same
  for (int i = 0; i < nCols_clp + 1; i++) (*csc_beg)[i] = A_csc_beg[i];
  for (int i = nCols_clp + 1; i < *nCols + 1; i++)
    (*csc_beg)[i] = (*csc_beg)[i - 1] + 1;

  // row idx changes
  for (int i = 0, k = 0; i < nCols_clp; i++) {
    // same order as in rhs
    // EQ or BOUND first
    for (int j = (*csc_beg)[i]; j < (*csc_beg)[i + 1]; j++) {
      if (constraint_type_clp[A_csc_idx[j]] == EQ ||
          constraint_type_clp[A_csc_idx[j]] == BOUND) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = A_csc_val[j];
        k++;
      }
    }
    // then LEQ or GEQ
    for (int j = (*csc_beg)[i]; j < (*csc_beg)[i + 1]; j++) {
      if (constraint_type_clp[A_csc_idx[j]] == LEQ) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = -A_csc_val[j];  // multiply -1
        k++;
      } else if (constraint_type_clp[A_csc_idx[j]] == GEQ) {
        (*csc_idx)[k] = (*constraint_new_idx)[A_csc_idx[j]];
        (*csc_val)[k] = A_csc_val[j];
        k++;
      }
    }
  }

  // slacks for BOUND
  for (int i = 0, j = nCols_clp; i < *nRows; i++) {
    if (constraint_type_clp[i] == BOUND) {
      (*csc_idx)[(*csc_beg)[j]] = (*constraint_new_idx)[i];
      (*csc_val)[(*csc_beg)[j]] = -1.0;
      j++;
    }
  }

  return retcode;
}

cupdlp_retcode problem_create(CUPDLPproblem **prob) {
  cupdlp_retcode retcode = RETCODE_OK;

  cupdlp_init_problem(*prob, 1);

  return retcode;
}

//cupdlp_retcode csc_create(CUPDLPcsc **csc_cpu) {
//  cupdlp_retcode retcode = RETCODE_OK;
//
//  cupdlp_init_csc_cpu(*csc_cpu, 1);
//
//  return retcode;
//}

cupdlp_retcode data_alloc(CUPDLPdata *data, cupdlp_int nRows, cupdlp_int nCols,
                          void *matrix, CUPDLP_MATRIX_FORMAT src_matrix_format,
                          CUPDLP_MATRIX_FORMAT dst_matrix_format) {
  cupdlp_retcode retcode = RETCODE_OK;

  data->nRows = nRows;
  data->nCols = nCols;
  data->matrix_format = dst_matrix_format;
  data->dense_matrix = cupdlp_NULL;
  data->csr_matrix = cupdlp_NULL;
  data->csc_matrix = cupdlp_NULL;
  data->device = CPU;

  switch (dst_matrix_format) {
    case DENSE:
      dense_create(&data->dense_matrix);
      dense_alloc_matrix(data->dense_matrix, nRows, nCols, matrix,
                                     src_matrix_format);
      break;
    case CSR:
      csr_create(&data->csr_matrix);
      csr_alloc_matrix(data->csr_matrix, nRows, nCols, matrix,
                                   src_matrix_format);
      break;
    case CSC:
      csc_create(&data->csc_matrix);
      csc_alloc_matrix(data->csc_matrix, nRows, nCols, matrix,
                                   src_matrix_format);
      break;
    case CSR_CSC:
      csc_create(&data->csc_matrix);
      csc_alloc_matrix(data->csc_matrix, nRows, nCols, matrix,
                                   src_matrix_format);
      csr_create(&data->csr_matrix);
      csr_alloc_matrix(data->csr_matrix, nRows, nCols, matrix,
                                   src_matrix_format);
      break;
    default:
      break;
  }
  // currently, only supprot that input matrix is CSC, and store both CSC and
  // CSR data->csc_matrix = matrix;

  return retcode;
}

cupdlp_retcode problem_alloc(
    CUPDLPproblem *prob, cupdlp_int nRows, cupdlp_int nCols, cupdlp_int nEqs,
    cupdlp_float *cost, cupdlp_float offset, cupdlp_float sign_origin,
    void *matrix, CUPDLP_MATRIX_FORMAT src_matrix_format,
    CUPDLP_MATRIX_FORMAT dst_matrix_format, cupdlp_float *rhs,
    cupdlp_float *lower, cupdlp_float *upper, cupdlp_float *alloc_matrix_time,
    cupdlp_float *copy_vec_time) {
  cupdlp_retcode retcode = RETCODE_OK;
  prob->nRows = nRows;
  prob->nCols = nCols;
  prob->nEqs = nEqs;
  prob->data = cupdlp_NULL;
  prob->cost = cupdlp_NULL;
  prob->offset = offset;
  prob->sign_origin = sign_origin;
  prob->rhs = cupdlp_NULL;
  prob->lower = cupdlp_NULL;
  prob->upper = cupdlp_NULL;

  cupdlp_float begin = getTimeStamp();

  cupdlp_init_data(prob->data, 1);
  cupdlp_init_vec_double(prob->cost, nCols);
  cupdlp_init_vec_double(prob->rhs, nRows);
  cupdlp_init_vec_double(prob->lower, nCols);
  cupdlp_init_vec_double(prob->upper, nCols);
  cupdlp_init_zero_vec_double(prob->hasLower, nCols);
  cupdlp_init_zero_vec_double(prob->hasUpper, nCols);

  data_alloc(prob->data, nRows, nCols, matrix, src_matrix_format,
                         dst_matrix_format);
  *alloc_matrix_time = getTimeStamp() - begin;

  prob->data->csc_matrix->MatElemNormInf = infNorm(
        ((CUPDLPcsc *)matrix)->colMatElem, ((CUPDLPcsc *)matrix)->nMatElem);

  begin = getTimeStamp();
  cupdlp_copy_vec(prob->cost, cost, cupdlp_float, nCols);
  cupdlp_copy_vec(prob->rhs, rhs, cupdlp_float, nRows);
  cupdlp_copy_vec(prob->lower, lower, cupdlp_float, nCols);
  cupdlp_copy_vec(prob->upper, upper, cupdlp_float, nCols);
  *copy_vec_time = getTimeStamp() - begin;

  // Keep
  cupdlp_haslb(prob->hasLower, prob->lower, -INFINITY, nCols);
  cupdlp_hasub(prob->hasUpper, prob->upper, +INFINITY, nCols);

  // TODO: cal dMaxCost, dMaxRhs, dMaxRowBound

  return retcode;
}

// ToDo: Why can linker not pick up infNorm, cupdlp_haslb and
// cupdlp_hasub from pdlp/cupdlp/cupdlp_linalg.c?
double infNorm(double *x, cupdlp_int n) {
  double norm = 0;
  for (HighsInt iX = 0; iX < n; iX++)
    norm = std::max(std::fabs(x[iX]), norm);
  return norm;
}
void cupdlp_haslb(cupdlp_float *haslb, const cupdlp_float *lb,
                  const cupdlp_float bound, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    haslb[i] = lb[i] > bound ? 1.0 : 0.0;
  }
}

void cupdlp_hasub(cupdlp_float *hasub, const cupdlp_float *ub,
                  const cupdlp_float bound, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    hasub[i] = ub[i] < bound ? 1.0 : 0.0;
  }
}


void reportParams(CUPDLPwork *w,
		  cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
		  cupdlp_bool *ifChangeFloatParam,
		  cupdlp_float *floatParam) {
  PDHG_PrintPDHGParam(w);  
}

void getUserParamsFromOptions(const HighsOptions& options,
			      cupdlp_bool *ifChangeIntParam,
			      cupdlp_int *intParam,
			      cupdlp_bool *ifChangeFloatParam,
			      cupdlp_float *floatParam) {
  for (cupdlp_int i = 0; i < N_INT_USER_PARAM; ++i) 
    ifChangeIntParam[i] = false;
  for (cupdlp_int i = 0; i < N_FLOAT_USER_PARAM; ++i) 
    ifChangeFloatParam[i] = false;
  // Assume all PDLP-related options in HiGHS cause changes
  ifChangeIntParam[N_ITER_LIM] = true;
  intParam[N_ITER_LIM] = options.pdlp_iteration_limit;
  //
  ifChangeIntParam[IF_SCALING] = true;
  intParam[IF_SCALING] = options.pdlp_scaling ? 1 : 0;
  //
  ifChangeFloatParam[D_PRIMAL_TOL] = true;
  floatParam[D_PRIMAL_TOL] = options.primal_feasibility_tolerance;
  //
  ifChangeFloatParam[D_DUAL_TOL] = true;
  floatParam[D_DUAL_TOL] = options.dual_feasibility_tolerance;
  //
  ifChangeFloatParam[D_GAP_TOL] = true;
  floatParam[D_GAP_TOL] = options.pdlp_d_gap_tol;
  //
  ifChangeFloatParam[D_TIME_LIM] = true;
  floatParam[D_TIME_LIM] = options.time_limit;
  //
  ifChangeIntParam[E_RESTART_METHOD] = true;
  intParam[E_RESTART_METHOD] = options.pdlp_e_restart_method;
}

