/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "HighsModelObject.h"
#include "HighsOptions.h"

/*
// Increment iteration count (here!) and (possibly) store the pivots for
// debugging NLA
void record_pivots(int columnIn, int columnOut, double alpha) {
  // NB This is where the iteration count is updated!
  if (columnIn >= 0) simplex_info_.iteration_count++;
#ifdef HiGHSDEV
  historyColumnIn.push_back(columnIn);
  historyColumnOut.push_back(columnOut);
  historyAlpha.push_back(alpha);
#endif
}
#ifdef HiGHSDEV
// Store and write out the pivots for debugging NLA
void writePivots(const char* suffix) {
  string filename = "z-" + solver_lp_->model_name_ + "-" + suffix;
  ofstream output(filename.c_str());
  int count = historyColumnIn.size();
  double current_run_highs_time = timer_->readRunHighsClock();
  output << solver_lp_->model_name_ << " " << count << "\t" <<
current_run_highs_time << endl; output << setprecision(12); for (int i = 0; i <
count; i++) { output << historyColumnIn[i] << "\t"; output <<
historyColumnOut[i] << "\t"; output << historyAlpha[i] << endl;
  }
  output.close();
}
#endif
*/
void clear_solver_lp_data(
    HighsModelObject &highs_model_object //!< Model object in which data for LP
                                         //!< to be solved is to be cleared
);

void clear_solver_lp(
    HighsModelObject &highs_model_object //!< Model object in which LP to be
                                         //!< solved is to be cleared
);

void options(
    HighsModelObject &highs_model_object, //!< Model object in which simplex
                                          //!< options are to be set
    const HighsOptions &opt               //!< HiGHS options
);

void update_solver_lp_status_flags(HighsModelObject &highs_model_object,
                                   LpAction action);

void report_basis(HighsModelObject &highs_model_object);

void report_solver_lp_status_flags(HighsModelObject &highs_model_object);
void compute_dual_objective_value(HighsModelObject &highs_model_object,
                                  int phase = 2);

void initialise_solver_lp_random_vectors(HighsModelObject &highs_model);

// TRANSPOSE:

void transpose_solver_lp(HighsModelObject &highs_model);

#ifdef HiGHSDEV
// Information on large costs
const double tlLargeCo = 1e5;
int numLargeCo;
vector<int> largeCostFlag;
double largeCostScale;
#endif

void scaleHighsModelInit(HighsModelObject &highs_model);

void scaleCosts(HighsModelObject &highs_model);

void scale_solver_lp(HighsModelObject &highs_model);

// PERMUTE:

void permute_solver_lp(HighsModelObject &highs_model);

// TIGHTEN:

void tighten_solver_lp(HighsModelObject &highs_model);

void initialise_basic_index(HighsModelObject &highs_model_object);

void allocate_work_and_base_arrays(HighsModelObject &highs_model_object);

void initialise_from_nonbasic(HighsModelObject &highs_model_object);

void replace_from_nonbasic(HighsModelObject &highs_model_object);

void initialise_with_logical_basis(HighsModelObject &highs_model_object);

void initialise_value_from_nonbasic(HighsModelObject &highs_model_object,
                                    int firstvar, int lastvar);

void initialise_value(HighsModelObject &highs_model_object);

void initialise_phase2_col_bound(HighsModelObject &highs_model_object,
                                 int firstcol, int lastcol);

void initialise_phase2_row_bound(HighsModelObject &highs_model_object,
                                 int firstrow, int lastrow);

void initialise_bound(HighsModelObject &highs_model_object, int phase = 2);

void initialise_phase2_col_cost(HighsModelObject &highs_model_object,
                                int firstcol, int lastcol);

void initialise_phase2_row_cost(HighsModelObject &highs_model_object,
                                int firstrow, int lastrow);

void initialise_cost(HighsModelObject &highs_model_object, int perturb = 0);

int get_nonbasicMove(HighsModelObject &highs_model_object, int var);

void populate_work_arrays(HighsModelObject &highs_model_object);

void replace_with_logical_basis(HighsModelObject &highs_model_object);

void replace_with_new_basis(HighsModelObject &highs_model_object,
                            const int *XbasicIndex);

void extend_with_logical_basis(HighsModelObject &highs_model_object,
                               int firstcol, int lastcol, int firstrow,
                               int lastrow);

void setup_num_basic_logicals(HighsModelObject &highs_model_object);

void setup_for_solve(HighsModelObject &highs_model_object);

bool nonbasic_flag_basic_index_ok(HighsModelObject &highs_model_object,
                                  int XnumCol, int XnumRow);

bool work_arrays_ok(HighsModelObject &highs_model_object, int phase);

bool one_nonbasic_move_vs_work_arrays_ok(HighsModelObject &highs_model_object,
                                         int var);

bool all_nonbasic_move_vs_work_arrays_ok(HighsModelObject &highs_model_object);

bool ok_to_solve(HighsModelObject &highs_model_object, int level, int phase);

void flip_bound(HighsModelObject &highs_model_object, int iCol);
/*
int handle_rank_deficiency(HighsModelObject &highs_model_object) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HFactor &factor = highs_model_object.factor_;
  HighsBasis &basis = highs_model_object.basis_;
  int rankDeficiency = factor.rankDeficiency;
  const int *noPvC = factor.getNoPvC();
  printf("Returned %d = factor.build();\n", rankDeficiency);
  fflush(stdout);
  vector<int> basicRows;
  const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
  basicRows.resize(numTot);
  //    printf("Before - basis.basicIndex_:"); for (int iRow=0;
iRow<solver_lp.numRow_; iRow++)
  //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
  for (int iRow = 0; iRow < solver_lp.numRow_; iRow++)
basicRows[basis.basicIndex_[iRow]] = iRow; for (int k = 0; k < rankDeficiency;
k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = solver_lp.numCol_ + factor.noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut,
simplex_info.workLower_[columnOut],
    //      simplex_info.workUpper_[columnOut]);
    if (basis.basicIndex_[rowOut] != columnOut) {
      printf("%d = basis.basicIndex_[rowOut] != noPvC[k] = %d\n",
basis.basicIndex_[rowOut], columnOut); fflush(stdout);
    }
    int sourceOut = setSourceOutFmBd(columnOut);
    updatePivots(columnIn, rowOut, sourceOut);
    updateMatrix(columnIn, columnOut);
  }
  //    printf("After  - basis.basicIndex_:"); for (int iRow=0;
iRow<solver_lp.numRow_; iRow++)
  //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
#ifdef HiGHSDEV
  factor.checkInvert();
#endif
  return 0;
}
*/
int compute_factor(HighsModelObject &highs_model_object);

void compute_primal(HighsModelObject &highs_model_object);

void compute_dual(HighsModelObject &highs_model_object);

void correct_dual(HighsModelObject &highs_model_object,
                  int *free_infeasibility_count);

void compute_dual_infeasible_in_dual(HighsModelObject &highs_model_object,
                                     int *dual_infeasibility_count);

void compute_dual_infeasible_in_primal(HighsModelObject &highs_model_object,
                                       int *dual_infeasibility_count);

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
int set_source_out_from_bound(HighsModelObject &highs_model_object,
                              const int column_out);

double
compute_primal_objective_function_value(HighsModelObject &highs_model_object);

// Record the shift in the cost of a particular column
double shift_cost(HighsModelObject &highs_model_object, int iCol,
                  double amount);

// Undo the shift in the cost of a particular column
double shift_back(HighsModelObject &highs_model_object, int iCol);

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void update_factor(HighsModelObject &highs_model_object, HVector *column,
                   HVector *row_ep, int *iRow, int *hint);

void update_pivots(HighsModelObject &highs_model_object, int columnIn,
                   int rowOut, int sourceOut);

void update_matrix(HighsModelObject &highs_model_object, int columnIn,
                   int columnOut);

#ifdef HiGHSDEV
void util_analyse_lp_solution(HighsModelObject &highs_model_object);
#endif

void report_iteration_count_dual_objective_value(
    HighsModelObject &highs_model_object, int i_v);

#endif // SIMPLEX_HSIMPLEX_H_
