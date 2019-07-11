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

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"

void setSimplexOptions(
	     HighsModelObject& highs_model_object  //!< Model object in which simplex
	                                           //!< options are to be set
	     );

SimplexSolutionStatus transition(
				 HighsModelObject& highs_model_object  //!< Model object
				 );

bool dual_infeasible(const double value, const double lower, const double upper,
                     const double dual, const double value_tolerance,
                     const double dual_tolerance);

// Methods not requiring HighsModelObject

void append_nonbasic_cols_to_basis(HighsLp& lp, HighsBasis& basis,
                                   int XnumNewCol);

void append_nonbasic_cols_to_basis(HighsLp& lp, SimplexBasis& simplex_basis,
                                   int XnumNewCol);

void append_basic_rows_to_basis(HighsLp& lp, HighsBasis& basis, int XnumNewRow);

bool highs_basis_ok(
		    //		    HighsLp& lp, HighsBasis& basis
		    );

bool nonbasic_flag_basic_index_ok(HighsLp& lp, SimplexBasis& simplex_basis);

#ifdef HiGHSDEV
void report_basis(HighsLp& lp, HighsBasis& basis);
void report_basis(HighsLp& lp, SimplexBasis& simplex_basis);
#endif

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
  string filename = "z-" + simplex_lp_->model_name_ + "-" + suffix;
  ofstream output(filename.c_str());
  int count = historyColumnIn.size();
  double current_run_highs_time = timer_->readRunHighsClock();
  output << simplex_lp_->model_name_ << " " << count << "\t" <<
current_run_highs_time << endl; output << setprecision(12); for (int i = 0; i <
count; i++) { output << historyColumnIn[i] << "\t"; output <<
historyColumnOut[i] << "\t"; output << historyAlpha[i] << endl;
  }
  output.close();
}
#endif
*/
void computeDualObjectiveValue(HighsModelObject& highs_model_object,
                               int phase = 2);

void computePrimalObjectiveValue(HighsModelObject& highs_model_object);

void initialiseSimplexLpRandomVectors(HighsModelObject& highs_model);

// SCALE:

void scaleHighsModelInit(HighsModelObject& highs_model);

void scaleCosts(HighsModelObject& highs_model);

void scaleFactorRanges(HighsModelObject& highs_model_object,
                       double& min_col_scale, double& max_col_scale,
                       double& min_row_scale, double& max_row_scale);

void scaleSimplexLp(HighsModelObject& highs_model);

// PERMUTE:

void permuteSimplexLp(HighsModelObject& highs_model);

void initialise_basic_index(HighsModelObject& highs_model_object);

void allocate_work_and_base_arrays(HighsModelObject& highs_model_object);

void initialise_from_nonbasic(HighsModelObject& highs_model_object);

void replace_from_nonbasic(HighsModelObject& highs_model_object);

void initialise_with_logical_basis(HighsModelObject& highs_model_object);

void initialise_value_from_nonbasic(HighsModelObject& highs_model_object,
                                    int firstvar, int lastvar);

void initialise_value(HighsModelObject& highs_model_object);

void initialise_phase2_col_bound(HighsModelObject& highs_model_object,
                                 int firstcol, int lastcol);

void initialise_phase2_row_bound(HighsModelObject& highs_model_object,
                                 int firstrow, int lastrow);

void initialise_bound(HighsModelObject& highs_model_object, int phase = 2);

void initialise_phase2_col_cost(HighsModelObject& highs_model_object,
                                int firstcol, int lastcol);

void initialise_phase2_row_cost(HighsModelObject& highs_model_object,
                                int firstrow, int lastrow);

void initialise_cost(HighsModelObject& highs_model_object, int perturb = 0);

int get_nonbasicMove(HighsModelObject& highs_model_object, int var);

void populate_work_arrays(HighsModelObject& highs_model_object);

void replace_with_logical_basis(HighsModelObject& highs_model_object);

void replace_with_new_basis(HighsModelObject& highs_model_object,
                            const int* XbasicIndex);

void setup_num_basic_logicals(HighsModelObject& highs_model_object);

#ifdef HiGHSDEV
void reportSimplexProfiling(HighsModelObject& highs_model_object);

#endif

/**
 * @brief Get the Hager condition number estimate for the basis matrix of a
 * model
 */
double computeBasisCondition(HighsModelObject& highs_model_object);

bool work_arrays_ok(HighsModelObject& highs_model_object, int phase);

bool one_nonbasic_move_vs_work_arrays_ok(HighsModelObject& highs_model_object,
                                         int var);

bool all_nonbasic_move_vs_work_arrays_ok(HighsModelObject& highs_model_object);

bool ok_to_solve(HighsModelObject& highs_model_object, int level, int phase);

void flip_bound(HighsModelObject& highs_model_object, int iCol);

int compute_factor(HighsModelObject& highs_model_object);

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void compute_primal(HighsModelObject& highs_model_object);

void computePrimalInfeasible(HighsModelObject& highs_model_object,
                             const bool report = false);

void computeDualInfeasible(HighsModelObject& highs_model_object,
                           const bool report = false);

void computeDualInfeasibleWithFlips(HighsModelObject& highs_model_object,
				    const bool report = false);

void compute_dual(HighsModelObject& highs_model_object);

void correct_dual(HighsModelObject& highs_model_object,
                  int* free_infeasibility_count);

// Record the shift in the cost of a particular column
void shift_cost(HighsModelObject& highs_model_object, int iCol, double amount);

// Undo the shift in the cost of a particular column
void shift_back(HighsModelObject& highs_model_object, int iCol);

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void update_factor(HighsModelObject& highs_model_object, HVector* column,
                   HVector* row_ep, int* iRow, int* hint);

void update_pivots(HighsModelObject& highs_model_object, int columnIn,
                   int rowOut, int sourceOut);

void update_matrix(HighsModelObject& highs_model_object, int columnIn,
                   int columnOut);

void logRebuild(HighsModelObject& highs_model_object, const bool primal,
                const int solve_phase);

std::string SimplexSolutionStatusToString(SimplexSolutionStatus status);

void reportSimplexLpStatus(
    HighsSimplexLpStatus&
        simplex_lp_status,  // !< Status of simplex LP to be reported
    const char* message = "");

void invalidateSimplexLpData(
    HighsSimplexLpStatus& simplex_lp_status  // !< Status of simplex LP whose
                                             // data are to be invalidated
);

void invalidateSimplexLp(
    HighsSimplexLpStatus&
        simplex_lp_status  // !< Status of simplex LP to be invalidated
);

void updateSimplexLpStatus(
    HighsSimplexLpStatus&
        simplex_lp_status,  // !< Status of simplex LP to be updated
    LpAction action         // !< Action prompting update
);

SimplexSolutionStatus solveUnconstrainedLp(HighsModelObject& highs_model_object);
#endif  // SIMPLEX_HSIMPLEX_H_
