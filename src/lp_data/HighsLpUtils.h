/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLpUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSLPUTILS_H_
#define LP_DATA_HIGHSLPUTILS_H_

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"

class HighsLp;

// Methods taking HighsLp as an argument
HighsStatus checkLp(
		    const HighsLp& lp
		    );

HighsStatus assessLp(
		     HighsLp& lp,
		     const HighsOptions& options,
		     const bool normalise
		     );

HighsStatus assessLpDimensions(
			       const HighsLp& lp
			       );

HighsStatus assess_costs(const int col_ix_os,
			 const int mask_num_col,
			 const bool interval,
			 const int from_col,
			 const int to_col,
			 const bool set,
			 const int num_set_entries,
			 const int* col_set,
			 const bool mask,
			 const int* col_mask,
			 const double* usr_col_cost,
			 const double infinite_cost
			 );

HighsStatus assess_bounds(
			  const char* type,
			  const int ix_os,
			  const int mask_num_ix,
			  const bool interval,
			  const int from_ix,
			  const int to_ix,
			  const bool set,
			  const int num_set_entries,
			  const int* ix_set,
			  const bool mask,
			  const int* ix_mask,
			  double* usr_lower,
			  double* usr_upper,
			  const double infinite_bound,
			  const bool normalise
			  );

HighsStatus assessMatrix(
			  const int vec_dim,
			  const int from_ix,
			  const int to_ix,
			  const int num_vec,
			  int& num_nz,
			  int* Xstart,
			  int* Xindex,
			  double* Xvalue,
			  const double small_matrix_value,
			  const double large_matrix_value,
			  bool normalise
			  );

HighsStatus add_lp_cols(
			HighsLp& lp,
			const int num_new_col,
			const double *XcolCost,
			const double *XcolLower,
			const double *XcolUpper,
			const int num_new_nz,
			const int *XAstart,
			const int *XAindex,
			const double *XAvalue,
			const HighsOptions& options
			);

HighsStatus append_lp_cols(
		   HighsLp& lp,
		   const int num_new_col,
		   const double *XcolCost,
		   const double *XcolLower,
		   const double *XcolUpper,
		   const int num_new_nz,
		   const int *XAstart,
		   const int *XAindex,
		   const double *XAvalue,
		   const HighsOptions& options,
		   const bool valid_matrix
		   );

HighsStatus append_cols_to_lp_vectors(
			       HighsLp &lp,
			       const int num_new_col,
			       const double *XcolCost,
			       const double *colLower,
			       const double *XcolUpper
			       );

HighsStatus append_cols_to_lp_matrix(
			   HighsLp &lp,
			   const int num_new_col,
			   const int num_new_nz,
			   const int *XAstart,
			   const int *XAindex,
			   const double *XAvalue
			   );

HighsStatus add_lp_rows(
		HighsLp& lp,
		const int num_new_row,
		const double *XrowLower,
		const double *XrowUpper,
		const int num_new_nz,
		const int *XARstart,
		const int *XARindex,
		const double *XARvalue,
		const HighsOptions& options
		);

HighsStatus append_lp_rows(
			   HighsLp& lp,
			   const int num_new_row,
			   const double *XrowLower,
			   const double *XrowUpper,
			   const int num_new_nz,
			   const int *XARstart,
			   const int *XARindex,
			   const double *XARvalue,
			   const HighsOptions& options,
			   bool valid_matrix
			   );

HighsStatus append_rows_to_lp_vectors(HighsLp &lp,
			       const int num_new_row,
			       const double *XrowLower,
			       const double *XrowUpper
			       );

HighsStatus append_rows_to_lp_matrix(HighsLp &lp,
			   const int num_new_row,
			   const int num_new_nz,
			   const int *XARstart,
			   const int *XARindex,
			   const double *XARvalue
			   );

HighsStatus delete_lp_cols(
			   HighsLp &lp,
			   const bool interval, const int from_col, const int to_col,
			   const bool set, const int num_set_entries, const int* col_set,
			   const bool mask, int* col_mask,
			   const bool valid_matrix);

HighsStatus delete_cols_from_lp_vectors(
			   HighsLp &lp,
			   int &new_num_col,
			   const bool interval, const int from_col, const int to_col,
			   const bool set, const int num_set_entries, const int* col_set,
			   const bool mask, const int* col_mask);

HighsStatus delete_cols_from_lp_matrix(
			   HighsLp &lp,
			   const bool interval, const int from_col, const int to_col,
			   const bool set, const int num_set_entries, const int* col_set,
			   const bool mask, int* col_mask);

HighsStatus delete_lp_rows(
			   HighsLp &lp,
			   const bool interval, const int from_row, const int to_row,
			   const bool set, const int num_set_entries, const int* row_set,
			   const bool mask, int* row_mask,
			   const bool valid_matrix);

HighsStatus delete_rows_from_lp_vectors(
			   HighsLp &lp,
			   int &new_num_row,
			   const bool interval, const int from_row, const int to_row,
			   const bool set, const int num_set_entries, const int* row_set,
			   const bool mask, const int* row_mask);

HighsStatus delete_rows_from_lp_matrix(
			   HighsLp &lp,
			   const bool interval, const int from_row, const int to_row,
			   const bool set, const int num_set_entries, const int* row_set,
			   const bool mask, int* row_mask);

HighsStatus change_lp_matrix_coefficient(
				  HighsLp &lp,
				  const int row,
				  const int col,
				  const double new_value
				  );

HighsStatus change_lp_costs(
			    HighsLp &lp,
			    const bool interval,
			    const int from_col,
			    const int to_col,
			    const bool set,
			    const int num_set_entries,
			    const int* col_set,
			    const bool mask,
			    const int* col_mask,
			    const double* usr_col_cost,
			    const double infinite_cost
			    );

HighsStatus change_lp_col_bounds(
				 HighsLp &lp,
				 const bool interval,
				 const int from_col,
				 const int to_col,
				 const bool set,
				 const int num_set_entries,
				 const int* col_set,
				 const bool mask,
				 const int* col_mask,
				 const double* usr_col_lower,
				 const double* usr_col_upper,
				 const double infinite_bound
				 );

HighsStatus change_lp_row_bounds(
				 HighsLp &lp,
				 const bool interval,
				 const int from_row,
				 const int to_row,
				 const bool set,
				 const int num_set_entries,
				 const int* row_set,
				 const bool mask,
				 const int* row_mask,
				 const double* usr_row_lower,
				 const double* usr_row_upper,
				 const double infinite_bound
				 );

HighsStatus change_bounds(
			  const char* type,
			  double* lower,
			  double* upper,
			  const int mask_num_ix,
			  const bool interval,
			  const int from_ix,
			  const int to_ix,
			  const bool set,
			  const int num_set_entries,
			  const int* ix_set,
			  const bool mask,
			  const int* ix_mask,
			  const double* usr_lower,
			  const double* usr_upper,
			  const double infinite_bound
			  );

/**
 * @brief Write out the LP as an MPS file
 */
bool writeLpAsMPS(
		  const char* filename,
		  const HighsLp& lp);		  

/**
 * @brief Report the data of an LP
 */
void reportLp(
	      const HighsLp &lp, //!< LP whose data are to be reported
	      const int report_level = 0 //!< 0 => scalar [dimensions]; 1=> vector [costs/bounds]; 2 => vector+matrix
	      );
/**
 * @brief Report the brief data of an LP 
 */
void reportLpBrief(
		   const HighsLp &lp //!< LP whose data are to be reported
		   );
/**
 * @brief Report the data of an LP
 */
void reportLpDimensions(
			const HighsLp &lp //!< LP whose data are to be reported
			);
/**
 * @brief Report the data of an LP
 */
void reportLpObjSense(
		      const HighsLp &lp //!< LP whose data are to be reported
		      );
/**
 * @brief Report the data of an LP
 */
void reportLpColVectors(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpRowVectors(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpColMatrix(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );

void reportMatrix(
	       const char* message,
	       const int num_col,
	       const int num_nz,
	       const int* start,
	       const int* index,
	       const double* value
	       );

// Get the costs for a contiguous set of columns
HighsStatus getLpCosts(
		const HighsLp& lp,
		const int from_col,
		const int to_col,
		double* XcolCost
		);

// Get the bounds for a contiguous set of columns
HighsStatus getLpColBounds(
		    const HighsLp& lp,
		    const int from_col,
		    const int to_col,
		    double* XcolLower,
		    double* XcolUpper
		    );

// Get the bounds for a contiguous set of rows
HighsStatus getLpRowBounds(
		    const HighsLp& lp,
		    const int from_row,
		    const int to_row,
		    double* XrowLower,
		    double* XrowUpper
		    );

HighsStatus getLpMatrixCoefficient(
			    const HighsLp& lp,
			    const int row,
			    const int col,
			    double *val
			    );



/*
  void reportLpSolution(
  HighsModelObject &highs_model //!< Model object whose LP solution is to be reported
  );
*/
#ifdef HiGHSDEV
// Analyse the data in an LP problem
void util_analyseLp(
		    const HighsLp &lp,
		    const char* message
		    );
#endif

HighsNewBasis getHighsBasis(const HighsLp& lp, const HighsBasis& basis);

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution);
HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution);
double calculateObjective(const HighsLp& lp, HighsSolution& solution);

HighsStatus assess_interval_set_mask(
				     const int max_ix, 
				     const bool interval,
				     const int from_ix,
				     const int to_ix,
				     const bool set,
				     const int num_set_entries,
				     const int* ix_set,
				     const bool mask,
				     const int* ix_mask,
				     int &from_k,
				     int &to_k
				     );

void update_out_in_ix(const int ix_dim, 
		      const bool interval,
		      const int from_ix,
		      const int to_ix,
		      const bool set,
		      const int num_set_entries,
		      const int* ix_set,
		      const bool mask,
		      const int* ix_mask,
		      int& out_from_ix,
		      int& out_to_ix,
		      int& in_from_ix,
		      int& in_to_ix,
		      int& current_set_entry);

bool isColDataNull(const double *usr_col_cost, const double *usr_col_lower,  const double *usr_col_upper);
bool isRowDataNull(const double *usr_row_lower,  const double *usr_row_upper);
bool isMatrixDataNull(const int *usr_matrix_start, const int *usr_matrix_index, const double *usr_matrix_value);

HighsLp transformIntoEqualityProblem(const HighsLp& lp);
HighsLp dualizeEqualityProblem(const HighsLp& lp);

#endif // LP_DATA_HIGHSLPUTILS_H_
