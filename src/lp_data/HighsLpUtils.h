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
		     bool normalise
		     );

HighsStatus assessLpDimensions(
			       const HighsLp& lp
			       );

HighsStatus assess_costs(int col_ix_os,
			 int mask_num_col,
			 bool interval,
			 int from_col,
			 int to_col,
			 bool set,
			 int num_set_entries,
			 const int* col_set,
			 bool mask,
			 const int* col_mask,
			 const double* usr_col_cost,
			 double infinite_cost
			 );

HighsStatus assess_bounds(
			  const char* type,
			  int ix_os,
			  int mask_num_ix,
			  bool interval,
			  int from_ix,
			  int to_ix,
			  bool set,
			  int num_set_entries,
			  const int* ix_set,
			  bool mask,
			  const int* ix_mask,
			  double* usr_lower,
			  double* usr_upper,
			  const double infinite_bound,
			  bool normalise
			  );

HighsStatus assessMatrix(
			  int Xvec_dim,
			  int Xfrom_ix,
			  int Xto_ix,
			  int Xnum_vec,
			  int Xnum_nz,
			  int* Xstart,
			  int* Xindex,
			  double* Xvalue,
			  double small_matrix_value,
			  double large_matrix_value,
			  bool normalise
			  );

HighsStatus assess_interval_set_mask(
				     int mask_num_ix, 
				     bool interval,
				     int from_ix,
				     int to_ix,
				     bool set,
				     int num_set_entries,
				     const int* ix_set,
				     bool mask,
				     const int* ix_mask,
				     int from_k,
				     int to_k
				     );

HighsStatus add_lp_cols(
			HighsLp& lp,
			int XnumNewCol,
			const double *XcolCost,
			const double *XcolLower,
			const double *XcolUpper,
			int XnumNewNZ,
			const int *XAstart,
			const int *XAindex,
			const double *XAvalue,
			const HighsOptions& options
			);

HighsStatus append_lp_cols(
		   HighsLp& lp,
		   int XnumNewCol,
		   const double *XcolCost,
		   const double *XcolLower,
		   const double *XcolUpper,
		   int XnumNewNZ,
		   const int *XAstart,
		   const int *XAindex,
		   const double *XAvalue,
		   const HighsOptions& options,
		   const bool valid_matrix
		   );

HighsStatus append_cols_to_lp_vectors(
			       HighsLp &lp,
			       int XnumNewCol,
			       const double *XcolCost,
			       const double *colLower,
			       const double *XcolUpper
			       );

HighsStatus append_cols_to_lp_matrix(
			   HighsLp &lp,
			   int XnumNewCol,
			   int XnumNewNZ,
			   const int *XAstart,
			   const int *XAindex,
			   const double *XAvalue
			   );

HighsStatus add_lp_rows(
		HighsLp& lp,
		int XnumNewRow,
		const double *XrowLower,
		const double *XrowUpper,
		int XnumNewNZ,
		const int *XARstart,
		const int *XARindex,
		const double *XARvalue,
		const HighsOptions& options
		);

HighsStatus append_lp_rows(
		   HighsLp& lp,
		   int XnumNewRow,
		   const double *XrowLower,
		   const double *XrowUpper,
		   int XnumNewNZ,
		   const int *XARstart,
		   const int *XARindex,
		   const double *XARvalue,
		   const HighsOptions& options
		   );

HighsStatus append_rows_to_lp_vectors(HighsLp &lp,
			       int XnumNewRow,
			       const double *XrowLower,
			       const double *XrowUpper
			       );

HighsStatus append_rows_to_lp_matrix(HighsLp &lp,
			   int XnumNewRow,
			   int XnumNewNZ,
			   const int *XARstart,
			   const int *XARindex,
			   const double *XARvalue
			   );

HighsStatus delete_lp_cols(
			   HighsLp &lp,
			   bool interval, int from_col, int to_col,
			   bool set, int num_set_entries, const int* col_set,
			   bool mask, const int* col_mask,
			   bool valid_matrix);

HighsStatus delete_cols_from_lp_vectors(
			   HighsLp &lp,
			   bool interval, int from_col, int to_col,
			   bool set, int num_set_entries, const int* col_set,
			   bool mask, const int* col_mask);

HighsStatus delete_cols_from_lp_matrix(
			   HighsLp &lp,
			   bool interval, int from_col, int to_col,
			   bool set, int num_set_entries, const int* col_set,
			   bool mask, const int* col_mask);

HighsStatus delete_lp_rows(
			   HighsLp &lp,
			   bool interval, int from_row, int to_row,
			   bool set, int num_set_entries, const int* row_set,
			   bool mask, const int* row_mask,
			   bool valid_matrix);

HighsStatus delete_rows_from_lp_vectors(
			   HighsLp &lp,
			   bool interval, int from_row, int to_row,
			   bool set, int num_set_entries, const int* row_set,
			   bool mask, const int* row_mask);

HighsStatus delete_rows_from_lp_matrix(
			   HighsLp &lp,
			   bool interval, int from_row, int to_row,
			   bool set, int num_set_entries, const int* row_set,
			   bool mask, const int* row_mask);

HighsStatus change_lp_matrix_coefficient(
				  HighsLp &lp,
				  int Xrow,
				  int Xcol,
				  const double XnewValue
				  );

HighsStatus change_lp_costs(
			    HighsLp &lp,
			    bool interval,
			    int from_col,
			    int to_col,
			    bool set,
			    int num_set_entries,
			    const int* col_set,
			    bool mask,
			    const int* col_mask,
			    const double* usr_col_cost,
			    const double infinite_cost
			    );

HighsStatus change_lp_col_bounds(
				 HighsLp &lp,
				 bool interval,
				 int from_col,
				 int to_col,
				 bool set,
				 int num_set_entries,
				 const int* col_set,
				 bool mask,
				 const int* col_mask,
				 const double* usr_col_lower,
				 const double* usr_col_upper,
				 const double infinite_bound
				 );

HighsStatus change_lp_row_bounds(
				 HighsLp &lp,
				 bool interval,
				 int from_row,
				 int to_row,
				 bool set,
				 int num_set_entries,
				 const int* row_set,
				 bool mask,
				 const int* row_mask,
				 const double* usr_row_lower,
				 const double* usr_row_upper,
				 const double infinite_bound
				 );

HighsStatus change_bounds(
			  const char* type,
			  double* lower,
			  double* upper,
			  int mask_num_ix,
			  bool interval,
			  int from_ix,
			  int to_ix,
			  bool set,
			  int num_set_entries,
			  const int* ix_set,
			  bool mask,
			  const int* ix_mask,
			  const double* usr_lower,
			  const double* usr_upper,
			  const double infinite_bound
			  );


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
void reportLpColVec(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpRowVec(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpColMtx(
		    const HighsLp &lp //!< LP whose data are to be reported
		    );

// Get the costs for a contiguous set of columns
HighsStatus getLpCosts(
		const HighsLp& lp,
		int Xfrom_col,
		int Xto_col,
		double* XcolCost
		);

// Get the bounds for a contiguous set of columns
HighsStatus getLpColBounds(
		    const HighsLp& lp,
		    int Xfrom_col,
		    int Xto_col,
		    double* XcolLower,
		    double* XcolUpper
		    );

// Get the bounds for a contiguous set of rows
HighsStatus getLpRowBounds(
		    const HighsLp& lp,
		    int Xfrom_row,
		    int Xto_row,
		    double* XrowLower,
		    double* XrowUpper
		    );

HighsStatus getLpMatrixCoefficient(
			    const HighsLp& lp,
			    int row,
			    int col,
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

HighsBasis_new getHighsBasis(const HighsLp& lp, const HighsBasis& basis);

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution);
HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution);

#endif // LP_DATA_HIGHSLPUTILS_H_
