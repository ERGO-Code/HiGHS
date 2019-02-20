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

HighsStatus checkLp(const HighsLp& lp);

// Methods taking HighsLp as an argument



// Get the costs for a contiguous set of columns
void getLpCosts(
		const HighsLp& lp,
		int firstcol,
		int lastcol,
		double* XcolCost
		);

// Get the bounds for a contiguous set of columns
void getLpColBounds(
		    const HighsLp& lp,
		    int firstcol,
		    int lastcol,
		    double* XcolLower,
		    double* XcolUpper
		    );

// Get the bounds for a contiguous set of rows
void getLpRowBounds(
		    const HighsLp& lp,
		    int firstrow,
		    int lastrow,
		    double* XrowLower,
		    double* XrowUpper
		    );

void getLpMatrixCoefficient(
			    const HighsLp& lp,
			    int row,
			    int col,
			    double *val
			    );

int add_lp_cols(
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
		const bool force = false
		);

int augment_lp_cols(
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
		const bool force = false
		);

int validate_col_bounds(
			int XnumCol,
			const double* XcolLower,
			const double* XcolUpper,
			double infinite_bound
			);

void add_cols_to_lp_vectors(
			    HighsLp &lp,
			    int XnumNewCol,
			    const double *XcolCost,
			    const double *colLower,
			    const double *XcolUpper
			    );

int filter_col_bounds(
		       HighsLp& lp,
		       int XfromCol,
		       int XtoCol,
		       double infinite_bound
		       );

int add_lp_rows(
		HighsLp& lp,
		int XnumNewRow,
		const double *XrowLower,
		const double *XrowUpper,
		int XnumNewNZ,
		const int *XARstart,
		const int *XARindex,
		const double *XARvalue,
		const HighsOptions& options,
		const bool force = false
		);

int augment_lp_rows(
		    HighsLp& lp,
		    int XnumNewRow,
		    const double *XrowLower,
		    const double *XrowUpper,
		    int XnumNewNZ,
		    const int *XARstart,
		    const int *XARindex,
		    const double *XARvalue,
		    const HighsOptions& options,
		    const bool force = false
		    );

int validate_row_bounds(
			int XnumRow,
			const double* XrowLower,
			const double* XrowUpper,
			double infinite_bound
			);

void add_rows_to_lp_vectors(HighsLp &lp,
			    int XnumNewRow,
			    const double *XrowLower,
			    const double *XrowUpper
			    );

int filter_row_bounds(
		       HighsLp& lp,
		       int XfromRow,
		       int XtoRow,
		       double infinite_bound
		       );

int validate_matrix_indices(
		     int XnumRow,
		     int XnumCol,
		     int XnumNZ,
		     const int* XAstart,
		     const int* XAindex);

void add_cols_to_lp_matrix(
			   HighsLp &lp,
			   int XnumNewCol,
			   int XnumNewNZ,
			   const int *XAstart,
			   const int *XAindex,
			   const double *XAvalue
			   );

void add_rows_to_lp_matrix(HighsLp &lp,
			   int XnumNewRow,
			   int XnumNewNZ,
			   const int *XARstart,
			   const int *XARindex,
			   const double *XARvalue
			   );

int filter_matrix_entries(
			  HighsLp& lp,
			  int XfromCol,
			  int XtoCol,
			  double small_matrix_value
			  );

int filter_row_matrix_entries(
			      int XnumCol,
			      int XnumRow,
			      int XnumNZ,
			      int* XARstart,
			      int* XARindex,
			      double* XARvalue,
			      double small_matrix_value
			      );


/**
 * @brief Report the data of an LP
 */
void reportLp(
	      const HighsLp &lp //!< LP whose data are to be reported
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
#endif // LP_DATA_HIGHSLPUTILS_H_
