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
#include "HighsLp.h"
#include "HighsStatus.h"

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

int validate_col_bounds(
			int XnumCol,
			const double* XcolLower,
			const double* XcolUpper
			);

void filter_col_bounds(
		       HighsLp& lp,
		       int XfromCol,
		       int XtoCol,
		       const double infinite_bound
		       );

int validate_row_bounds(
			int XnumRow,
			const double* XrowLower,
			const double* XrowUpper
			);

void filter_row_bounds(
		       HighsLp& lp,
		       int XfromRow,
		       int XtoRow,
		       const double infinite_bound
		       );

int validate_matrix(
		     int XnumRow,
		     int XnumCol,
		     int XnumNZ,
		     const int* XAstart,
		     const int* XAindex,
		     const double* XAvalue
		     );

void filter_matrix_values(
			  HighsLp& lp,
			  int XfromCol,
			  int XtoCol,
			  double small_matrix_value
			  );

void filter_row_matrix_values(
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
void util_analyseLp(const HighsLp &lp, const char* message);
#endif
#endif // LP_DATA_HIGHSLPUTILS_H_
