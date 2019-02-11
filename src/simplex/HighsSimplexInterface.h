/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexInterface.h
 * @brief Dual simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSIMPLEXINTERFACE_H_
#define SIMPLEX_HIGHSSIMPLEXINTERFACE_H_

#include <vector>

#include "HConfig.h"
#include "HVector.h"
#include "HighsModelObject.h"
#include "HSimplex.h"

//class HFactor;

/**
 * @brief Dual simplex solver for HiGHS
 */
class HighsSimplexInterface {
 public:
 HighsSimplexInterface(HighsModelObject& highs_model_object) : highs_model_object(highs_model_object) {}
  
  HighsModelObject &highs_model_object;
  
  /**
   * @brief Load a model from arrays 
   */
  void load_from_arrays(
			int XnumCol,
			int XobjSense,
			const double* XcolCost,
			const double* XcolLower,
			const double* XcolUpper,
			int XnumRow,
			const double* XrowLower,
			const double* XrowUpper,
			int XnumNz,
			const int* XAstart,
			const int* XAindex,
			const double* XAvalue
			);

  // Methods for brief reports
  void report_simplex_outcome(const char* message);
  void report_simplex_solution_status();

  /**
   * @brief Get the column and row (primal) values and dual (values)
   */
  void get_primal_dual_values(
			      vector<double> &XcolValue,
			      vector<double> &XcolDual,
			      vector<double> &XrowValue,
			      vector<double> &XrowDual
			      );

  void get_nonbasicMove(vector<int> &XnonbasicMove);

  void get_basicIndex_nonbasicFlag(
				   vector<int> &XbasicIndex,
				   vector<int> &XnonbasicFlag
				   );

  // Utility to get the indices of the basic variables for SCIP
  int get_basic_indices(
			int *bind
			);
  // Utilities to convert model basic/nonbasic status to/from SCIP-like status
  // Convert model basic/nonbasic status from SCIP-like status
  // Postive  return value k implies invalid basis status for column k-1
  // Negative return value k implies invalid basis status for row   -k-1
  int convert_baseStat_to_working(
				  const int* cstat,
				  const int* rstat
				  );
  int convert_Working_to_BaseStat(
				  int* cstat,
				  int* rstat
				  );



  /**
   * @brief Get the LP objective function value from column values
   */
  double get_lp_objective_value(
			      vector<double> &XcolValue
			      );

#ifdef HiGHSDEV
  /**
   * @brief Test load_from_arrays
   */
  void check_load_from_arrays();

  /**
   * @brief Check that what's passed from postsolve is valid
   */
  void check_load_from_postsolve();
#endif

  /**
   * @brief Add a contiguous set of columns to the model data---making them nonbasic
   */
  void util_add_cols(
		     int ncols,
		     const double *XcolCost,
		     const double *colLower,
		     const double *XcolUpper,
		     int nnonz,
		     const int *XAstart,
		     const int *XAindex,
		     const double *XAvalue
		     );

  void util_delete_cols(
			int firstcol,
			int lastcol
			);

  void util_delete_col_set(
			  vector<int>& dstat
			  );

  void util_extract_cols(
			 int firstcol, int lastcol,
			 double* XcolLower,
			 double* XcolUpper,
			 int* nnonz,
			 int* XAstart,
			 int* XAindex,
			 double* XAvalue
			 );


  /**
   * @brief Add a contiguous set of rows to the model data---making them basic
   */
  void util_add_rows(
		     int nrows,
		     const double *XrowLower,
		     const double *XrowUpper,
		     int nnonz,
		     const int *XARstart,
		     const int *XARindex,
		     const double *XARvalue
		     );
  void util_delete_rows(
			int firstrow,
			int lastrow
			);

  void util_delete_row_set(
			  vector<int>& dstat
			  );

  void util_extract_rows(
			 int firstrow,
			 int lastrow,
			 double* XrowLower,
			 double* XrowUpper,
			 int* nnonz,
			 int* XARstart,
			 int* XARindex,
			 double* XARvalue
			 );

  void util_change_coefficient(
			       int row,
			       int col,
			       const double newval
			       );

  // Shift the objective
  void shift_objective_value(
			     double shift
			     );

  // Utilities to get/change costs and bounds
  // Change the objective sense
  int change_ObjSense(
		      int Xsense
		      );

// Change the costs for all columns
  int change_costs_all(
		       const double* XcolCost
		       );

// Change the costs for a set of columns
  int change_costs_set(
		       int ncols,
		       const int* XcolCostIndex,
                       const double* XcolCostValues
		       );

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
  int change_col_bounds_all(
			    const double* XcolLower,
			    const double* XcolUpper
			    );

// Change the bounds for a set of columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
  int change_col_bounds_set(
			    int ncols,
			    const int* XcolBoundIndex,
			    const double* XcolLowerValues,
			    const double* XcolUpperValues
			    );

// Change the bounds for all rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
  int change_row_bounds_all(
			    const double* XrowLower,
			    const double* XrowUpper
			    );

// Change the bounds for a set of rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
  int change_row_bounds_set(
			    int nrows,
			    const int* XrowBoundIndex,
			    const double* XrowLowerValues,
			    const double* XrowUpperValues
			    );
  
  int write_to_mps(
		   const char* filename
		   );
#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void change_update_method(
			    int updateMethod
			    );
#endif
};

#endif /* SIMPLEX_HIGHSSIMPLEXINTERFACE_H_ */
