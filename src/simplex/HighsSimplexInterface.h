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
 * @brief Return or report data from simplex solves and interface that data with changes to models
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSIMPLEXINTERFACE_H_
#define SIMPLEX_HIGHSSIMPLEXINTERFACE_H_

#include <vector>

#include "HConfig.h"
#include "simplex/HVector.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HSimplex.h"

//class HFactor;

/**
 * @brief Return or report data from simplex solves and interface that data with changes to models
 */
class HighsSimplexInterface {
 public:
 HighsSimplexInterface(HighsModelObject& highs_model_object) : highs_model_object(highs_model_object) {}
  
  HighsModelObject &highs_model_object;
  
  /**
   * @brief Report the outcome of a simplex solve, printing a message first to contextualise the call
   */
  void report_simplex_outcome(
			      const char* message
			      );
  /**
   * @brief Compute the LP objective function value from column values
   */
  double get_lp_objective_value(
				vector<double> &XcolValue
				);

  /**
   * @brief Get vectors of column and row (primal) values and dual (values)
   */
  void get_primal_dual_values(
			      vector<double> &XcolValue, //!> Column primal activities
			      vector<double> &XcolDual,  //!> Column dual activities
			      vector<double> &XrowValue, //!> Row primal activities
			      vector<double> &XrowDual   //!> Row dual activities
			      );

  /**
   * @brief Get the basicIndex and nonbasicFlag vectors - Used?
   */
  void get_basicIndex_nonbasicFlag(
				   vector<int> &XbasicIndex,  //!> Indices of basic variables
				   vector<int> &XnonbasicFlag //!> Flag to indicate which variables are nonbasic
				   );

  /**
   * @brief Get the indices of the basic variables for SCIP
   */
  int get_basic_indices(
			int *bind //!> Indices of basic variables
			);

  /**
   * @brief Convert a SCIP baseStat for columns and rows to HiGHS basis
   * Postive  return value k implies invalid basis status for column k-1
   * Negative return value k implies invalid basis status for row   -k-1
   */
  int convert_baseStat_to_working(
				  const int* cstat, //!> Column baseStat
				  const int* rstat  //!> Row baseStat
				  );

  /**
   * @brief Convert a HiGHS basis to SCIP baseStat for columns and rows
   * Postive  return value k implies invalid basis status for column k-1
   * Negative return value k implies invalid basis status for row   -k-1
   */
  int convert_Working_to_BaseStat(
				  int* cstat, //!> Column baseStat
				  int* rstat  //!> Row baseStat
				  );

#ifdef HiGHSDEV
  /**
   * @brief Check that what's passed from postsolve is valid - Used?
   */
  void check_load_from_postsolve();
#endif

  /**
   * @brief Add a contiguous set of columns to the model data---making them nonbasic
   */
  int util_add_cols(
		    int XnumCol,
		    const double *XcolCost,
		    const double *XcolLower,
		    const double *XcolUpper,
		    int XnumNZ,
		    const int *XAstart,
		    const int *XAindex,
		    const double *XAvalue,
		    const bool force = false
		    );
  
  void util_delete_cols(
			int XfromCol,
			int XtoCol
			);
  
  void util_delete_col_set(
			   vector<int>& dstat
			   );
  
  void util_extract_cols(
			 int XfromCol, int XtoCol,
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
  int util_add_rows(
		     int XnumNewRow,
		     const double *XrowLower,
		     const double *XrowUpper,
		     int XnumNewNZ,
		     const int *XARstart,
		     const int *XARindex,
		     const double *XARvalue,
		     const bool force = false
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
			       int Xrow,
			       int Xcol,
			       const double XnewValue
			       );

  // Shift the objective
  void shift_objective_value(
			     double Xshift
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
		       int XnumColInSet,
		       const int* XcolCostIndex,
                       const double* XcolCostValue
		       );

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
  int change_col_bounds_all(
			    const double* XcolLower,
			    const double* XcolUpper,
			    bool force = false
			    );

// Change the bounds for a set of columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
  int change_col_bounds_set(
			    int ncols,
			    const int* XcolBoundIndex,
			    const double* XcolLowerValues,
			    const double* XcolUpperValues,
			    bool force = false
			    );

// Change the bounds for all rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
  int change_row_bounds_all(
			    const double* XrowLower,
			    const double* XrowUpper,
			    bool force = false
			    );

// Change the bounds for a set of rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
  int change_row_bounds_set(
			    int nrows,
			    const int* XrowBoundIndex,
			    const double* XrowLowerValues,
			    const double* XrowUpperValues,
			    bool force = false
			    );
  
#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void change_update_method(
			    int updateMethod
			    );
#endif
};

#endif /* SIMPLEX_HIGHSSIMPLEXINTERFACE_H_ */
