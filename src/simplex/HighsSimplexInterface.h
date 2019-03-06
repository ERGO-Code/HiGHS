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
   * @brief Add a contiguous set of columns to the model data---making them nonbasic
   */
  HighsStatus util_add_cols(
			    int XnumCol,
			    const double *XcolCost,
			    const double *XcolLower,
			    const double *XcolUpper,
			    int XnumNZ,
			    const int *XAstart,
			    const int *XAindex,
			    const double *XAvalue
			    );
  
  HighsStatus delete_cols(
			  int from_col,
			  int to_col
			  );
  

  HighsStatus delete_cols(
			  int num_set_entries,
			  const int* col_set
			  );
  
  HighsStatus delete_cols(
			  const int* col_mask
			  );
  
  HighsStatus delete_cols_general(
				  bool interval,
				  int from_col,
				  int to_col,
				  bool set,
				  int num_set_entries,
				  const int* col_set,
				  bool mask,
				  const int* col_mask
				  );
  
  HighsStatus getCols(
		      const int from_col,
		      const int to_col,
		      int num_col,
		      double *col_cost,
		      double *col_lower,
		      double *col_upper,
		      int num_nz,
		      int *col_matrix_start,
		      int *col_matrix_index,
		      double *col_matrix_value
		      );
  
  
  HighsStatus getCols(
		      const int num_set_entries,
		      const int* col_set,
		      int num_col,
		      double *col_cost,
		      double *col_lower,
		      double *col_upper,
		      int num_nz,
		      int *col_matrix_start,
		      int *col_matrix_index,
		      double *col_matrix_value
		      );
  
  HighsStatus getCols(
		      const int* col_mask,
		      int num_col,
		      double *col_cost,
		      double *col_lower,
		      double *col_upper,
		      int num_nz,
		      int *col_matrix_start,
		      int *col_matrix_index,
		      double *col_matrix_value
		      );
  
  HighsStatus getColsGeneral(
			     const bool interval,
			     const int from_col,
			     const int to_col,
			     const bool set,
			     const int num_set_entries,
			     const int* col_set,
			     const bool mask,
			     const int* col_mask,
			     int num_col,
			     double *col_cost,
			     double *col_lower,
			     double *col_upper,
			     int num_nz,
			     int *col_matrix_start,
			     int *col_matrix_index,
			     double *col_matrix_value
			     );
  
  HighsStatus getRows(
		      const int from_row,
		      const int to_row,
		      int num_row,
		      double *row_lower,
		      double *row_upper,
		      int num_nz,
		      int *row_matrix_start,
		      int *row_matrix_index,
		      double *row_matrix_value
		      );
  
  
  HighsStatus getRows(
		      const int num_set_entries,
		      const int* row_set,
		      int num_row,
		      double *row_lower,
		      double *row_upper,
		      int num_nz,
		      int *row_matrix_start,
		      int *row_matrix_index,
		      double *row_matrix_value
		      );
  
  HighsStatus getRows(
		      const int* row_mask,
		      int num_row,
		      double *row_lower,
		      double *row_upper,
		      int num_nz,
		      int *row_matrix_start,
		      int *row_matrix_index,
		      double *row_matrix_value
		      );
  
 HighsStatus getRowsGeneral(
			    const bool interval,
			    const int from_row,
			    const int to_row,
			    const bool set,
			    const int num_set_entries,
			    const int* row_set,
			    const bool mask,
			    const int* row_mask,
			    int num_row,
			    double *row_lower,
			    double *row_upper,
			    int num_nz,
			    int *row_matrix_start,
			    int *row_matrix_index,
			    double *row_matrix_value
			    );


  /**
   * @brief Add a contiguous set of rows to the model data---making them basic
   */
  HighsStatus util_add_rows(
			    int XnumNewRow,
			    const double *XrowLower,
			    const double *XrowUpper,
			    int XnumNewNZ,
			    const int *XARstart,
			    const int *XARindex,
			    const double *XARvalue
			    );

  HighsStatus delete_rows(
			  int from_row,
			  int to_row
			  );
  

  HighsStatus delete_rows(
			  int num_set_entries,
			  const int* row_set
			  );
  
  HighsStatus delete_rows(
			  const int* row_mask
			  );
  
  HighsStatus delete_rows_general(
				  bool interval,
				  int from_row,
				  int to_row,
				  bool set,
				  int num_set_entries,
				  const int* row_set,
				  bool mask,
				  const int* row_mask
				  );
  
  HighsStatus util_change_coefficient(
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
  HighsStatus change_ObjSense(
			      int Xsense
			      );

// Change the costs for an interval of columns
  HighsStatus change_costs(
			  int from_col,
			  int to_col,
			  const double* usr_col_cost
			  );
  

// Change the costs from an ordered set of indices
  HighsStatus change_costs(
			   int num_set_entries,
			   const int* col_set,
			   const double* usr_col_cost
			   );
  
// Change the costs with a mask
  HighsStatus change_costs(
			   const int* col_mask,
			   const double* usr_col_cost
			   );
  
  HighsStatus change_costs_general(
				  bool interval,
				  int from_col,
				  int to_col,
				  bool set,
				  int num_set_entries,
				  const int* col_set,
				  bool mask,
				  const int* col_mask,
				  const double* usr_col_cost
				  );
  
// Change the bounds for an interval of columns
  HighsStatus change_col_bounds(
			  int from_col,
			  int to_col,
			  const double* usr_col_lower,
			  const double* usr_col_upper
			  );
  

// Change the bounds from an ordered set of indices
  HighsStatus change_col_bounds(
			   int num_set_entries,
			   const int* col_set,
			   const double* usr_col_lower,
			   const double* usr_col_upper
			   );
  
// Change the bounds with a mask
  HighsStatus change_col_bounds(
			   const int* col_mask,
			   const double* usr_col_lower,
			   const double* usr_col_upper
			   );
  
  HighsStatus change_col_bounds_general(
				  bool interval,
				  int from_col,
				  int to_col,
				  bool set,
				  int num_set_entries,
				  const int* col_set,
				  bool mask,
				  const int* col_mask,
				  const double* usr_col_lower,
				  const double* usr_col_upper
				  );

// Change the bounds for an interval of rows
  HighsStatus change_row_bounds(
			  int from_row,
			  int to_row,
			  const double* usr_row_lower,
			  const double* usr_row_upper
			  );
  

// Change the bounds from an ordered set of indices
  HighsStatus change_row_bounds(
			   int num_set_entries,
			   const int* row_set,
			   const double* usr_row_lower,
			   const double* usr_row_upper
			   );
  
// Change the bounds with a mask
  HighsStatus change_row_bounds(
			   const int* row_mask,
			   const double* usr_row_lower,
			   const double* usr_row_upper
			   );
  
  HighsStatus change_row_bounds_general(
				  bool interval,
				  int from_row,
				  int to_row,
				  bool set,
				  int num_set_entries,
				  const int* row_set,
				  bool mask,
				  const int* row_mask,
				  const double* usr_row_lower,
				  const double* usr_row_upper
				  );

#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void change_update_method(
			    int updateMethod
			    );
#endif

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


};

#endif /* SIMPLEX_HIGHSSIMPLEXINTERFACE_H_ */


