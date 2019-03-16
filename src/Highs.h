#ifndef HIGHS_H_
#define HIGHS_H_

#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsTimer.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"
#include "lp_data/HighsModelBuilder.h"

// Class to set parameters and run HiGHS
class Highs
{
public:
  // see if an empty lp should have Astart[0] = 0
  Highs() {}
  Highs(HighsOptions& options) { 
    options_ = options;
  }

  // The public method run() calls runSolver to solve problem before
  // or after presolve (or crash later) depending on the specified options.
  HighsStatus initializeLp(
			   const HighsLp &lp
			   );
  HighsStatus run();

  const HighsLp& getLp() const;
  const HighsSolution& getSolution() const;
  // todo: rename to HighsBasis when the current HighsBasis
  // becomes SimplexBasis
  const HighsBasis_new& getBasis() const; 

  double getRowValue(
		     const int row
		     ) const;
  double getObjectiveValue() const;
  const int getIterationCount() const;
  // todo: getRangingInformation(..)

  // In the solution passed as a parameter below can have one or many of
  // col_value, col_dual and row_dual set. If any of them are not set the
  // solution in Highs does not get updated.
  HighsStatus setSolution(
			  const HighsSolution& solution
			  );

  HighsStatus setBasis(
		       const HighsBasis_new& basis
		       ) {
    basis_ = basis;
    return HighsStatus::OK;
  }

   /**
   * @brief Adds a row to the model
   */
  bool addRow(const double lower_bound, //!< lower bound of the row
              const double upper_bound, //!< upper bound of the row
              const int num_new_nz, //!< number of nonzeros in the row
              const int *columns,  //!< array of size num_new_nz with column indices
              const double *values //!< array of size num_new_nz with coefficients
              );

   /**
   * @brief Adds multiple rows to the model
   */
  bool addRows(const int num_new_rows, //!< number of new rows
	       const double *lower_bounds,  //!< array of size num_new_rows with lower bounds
	       const double *upper_bounds, //!< array of size num_new_rows with upper bounds
	       const int *row_starts, //!< array of size num_new_rows+1 with start indices of the rows
	       const int num_new_nz, //!< number of total new nonzeros
	       const int *columns,  //!< array of size num_new_nz with column indices for all rows
	       const double *values //!< array of size num_new_nz with coefficients for all rows
	       );

  bool addCol(const double cost,
              const double lower_bound, 
              const double upper_bound,
              const int num_new_nz,
              const int *rows, 
              const double *values
	      );

  bool addCols(const int num_new_rows, 
	       const double* column_costs,
	       const double *lower_bounds, 
	       const double *upper_bounds,
	       const int *col_starts,
	       const int num_new_nz,
	       const int *rows, 
	       const double *values
	       );

  bool changeObjectiveSense(
			    int sense
			    );

  bool changeColCost(
		     int index, 
		     double coef
		     );

  bool changeColsCost(
		      int from_col,
		      int to_col, 
		      double* coef
		      );

  bool changeColsCost(
		      int n, 
		      int* index, 
		      double* coef
		      );

  bool changeColsCost(
		      int* mask, 
		      double* coef
		      );

  bool changeColBounds(
		       int index, 
		       double lower, 
		       double higher
		       );

  bool changeColsBounds(
			int n, 
			int* index, 
			double* lower, 
			double* higher
			);

  bool changeColsBounds(
			int* mask, 
			double* lower, 
			double* higher
			);

  bool changeRowBounds(
		       int index, 
		       double lower, 
		       double higher
		       );
  
  bool changeRowsBounds(
			int n, 
			int* index, 
			double* lower, 
			double* higher
			);

  bool changeRowsBounds(
			int* mask, 
			double* lower, 
			double* higher
			);

  bool getCols(
	       const int from_col,
	       const int to_col,
	       int& num_col,
	       double *col_costs,
	       double *col_lower,
	       double *col_upper,
	       int& num_nz,
	       int *col_matrix_start,
	       int *col_matrix_index,
	       double *col_matrix_value
	       );

  bool getCols(
	       const int n, 
	       const int* indices,
	       int& num_col,
	       double *col_costs,
	       double *col_lower,
	       double *col_upper,
	       int& num_nz,
	       int *col_matrix_start,
	       int *col_matrix_index,
	       double *col_matrix_value
	       );


  bool getCols(
	       const int* col_mask,
	       int& num_col,
	       double *col_costs,
	       double *col_lower,
	       double *col_upper,
	       int& num_nz,
	       int *col_matrix_start,
	       int *col_matrix_index,
	       double *col_matrix_value
	       );

  bool getRows(
	       const int from_row,
	       const int to_row,
	       int& num_row,
	       double *row_lower,
	       double *row_upper,
	       int& num_nz,
	       int *row_matrix_start,
	       int *row_matrix_index,
	       double *row_matrix_value
	       );


  bool getRows(
	       const int n, 
	       const int* indices,
	       int& num_row,
	       double *row_lower,
	       double *row_upper,
	       int& num_nz,
	       int *row_matrix_start,
	       int *row_matrix_index,
	       double *row_matrix_value
	       );


  bool getRows(
	       const int* row_mask,
	       int& num_row,
	       double *row_lower,
	       double *row_upper,
	       int& num_nz,
	       int *row_matrix_start,
	       int *row_matrix_index,
	       double *row_matrix_value
	       );

  bool deleteCols(
		  const int from_col,
		  const int to_col
		  );

  bool deleteCols(
		  const int n, 
		  const int* indices
		  );

  bool deleteCols(
		  int* mask
		  );

  bool deleteRows(
		  const int from_row, 
		  const int to_row
		  );

  bool deleteRows(
		  const int n, 
		  const int* indices
		  );

  bool deleteRows(
		  int* mask
		  );

  // change coeff (int row, int col) | ...
  // ipx (not implemented)


  // todo: Set warm/hot start methods

#ifdef OSI_FOUND
  friend class OsiHiGHSSolverInterface;
#endif

  HighsOptions options_;

private:
  HighsSolution solution_;
  HighsBasis_new basis_;
  HighsLp lp_;

  HighsTimer timer_;

  // Each HighsModelObject holds a const ref to its lp_. There is potentially
  // several hmos_ to allow for the solution of several different modified
  // versions of the original LP for instance different levels of presolve.
  std::vector<HighsModelObject> hmos_;

  bool simplex_has_run_;

  HighsPresolveStatus runPresolve(PresolveInfo &presolve_info);
  HighsPostsolveStatus runPostsolve(PresolveInfo &presolve_info);
  HighsStatus runSolver(HighsModelObject &model);

  HighsStatus runBnb();
  HighsStatus solveRootNode();
  HighsStatus solveNode();

};

#endif
