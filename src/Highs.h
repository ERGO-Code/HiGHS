/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.h
 * @brief The HiGHS class
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef HIGHS_H_
#define HIGHS_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsModelBuilder.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"
#include "mip/SolveMip.h"
#include "util/HighsTimer.h"

/**
 * @brief Class to set parameters and run HiGHS
 */
class Highs {
 public:
  // see if an empty lp should have Astart[0] = 0
  Highs();
  Highs(HighsOptions& options) { options_ = options; }

  HighsStatus setHighsOptionValue(const std::string& option,
                                  const std::string& value) {
    OptionStatus status = setOptionValue(options_, option, value);
    if (status != OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }

  /**
   * @brief Clears the vector of HighsModelObjects (hmos), creates a
   * HighsModelObject for this LP and makes it the first of the vector
   * of HighsModelObjects
   */
  HighsStatus initializeLp(
      const HighsLp& lp  //!< The HighsLp instance for this LP
  );

  /**
   * @brief reads a model from a file and initializes the Highs object
   */
  HighsStatus initializeFromFile(const std::string filename  //!< the filename
  );

  /**
   * @brief writes the current model to a file
   */
  HighsStatus writeToFile(const std::string filename  //!< the filename
  );

  /**
   * @brief Calls runSolver to solve the LP according to the
   * specified options
   */
  HighsStatus run();

  /**
   * @brief Returns the HighsLp instance for the LP of the (first?)
   * HighsModelObject
   */
  const HighsLp& getLp() const;

  /**
   * @brief Returns the HighsSolution instance for the LP of the
   * (first?)  HighsModelObject
   */
  const HighsSolution& getSolution() const;

  /**
   * @brief Returns the HighsBasis instance for the LP of the
   * (first?) HighsModelObject TODO: rename to HighsBasis when the
   * current HighsBasis becomes SimplexBasis
   */
  const HighsBasis& getBasis() const;

  /**
   * @brief Returns the (dual) objective function value for the LP of
   * the (first?) HighsModelObject
   */
  double getObjectiveValue() const;

  /**
   * @brief Returns the number of simplex iterations for the LP of the
   * (first?) HighsModelObject
   */
  int getIterationCount() const;
  // todo: getRangingInformation(..)

  /**
   * @brief Uses the HighsSolution passed to set the solution for the
   * LP of the (first?) HighsModelObject, according to ?? TODO Understand this
   * ??
   */
  // In the solution passed as a parameter below can have one or many of
  // col_value, col_dual and row_dual set. If any of them are not set the
  // solution in Highs does not get updated.
  HighsStatus setSolution(const HighsSolution& solution);

  /**
   * @brief Uses the HighsBasis passed to set the basis for the
   * LP of the (first?) HighsModelObject
   */
  HighsStatus setBasis(const HighsBasis& basis);

  /**
   * @brief Reports the solution and basis status for the LP of the
   * (first?) HighsModelObject
   */
  void reportSolution();

  /**
   * @brief Adds a row to the model
   */
  bool addRow(
      const double lower,    //!< Lower bound of the row
      const double upper,    //!< Upper bound of the row
      const int num_new_nz,  //!< Number of nonzeros in the row
      const int* indices,    //!< Array of size num_new_nz with column indices
      const double* values   //!< Array of size num_new_nz with column values
  );

  /**
   * @brief Adds multiple rows to the model
   */
  bool addRows(
      const int num_new_row,  //!< Number of new rows
      const double* lower,    //!< Array of size num_new_row with lower bounds
      const double* upper,    //!< Array of size num_new_row with upper bounds
      const int num_new_nz,   //!< Number of new nonzeros
      const int*
          starts,  //!< Array of size num_new_row with start indices of the rows
      const int* indices,  //!< Array of size num_new_nz with column indices for
                           //!< all rows
      const double*
          values  //!< Array of size num_new_nz with column values for all rows
  );

  /**
   * @brief Adds a column to the model
   */
  bool addCol(
      const double cost,     //!< Cost of the column
      const double lower,    //!< Lower bound of the column
      const double upper,    //!< Upper bound of the column
      const int num_new_nz,  //!< Number of nonzeros in the column
      const int* indices,    //!< Array of size num_new_nz with row indices
      const double* values   //!< Array of size num_new_nz with row values
  );

  /**
   * @brief Adds multiple columns to the model
   */
  bool addCols(
      const int num_new_col,  //!< Number of new columns
      const double* costs,    //!< Array of size num_new_col with costs
      const double* lower,    //!< Array of size num_new_col with lower bounds
      const double* upper,    //!< Array of size num_new_col with upper bounds
      const int num_new_nz,   //!< Number of new nonzeros
      const int* starts,   //!< Array of size num_new_row with start indices of
                           //!< the columns
      const int* indices,  //!< Array of size num_new_nz with row indices for
                           //!< all columns
      const double*
          values  //!< Array of size num_new_nz with row values for all columns
  );

  /**
   * @brief Change the objective sense of the model
   */
  bool changeObjectiveSense(const int sense  //!< New objective sense
  );

  /**
   * @brief Change the cost of a column
   */
  bool changeColCost(
      const int col,     //!< The index of the column whose cost is to change
      const double cost  //!< The new cost
  );

  /**
   * @brief Change the cost of multiple columns given by a set of indices
   */
  bool changeColsCost(
      const int num_set_entries,  //!< The number of indides in the set
      const int* set,     //!< Array of size num_set_entries with indices of
                          //!< columns whose costs change
      const double* cost  //!< Array of size num_set_entries with new costs
  );

  /**
   * @brief Change the cost of multiple columns given by a mask
   */
  bool changeColsCost(
      const int* mask,    //!< Full length array with 1 => change; 0 => not
      const double* cost  //!< Full length array of new costs
  );

  /**
   * @brief Change the bounds of a column
   */
  bool changeColBounds(
      const int col,  //!< The index of the column whose bounds are to change
      const double lower,  //!< The new lower bound
      const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple columns given by an interval
   */
  bool changeColsBounds(
      const int
          from_col,      //!< The index of the first column whose bounds change
      const int to_col,  //!< One more than the index of the last column whose
                         //!< bounds change
      const double*
          lower,  //!< Array of size to_col-from_col with new lower bounds
      const double*
          upper  //!< Array of size to_col-from_col with new upper bounds
  );

  /**
   * @brief Change the bounds of multiple columns given by a set of indices
   */
  bool changeColsBounds(
      const int num_set_entries,  //!< The number of indides in the set
      const int* set,  //!< Array of size num_set_entries with indices of
                       //!< columns whose bounds change
      const double*
          lower,  //!< Array of size num_set_entries with new lower bounds
      const double*
          upper  //!< Array of size num_set_entries with new upper bounds
  );

  /**
   * @brief Change the cost of multiple columns given by a mask
   */
  bool changeColsBounds(
      const int* mask,      //!< Full length array with 1 => change; 0 => not
      const double* lower,  //!< Full length array of new lower bounds
      const double* upper   //!< Full length array of new upper bounds
  );

  /**
   * @brief Change the bounds of a row
   */
  bool changeRowBounds(
      const int row,       //!< The index of the row whose bounds are to change
      const double lower,  //!< The new lower bound
      const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple rows given by a set of indices
   */
  bool changeRowsBounds(
      const int num_set_entries,  //!< The number of indides in the set
      const int* set,  //!< Array of size num_set_entries with indices of rows
                       //!< whose bounds change
      const double*
          lower,  //!< Array of size num_set_entries with new lower bounds
      const double*
          upper  //!< Array of size num_set_entries with new upper bounds
  );

  /**
   * @brief Change the cost of multiple rows given by a mask
   */
  bool changeRowsBounds(
      const int* mask,      //!< Full length array with 1 => change; 0 => not
      const double* lower,  //!< Full length array of new lower bounds
      const double* upper   //!< Full length array of new upper bounds
  );

  /**
   * @brief Get multiple columns from the model given by an interval
   */
  bool getCols(const int from_col,  //!< The index of the first column to get
                                    //!< from the model
               const int to_col,  //!< One more than the last column to get from
                                  //!< the model
               int& num_col,      //!< Number of columns got from the model
               double* costs,     //!< Array of size num_col with costs
               double* lower,     //!< Array of size num_col with lower bounds
               double* upper,     //!< Array of size num_col with upper bounds
               int& num_nz,       //!< Number of nonzeros got from the model
               int* matrix_start,  //!< Array of size num_col with start indices
                                   //!< of the columns
               int* matrix_index,  //!< Array of size num_nz with row indices
                                   //!< for the columns
               double* matrix_value  //!< Array of size num_nz with row values
                                     //!< for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a set
   */
  bool getCols(const int num_set_entries,  //!< The number of indides in the set
               const int* set,  //!< Array of size num_set_entries with indices
                                //!< of columns to get
               int& num_col,    //!< Number of columns got from the model
               double* costs,   //!< Array of size num_col with costs
               double* lower,   //!< Array of size num_col with lower bounds
               double* upper,   //!< Array of size num_col with upper bounds
               int& num_nz,     //!< Number of nonzeros got from the model
               int* matrix_start,  //!< Array of size num_col with start indices
                                   //!< of the columns
               int* matrix_index,  //!< Array of size num_nz with row indices
                                   //!< for the columns
               double* matrix_value  //!< Array of size num_nz with row values
                                     //!< for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a mask
   */
  bool getCols(const int* mask,  //!< Full length array with 1 => get; 0 => not
               int& num_col,     //!< Number of columns got from the model
               double* costs,    //!< Array of size num_col with costs
               double* lower,    //!< Array of size num_col with lower bounds
               double* upper,    //!< Array of size num_col with upper bounds
               int& num_nz,      //!< Number of nonzeros got from the model
               int* matrix_start,    //!<  Array of size num_col with start
                                     //!<  indices of the columns
               int* matrix_index,    //!<  Array of size num_nz with row indices
                                     //!<  for the columns
               double* matrix_value  //!<  Array of size num_nz with row values
                                     //!<  for the columns
  );

  /**
   * @brief Get multiple rows from the model given by an interval
   */
  bool getRows(
      const int from_row,  //!< The index of the first row to get from the model
      const int to_row,    //!< One more than the last row get from the model
      int& num_row,        //!< Number of rows got from the model
      double* lower,       //!< Array of size num_row with lower bounds
      double* upper,       //!< Array of size num_row with upper bounds
      int& num_nz,         //!< Number of nonzeros got from the model
      int* matrix_start,   //!< Array of size num_row with start indices of the
                           //!< rows
      int* matrix_index,   //!< Array of size num_nz with column indices for the
                           //!< rows
      double* matrix_value  //!< Array of size num_nz with column values for the
                            //!< rows
  );

  /**
   * @brief Get multiple rows from the model given by a set
   */
  bool getRows(const int num_set_entries,  //!< The number of indides in the set
               const int* set,  //!< Array of size num_set_entries with indices
                                //!< of rows to get
               int& num_row,    //!< Number of rows got from the model
               double* lower,   //!< Array of size num_row with lower bounds
               double* upper,   //!< Array of size num_row with upper bounds
               int& num_nz,     //!< Number of nonzeros got from the model
               int* matrix_start,  //!< Array of size num_row with start indices
                                   //!< of the rows
               int* matrix_index,  //!< Array of size num_nz with column indices
                                   //!< for the rows
               double* matrix_value  //!< Array of size num_nz with column
                                     //!< values for the rows
  );

  /**
   * @brief Get multiple rows from the model given by a mask
   */
  bool getRows(const int* mask,  //!< Full length array with 1 => get; 0 => not
               int& num_row,     //!< Number of rows got from the model
               double* lower,    //!< Array of size num_row with lower bounds
               double* upper,    //!< Array of size num_row with upper bounds
               int& num_nz,      //!< Number of nonzeros got from the model
               int* matrix_start,  //!< Array of size num_row with start indices
                                   //!< of the rows
               int* matrix_index,  //!< Array of size num_nz with column indices
                                   //!< for the rows
               double* matrix_value  //!< Array of size num_nz with column
                                     //!< values for the rows
  );

  /**
   * @brief Delete multiple columns from the model given by an interval
   */
  bool deleteCols(const int from_col,  //!< The index of the first column to
                                       //!< delete from the model
                  const int to_col  //!< One more than the last column to delete
                                    //!< from the model
  );

  /**
   * @brief Delete multiple columns from the model given by a set
   */
  bool deleteCols(
      const int num_set_entries,  //!< The number of indides in the set
      const int* set  //!< Array of size num_set_entries with indices of columns
                      //!< to delete
  );

  /**
   * @brief Delete multiple columns from the model given by a mask
   */
  bool deleteCols(int* mask  //!< Full length array with 1 => delete; 0 => not
  );

  /**
   * @brief Delete multiple rows from the model given by an interval
   */
  bool deleteRows(
      const int
          from_row,     //!< The index of the first row to delete from the model
      const int to_row  //!< One more than the last row delete from the model
  );

  /**
   * @brief Delete multiple rows from the model given by a set
   */
  bool deleteRows(
      const int num_set_entries,  //!< The number of indides in the set
      const int* set  //!< Array of size num_set_entries with indices of columns
                      //!< to delete
  );

  /**
   * @brief Delete multiple rows from the model given by a mask
   */
  bool deleteRows(int* mask  //!< Full length array with 1 => delete; 0 => not
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
  HighsBasis basis_;
  HighsLp lp_;

  HighsTimer timer_;

  // Each HighsModelObject holds a const ref to its lp_. There are potentially
  // several hmos_ to allow for the solution of several different modified
  // versions of the original LP for instance different levels of presolve.
  std::vector<HighsModelObject> hmos_;

  bool simplex_has_run_;

  HighsStatus callRunSolver(HighsModelObject& model, int& iteration_count,
                            const string message);
  HighsStatus runSolver(HighsModelObject& model);

  HighsPresolveStatus runPresolve(PresolveInfo& presolve_info);
  HighsPostsolveStatus runPostsolve(PresolveInfo& presolve_info);

  HighsStatus runBnb();
  HighsStatus solveRootNode(Node& root);
  HighsStatus solveNode(Node& node);
};

#endif
