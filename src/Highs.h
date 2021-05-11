/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file Highs.h
 * @brief The HiGHS class
 */
#ifndef HIGHS_H_
#define HIGHS_H_

#include <sstream>

#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsRanging.h"
#include "lp_data/HighsSolutionDebug.h"
#include "presolve/PresolveComponent.h"

/**
 * @brief Class to set parameters and run HiGHS
 */
class Highs {
 public:
  Highs();
  virtual ~Highs() {}

  /**
   * Methods for model input
   */

  /**
   * @brief Every model loading module eventually uses passModel to
   * communicate the model to HiGHS. It clears the vector of
   * HighsModelObjects (hmos), creates a HighsModelObject for this LP
   * and makes it the first of the vector of HighsModelObjects
   */
  HighsStatus passModel(HighsLp lp  //!< The HighsLp instance for this LP
  );

  HighsStatus passModel(const HighsInt num_col, const HighsInt num_row,
                        const HighsInt num_nz, const bool rowwise,
                        const double* costs, const double* col_lower,
                        const double* col_upper, const double* row_lower,
                        const double* row_upper, const HighsInt* astart,
                        const HighsInt* aindex, const double* avalue,
                        const HighsInt* integrality = NULL);

  /**
   * @brief reads in a model and initializes the HighsModelObject
   */
  HighsStatus readModel(const std::string filename  //!< the filename
  );

  /**
   * @brief reads in a basis
   */
  HighsStatus readBasis(const std::string filename  //!< the filename
  );

  /**
   * @brief Clears the current model
   */
  HighsStatus clearModel();

  /**
   * @brief Solves the model according to the specified options
   */
  HighsStatus run();

  /**
   * @brief writes the current solution to a file
   */
  HighsStatus writeSolution(const std::string filename,  //!< the filename
                            const bool pretty = false)
      const;  //!< Write in pretty (human-readable) format

  /**
   * Methods for HiGHS option input/output
   */

  /**
   * @brief Sets an option to the bool/int/double/string  value if it's
   * legal and, for bool/int/double, only if it's of the correct type
   */

  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const bool value            //!< The option value
  );

  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const HighsInt value        //!< The option value
  );

#ifdef HIGHSINT64
  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const int value             //!< The option value
  ) {
    return setOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const double value          //!< The option value
  );

  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const std::string value     //!< The option value
  );

  HighsStatus setOptionValue(const std::string& option,  //!< The option name
                             const char* value           //!< The option value
  );

  HighsStatus readOptions(const std::string filename  //!< The filename
  );

  HighsStatus passOptions(const HighsOptions& options  //!< The options
  );

  const HighsOptions& getOptions();

  /**
   * @brief Gets an option value as bool/int/double/string and, for
   * bool/int/double, only if it's of the correct type.
   */
  HighsStatus getOptionValue(const std::string& option,  //!< The option name
                             bool& value                 //!< The option value
  );

  HighsStatus getOptionValue(const std::string& option,  //!< The option name
                             HighsInt& value             //!< The option value
  );

  HighsStatus getOptionValue(const std::string& option,  //!< The option name
                             double& value               //!< The option value
  );

  HighsStatus getOptionValue(const std::string& option,  //!< The option name
                             std::string& value          //!< The option value
  );

  /**
   * @brief Get the type expected by an option
   */
  HighsStatus getOptionType(const std::string& option,  //!< The option name
                            HighsOptionType& type       //!< The option type
  );

  const HighsOptions& getOptions() const;

  HighsStatus resetOptions();

  HighsStatus writeOptions(const std::string filename,  //!< The filename
                           const bool report_only_non_default_values = true);

  /**
   * @brief Gets an option value as int/double, and only if it's of the correct
   * type.
   */

  const HighsInfo& getInfo() const;

  HighsStatus getInfoValue(const std::string& info,  //!< The info name
                           HighsInt& value           //!< The info value
  );

  HighsStatus getInfoValue(const std::string& info,  //!< The info name
                           double& value) const;     //!< The info value

  HighsStatus writeInfo(const std::string filename  //!< The filename
  );
  /**
   * Methods for model output
   */

  /**
   * @brief Returns the HighsLp instance of the model in HiGHS
   */
  const HighsLp& getModel() const { return lp_; }

  /**
   * @brief Returns the HighsSolution
   */
  const HighsSolution& getSolution() const { return solution_; }

  /**
   * @brief Returns the HighsBasis
   */
  const HighsBasis& getBasis() const { return basis_; }

  /**
   * @brief Returns the current model status
   */
  const HighsModelStatus& getModelStatus(const bool scaled_model = false) const;

  /**
   * @brief Indicates whether a dual unbounded ray exdists, and gets
   * it if it does and dual_ray is not NULL
   */
  HighsStatus getDualRay(bool& has_dual_ray, double* dual_ray_value = NULL);

  /**
   * @brief Indicates whether a primal unbounded ray exdists, and gets
   * it if it does and primal_ray is not NULL
   */
  HighsStatus getPrimalRay(bool& has_primal_ray,
                           double* primal_ray_value = NULL);

  /**
   * @brief Gets the ranging information for the current LP
   */
  HighsStatus getRanging(HighsRanging& ranging);

  /**
   * Methods for operations with the invertible representation of the
   * current basis matrix
   */

  /**
   * @brief Gets the basic variables in the order corresponding to
   * calls to getBasisInverseRow, getBasisInverseCol, getBasisSolve,
   * getBasisTransposeSolve, getReducedRow and getReducedColumn. As
   * required by SCIP, non-negative entries are indices of columns,
   * and negative entries are -(row_index+1).
   */
  HighsStatus getBasicVariables(HighsInt* basic_variables  //!< Basic variables
  );
  /**
   * @brief Gets a row of \f$B^{-1}\f$ for basis matrix \f$B\f$
   */
  HighsStatus getBasisInverseRow(
      const HighsInt row,           //!< Index of row required
      double* row_vector,           //!< Row required
      HighsInt* row_num_nz = NULL,  //!< Number of nonzeros
      HighsInt* row_indices = NULL  //!< Indices of nonzeros
  );

  /**
   * @brief Gets a column of \f$B^{-1}\f$ for basis matrix \f$B\f$
   */
  HighsStatus getBasisInverseCol(
      const HighsInt col,           //!< Index of column required
      double* col_vector,           //!< Column required
      HighsInt* col_num_nz = NULL,  //!< Number of nonzeros
      HighsInt* col_indices = NULL  //!< Indices of nonzeros
  );

  /**
   * @brief Forms \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$
   */
  HighsStatus getBasisSolve(
      const double* rhs,                 //!< RHS \f$\mathbf{b}\f$
      double* solution_vector,           //!< Solution  \f$\mathbf{x}\f$
      HighsInt* solution_num_nz = NULL,  //!< Number of nonzeros
      HighsInt* solution_indices = NULL  //!< Indices of nonzeros
  );

  /**
   * @brief Forms \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$
   */
  HighsStatus getBasisTransposeSolve(
      const double* rhs,                 //!< RHS \f$\mathbf{b}\f$
      double* solution_vector,           //!< Solution  \f$\mathbf{x}\f$
      HighsInt* solution_nz = NULL,      //!< Number of nonzeros
      HighsInt* solution_indices = NULL  //!< Indices of nonzeros
  );

  /**
   * @brief Forms a row of \f$B^{-1}A\f$
   */
  HighsStatus getReducedRow(
      const HighsInt row,            //!< Index of row required
      double* row_vector,            //!< Row required
      HighsInt* row_num_nz = NULL,   //!< Number of nonzeros
      HighsInt* row_indices = NULL,  //!< Indices of nonzeros
      const double* pass_basis_inverse_row_vector =
          NULL  //!< Necessary row of \f$B^{-1}\f$
  );

  /**
   * @brief Forms a column of \f$B^{-1}A\f$
   */
  HighsStatus getReducedColumn(
      const HighsInt col,           //!< Index of column required
      double* col_vector,           //!< Column required
      HighsInt* col_num_nz = NULL,  //!< Number of nonzeros
      HighsInt* col_indices = NULL  //!< Indices of nonzeros
  );

  /**
   * @brief Get the number of columns in the LP of the (first?)
   * HighsModelObject
   */
  HighsInt getNumCols() const { return lp_.numCol_; }

  /**
   * @brief Get the number of rows in the LP of the (first?)
   * HighsModelObject
   */
  HighsInt getNumRows() const { return lp_.numRow_; }

  /**
   * @brief Get the number of entries in the LP of the (first?)
   * HighsModelObject
   */
  HighsInt getNumEntries() {
    if (lp_.numCol_) return lp_.Astart_[lp_.numCol_];
    return 0;
  }

  /**
   * @brief Get the objective sense of the model
   */
  bool getObjectiveSense(ObjSense& sense);

  /**
   * @brief Get multiple columns from the model given by an interval
   */
  bool getCols(const HighsInt from_col,  //!< The index of the first column to
                                         //!< get from the model
               const HighsInt to_col,  //!< One more than the last column to get
                                       //!< from the model
               HighsInt& num_col,      //!< Number of columns got from the model
               double* costs,          //!< Array of size num_col with costs
               double* lower,     //!< Array of size num_col with lower bounds
               double* upper,     //!< Array of size num_col with upper bounds
               HighsInt& num_nz,  //!< Number of nonzeros got from the model
               HighsInt* matrix_start,  //!< Array of size num_col with start
                                        //!< indices of the columns
               HighsInt* matrix_index,  //!< Array of size num_nz with row
                                        //!< indices for the columns
               double* matrix_value  //!< Array of size num_nz with row values
                                     //!< for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a set
   */
  bool getCols(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,     //!< Array of size num_set_entries with indices
                               //!< of columns to get
      HighsInt& num_col,       //!< Number of columns got from the model
      double* costs,           //!< Array of size num_col with costs
      double* lower,           //!< Array of size num_col with lower bounds
      double* upper,           //!< Array of size num_col with upper bounds
      HighsInt& num_nz,        //!< Number of nonzeros got from the model
      HighsInt* matrix_start,  //!< Array of size num_col with start indices
                               //!< of the columns
      HighsInt* matrix_index,  //!< Array of size num_nz with row indices
                               //!< for the columns
      double* matrix_value     //!< Array of size num_nz with row values
                               //!< for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a mask
   */
  bool getCols(
      const HighsInt* mask,    //!< Full length array with 1 => get; 0 => not
      HighsInt& num_col,       //!< Number of columns got from the model
      double* costs,           //!< Array of size num_col with costs
      double* lower,           //!< Array of size num_col with lower bounds
      double* upper,           //!< Array of size num_col with upper bounds
      HighsInt& num_nz,        //!< Number of nonzeros got from the model
      HighsInt* matrix_start,  //!<  Array of size num_col with start
                               //!<  indices of the columns
      HighsInt* matrix_index,  //!<  Array of size num_nz with row indices
                               //!<  for the columns
      double* matrix_value     //!<  Array of size num_nz with row values
                               //!<  for the columns
  );

  /**
   * @brief Get multiple rows from the model given by an interval
   */
  bool getRows(
      const HighsInt
          from_row,  //!< The index of the first row to get from the model
      const HighsInt to_row,  //!< One more than the last row get from the model
      HighsInt& num_row,      //!< Number of rows got from the model
      double* lower,          //!< Array of size num_row with lower bounds
      double* upper,          //!< Array of size num_row with upper bounds
      HighsInt& num_nz,       //!< Number of nonzeros got from the model
      HighsInt* matrix_start,  //!< Array of size num_row with start indices of
                               //!< the rows
      HighsInt* matrix_index,  //!< Array of size num_nz with column indices for
                               //!< the rows
      double* matrix_value  //!< Array of size num_nz with column values for the
                            //!< rows
  );

  /**
   * @brief Get multiple rows from the model given by a set
   */
  bool getRows(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,     //!< Array of size num_set_entries with indices
                               //!< of rows to get
      HighsInt& num_row,       //!< Number of rows got from the model
      double* lower,           //!< Array of size num_row with lower bounds
      double* upper,           //!< Array of size num_row with upper bounds
      HighsInt& num_nz,        //!< Number of nonzeros got from the model
      HighsInt* matrix_start,  //!< Array of size num_row with start indices
                               //!< of the rows
      HighsInt* matrix_index,  //!< Array of size num_nz with column indices
                               //!< for the rows
      double* matrix_value     //!< Array of size num_nz with column
                               //!< values for the rows
  );

  /**
   * @brief Get multiple rows from the model given by a mask
   */
  bool getRows(
      const HighsInt* mask,    //!< Full length array with 1 => get; 0 => not
      HighsInt& num_row,       //!< Number of rows got from the model
      double* lower,           //!< Array of size num_row with lower bounds
      double* upper,           //!< Array of size num_row with upper bounds
      HighsInt& num_nz,        //!< Number of nonzeros got from the model
      HighsInt* matrix_start,  //!< Array of size num_row with start indices
                               //!< of the rows
      HighsInt* matrix_index,  //!< Array of size num_nz with column indices
                               //!< for the rows
      double* matrix_value     //!< Array of size num_nz with column
                               //!< values for the rows
  );

  /**
   * @brief Get a matrix coefficient
   */
  bool getCoeff(const HighsInt row,  //!< Row of coefficient to be got
                const HighsInt col,  //!< Column of coefficient to be got
                double& value        //!< Coefficient
  );

  /**
   * @brief writes out current model
   */
  HighsStatus writeModel(const std::string filename  //!< the filename
  );

  /**
   * @brief writes out current basis
   */
  HighsStatus writeBasis(const std::string filename  //!< the filename
  );

  /**
   * Methods for model modification
   */

  /**
   * @brief Change the objective sense of the model
   */
  bool changeObjectiveSense(const ObjSense sense  //!< New objective sense
  );

  /**
   * @brief Change the integrality of a column
   */
  bool changeColIntegrality(
      const HighsInt
          col,  //!< The index of the column whose integrality is to change
      const HighsVarType integrality  //!< The new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by an interval
   */
  bool changeColsIntegrality(
      const HighsInt from_col,  //!< The index of the first column whose
                                //!< integrality changes
      const HighsInt to_col,    //!< One more than the index of the last column
                                //!< whose integrality changes
      const HighsVarType*
          integrality  //!< Array of size num_set_entries with new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by a set of indices
   */
  bool changeColsIntegrality(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,  //!< Array of size num_set_entries with indices of
                            //!< columns whose integralitys change
      const HighsVarType*
          integrality  //!< Array of size num_set_entries with new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by a mask
   */
  bool changeColsIntegrality(
      const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
      const HighsVarType* integrality  //!< Full length array of new integrality
  );

  /**
   * @brief Change the cost of a column
   */
  bool changeColCost(
      const HighsInt col,  //!< The index of the column whose cost is to change
      const double cost    //!< The new cost
  );

  /**
   * @brief Change the cost of multiple columns given by an interval
   */
  bool changeColsCost(
      const HighsInt
          from_col,  //!< The index of the first column whose cost changes
      const HighsInt to_col,  //!< One more than the index of the last column
                              //!< whose cost changes
      const double* cost      //!< Array of size num_set_entries with new costs
  );

  /**
   * @brief Change the cost of multiple columns given by a set of indices
   */
  bool changeColsCost(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,  //!< Array of size num_set_entries with indices of
                            //!< columns whose costs change
      const double* cost    //!< Array of size num_set_entries with new costs
  );

  /**
   * @brief Change the cost of multiple columns given by a mask
   */
  bool changeColsCost(
      const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
      const double* cost     //!< Full length array of new costs
  );

  /**
   * @brief Change the bounds of a column
   */
  bool changeColBounds(const HighsInt col,  //!< The index of the column whose
                                            //!< bounds are to change
                       const double lower,  //!< The new lower bound
                       const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple columns given by an interval
   */
  bool changeColsBounds(
      const HighsInt
          from_col,  //!< The index of the first column whose bounds change
      const HighsInt to_col,  //!< One more than the index of the last column
                              //!< whose bounds change
      const double*
          lower,  //!< Array of size to_col-from_col with new lower bounds
      const double*
          upper  //!< Array of size to_col-from_col with new upper bounds
  );

  /**
   * @brief Change the bounds of multiple columns given by a set of indices
   */
  bool changeColsBounds(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,  //!< Array of size num_set_entries with indices of
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
      const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
      const double* lower,   //!< Full length array of new lower bounds
      const double* upper    //!< Full length array of new upper bounds
  );

  /**
   * @brief Change the bounds of a row
   */
  bool changeRowBounds(
      const HighsInt row,  //!< The index of the row whose bounds are to change
      const double lower,  //!< The new lower bound
      const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple rows given by an interval
   */
  bool changeRowsBounds(const HighsInt from_row, const HighsInt to_row,
                        const double* lower, const double* upper);

  /**
   * @brief Change the bounds of multiple rows given by a set of indices
   */
  bool changeRowsBounds(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set,  //!< Array of size num_set_entries with indices of
                            //!< rows whose bounds change
      const double*
          lower,  //!< Array of size num_set_entries with new lower bounds
      const double*
          upper  //!< Array of size num_set_entries with new upper bounds
  );

  /**
   * @brief Change the cost of multiple rows given by a mask
   */
  bool changeRowsBounds(
      const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
      const double* lower,   //!< Full length array of new lower bounds
      const double* upper    //!< Full length array of new upper bounds
  );

  /**
   * @brief Change a matrix coefficient
   */
  bool changeCoeff(const HighsInt row,  //!< Row of coefficient to be changed
                   const HighsInt col,  //!< Column of coefficient to be changed
                   const double value   //!< Coefficient
  );
  /**
   * @brief Adds a row to the model
   */
  bool addRow(
      const double lower,         //!< Lower bound of the row
      const double upper,         //!< Upper bound of the row
      const HighsInt num_new_nz,  //!< Number of nonzeros in the row
      const HighsInt*
          indices,          //!< Array of size num_new_nz with column indices
      const double* values  //!< Array of size num_new_nz with column values
  );

  /**
   * @brief Adds multiple rows to the model
   */
  bool addRows(
      const HighsInt num_new_row,  //!< Number of new rows
      const double* lower,  //!< Array of size num_new_row with lower bounds
      const double* upper,  //!< Array of size num_new_row with upper bounds
      const HighsInt num_new_nz,  //!< Number of new nonzeros
      const HighsInt*
          starts,  //!< Array of size num_new_row with start indices of the rows
      const HighsInt* indices,  //!< Array of size num_new_nz with column
                                //!< indices for all rows
      const double*
          values  //!< Array of size num_new_nz with column values for all rows
  );

  /**
   * @brief Adds a column to the model
   */
  bool addCol(
      const double cost,          //!< Cost of the column
      const double lower,         //!< Lower bound of the column
      const double upper,         //!< Upper bound of the column
      const HighsInt num_new_nz,  //!< Number of nonzeros in the column
      const HighsInt* indices,    //!< Array of size num_new_nz with row indices
      const double* values        //!< Array of size num_new_nz with row values
  );

  /**
   * @brief Adds multiple columns to the model
   */
  bool addCols(
      const HighsInt num_new_col,  //!< Number of new columns
      const double* costs,         //!< Array of size num_new_col with costs
      const double* lower,  //!< Array of size num_new_col with lower bounds
      const double* upper,  //!< Array of size num_new_col with upper bounds
      const HighsInt num_new_nz,  //!< Number of new nonzeros
      const HighsInt* starts,  //!< Array of size num_new_row with start indices
                               //!< of the columns
      const HighsInt* indices,  //!< Array of size num_new_nz with row indices
                                //!< for all columns
      const double*
          values  //!< Array of size num_new_nz with row values for all columns
  );

  /**
   * @brief Delete multiple columns from the model given by an interval
   */
  bool deleteCols(const HighsInt from_col,  //!< The index of the first column
                                            //!< to delete from the model
                  const HighsInt to_col  //!< One more than the last column to
                                         //!< delete from the model
  );

  /**
   * @brief Delete multiple columns from the model given by a set
   */
  bool deleteCols(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set  //!< Array of size num_set_entries with indices of
                           //!< columns to delete
  );

  /**
   * @brief Delete multiple columns from the model given by a mask
   */
  bool deleteCols(
      HighsInt* mask  //!< Full length array with 1 => delete; !0 => not. The
                      //!< new index of any column
                      //! not deleted is returned in place of the value 0.
  );

  /**
   * @brief Delete multiple rows from the model given by an interval
   */
  bool deleteRows(const HighsInt from_row,  //!< The index of the first row to
                                            //!< delete from the model
                  const HighsInt to_row  //!< One more than the last row delete
                                         //!< from the model
  );

  /**
   * @brief Delete multiple rows from the model given by a set
   */
  bool deleteRows(
      const HighsInt num_set_entries,  //!< The number of indides in the set
      const HighsInt* set  //!< Array of size num_set_entries with indices of
                           //!< columns to delete
  );

  /**
   * @brief Delete multiple rows from the model given by a mask
   */
  bool deleteRows(HighsInt* mask  //!< Full length array with 1 => delete; 0 =>
                                  //!< not. The new index of any row not deleted
                                  //!< is returned in place of the value 0.
  );

  /**
   * @brief Scale a matrix column (and cost) by a constant - flipping bounds if
   * the constant is negative
   */
  bool scaleCol(const HighsInt col,    //!< Column to change
                const double scaleval  //!< Scaling value
  );

  /**
   * @brief Scale a matrix row by a constant - flipping bounds if the constant
   * is negative
   */
  bool scaleRow(const HighsInt row,    //!< Row to change
                const double scaleval  //!< Scaling value
  );

  /**
   * Other methods for specialist applications
   */

  /**
   * Methods for setting the basis and solution
   */

  /**
   * @brief Uses the HighsSolution passed to set the solution for the
   * LP of the (first?) HighsModelObject, according to ?? TODO Understand this
   * ??
   */
  // In the solution passed as a parameter below can have one or many of
  // col_value, col_dual and row_dual set. If any of them are not set the
  // solution in Highs does not get updated.
  HighsStatus setSolution(
      const HighsSolution& solution  //!< Solution to be used
  );

  /**
   * @brief Uses the HighsBasis passed to set the basis for the
   * LP of the (first?) HighsModelObject
   */
  HighsStatus setBasis(const HighsBasis& basis  //!< Basis to be used
  );

  /**
   * @brief Clears the HighsBasis for the LP of the HighsModelObject
   */
  HighsStatus setBasis();

  /**
   * @brief Gets the value of infinity used by HiGHS
   */
  double getInfinity();

  /**
   * @brief Gets the run time of HiGHS
   */
  double getRunTime();
  /**
   * @brief Clear data associated with solving the model: basis, solution and
   * internal data etc
   */
  HighsStatus clearSolver();

#ifdef HiGHSDEV
  /**
   * @brief Report the model status, solution and basis vector sizes and basis
   * validity
   */
  void reportModelStatusSolutionBasis(const std::string message,
                                      const HighsInt hmo_ix = -1);
#endif

  std::string modelStatusToString(const HighsModelStatus model_status) const;

  std::string solutionStatusToString(const HighsInt solution_status);

  void setMatrixOrientation(const MatrixOrientation& desired_orientation =
                                MatrixOrientation::kColwise);

#ifdef OSI_FOUND
  friend class OsiHiGHSSolverInterface;
#endif
  // Start of deprecated methods

  const HighsLp& getLp() const { return getModel(); }

  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const bool value            //!< The option value
  );

  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const HighsInt value        //!< The option value
  );

#ifdef HIGHSINT64
  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const int value             //!< The option value
  ) {
    return setHighsOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const double value          //!< The option value
  );

  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const std::string value     //!< The option value
  );

  HighsStatus setHighsOptionValue(
      const std::string& option,  //!< The option name
      const char* value           //!< The option value
  );

  HighsStatus readHighsOptions(const std::string filename  //!< The filename
  );

  HighsStatus passHighsOptions(const HighsOptions& options  //!< The options
  );

  HighsStatus getHighsOptionValue(
      const std::string& option,  //!< The option name
      bool& value                 //!< The option value
  );

  HighsStatus getHighsOptionValue(
      const std::string& option,  //!< The option name
      HighsInt& value             //!< The option value
  );

  HighsStatus getHighsOptionValue(
      const std::string& option,  //!< The option name
      double& value               //!< The option value
  );

  HighsStatus getHighsOptionValue(
      const std::string& option,  //!< The option name
      std::string& value          //!< The option value
  );

  HighsStatus getHighsOptionType(
      const std::string& option,  //!< The option name
      HighsOptionType& type       //!< The option type
  );

  const HighsOptions& getHighsOptions() const;

  HighsStatus resetHighsOptions();

  HighsStatus writeHighsOptions(
      const std::string filename,  //!< The filename
      const bool report_only_non_default_values = true);

  double getObjectiveValue() { return info_.objective_function_value; }

  HighsInt getSimplexIterationCount() { return info_.simplex_iteration_count; }

  HighsStatus setHighsLogfile(FILE* logfile = NULL);

  HighsStatus setHighsOutput(FILE* output = NULL);

  const HighsInfo& getHighsInfo() const;

  HighsStatus getHighsInfoValue(const std::string& info,  //!< The info name
                                HighsInt& value           //!< The info value
  );

  HighsStatus getHighsInfoValue(const std::string& info,  //!< The info name
                                double& value) const;     //!< The info value

  HighsStatus writeHighsInfo(const std::string filename  //!< The filename
  );

  double getHighsInfinity();

  double getHighsRunTime();
  // End of deprecated methods
 private:
  HighsSolution solution_;
  HighsBasis basis_;
  HighsLp lp_;

  HighsTimer timer_;

  HighsOptions options_;
  HighsIterationCounts iteration_counts_;
  HighsInfo info_;

  HighsModelStatus model_status_ = HighsModelStatus::kNotset;
  HighsModelStatus scaled_model_status_ = HighsModelStatus::kNotset;

  // Each HighsModelObject holds a const ref to its lp_. There are at most two
  // entries in hmos_: the original LP and the LP reduced by presolve
  std::vector<HighsModelObject> hmos_;

  // Record of maximum number of OMP threads. If OMP is available then
  // it's set to the correct positive number in Highs::run()
  HighsInt omp_max_threads = 0;

  // This is strictly for debugging. It's used to check whether
  // returnFromRun() was called after the previous call to
  // Highs::run() and, assuming that this is always done, it checks
  // whether Highs::run() is called recursively.
  bool called_return_from_run = true;

  HighsStatus callSolveLp(const HighsInt model_index, const string message);
  HighsStatus callSolveMip();

  PresolveComponent presolve_;
  HighsPresolveStatus runPresolve();
  HighsPostsolveStatus runPostsolve();

  HighsStatus openWriteFile(const string filename, const string method_name,
                            FILE*& file, bool& html) const;

  HighsStatus getUseModelStatus(
      HighsModelStatus& use_model_status,
      const double unscaled_primal_feasibility_tolerance,
      const double unscaled_dual_feasibility_tolerance,
      const bool rerun_from_logical_basis = false);

  bool unscaledOptimal(const double unscaled_primal_feasibility_tolerance,
                       const double unscaled_dual_feasibility_tolerance,
                       const bool report = false);

  bool haveHmo(const string method_name) const;

  void newHighsBasis();
  void forceHighsSolutionBasisSize();
  void setHighsModelStatusAndInfo(const HighsModelStatus model_status);
  void setHighsModelStatusBasisSolutionAndInfo();

  HighsStatus reset();

  void clearModelStatus();
  void clearSolution();
  void clearBasis();
  void clearInfo();
  void noSolution();

  void underDevelopmentLogMessage(const string method_name);
  HighsStatus returnFromRun(const HighsStatus return_status);
  HighsStatus returnFromHighs(const HighsStatus return_status);

  // Interface methods
  HighsStatus addColsInterface(HighsInt XnumNewCol, const double* XcolCost,
                               const double* XcolLower, const double* XcolUpper,
                               HighsInt XnumNewNZ, const HighsInt* XAstart,
                               const HighsInt* XAindex, const double* XAvalue);

  HighsStatus addRowsInterface(HighsInt XnumNewRow, const double* XrowLower,
                               const double* XrowUpper, HighsInt XnumNewNZ,
                               const HighsInt* XARstart,
                               const HighsInt* XARindex,
                               const double* XARvalue);

  HighsStatus deleteColsInterface(HighsIndexCollection& index_collection);

  HighsStatus deleteRowsInterface(HighsIndexCollection& index_collection);

  HighsStatus getColsInterface(const HighsIndexCollection& index_collection,
                               HighsInt& num_col, double* col_cost,
                               double* col_lower, double* col_upper,
                               HighsInt& num_nz, HighsInt* col_matrix_start,
                               HighsInt* col_matrix_index,
                               double* col_matrix_value);

  HighsStatus getRowsInterface(const HighsIndexCollection& index_collection,
                               HighsInt& num_row, double* row_lower,
                               double* row_upper, HighsInt& num_nz,
                               HighsInt* row_matrix_start,
                               HighsInt* row_matrix_index,
                               double* row_matrix_value);

  HighsStatus getCoefficientInterface(const HighsInt Xrow, const HighsInt Xcol,
                                      double& value);

  HighsStatus changeObjectiveSenseInterface(const ObjSense Xsense);
  HighsStatus changeIntegralityInterface(HighsIndexCollection& index_collection,
                                         const HighsVarType* usr_inegrality);
  HighsStatus changeCostsInterface(HighsIndexCollection& index_collection,
                                   const double* usr_col_cost);
  HighsStatus changeColBoundsInterface(HighsIndexCollection& index_collection,
                                       const double* usr_col_lower,
                                       const double* usr_col_upper);
  HighsStatus changeRowBoundsInterface(HighsIndexCollection& index_collection,
                                       const double* usr_row_lower,
                                       const double* usr_row_upper);
  HighsStatus changeCoefficientInterface(const HighsInt Xrow,
                                         const HighsInt Xcol,
                                         const double XnewValue);
  HighsStatus scaleColInterface(const HighsInt col, const double scaleval);
  HighsStatus scaleRowInterface(const HighsInt row, const double scaleval);
  HighsStatus setNonbasicStatusInterface(
      const HighsIndexCollection& index_collection, const bool columns);
  HighsStatus getBasicVariablesInterface(HighsInt* basic_variables);
  HighsStatus basisSolveInterface(const vector<double>& rhs,
                                  double* solution_vector,
                                  HighsInt* solution_num_nz,
                                  HighsInt* solution_indices, bool transpose);
  void clearBasisInterface();

  HighsStatus getDualRayInterface(bool& has_dual_ray, double* dual_ray_value);

  HighsStatus getPrimalRayInterface(bool& has_primal_ray,
                                    double* primal_ray_value);

  friend class HighsMipSolver;
};

#endif
