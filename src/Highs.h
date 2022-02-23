/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file Highs.h
 * @brief The HiGHS class
 */
#ifndef HIGHS_H_
#define HIGHS_H_

#include <sstream>

#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsRanging.h"
#include "lp_data/HighsSolutionDebug.h"
#include "model/HighsModel.h"
#include "presolve/PresolveComponent.h"

/**
 * @brief Class to set parameters and run HiGHS
 */
class Highs {
 public:
  Highs();
  virtual ~Highs() {
    FILE* log_file_stream = options_.log_options.log_file_stream;
    if (log_file_stream != nullptr) {
      assert(log_file_stream != stdout);
      fclose(log_file_stream);
    }
  }

  /**
   * @brief Reset the options and then call clearModel()
   */
  HighsStatus clear();

  /**
   * @brief Clear the incumbent model and then call clearSolver()
   */
  HighsStatus clearModel();

  /**
   * @brief Clear all solution data associated with the model
   */
  HighsStatus clearSolver();

  /**
   * Methods for model input
   */

  /**
   * Every model loading module eventually uses
   * passModel(HighsModel model) to communicate the model to HiGHS.
   */

  /**
   * @brief Pass a HighsModel instance to Highs
   */
  HighsStatus passModel(HighsModel model);

  /**
   * @brief Pass a HighsLp instance to Highs
   */
  HighsStatus passModel(HighsLp lp);

  /**
   * @brief Pass a QP (possibly with integrality data) via pointers to vectors of data
   */
  HighsStatus passModel(const HighsInt num_col, const HighsInt num_row, const HighsInt num_nz,
			const HighsInt q_num_nz, const HighsInt a_format, const HighsInt q_format,
			const HighsInt sense, const double offset, const double* costs,
			const double* col_lower, const double* col_upper, const double* row_lower,
			const double* row_upper, const HighsInt* a_start, const HighsInt* a_index,
			const double* a_value, const HighsInt* q_start, const HighsInt* q_index,
			const double* q_value, const HighsInt* integrality = nullptr);

  /**
   * @brief Pass an LP (possibly with integrality data) via pointers to vectors of data
   */
  HighsStatus passModel(const HighsInt num_col, const HighsInt num_row,
                        const HighsInt num_nz, const HighsInt a_format,
                        const HighsInt sense, const double offset,
                        const double* costs, const double* col_lower,
                        const double* col_upper, const double* row_lower,
                        const double* row_upper, const HighsInt* a_start,
                        const HighsInt* a_index, const double* a_value,
                        const HighsInt* integrality = nullptr);

  /**
   * @brief Pass a HighsHessian instance for the incumbent model
   */
  HighsStatus passHessian(HighsHessian hessian_);

  /**
   * @brief Pass the Hessian for the incumbent model via pointers to vectors of data
   */
  HighsStatus passHessian(const HighsInt dim, const HighsInt num_nz,
                          const HighsInt format, const HighsInt* start,
                          const HighsInt* index, const double* value);

  /**
   * @brief Read in a model
   */
  HighsStatus readModel(const std::string filename);

  /**
   * @brief Read in a basis
   */
  HighsStatus readBasis(const std::string filename);

  /**
   * @brief Presolve the incumbent model
   */
  HighsStatus presolve();

  /**
   * @brief Solve the incumbent model according to the specified options
   */
  HighsStatus run();

  /**
   * @brief Postsolve the incumbent model
   */
  HighsStatus postsolve(const HighsSolution& solution, const HighsBasis& basis);

  /**
   * @brief Write the current solution to a file in a given style
   */
  HighsStatus writeSolution(const std::string filename,
                            const HighsInt style);

  /**
   * @brief Read a HiGHS solution file in a given style
   */
  HighsStatus readSolution(const std::string filename,
                           const HighsInt style);

  /**
   * @brief Check the feasibility of the current solution. Of value
   * after calling Highs::readSolution
   */
  HighsStatus checkSolutionFeasibility();

  /**
   * Methods for HiGHS option input/output
   */

  /**
   * @brief Set an option to the bool/HighsInt/double/string value if it's
   * legal and, for bool/HighsInt/double, only if it's of the correct type
   */

  HighsStatus setOptionValue(const std::string& option,
                             const bool value          
  );

  HighsStatus setOptionValue(const std::string& option,
                             const HighsInt value      
  );

#ifdef HIGHSINT64
  HighsStatus setOptionValue(const std::string& option,
                             const int value           
  ) {
    return setOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setOptionValue(const std::string& option,
                             const double value        
  );

  HighsStatus setOptionValue(const std::string& option,
                             const std::string value   
  );

  HighsStatus setOptionValue(const std::string& option,
                             const char* value         
  );

  /**
   * @brief Read option values from a file
   */
  HighsStatus readOptions(const std::string filename);

  /**
   * @brief Pass a HighsOptions instance to Highs
   */
  HighsStatus passOptions(const HighsOptions& options);

  /**
   * @brief Get a const reference to the internal option values
   */
  const HighsOptions& getOptions() const { return options_; }

  /**
   * @brief Gets an option value as bool/HighsInt/double/string and, for
   * bool/int/double, only if it's of the correct type.
   */
  HighsStatus getOptionValue(const std::string& option,
                             bool& value) const;

  HighsStatus getOptionValue(const std::string& option,
                             HighsInt& value) const;

  HighsStatus getOptionValue(const std::string& option,
                             double& value) const;

  HighsStatus getOptionValue(const std::string& option,
                             std::string& value) const;

  /**
   * @brief Get the type expected by an option
   */
  HighsStatus getOptionType(const std::string& option,
                            HighsOptionType& type) const;

  /**
   * @brief Reset the options to the default values
   */
  HighsStatus resetOptions();

  /**
   * @brief Write (deviations from default values of) the options to a
   * file, with the extension ".html" producing HTML, otherwise using
   * the standard format used to read options from a file.
   */
  HighsStatus writeOptions(const std::string filename,  //!< The filename
                           const bool report_only_deviations = false) const;

  /**
   * @brief Get a const reference to the internal info values
   * type.
   */
  const HighsInfo& getInfo() const { return info_; }

  /**
   * @brief Get an info value as HighsInt/int64_t/double, and only if
   * it's of the correct type.
   */

  HighsStatus getInfoValue(const std::string& info,
                           HighsInt& value) const;

#ifndef HIGHSINT64
  HighsStatus getInfoValue(const std::string& info,
                           int64_t& value) const;
#endif

  HighsStatus getInfoValue(const std::string& info,
                           double& value) const;

  /**
   * @brief Write info values to a file, with the extension ".html"
   * producing HTML, otherwise using the standard format used to read
   * options from a file.
   */
  HighsStatus writeInfo(const std::string filename) const;

  /**
   * Methods for model output
   */

  /**
   * @brief Returns the presolved HighsModel instance in HiGHS
   */
  const HighsLp& getPresolvedLp() const { return presolved_model_.lp_; }

  /**
   * @brief Returns the presolved HighsModel instance in HiGHS
   */
  const HighsModel& getPresolvedModel() const { return presolved_model_; }

  /**
   * @brief Returns the HighsLp instance in the HiGHS model
   */
  const HighsLp& getLp() const { return model_.lp_; }

  /**
   * @brief Returns the model in HiGHS
   */
  const HighsModel& getModel() const { return model_; }

  /**
   * @brief Returns the HighsSolution
   */
  const HighsSolution& getSolution() const { return solution_; }

  /**
   * @brief Returns the HighsBasis
   */
  const HighsBasis& getBasis() const { return basis_; }

  /**
   * @brief Gets the hot start basis data from the most recent simplex
   * solve. Advanced method: for HiGHS MIP solver
   */
  const HotStart& getHotStart() const { return ekk_instance_.hot_start_; }

  /**
   * @brief Returns the current model status
   */
  const HighsModelStatus& getModelStatus(
      const bool scaled_model = false) const {
    return scaled_model ? scaled_model_status_ : model_status_;
  }

  /**
   * @brief Indicates whether a dual unbounded ray exdists, and gets
   * it if it does and dual_ray is not nullptr
   */
  HighsStatus getDualRay(bool& has_dual_ray, double* dual_ray_value = nullptr);

  /**
   * @brief Indicates whether a primal unbounded ray exdists, and gets
   * it if it does and primal_ray is not nullptr
   */
  HighsStatus getPrimalRay(bool& has_primal_ray,
                           double* primal_ray_value = nullptr);

  /**
   * @brief Gets the ranging information for the current LP, possibly
   * returning it, as well as holding it internally
   */
  HighsStatus getRanging();
  HighsStatus getRanging(HighsRanging& ranging);

  /**
   * @brief Gets the current model objective value
   */
  double getObjectiveValue() { return info_.objective_function_value; }

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
      const HighsInt row,              //!< Index of row required
      double* row_vector,              //!< Row required
      HighsInt* row_num_nz = nullptr,  //!< Number of nonzeros
      HighsInt* row_indices = nullptr  //!< Indices of nonzeros
  );

  /**
   * @brief Gets a column of \f$B^{-1}\f$ for basis matrix \f$B\f$
   */
  HighsStatus getBasisInverseCol(
      const HighsInt col,              //!< Index of column required
      double* col_vector,              //!< Column required
      HighsInt* col_num_nz = nullptr,  //!< Number of nonzeros
      HighsInt* col_indices = nullptr  //!< Indices of nonzeros
  );

  /**
   * @brief Forms \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$
   */
  HighsStatus getBasisSolve(
      const double* rhs,                    //!< RHS \f$\mathbf{b}\f$
      double* solution_vector,              //!< Solution  \f$\mathbf{x}\f$
      HighsInt* solution_num_nz = nullptr,  //!< Number of nonzeros
      HighsInt* solution_indices = nullptr  //!< Indices of nonzeros
  );

  /**
   * @brief Forms \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
   * \f$\mathbf{b}\f$
   */
  HighsStatus getBasisTransposeSolve(
      const double* rhs,                    //!< RHS \f$\mathbf{b}\f$
      double* solution_vector,              //!< Solution  \f$\mathbf{x}\f$
      HighsInt* solution_nz = nullptr,      //!< Number of nonzeros
      HighsInt* solution_indices = nullptr  //!< Indices of nonzeros
  );

  /**
   * @brief Forms a row of \f$B^{-1}A\f$
   */
  HighsStatus getReducedRow(
      const HighsInt row,               //!< Index of row required
      double* row_vector,               //!< Row required
      HighsInt* row_num_nz = nullptr,   //!< Number of nonzeros
      HighsInt* row_indices = nullptr,  //!< Indices of nonzeros
      const double* pass_basis_inverse_row_vector =
          nullptr  //!< Necessary row of \f$B^{-1}\f$
  );

  /**
   * @brief Forms a column of \f$B^{-1}A\f$
   */
  HighsStatus getReducedColumn(
      const HighsInt col,              //!< Index of column required
      double* col_vector,              //!< Column required
      HighsInt* col_num_nz = nullptr,  //!< Number of nonzeros
      HighsInt* col_indices = nullptr  //!< Indices of nonzeros
  );

  /**
   * @brief Get the number of columns in the incumbent model
   */
  HighsInt getNumCol() const { return model_.lp_.num_col_; }

  /**
   * @brief Get the number of rows in the incumbent model
   */
  HighsInt getNumRow() const { return model_.lp_.num_row_; }

  /**
   * @brief Get the number of (constraint matrix) nonzeros in the incumbent
   * model
   */
  HighsInt getNumNz() const { return model_.lp_.a_matrix_.numNz(); }

  /**
   * @brief Get the number of Hessian matrix nonzeros in the incumbent model
   */
  HighsInt getHessianNumNz() { return model_.hessian_.numNz(); }

  /**
   * @brief Get the objective sense of the model
   */
  HighsStatus getObjectiveSense(ObjSense& sense);

  /**
   * @brief Get the objective offset of the model
   */
  HighsStatus getObjectiveOffset(double& offset);

  /**
   * @brief Get multiple columns from the model given by an interval
   */
  HighsStatus getCols(const HighsInt from_col, //!< The index of the first column to get from the model
		      const HighsInt to_col,   //!< The index of the last column to get from the model
		      HighsInt& num_col,       //!< Number of columns got from the model
		      double* costs,           //!< Array of size num_col with costs
		      double* lower,           //!< Array of size num_col with lower bounds
		      double* upper,           //!< Array of size num_col with upper bounds
		      HighsInt& num_nz,        //!< Number of nonzeros got from the model
		      HighsInt* start,         //!< Array of size num_col with start indices of the columns
		      HighsInt* index,         //!< Array of size num_nz with row indices for the columns
		      double* value            //!< Array of size num_nz with row values for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a set
   */
  HighsStatus getCols(const HighsInt num_set_entries, //!< The number of indides in the set
		      const HighsInt* set,            //!< Array of size num_set_entries with indices of columns to get
		      HighsInt& num_col,              //!< Number of columns got from the model
		      double* costs,                  //!< Array of size num_col with costs
		      double* lower,                  //!< Array of size num_col with lower bounds
		      double* upper,                  //!< Array of size num_col with upper bounds
		      HighsInt& num_nz,               //!< Number of nonzeros got from the model
		      HighsInt* start,                //!< Array of size num_col with start indices of the columns
		      HighsInt* index,                //!< Array of size num_nz with row indices for the columns
		      double* value                   //!< Array of size num_nz with row values for the columns
  );

  /**
   * @brief Get multiple columns from the model given by a mask
   */
  HighsStatus getCols(const HighsInt* mask, //!< Full length array with 1 => get; 0 => not
		      HighsInt& num_col,    //!< Number of columns got from the model
		      double* costs,        //!< Array of size num_col with costs
		      double* lower,        //!< Array of size num_col with lower bounds
		      double* upper,        //!< Array of size num_col with upper bounds
		      HighsInt& num_nz,     //!< Number of nonzeros got from the model
		      HighsInt* start,      //!<  Array of size num_col with start indices of the columns
		      HighsInt* index,      //!<  Array of size num_nz with row indices for the columns
		      double* value         //!<  Array of size num_nz with row values for the columns
  );

  /**
   * @brief Get multiple rows from the model given by an interval
   */
  HighsStatus getRows(const HighsInt from_row, //!< The index of the first row to get from the model
		      const HighsInt to_row,   //!< The index of the last row get from the model
		      HighsInt& num_row,       //!< Number of rows got from the model
		      double* lower,           //!< Array of size num_row with lower bounds
		      double* upper,           //!< Array of size num_row with upper bounds
		      HighsInt& num_nz,        //!< Number of nonzeros got from the model
		      HighsInt* start,         //!< Array of size num_row with start indices of the rows
		      HighsInt* index,         //!< Array of size num_nz with column indices for the rows
		      double* value            //!< Array of size num_nz with column values for the rows
  );

  /**
   * @brief Get multiple rows from the model given by a set
   */
  HighsStatus getRows(const HighsInt num_set_entries, //!< The number of indides in the set
		      const HighsInt* set,            //!< Array of size num_set_entries with indices of rows to get
		      HighsInt& num_row,              //!< Number of rows got from the model
		      double* lower,                  //!< Array of size num_row with lower bounds
		      double* upper,                  //!< Array of size num_row with upper bounds
		      HighsInt& num_nz,               //!< Number of nonzeros got from the model
		      HighsInt* start,                //!< Array of size num_row with start indices of the rows
		      HighsInt* index,                //!< Array of size num_nz with column indices for the rows
		      double* value                   //!< Array of size num_nz with column values for the rows
  );

  /**
   * @brief Get multiple rows from the model given by a mask
   */
  HighsStatus getRows(const HighsInt* mask, //!< Full length array with 1 => get; 0 => not
		      HighsInt& num_row,    //!< Number of rows got from the model
		      double* lower,        //!< Array of size num_row with lower bounds
		      double* upper,        //!< Array of size num_row with upper bounds
		      HighsInt& num_nz,     //!< Number of nonzeros got from the model
		      HighsInt* start,      //!< Array of size num_row with start indices of the rows
		      HighsInt* index,      //!< Array of size num_nz with column indices for the rows
		      double* value         //!< Array of size num_nz with column values for the rows
  );

  /**
   * @brief Get a matrix coefficient
   */
  HighsStatus getCoeff(const HighsInt row,  //!< Row of coefficient to be got
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
  HighsStatus changeObjectiveSense(const ObjSense sense  //!< New objective sense
  );

  /**
   * @brief Change the objective offset of the model
   */
  HighsStatus changeObjectiveOffset(const double offset  //!< New objective offset
  );

  /**
   * @brief Change the integrality of a column
   */
  HighsStatus changeColIntegrality(const HighsInt col,            //!< The index of the column whose integrality is to change
				   const HighsVarType integrality //!< The new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by an interval
   */
  HighsStatus changeColsIntegrality(const HighsInt from_col, //!< The index of the first column whose integrality changes
				    const HighsInt to_col,   //!< The index of the last column whose integrality changes
				    const HighsVarType* integrality //!< Array of size num_set_entries with new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by a set of indices
   */
  HighsStatus changeColsIntegrality(const HighsInt num_set_entries,  //!< The number of indides in the set
				    const HighsInt* set,  //!< Array of size num_set_entries with indices of
				    //!< columns whose integralitys change
				    const HighsVarType*
				    integrality  //!< Array of size num_set_entries with new integrality
  );

  /**
   * @brief Change the integrality of multiple columns given by a mask
   */
  HighsStatus changeColsIntegrality(const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
				    const HighsVarType* integrality  //!< Full length array of new integrality
  );

  /**
   * @brief Change the cost of a column
   */
  HighsStatus changeColCost(const HighsInt col,  //!< The index of the column whose cost is to change
			    const double cost    //!< The new cost
  );

  /**
   * @brief Change the cost of multiple columns given by an interval
   */
  HighsStatus changeColsCost(const HighsInt from_col,  //!< The index of the first column whose cost changes
			     const HighsInt to_col,  //!< The index of the last column whose cost changes
			     const double* cost      //!< Array of size num_set_entries with new costs
  );

  /**
   * @brief Change the cost of multiple columns given by a set of indices
   */
  HighsStatus changeColsCost(const HighsInt num_set_entries,  //!< The number of indides in the set
			     const HighsInt* set,  //!< Array of size num_set_entries with indices of
			     //!< columns whose costs change
			     const double* cost    //!< Array of size num_set_entries with new costs
  );

  /**
   * @brief Change the cost of multiple columns given by a mask
   */
  HighsStatus changeColsCost(const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
			     const double* cost     //!< Full length array of new costs
  );

  /**
   * @brief Change the bounds of a column
   */
  HighsStatus changeColBounds(const HighsInt col,  //!< The index of the column whose bounds are to change
			      const double lower,  //!< The new lower bound
			      const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple columns given by an interval
   */
  HighsStatus changeColsBounds(const HighsInt from_col,  //!< The index of the first column whose bounds change
			       const HighsInt to_col,  //!< The index of the last column whose bounds change
			       const double* lower,  //!< Array of size to_col-from_col with new lower bounds
			       const double* upper  //!< Array of size to_col-from_col with new upper bounds
  );

  /**
   * @brief Change the bounds of multiple columns given by a set of indices
   */
  HighsStatus changeColsBounds(const HighsInt num_set_entries,  //!< The number of indides in the set
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
  HighsStatus changeColsBounds(const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
			       const double* lower,   //!< Full length array of new lower bounds
			       const double* upper    //!< Full length array of new upper bounds
  );

  /**
   * @brief Change the bounds of a row
   */
  HighsStatus changeRowBounds(const HighsInt row,  //!< The index of the row whose bounds are to change
			      const double lower,  //!< The new lower bound
			      const double upper   //!< The new upper bound
  );

  /**
   * @brief Change the bounds of multiple rows given by an interval
   */
  HighsStatus changeRowsBounds(const HighsInt from_row,  //!< The index of the first column whose bounds change
			       const HighsInt to_row,  //!< The index of the last column whose bounds change
                               const double* lower,//!< Array of size to_col-from_col with new lower bounds
			       const double* upper //!< Array of size to_col-from_col with new upper bounds
  );

  /**
   * @brief Change the bounds of multiple rows given by a set of indices
   */
  HighsStatus changeRowsBounds(const HighsInt num_set_entries,  //!< The number of indides in the set
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
  HighsStatus changeRowsBounds(const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
			       const double* lower,   //!< Full length array of new lower bounds
			       const double* upper    //!< Full length array of new upper bounds
			       );

  /**
   * @brief Change a matrix coefficient
   */
  HighsStatus changeCoeff(const HighsInt row,  //!< Row of coefficient to be changed
			  const HighsInt col,  //!< Column of coefficient to be changed
			  const double value   //!< Coefficient
  );
  /**
   * @brief Adds a row to the model
   */
  HighsStatus addRow(const double lower,         //!< Lower bound of the row
		     const double upper,         //!< Upper bound of the row
		     const HighsInt num_new_nz,  //!< Number of nonzeros in the row
		     const HighsInt*
		     indices,          //!< Array of size num_new_nz with column indices
		     const double* values  //!< Array of size num_new_nz with column values
  );

  /**
   * @brief Adds multiple rows to the model
   */
  HighsStatus addRows(const HighsInt num_new_row,  //!< Number of new rows
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
  HighsStatus addCol(const double cost,          //!< Cost of the column
		     const double lower,         //!< Lower bound of the column
		     const double upper,         //!< Upper bound of the column
		     const HighsInt num_new_nz,  //!< Number of nonzeros in the column
		     const HighsInt* indices,    //!< Array of size num_new_nz with row indices
		     const double* values        //!< Array of size num_new_nz with row values
  );

  /**
   * @brief Adds multiple columns to the model
   */
  HighsStatus addCols(const HighsInt num_new_col,  //!< Number of new columns
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
  HighsStatus deleteCols(const HighsInt from_col,  //!< The index of the first column to delete from the model
			 const HighsInt to_col     //!< The index of the last column to delete from the model
  );

  /**
   * @brief Delete multiple columns from the model given by a set
   */
  HighsStatus deleteCols(const HighsInt num_set_entries,  //!< The number of indides in the set
			 const HighsInt* set  //!< Array of size num_set_entries with indices of
			 //!< columns to delete
  );

  /**
   * @brief Delete multiple columns from the model given by a mask
   */
  HighsStatus deleteCols(HighsInt* mask  //!< Full length array with 1 => delete; !0 => not. The
			 //!< new index of any column
			 //! not deleted is returned in place of the value 0.
  );

  /**
   * @brief Delete multiple rows from the model given by an interval
   */
  HighsStatus deleteRows(const HighsInt from_row,  //!< The index of the first row to delete from the model
			 const HighsInt to_row     //!< The index of the last row to delete from the model
  );

  /**
   * @brief Delete multiple rows from the model given by a set
   */
  HighsStatus deleteRows(const HighsInt num_set_entries,  //!< The number of indides in the set
			 const HighsInt* set  //!< Array of size num_set_entries with indices of
			 //!< columns to delete
  );

  /**
   * @brief Delete multiple rows from the model given by a mask
   */
  HighsStatus deleteRows(HighsInt* mask  //!< Full length array with 1 => delete; 0 =>
			 //!< not. The new index of any row not deleted
			 //!< is returned in place of the value 0.
  );

  /**
   * @brief Scale a matrix column (and cost) by a constant - flipping bounds if
   * the constant is negative
   */
  HighsStatus scaleCol(const HighsInt col,    //!< Column to change
                       const double scale_value  //!< Scaling value
  );

  /**
   * @brief Scale a matrix row by a constant - flipping bounds if the constant
   * is negative
   */
  HighsStatus scaleRow(const HighsInt row,    //!< Row to change
                       const double scale_value  //!< Scaling value
  );

  /**
   * Other methods for specialist applications
   */

  /**
   * Methods for setting basis_ and solution_
   */

  /**
   * @brief Uses the HighsSolution passed to set solution_
   */
  // In the solution passed as a parameter below can have one or many of
  // col_value, col_dual and row_dual set. If any of them are not set the
  // solution in Highs does not get updated.
  HighsStatus setSolution(const HighsSolution& solution);

  /**
   * @brief Sets the callback method and user data to use for logging
   */
  HighsStatus setLogCallback(void (*log_callback)(HighsLogType, const char*,
                                                  void*),
                             void* log_callback_data = nullptr);

  /**
   * @brief Uses the HighsBasis passed to set basis_
   */
  HighsStatus setBasis(const HighsBasis& basis, const std::string origin = "");

  /**
   * @brief Clears basis_
   */
  HighsStatus setBasis();

  /**
   * @brief Sets up for simpelx using the supplied hot start
   * data. Advanced method: for HiGHS MIP solver
   */
  HighsStatus setHotStart(const HotStart& hot_start);

  /**
   * @brief Freezes the current basis and standard NLA, returning a
   * value to be used to recover this basis and standard NLA at
   * minimal cost. Advanced method: for HiGHS MIP solver
   */
  HighsStatus freezeBasis(HighsInt& frozen_basis_id);

  /**
   * @brief Unfreeze a frozen basis and standard NLA (if
   * possible). Advanced method: for HiGHS MIP solver
   */
  HighsStatus unfreezeBasis(const HighsInt frozen_basis_id);

  /**
   * @brief Checks that all frozen basis data has been
   * cleared. Advanced method: for HiGHS MIP solver
   */
  HighsStatus frozenBasisAllDataClear() {
    return ekk_instance_.frozenBasisAllDataClear();
  }

  /**
   * @Brief Put a copy of the current iterate - basis; invertible
   * representation and dual edge weights - into storage within
   * HSimplexNla. Advanced method: for HiGHS MIP solver
   */
  HighsStatus putIterate();

  /**
   * @Brief Get a copy of the iterate stored within HSimplexNla and
   * overwrite the current iterate. Advanced method: for HiGHS MIP
   * solver
   */
  HighsStatus getIterate();

  /**
   * @brief Gets the value of infinity used by HiGHS
   */
  double getInfinity() { return kHighsInf; }

  /**
   * @brief Gets the run time of HiGHS
   */
  double getRunTime() { return timer_.readRunHighsClock(); }

  /**
   * @brief Gets the dual edge weights (steepest/devex) in the order of the
   * basic indices or nullptr when they are not available.
   */
  const double* getDualEdgeWeights() const {
    return ekk_instance_.dual_edge_weight_.empty()
               ? nullptr
               : ekk_instance_.dual_edge_weight_.data();
  }

  /**
   * @brief Runs ipx crossover and if successful loads basis into Highs::basis_
   */
  HighsStatus crossover();
  HighsStatus crossover(HighsSolution& solution);

  /**
   * @brief Opens a named log file
   */
  HighsStatus openLogFile(const std::string log_file = "");

  std::string modelStatusToString(const HighsModelStatus model_status) const;

  std::string solutionStatusToString(const HighsInt solution_status) const;

  std::string basisStatusToString(const HighsBasisStatus basis_status) const;

  std::string basisValidityToString(const HighsInt basis_validity) const;

  HighsStatus setMatrixFormat(const MatrixFormat desired_format) {
    this->model_.lp_.setFormat(desired_format);
    return HighsStatus::kOk;
  }

#ifdef OSI_FOUND
  friend class OsiHiGHSSolverInterface;
#endif
  // Start of deprecated methods

  HighsInt getNumCols() const {
    deprecationMessage("getNumCols", "getNumCol");
    return getNumCol();
  }
  HighsInt getNumRows() const {
    deprecationMessage("getNumRows", "getNumRow");
    return getNumRow();
  }
  HighsInt getNumEntries() {
    deprecationMessage("getNumEntries", "getNumNz");
    return getNumNz();
  }

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const bool value    );

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const HighsInt value);

#ifdef HIGHSINT64
  HighsStatus setHighsOptionValue(
      const std::string& option,
      const int value             //!< The option value
  ) {
    deprecationMessage("setHighsOptionValue", "setOptionValue");
    return setOptionValue(option, HighsInt{value});
  }
#endif

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const double value  );

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const std::string value     //!< The option value
  );

  HighsStatus setHighsOptionValue(
      const std::string& option,
      const char* value   );

  HighsStatus readHighsOptions(const std::string filename  //!< The filename
  );

  HighsStatus passHighsOptions(const HighsOptions& options  //!< The options
  );

  HighsStatus getHighsOptionValue(
      const std::string& option,
      bool& value         );

  HighsStatus getHighsOptionValue(
      const std::string& option,
      HighsInt& value     );

  HighsStatus getHighsOptionValue(
      const std::string& option,
      double& value       );

  HighsStatus getHighsOptionValue(
      const std::string& option,
      std::string& value  );

  HighsStatus getHighsOptionType(
      const std::string& option,
      HighsOptionType& type       //!< The option type
  );

  const HighsOptions& getHighsOptions() const;

  HighsStatus resetHighsOptions();

  HighsStatus writeHighsOptions(
      const std::string filename,  //!< The filename
      const bool report_only_non_default_values = true);

  HighsInt getSimplexIterationCount() {
    deprecationMessage("getSimplexIterationCount", "None");
    return info_.simplex_iteration_count;
  }

  HighsStatus setHighsLogfile(FILE* logfile = nullptr);

  HighsStatus setHighsOutput(FILE* output = nullptr);

  const HighsInfo& getHighsInfo() const;

  HighsStatus getHighsInfoValue(const std::string& info,
                                HighsInt& value      );

  HighsStatus getHighsInfoValue(const std::string& info,
                                double& value) const;     //!< The info value

  HighsStatus writeHighsInfo(const std::string filename  //!< The filename
  );

  double getHighsInfinity();

  double getHighsRunTime();

  HighsStatus writeSolution(const std::string filename,
                            const bool pretty = false) const;

  void deprecationMessage(const std::string method_name,
                          const std::string alt_method_name) const;

  // End of deprecated methods
 private:
  HighsSolution solution_;
  HighsBasis basis_;
  HighsModel model_;
  HighsModel presolved_model_;
  HighsTimer timer_;

  HighsOptions options_;
  HighsInfo info_;
  HighsRanging ranging_;

  HighsPresolveStatus model_presolve_status_ =
      HighsPresolveStatus::kNotPresolved;
  HighsModelStatus model_status_ = HighsModelStatus::kNotset;
  HighsModelStatus scaled_model_status_ = HighsModelStatus::kNotset;

  HEkk ekk_instance_;

  HighsInt max_threads = 0;
  // This is strictly for debugging. It's used to check whether
  // returnFromRun() was called after the previous call to
  // Highs::run() and, assuming that this is always done, it checks
  // whether Highs::run() is called recursively.
  bool called_return_from_run = true;
  HighsInt debug_run_call_num_ = 0;

  void exactResizeModel() {
    this->model_.lp_.exactResize();
    this->model_.hessian_.exactResize();
  }

  HighsStatus callSolveLp(HighsLp& lp, const string message);
  HighsStatus callSolveQp();
  HighsStatus callSolveMip();

  PresolveComponent presolve_;
  HighsPresolveStatus runPresolve();
  HighsPostsolveStatus runPostsolve();

  HighsStatus openWriteFile(const string filename, const string method_name,
                            FILE*& file, bool& html) const;

  void reportModel();
  void newHighsBasis();
  void forceHighsSolutionBasisSize();
  //
  // For cases where there is no solution data for the model, but its
  // status is proved otherwise. Sets the model status, then clears any solution
  // and basis data
  void setHighsModelStatusAndClearSolutionAndBasis(
      const HighsModelStatus model_status);
  //
  // Sets unscaled and scaled model status, basis, solution and info
  // from the highs_model_object
  void setBasisValidity();
  //
  // Clears the presolved model and its status
  void clearPresolve();
  //
  // Methods to clear solver data for users in Highs class members
  // before (possibly) updating them with data from trying to solve
  // the inumcumbent model.
  //
  // Clears all solver data in Highs class members by calling
  // clearModelStatus(), clearSolution(), clearBasis(),
  // clearInfo() and clearEkk()
  void clearUserSolverData();
  //
  // Clears the model status, solution_ and info_
  void clearModelStatusSolutionAndInfo();
  //
  // Sets unscaled and scaled model status to HighsModelStatus::kNotset
  void clearModelStatus();
  //
  // Sets primal and dual solution status to
  // kSolutionStatusNone, and clears solution_ vectors
  void clearSolution();
  //
  // Invalidates basis and clears basis_ vectors
  void clearBasis();
  //
  // Invalidates info_ and resets the values of its members
  void clearInfo();
  //
  // Invalidates ranging_ and clears its vectors
  void clearRanging();

  // Invalidates ekk_instance_
  void clearEkk();

  HighsStatus returnFromRun(const HighsStatus return_status);
  HighsStatus returnFromHighs(const HighsStatus return_status);
  void reportSolvedLpQpStats();

  void underDevelopmentLogMessage(const std::string method_name);

  // Interface methods
  HighsStatus addColsInterface(HighsInt ext_num_new_col, const double* ext_col_cost,
                               const double* ext_col_lower, const double* ext_col_upper,
                               HighsInt ext_num_new_nz, const HighsInt* ext_a_start,
                               const HighsInt* ext_a_index, const double* ext_a_value);

  HighsStatus addRowsInterface(HighsInt ext_num_new_row, const double* ext_row_lower,
                               const double* ext_row_upper, HighsInt ext_num_new_nz,
                               const HighsInt* ext_ar_start,
                               const HighsInt* ext_ar_index,
                               const double* ext_ar_value);

  void deleteColsInterface(HighsIndexCollection& index_collection);

  void deleteRowsInterface(HighsIndexCollection& index_collection);

  void getColsInterface(const HighsIndexCollection& index_collection,
                        HighsInt& num_col, double* col_cost, double* col_lower,
                        double* col_upper, HighsInt& num_nz,
                        HighsInt* col_matrix_start, HighsInt* col_matrix_index,
                        double* col_matrix_value);

  void getRowsInterface(const HighsIndexCollection& index_collection,
                        HighsInt& num_row, double* row_lower, double* row_upper,
                        HighsInt& num_nz, HighsInt* row_matrix_start,
                        HighsInt* row_matrix_index, double* row_matrix_value);

  void getCoefficientInterface(const HighsInt ext_row, const HighsInt ext_col,
                               double& value);

  HighsStatus changeObjectiveSenseInterface(const ObjSense ext_sense);
  HighsStatus changeObjectiveOffsetInterface(const double ext_offset);
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
  void changeCoefficientInterface(const HighsInt ext_row, const HighsInt ext_col,
                                  const double ext_new_value);
  HighsStatus scaleColInterface(const HighsInt col, const double scale_value);
  HighsStatus scaleRowInterface(const HighsInt row, const double scale_value);

  void setNonbasicStatusInterface(const HighsIndexCollection& index_collection,
                                  const bool columns);
  void appendNonbasicColsToBasisInterface(const HighsInt ext_num_new_col);
  void appendBasicRowsToBasisInterface(const HighsInt ext_num_new_row);

  HighsStatus getBasicVariablesInterface(HighsInt* basic_variables);
  HighsStatus basisSolveInterface(const vector<double>& rhs,
                                  double* solution_vector,
                                  HighsInt* solution_num_nz,
                                  HighsInt* solution_indices, bool transpose);

  HighsStatus setHotStartInterface(const HotStart& hot_start);

  void zeroIterationCounts();

  HighsStatus getDualRayInterface(bool& has_dual_ray, double* dual_ray_value);

  HighsStatus getPrimalRayInterface(bool& has_primal_ray,
                                    double* primal_ray_value);
  bool aFormatOk(const HighsInt num_nz, const HighsInt format);
  bool qFormatOk(const HighsInt num_nz, const HighsInt format);
  void clearZeroHessian();
  HighsStatus checkOptimality(const std::string solver_type,
                              HighsStatus return_status);
  HighsStatus invertRequirementError(std::string method_name);
};

#endif
