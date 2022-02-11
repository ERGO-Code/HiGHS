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
#ifndef HIGHS_C_API
#define HIGHS_C_API

#include "util/HighsInt.h"
const HighsInt HighsStatuskError = -1;
const HighsInt HighsStatuskOk = 0;
const HighsInt HighsStatuskWarning = 1;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Formulate and solve a linear program using HiGHS.
 *
 * @param numcol    the number of columns
 * @param numrow    the number of rows
 * @param numnz     the number of nonzeros in the constraint matrix
 * @param a_format  the format of the constraint matrix
 * @param sense     the optimization sense. Use -1 for maximization, 1 otherwise
 * @param offset    the objective constant
 * @param colcost   array of length [numcol] with the column costs
 * @param collower  array of length [numcol] with the column lower bounds
 * @param colupper  array of length [numcol] with the column upper bounds
 * @param rowlower  array of length [numrow] with the row lower bounds
 * @param rowupper  array of length [numrow] with the row upper bounds
 * @param astart    the constraint matrix is provided to HiGHS in compressed
 *                  sparse column form (if `a_format` is `kColwise`, otherwise
 *                  compressed sparse row form). The sparse matrix consists of
 *                  three arrays, `astart`, `aindex`, and `avalue`. `astart` is
 *                  an array of length [numcol] containing the starting index of
 *                  each column in `aindex`. If `a_format` is `kRowwise` the
 *                  array is of length [numrow] corresponding to each row.
 * @param aindex    array of length [numnz] with indices of matrix entries
 * @param avalue    array of length [numnz] with values of matrix entries
 *
 * @param colvalue       array of length [numcol], filled with the primal column
 *                       solution
 * @param coldual        array of length [numcol], filled with the dual column
 *                       solution
 * @param rowvalue       array of length [numrow], filled with the primal row
 *                       solution
 * @param rowdual        array of length [numrow], filled with the dual row
 *                       solution
 * @param colbasisstatus array of length [numcol], filled with the basis status
 *                       of the columns
 * @param rowbasisstatus array of length [numrow], filled with the basis status
 *                       of the rows
 * @param modelstatus    termination status of the model after the solve
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_lpCall(const HighsInt numcol, const HighsInt numrow,
                      const HighsInt numnz, const HighsInt a_format,
                      const HighsInt sense, const double offset,
                      const double* colcost, const double* collower,
                      const double* colupper, const double* rowlower,
                      const double* rowupper, const HighsInt* astart,
                      const HighsInt* aindex, const double* avalue,
                      double* colvalue, double* coldual, double* rowvalue,
                      double* rowdual, HighsInt* colbasisstatus,
                      HighsInt* rowbasisstatus, HighsInt* modelstatus);

/**
 * Formulate and solve a mixed-integer linear program using HiGHS.
 *
 * The signature of this method is identical to `Highs_lpCall`, except that it
 * has an additional `integrality` argument, and that it is missing the
 * `coldual`, `rowdual`, `colbasisstatus` and `rowbasisstatus` arguments.
 *
 * @param integrality   array of length [numcol] indicating whether the columns
 *                      are continuous (0) or integer (1)
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_mipCall(const HighsInt numcol, const HighsInt numrow,
                       const HighsInt numnz, const HighsInt a_format,
                       const HighsInt sense, const double offset,
                       const double* colcost, const double* collower,
                       const double* colupper, const double* rowlower,
                       const double* rowupper, const HighsInt* astart,
                       const HighsInt* aindex, const double* avalue,
                       const HighsInt* integrality, double* colvalue,
                       double* rowvalue, HighsInt* modelstatus);

/**
 * Formulate and solve a quadratic program using HiGHS.
 *
 * The signature of this method is identical to `Highs_lpCall`, except that it
 * has additional arguments for specifying the Hessian matrix.

 * @param q_numnz   the number of nonzeros in the Hessian matrix
 * @param q_format  the format of the Hessian matrix. If q_numnz > 0, this
                    must be 1
 * @param qstart    the Hessian matrix is provided in the same format as the
 *                  constraint matrix, using `qstart`, `qindex`, and `qvalue` in
 *                  the place of `astart`, `aindex`, and `avalue`
 * @param qindex    array of length [q_numnz] with indices of matrix entries
 * @param qvalue    array of length [q_numnz] with values of matrix entries
  *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_qpCall(
    const HighsInt numcol, const HighsInt numrow, const HighsInt numnz,
    const HighsInt q_numnz, const HighsInt a_format, const HighsInt q_format,
    const HighsInt sense, const double offset, const double* colcost,
    const double* collower, const double* colupper, const double* rowlower,
    const double* rowupper, const HighsInt* astart, const HighsInt* aindex,
    const double* avalue, const HighsInt* qstart, const HighsInt* qindex,
    const double* qvalue, double* colvalue, double* coldual, double* rowvalue,
    double* rowdual, HighsInt* colbasisstatus, HighsInt* rowbasisstatus,
    HighsInt* modelstatus);

/**
 * Create a HiGHS model object and return the reference.
 *
 * Call `Highs_destroy` on the returned reference to clean up allocated memory.
 *
 * @returns A pointer to the HiGHS model object
 */
void* Highs_create(void);

/**
 * Destroy the model `highs` created by `Highs_create` and free all
 * corresponding memory. Future calls using `highs` are not allowed.
 *
 * To empty a model without invalidating `highs`, see `Highs_clearModel`.
 *
 * @param highs a pointer to the HiGHS model object
 */
void Highs_destroy(void* highs);

/**
 * Read a model from `filename` into `highs`.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the filename to read.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_readModel(void* highs, const char* filename);

/**
 * Write the model `highs` to `filename`.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the filename to write.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_writeModel(void* highs, const char* filename);

/**
 * Remove all variables and constraints from the model `highs`, but do not
 * invalidate the pointer `highs`. Future calls (for example, adding new
 * variables and constraints) are allowed.
 *
 * See `Highs_destroy` to clear the model and free all associated memory.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_clearModel(void* highs);

/**
 * Optimize a model. The algorithm used by HiGHS depends on the options that
 * have been set.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_run(void* highs);

/**
 * Write the solution information (including basis status if available) to a
 * file.
 *
 * See also: `Highs_writeSolutionPretty`.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the name of the file to write the results to
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_writeSolution(const void* highs, const char* filename);

/**
 * Write the solution information (including basis status if available) to a
 * file in a human-readable format.
 *
 * The method identical to `Highs_writeSolution`, except that the printout is in
 * a human-readiable format.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the name of the file to write the results to
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_writeSolutionPretty(const void* highs, const char* filename);

/**
 * Pass a linear program (LP) to HiGHS in a single function call.
 *
 * The signature of this function is identical to `Highs_passModel`, without the
 * arguments for passing the Hessian matrix of a quadratic program and the
 * integrality vector.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_passLp(void* highs, const HighsInt numcol, const HighsInt numrow,
                      const HighsInt numnz, const HighsInt a_format,
                      const HighsInt sense, const double offset,
                      const double* colcost, const double* collower,
                      const double* colupper, const double* rowlower,
                      const double* rowupper, const HighsInt* astart,
                      const HighsInt* aindex, const double* avalue);

/**
 * Pass a mixed-integer linear program (MILP) to HiGHS in a single function
 * call.
 *
 * The signature of function is identical to `Highs_passModel`, without the
 * arguments for passing the Hessian matrix of a quadratic program.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_passMip(void* highs, const HighsInt numcol,
                       const HighsInt numrow, const HighsInt numnz,
                       const HighsInt a_format, const HighsInt sense,
                       const double offset, const double* colcost,
                       const double* collower, const double* colupper,
                       const double* rowlower, const double* rowupper,
                       const HighsInt* astart, const HighsInt* aindex,
                       const double* avalue, const HighsInt* integrality);

/**
 * Pass a model to HiGHS in a single function call. This is faster than
 * constructing the model using `Highs_addRow` and `Highs_addCol`.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param numcol    the number of columns
 * @param numrow    the number of rows
 * @param numnz     the number of elements in the constraint matrix
 * @param q_num_nz  the number of elements in the Hessian matrix
 * @param a_format  the format of the constraint matrix to use. See
 *                  HighsMatrixFormat for more details.
 * @param q_format  the format of the Hessian matrix to use. See
 *                  HighsMatrixFormat for more details.
 * @param sense     the optimization sense (-1 for maximization, 1 for
 *                  minimization)
 * @param offset    the constant term in the objective function
 * @param colcost   array of length [numcol] with the objective coefficients
 * @param collower  array of length [numcol] with the lower column bounds
 * @param colupper  array of length [numcol] with the upper column bounds
 * @param rowlower  array of length [numrow] with the upper row bounds
 * @param rowupper  array of length [numrow] with the upper row bounds
 * @param astart    the constraint matrix is provided to HiGHS in compressed
 *                  sparse column form (if `a_format` is `kColwise`, otherwise
 *                  compressed sparse row form). The sparse matrix consists of
 *                  three arrays, `astart`, `aindex`, and `avalue`. `astart` is
 *                  an array of length [numcol] containing the starting index of
 *                  each column in `aindex`. If `a_format` is `kRowwise` the
 *                  array is of length [numrow] corresponding to each row.
 * @param aindex    array of length [numnz] with indices of matrix entries
 * @param avalue    array of length [numnz] with values of matrix entries
 * @param qstart    the Hessian matrix is provided in the same format as the
 *                  constraint matrix, using `qstart`, `qindex`, and `qvalue` in
 *                  the place of `astart`, `aindex`, and `avalue`. If the model
 *                  is linear, pass NULL.
 * @param qindex    array of length [q_numnz] with indices of matrix entries. If
 *                  the model is linear, pass NULL.
 * @param qvalue    array of length [q_numnz] with values of matrix entries. If
 *                  the model is linear, pass NULL.
 * @param integrality an array of length [numcol] indicating whether the column
 *                    is continuous (0) or integer (1). If the model is linear,
 *                    pass NULL.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_passModel(
    void* highs, const HighsInt numcol, const HighsInt numrow,
    const HighsInt numnz, const HighsInt q_num_nz, const HighsInt a_format,
    const HighsInt q_format, const HighsInt sense, const double offset,
    const double* colcost, const double* collower, const double* colupper,
    const double* rowlower, const double* rowupper, const HighsInt* astart,
    const HighsInt* aindex, const double* avalue, const HighsInt* qstart,
    const HighsInt* qindex, const double* qvalue, const HighsInt* integrality);

// TODO(odow): document this function
/**
 * Set the Hessian matrix for a quadratic objective.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param dim
 * @param num_nz    the number of non-zero elements in the Hessian matrix
 * @param format
 * @param start
 * @param index
 * @param value
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_passHessian(void* highs, const HighsInt dim,
                           const HighsInt num_nz, const HighsInt format,
                           const HighsInt* start, const HighsInt* index,
                           const double* value);

/**
 * Set a boolean-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setBoolOptionValue(void* highs, const char* option,
                                  const HighsInt value);

/**
 * Set an int-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setIntOptionValue(void* highs, const char* option,
                                 const HighsInt value);

/**
 * Set a double-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setDoubleOptionValue(void* highs, const char* option,
                                    const double value);

/**
 * Set a string-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setStringOptionValue(void* highs, const char* option,
                                    const char* value);

/**
 * Get a boolean-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     storage for the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBoolOptionValue(const void* highs, const char* option,
                                  HighsInt* value);

/**
 * Get an int-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     storage for the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getIntOptionValue(const void* highs, const char* option,
                                 HighsInt* value);

/**
 * Get a double-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     storage for the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getDoubleOptionValue(const void* highs, const char* option,
                                    double* value);

/**
 * Get a string-valued option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param value     pointer to allocated memory to store the value of the option
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getStringOptionValue(const void* highs, const char* option,
                                    char* value);

/**
 * Get the type expected by an option.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param option    the name of the option
 * @param type      int in which to store the option type. See `HighsOptionType`
 *                  for more details.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getOptionType(const void* highs, const char* option,
                             HighsInt* type);

/**
 * Reset all options to their default value.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_resetOptions(void* highs);

/**
 * Write the current options to file.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the filename to write the options to
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_writeOptions(const void* highs, const char* filename);

/**
 * Write the value of non-default options to file.
 *
 * This is similar to `Highs_writeOptions`, except only non-default options are
 * written to `filename`.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param filename  the filename to write the options to
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_writeOptionsDeviations(const void* highs, const char* filename);

/**
 * Get an int-valued info value.
 *
 * @param highs a pointer to the HiGHS model object
 * @param info  the name of the info
 * @param value a reference to an integer that the result will be stored in
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getIntInfoValue(const void* highs, const char* info,
                               HighsInt* value);

/**
 * Get a double-valued info value.
 *
 * @param highs a pointer to the HiGHS model object
 * @param info  the name of the info
 * @param value a reference to an double that the result will be stored in
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getDoubleInfoValue(const void* highs, const char* info,
                                  double* value);

/**
 * Get the primal and dual solution from an optimized model.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param colvalue  array of length [numcol], filled with primal column values
 * @param coldual   array of length [numcol], filled with dual column values
 * @param rowvalue  array of length [numrow], filled with primal row values
 * @param rowdual   array of length [numrow], filled with dual row values
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getSolution(const void* highs, double* colvalue, double* coldual,
                           double* rowvalue, double* rowdual);

/**
 * Given a linear program with a basic feasible solution, get the column and row
 * basis statuses.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param colstatus array of length [numcol], to be filled with the column basis
 *                  statuses
 * @param rowstatus array of length [numrow], to be filled with the row basis
 *                  statuses
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasis(const void* highs, HighsInt* colstatus,
                        HighsInt* rowstatus);

/**
 * Return the optimization status of the model.
 *
 * See `HighsModelStatus` for a list of possible returns.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns an integer corresponding to the `HighsModelStatus` enum
 */
HighsInt Highs_getModelStatus(const void* highs);

/**
 * Return the optimization status of the scaled model.
 *
 * See `HighsModelStatus` for a list of possible returns.
 *
 * Note that the status of the scaled model may not match the status of the
 * unscaled model as returned by `Highs_getModelStatus`.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns an integer corresponding to the `HighsModelStatus` enum
 */
HighsInt Highs_getScaledModelStatus(const void* highs);

/**
 * Get an unbounded dual ray that is a certificate of primal infeasibility.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param has_dual_ray      a pointer to an int to store 1 if the dual ray
 *                          exists
 * @param dual_ray_value    an array of length [numrow] filled with the
 *                          unbounded ray
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getDualRay(const void* highs, HighsInt* has_dual_ray,
                          double* dual_ray_value);

/**
 * Get an unbounded primal ray that is a certificate of dual infeasibility.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param has_primal_ray    a pointer to an int to store 1 if the primal ray
 *                          exists
 * @param primal_ray_value  an array of length [numcol] filled with the
 *                          unbounded ray
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getPrimalRay(const void* highs, HighsInt* has_primal_ray,
                            double* primal_ray_value);

/**
 * Get the primal objective function value.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the primal objective function value
 */
double Highs_getObjectiveValue(const void* highs);

/**
 * Get the indices of the rows and columns that make up the basis matrix of a
 * basic feasible solution.
 *
 * Non-negative entries are indices of columns, and negative entries are
 * `-row_index - 1`. For example, `{1, -1}` would be the second column and first
 * row.
 *
 * The order of these rows and columns is important for calls to the functions:
 *  - `Highs_getBasisInverseRow`
 *  - `Highs_getBasisInverseCol`
 *  - `Highs_getBasisSolve`
 *  - `Highs_getBasisTransposeSolve`
 *  - `Highs_getReducedRow`
 *  - `Highs_getReducedColumn`
 *
 * @param highs             a pointer to the HiGHS model object
 * @param basic_variables   array of size [num_rows], filled with the indices of
 *                          the basic variables
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasicVariables(const void* highs, HighsInt* basic_variables);

/**
 * Get a row of the inverse basis matrix \f$B^{-1}\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `row_vector` and `row_indices` must have an allocated length of
 * [numrow]. However, check `row_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param row           index of the row to compute
 * @param row_vector    values of the non-zero elements
 * @param row_num_nz    the number of non-zeros in the row
 * @param row_indices   indices of the non-zero elements
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasisInverseRow(const void* highs, const HighsInt row,
                                  double* row_vector, HighsInt* row_num_nz,
                                  HighsInt* row_indices);

/**
 * Get a column of the inverse basis matrix \f$B^{-1}\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `col_vector` and `col_indices` must have an allocated length of
 * [numrow]. However, check `col_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param col           index of the column to compute
 * @param col_vector    values of the non-zero elements
 * @param col_num_nz    the number of non-zeros in the column
 * @param col_indices   indices of the non-zero elements

 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasisInverseCol(const void* highs, const HighsInt col,
                                  double* col_vector, HighsInt* col_num_nz,
                                  HighsInt* col_indices);

/**
 * Compute \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `solution_vector` and `solution_indices` must have an allocated
 * length of [numrow]. However, check `solution_num_nz` to see how many
 * non-zero elements are actually stored.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param rhs               the right-hand side vector `b`
 * @param solution_vector   values of the non-zero elements
 * @param solution_num_nz   the number of non-zeros in the solution
 * @param solution_indices  indices of the non-zero elements
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasisSolve(const void* highs, const double* rhs,
                             double* solution_vector, HighsInt* solution_num_nz,
                             HighsInt* solution_indices);

/**
 * Compute \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `solution_vector` and `solution_indices` must have an allocated
 * length of [numrow]. However, check `solution_num_nz` to see how many
 * non-zero elements are actually stored.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param rhs               the right-hand side vector `b`
 * @param solution_vector   values of the non-zero elements
 * @param solution_num_nz   the number of non-zeros in the solution
 * @param solution_indices  indices of the non-zero elements
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getBasisTransposeSolve(const void* highs, const double* rhs,
                                      double* solution_vector,
                                      HighsInt* solution_nz,
                                      HighsInt* solution_indices);

/**
 * Compute a row of \f$B^{-1}A\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `row_vector` and `row_indices` must have an allocated length of
 * [numrow]. However, check `row_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param row           index of the row to compute
 * @param row_vector    values of the non-zero elements
 * @param row_num_nz    the number of non-zeros in the row
 * @param row_indices   indices of the non-zero elements
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getReducedRow(const void* highs, const HighsInt row,
                             double* row_vector, HighsInt* row_num_nz,
                             HighsInt* row_indices);

/**
 * Compute a column of \f$B^{-1}A\f$.
 *
 * See `Highs_getBasicVariables` for a description of the `B` matrix.
 *
 * The arrays `col_vector` and `col_indices` must have an allocated length of
 * [numrow]. However, check `col_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param col           index of the column to compute
 * @param col_vector    values of the non-zero elements
 * @param col_num_nz    the number of non-zeros in the column
 * @param col_indices   indices of the non-zero elements

 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getReducedColumn(const void* highs, const HighsInt col,
                                double* col_vector, HighsInt* col_num_nz,
                                HighsInt* col_indices);

/**
 * Set a basic feasible solution by passing the column and row basis statuses to
 * the model.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param colstatus an array of length [numcol] with the column basis status
 * @param rowstatus an array of length [numrow] with the row basis status
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setBasis(void* highs, const HighsInt* colstatus,
                        const HighsInt* rowstatus);

/**
 * Set a logical basis in the model.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_setLogicalBasis(void* highs);

/**
 * Return the cumulative wall-clock time spent in `Highs_run`.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the cumulative wall-clock time spent in `Highs_run`
 */
double Highs_getRunTime(const void* highs);

/**
 * Add a new row (a linear constraint) to the model.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param lower         lower bound of the row
 * @param upper         upper bound of the row
 * @param num_new_nz    number of non-zeros in the row
 * @param indices       array of size [num_new_nz] with column indices
 * @param values        array of size [num_new_nz] with column values
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_addRow(void* highs, const double lower, const double upper,
                      const HighsInt num_new_nz, const HighsInt* indices,
                      const double* values);

/**
 * Add multiple rows (linear constraints) to the model.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param num_new_row   the number of new rows to add
 * @param lower         array of size [num_new_row] with the lower bounds of the
 *                      rows
 * @param upper         array of size [num_new_row] with the upper bounds of the
 *                      rows
 * @param num_new_nz    number of non-zeros in the rows
 * @param starts        the constraint coefficients are given as a matrix in
 *                      compressed sparse row form by the arrays `starts`,
 *                      `indices`, and `values`. `starts` is an array of size
 *                      [num_new_rows] with the start index of each row in
 *                      indices and values.
 * @param indices       array of size [num_new_nz] with column indices
 * @param values        array of size [num_new_nz] with column values
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_addRows(void* highs, const HighsInt num_new_row,
                       const double* lower, const double* upper,
                       const HighsInt num_new_nz, const HighsInt* starts,
                       const HighsInt* indices, const double* values);

/**
 * Add a new column (variable) to the model.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param cost          objective coefficient of the column
 * @param lower         lower bound of the column
 * @param upper         upper bound of the column
 * @param num_new_nz    number of non-zeros in the column
 * @param indices       array of size [num_new_nz] with the row indices
 * @param values        array of size [num_new_nz] with row values
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_addCol(void* highs, const double cost, const double lower,
                      const double upper, const HighsInt num_new_nz,
                      const HighsInt* indices, const double* values);

/**
 * Add multiple columns (linear constraints) to the model.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param num_new_col   number of new columns to add
 * @param costs         array of size [num_new_col] with objective coefficients
 * @param lower         array of size [num_new_col] with lower bounds
 * @param upper         array of size [num_new_col] with upper bounds
 * @param num_new_nz    number of new nonzeros in the constraint matrix
 * @param starts        the constraint coefficients are given as a matrix in
 *                      compressed sparse column form by the arrays `starts`,
 *                      `indices`, and `values`. `starts` is an array of size
 *                      [num_new_cols] with the start index of each row in
 *                      indices and values.
 * @param indices       array of size [num_new_nz] with row indices
 * @param values        array of size [num_new_nz] with row values
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_addCols(void* highs, const HighsInt num_new_col,
                       const double* costs, const double* lower,
                       const double* upper, const HighsInt num_new_nz,
                       const HighsInt* starts, const HighsInt* indices,
                       const double* values);

/**
 * Change the objective sense of the model.
 *
 * @param highs a pointer to the HiGHS model object
 * @param sense the new optimization sense. Use -1 for maximization and 1 for
 *              minimization.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeObjectiveSense(void* highs, const HighsInt sense);

/**
 * Change the objective offset of the model.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param offset    the new objective offset
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeObjectiveOffset(void* highs, const double offset);

/**
 * Change the integrality of a column.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param col           the column index to change
 * @param integrality   the new integrality of the column. Use 0 for continuous
 *                      and 1 for integer.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColIntegrality(void* highs, const HighsInt col,
                                    const HighsInt integrality);

/**
 * Change the integrality of multiple adjacent columns.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param from_col      the index of the first column whose integrality changes
 * @param to_col        the index of the last column whose integrality
 *                      changes
 * @param integrality   an array of length [to_col - from_col + 1] with the new
 *                      integralities of the columns
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsIntegralityByRange(void* highs,
                                            const HighsInt from_col,
                                            const HighsInt to_col,
                                            const HighsInt* integrality);

/**
 * Change the integrality of multiple columns given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of columns to change
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the columns to change
 * @param integrality       an array of length [num_set_entries] with the new
 *                          integralities of the columns.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsIntegralityBySet(void* highs,
                                          const HighsInt num_set_entries,
                                          const HighsInt* set,
                                          const HighsInt* integrality);

/**
 * Change the integrality of multiple columns given by a mask.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param mask          an array of length [numcol] with 1 if the column
 *                      integrality should be changed and 0 otherwise
 * @param integrality   an array of length [numcol] with the new
 *                      integralities of the columns.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsIntegralityByMask(void* highs, const HighsInt* mask,
                                           const HighsInt* integrality);

/**
 * Change the objective coefficient of a column.
 *
 * @param highs a pointer to the HiGHS model object
 * @param col   the index of the column fo change
 * @param cost  the new objective coefficient
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColCost(void* highs, const HighsInt col,
                             const double cost);

/**
 * Change the cost coefficients of multiple adjacent columns.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param from_col  the index of the first column whose cost changes
 * @param to_col    the index of the last column whose cost changes
 * @param cost      an array of length [to_col - from_col + 1] with the new
 *                  objective coefficients
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsCostByRange(void* highs, const HighsInt from_col,
                                     const HighsInt to_col, const double* cost);

/**
 * Change the cost of multiple columns given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of columns to change
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the columns to change
 * @param cost              an array of length [num_set_entries] with the new
 *                          costs of the columns.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsCostBySet(void* highs, const HighsInt num_set_entries,
                                   const HighsInt* set, const double* cost);

/**
 * Change the cost of multiple columns given by a mask.
 *
 * @param highs a pointer to the HiGHS model object
 * @param mask  an array of length [numcol] with 1 if the column
 *              cost should be changed and 0 otherwise
 * @param cost  an array of length [numcol] with the new costs
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsCostByMask(void* highs, const HighsInt* mask,
                                    const double* cost);

/**
 * Change the variable bounds of a column.
 *
 * @param highs a pointer to the HiGHS model object
 * @param col   the index of the column whose bounds are to change
 * @param lower the new lower bound
 * @param upper the new upper bound
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColBounds(void* highs, const HighsInt col,
                               const double lower, const double upper);

/**
 * Change the variable bounds of multiple adjacent columns.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param from_col  the index of the first column whose bound changes
 * @param to_col    the index of the last column whose bound changes
 * @param lower     an array of length [to_col - from_col + 1] with the new
 *                  lower bounds
 * @param upper     an array of length [to_col - from_col + 1] with the new
 *                  upper bounds
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsBoundsByRange(void* highs, const HighsInt from_col,
                                       const HighsInt to_col,
                                       const double* lower,
                                       const double* upper);

/**
 * Change the bounds of multiple columns given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of columns to change
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the columns to change
 * @param lower             an array of length [num_set_entries] with the new
 *                          lower bounds
 * @param upper             an array of length [num_set_entries] with the new
 *                          upper bounds
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper);

/**
 * Change the variable bounds of multiple columns given by a mask.
 *
 * @param highs a pointer to the HiGHS model object
 * @param mask  an array of length [numcol] with 1 if the column
 *              bounds should be changed and 0 otherwise
 * @param lower an array of length [numcol] with the new lower bounds
 * @param upper an array of length [numcol] with the new upper bounds
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeColsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower, const double* upper);

/**
 * Change the bounds of a row.
 *
 * @param highs a pointer to the HiGHS model object
 * @param row   the index of the row whose bounds are to change
 * @param lower the new lower bound
 * @param upper the new upper bound
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeRowBounds(void* highs, const HighsInt row,
                               const double lower, const double upper);

/**
 * Change the bounds of multiple rows given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of rows to change
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the rows to change
 * @param lower             an array of length [num_set_entries] with the new
 *                          lower bounds
 * @param upper             an array of length [num_set_entries] with the new
 *                          upper bounds
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeRowsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper);

/**
 * Change the bounds of multiple rows given by a mask.
 *
 * @param highs a pointer to the HiGHS model object
 * @param mask  an array of length [numrow] with 1 if the row
 *              bounds should be changed and 0 otherwise
 * @param lower an array of length [numrow] with the new lower bounds
 * @param upper an array of length [numrow] with the new upper bounds
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeRowsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower, const double* upper);

/**
 * Change a coefficient in the constraint matrix.
 *
 * @param highs a pointer to the HiGHS model object
 * @param row the index of the row to change
 * @param col the index of the col to change
 * @param value the new constraint coefficient
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_changeCoeff(void* highs, const HighsInt row, const HighsInt col,
                           const double value);

/**
 * Get the objective sense.
 *
 * @param highs a pointer to the HiGHS model object
 * @param sense stores the current objective sense
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getObjectiveSense(const void* highs, HighsInt* sense);

/**
 * Get the objective offset.
 *
 * @param highs a pointer to the HiGHS model object
 * @param offset stores the current objective offset
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getObjectiveOffset(const void* highs, double* offset);

/**
 * Get data associated with multiple adjacent columns from the model.
 *
 * To query the constraint coefficients, this function should be called twice:
 *  - First, call this function with `matrix_start`, `matrix_index`, and
 *    `matrix_value` as `NULL`. This call will populate `num_nz` with the
 *    number of nonzero elements in the corresponding section of the constraint
 *    matrix.
 *  - Second, allocate new `matrix_index` and `matrix_value` arrays of length
 *    `num_nz` and call this function again to populate the new arrays with
 *    their contents.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param from_col      the first column for which to query data for
 * @param to_col        the last column (inclusive) for which to query data for
 * @param num_col       an integer populated with the number of columns got from
 *                      the model (this should equal `to_col - from_col + 1`)
 * @param costs         array of size [to_col - from_col + 1] for the column
 *                      cost coefficients
 * @param lower         array of size [to_col - from_col + 1] for the column
 *                      lower bounds
 * @param upper         array of size [to_col - from_col + 1] for the column
 *                      upper bounds
 * @param num_nz        an integer populated with the number of non-zero
 *                      elements in the constraint matrix
 * @param matrix_start  array of size [to_col - from_col + 1] with the start
 *                      indices of each
 *                      column in `matrix_index` and `matrix_value`
 * @param matrix_index  array of size [num_nz] with the row indices of each
 *                      element in the constraint matrix
 * @param matrix_value  array of size [num_nz] with the non-zero elements of the
 *                      constraint matrix.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getColsByRange(const void* highs, const HighsInt from_col,
                              const HighsInt to_col, HighsInt* num_col,
                              double* costs, double* lower, double* upper,
                              HighsInt* num_nz, HighsInt* matrix_start,
                              HighsInt* matrix_index, double* matrix_value);

/**
 * Get data associated with multiple columns given by an array.
 *
 * This function is identical to `Highs_getColsByRange`, except for how the
 * columns are specified.
 *
 * @param num_set_indices   the number of indices in the set
 * @param set               array of size [num_set_entries] with the column
 *                          indices to get
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getColsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_col,
                            double* costs, double* lower, double* upper,
                            HighsInt* num_nz, HighsInt* matrix_start,
                            HighsInt* matrix_index, double* matrix_value);

/**
 * Get data associated with multiple columns given by a mask.
 *
 * This function is identical to `Highs_getColsByRange`, except for how the
 * columns are specified.
 *
 * @param mask  array of length [num_col] containing a 1 to get the column and 0
 *              otherwise
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getColsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_col, double* costs, double* lower,
                             double* upper, HighsInt* num_nz,
                             HighsInt* matrix_start, HighsInt* matrix_index,
                             double* matrix_value);

/**
 * Get data associated with multiple adjacent rows from the model.
 *
 * To query the constraint coefficients, this function should be called twice:
 *  - First, call this function with `matrix_start`, `matrix_index`, and
 *    `matrix_value` as `NULL`. This call will populate `num_nz` with the
 *    number of nonzero elements in the corresponding section of the constraint
 *    matrix.
 *  - Second, allocate new `matrix_index` and `matrix_value` arrays of length
 *    `num_nz` and call this function again to populate the new arrays with
 *    their contents.
 *
 * @param highs         a pointer to the HiGHS model object
 * @param from_row      the first row for which to query data for
 * @param to_row        the last row (inclusive) for which to query data for
 * @param num_row       an integer populated with the number of row got from the
 *                      model
 * @param lower         array of size [to_row - from_row + 1] for the row lower
 *                      bounds
 * @param upper         array of size [to_row - from_row + 1] for the row upper
 *                      bounds
 * @param num_nz        an integer populated with the number of non-zero
 *                      elements in the constraint matrix
 * @param matrix_start  array of size [to_row - from_row + 1] with the start
 *                      indices of each row in `matrix_index` and `matrix_value`
 * @param matrix_index  array of size [num_nz] with the column indices of each
 *                      element in the constraint matrix
 * @param matrix_value  array of size [num_nz] with the non-zero elements of the
 *                      constraint matrix.
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getRowsByRange(const void* highs, const HighsInt from_row,
                              const HighsInt to_row, HighsInt* num_row,
                              double* lower, double* upper, HighsInt* num_nz,
                              HighsInt* matrix_start, HighsInt* matrix_index,
                              double* matrix_value);

/**
 * Get data associated with multiple rows given by an array.
 *
 * This function is identical to `Highs_getRowsByRange`, except for how the
 * rows are specified.
 *
 * @param num_set_indices   the number of indices in the set
 * @param set               array of size [num_set_entries] with the row indices
 *                          to get
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getRowsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_row,
                            double* lower, double* upper, HighsInt* num_nz,
                            HighsInt* matrix_start, HighsInt* matrix_index,
                            double* matrix_value);

/**
 * Get data associated with multiple rows given by a mask.
 *
 * This function is identical to `Highs_getRowsByRange`, except for how the
 * rows are specified.
 *
 * @param mask  array of length [num_row] containing a 1 to get the row and 0
 *              otherwise
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getRowsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_row, double* lower, double* upper,
                             HighsInt* num_nz, HighsInt* matrix_start,
                             HighsInt* matrix_index, double* matrix_value);

/**
 * Delete multiple adjacent columns.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param from_col  the index of the first column to delete
 * @param to_col    the index of the last column to delete
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteColsByRange(void* highs, const HighsInt from_col,
                                 const HighsInt to_col);

/**
 * Delete multiple columns given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of columns to delete
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the columns to delete
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteColsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set);

/**
 * Delete multiple columns given by a mask.
 *
 * @param highs a pointer to the HiGHS model object
 * @param mask  an array of length [numcol] with 1 if the column
 *              should be deleted and 0 otherwise
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteColsByMask(void* highs, HighsInt* mask);

/**
 * Delete multiple adjacent rows.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param from_row  the index of the first row to delete
 * @param to_row    the index of the last row to delete
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteRowsByRange(void* highs, const int from_row,
                                 const HighsInt to_row);

/**
 * Delete multiple rows given by an array of indices.
 *
 * @param highs             a pointer to the HiGHS model object
 * @param num_set_entries   the number of rows to delete
 * @param set               an array of size [num_set_entries] with the indices
 *                          of the rows to delete
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteRowsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set);

/**
 * Delete multiple rows given by a mask.
 *
 * @param highs a pointer to the HiGHS model object
 * @param mask  an array of length [numrow] with 1 if the row should be deleted
 *              and 0 otherwise
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_deleteRowsByMask(void* highs, HighsInt* mask);

/**
 * Scale a column by a constant.
 *
 * Scaling a column modifies the elements in the constraint matrix, the variable
 * bounds, and the objective coefficient.
 *
 * If scaleval < 0, the variable bounds flipped.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param col       the index of the column to scale
 * @param scaleval  the value by which to scale the column
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_scaleCol(void* highs, const HighsInt col, const double scaleval);

/**
 * Scale a row by a constant.
 *
 * If scaleval < 0, the row bounds are flipped.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param row       the index of the row to scale
 * @param scaleval  the value by which to scale the row
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_scaleRow(void* highs, const HighsInt row, const double scaleval);

/**
 * Return the value of infinity used by HiGHS.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the value of infinity used by HiGHS
 */
double Highs_getInfinity(const void* highs);

/**
 * Return the number of columns in the model.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the number of columns in the model
 */
HighsInt Highs_getNumCol(const void* highs);

/**
 * Return the number of rows in the model.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the number of rows in the model.
 */
HighsInt Highs_getNumRow(const void* highs);

/**
 * Return the number of nonzeros in the constraint matrix of the model.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the number of nonzeros in the constraint matrix of the model.
 */
HighsInt Highs_getNumNz(const void* highs);

/**
 * Return the number of nonzeroes in the Hessian matrix of the model.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns the number of nonzeroes in the Hessian matrix of the model.
 */
HighsInt Highs_getHessianNumNz(const void* highs);

/**
 * Get the data from a HiGHS model.
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getNumCol`
 *  - `Highs_getNumRow`
 *  - `Highs_getNumNz`
 *  - `Highs_getHessianNumNz`
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getModel(const void* highs, const HighsInt a_format,
                        const HighsInt q_format, HighsInt* numcol,
                        HighsInt* numrow, HighsInt* numnz,
                        HighsInt* hessian_num_nz, HighsInt* sense,
                        double* offset, double* colcost, double* collower,
                        double* colupper, double* rowlower, double* rowupper,
                        HighsInt* astart, HighsInt* aindex, double* avalue,
                        HighsInt* qstart, HighsInt* qindex, double* qvalue,
                        HighsInt* integrality);

// Fails on Windows and MacOS since string_model_status is destroyed
// after the method returns, so what's returned is a pointer to
// something that no longer exists.
//
// /**
//  * @brief Returns a pointer to a character representation of a model
//  * status
//  */
// const char* Highs_modelStatusToChar(
//     void* highs,
//     HighsInt int_model_status  //!< Status to interpret
// );
//
// Fails on Windows and MacOS since string_solution_status is
// destroyed after the method returns, so what's returned is a pointer
// to something that no longer exists.
//
// /**
//  * @brief Returns a pointer to a character representation of a
//  * solution status
//  */
// const char* Highs_solutionStatusToChar(
//     void* highs,
//     HighsInt int_solution_status  //!< Status to interpret
// );

/**
 * Given a model solved with an interior point method, run crossover to compute
 * a basic feasible solution.
 *
 * @param highs a pointer to the HiGHS model object
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_crossover(void* highs);

/**
 * Set a primal-dual solution as a starting point, then run crossover to compute
 * a basic feasible solution.
 *
 * @param highs     a pointer to the HiGHS model object
 * @param n         the number of variables
 * @param m         the number of rows
 * @param col_value array of length [m] with optimal primal solution for each
 *                  column
 * @param col_dual  array of length [n] with optimal dual solution for each
 *                  column
 * @param row_dual  array of length [m] with optimal dual solution for each row
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_crossover_set(void* highs, const int n, const int m,
                             double* col_value, double* col_dual,
                             double* row_dual);

// *********************
// * Deprecated methods*
// *********************

HighsInt Highs_call(const HighsInt numcol, const HighsInt numrow,
                    const HighsInt numnz, const double* colcost,
                    const double* collower, const double* colupper,
                    const double* rowlower, const double* rowupper,
                    const HighsInt* astart, const HighsInt* aindex,
                    const double* avalue, double* colvalue, double* coldual,
                    double* rowvalue, double* rowdual, HighsInt* colbasisstatus,
                    HighsInt* rowbasisstatus, HighsInt* modelstatus);

HighsInt Highs_runQuiet(void* highs);

HighsInt Highs_setHighsLogfile(void* highs, const void* logfile);

HighsInt Highs_setHighsOutput(void* highs, const void* outputfile);

HighsInt Highs_getIterationCount(const void* highs);

HighsInt Highs_getSimplexIterationCount(const void* highs);

HighsInt Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                       const HighsInt value);

// void Highs_getLp(
//     void *highs,       //!< HiGHS object reference
//     int* numcol,        //!< number of columns
//     int* numrow,        //!< number of rows
//     int* numnz,         //!< number of entries in the constraint matrix
//     double *colcost,   //!< array of length [numcol] with column costs
//     double *collower,  //!< array of length [numcol] with lower column bounds
//     double *colupper,  //!< array of length [numcol] with upper column bounds
//     double *rowlower,  //!< array of length [numrow] with lower row bounds
//     double *rowupper,  //!< array of length [numrow] with upper row bounds
//     int *astart,       //!< array of length [numcol] with column start
//     indices int *
//         aindex,  //!< array of length [numnz] with row indices of matrix
//         entries
//     double *avalue  //!< array of length [numnz] with value of matrix entries
// );
HighsInt Highs_setHighsIntOptionValue(void* highs, const char* option,
                                      const HighsInt value);

HighsInt Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                         const double value);

HighsInt Highs_setHighsStringOptionValue(void* highs, const char* option,
                                         const char* value);

HighsInt Highs_setHighsOptionValue(void* highs, const char* option,
                                   const char* value);

HighsInt Highs_getHighsBoolOptionValue(const void* highs, const char* option,
                                       HighsInt* value);

HighsInt Highs_getHighsIntOptionValue(const void* highs, const char* option,
                                      HighsInt* value);

HighsInt Highs_getHighsDoubleOptionValue(const void* highs, const char* option,
                                         double* value);

HighsInt Highs_getHighsStringOptionValue(const void* highs, const char* option,
                                         char* value);

HighsInt Highs_getHighsOptionType(const void* highs, const char* option,
                                  HighsInt* type);

HighsInt Highs_resetHighsOptions(void* highs);

HighsInt Highs_getHighsIntInfoValue(const void* highs, const char* info,
                                    HighsInt* value);

HighsInt Highs_getHighsDoubleInfoValue(const void* highs, const char* info,
                                       double* value);

HighsInt Highs_getNumCols(const void* highs);

HighsInt Highs_getNumRows(const void* highs);

double Highs_getHighsInfinity(const void* highs);

double Highs_getHighsRunTime(const void* highs);

// const char* Highs_highsModelStatusToChar(void* highs,
//                                          HighsInt int_model_status);

HighsInt Highs_setOptionValue(void* highs, const char* option,
                              const char* value);

#ifdef __cplusplus
}
#endif

#endif
