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

/*
 * @brief solves an LP using HiGHS
 */
HighsInt Highs_lpCall(
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt*
        astart,  //!< array of length [numcol] with column start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with row indices of matrix entries
    const double*
        avalue,        //!< array of length [numnz] with value of matrix entries
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* coldual,   //!< array of length [numcol], filled with column duals
    double* rowvalue,  //!< array of length [numrow], filled with row values
    double* rowdual,   //!< array of length [numrow], filled with row duals
    HighsInt* colbasisstatus,  //!< array of length [numcol], filled with column
                               //!< basis status
    HighsInt* rowbasisstatus,  //!< array of length [numrow], filled with row
                               //!< basis status
    HighsInt* modelstatus      //!< status of the model will be saved here
);

/*
 * @brief solves a MIP using HiGHS
 */
HighsInt Highs_mipCall(
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt*
        astart,  //!< array of length [numcol] with column start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with row indices of matrix entries
    const double*
        avalue,  //!< array of length [numnz] with value of matrix entries
    const HighsInt*
        integrality,   //!< array of length [numcol] indicating whether
                       //!< variables are continuous (0) or integer (1)
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* rowvalue,  //!< array of length [numrow], filled with row values
    HighsInt* modelstatus  //!< status of the model will be saved here
);

/*
 * @brief solves a QP using HiGHS
 */
HighsInt Highs_qpCall(
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt q_numnz,   //!< number of entries in the Hessian matrix
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt q_format,  //!< format of the Hessian matrix (if numnz>0 it
                              //!< must be 1; Otherwise N/A)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt*
        astart,  //!< array of length [numcol] with column start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with row indices of matrix entries
    const double*
        avalue,  //!< array of length [numnz] with value of matrix entries
    const HighsInt* qstart,  //!< array of length [numcol] with Hessian column
                             //!< start indices
    const HighsInt* qindex,  //!< array of length [numnz] with row indices of
                             //!< Hessian matrix entries
    const double* qvalue,    //!< array of length [numnz] with value of Hessian
                             //!< matrix entries
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* coldual,   //!< array of length [numcol], filled with column duals
    double* rowvalue,  //!< array of length [numrow], filled with row values
    double* rowdual,   //!< array of length [numrow], filled with row duals
    HighsInt* colbasisstatus,  //!< array of length [numcol], filled with column
                               //!< basis status
    HighsInt* rowbasisstatus,  //!< array of length [numrow], filled with row
                               //!< basis status
    HighsInt* modelstatus      //!< status of the model will be saved here
);

/*
 * @brief creates a HiGHS object and returns the reference
 */
void* Highs_create(void);

/*
 * @brief destroys a HiGHS object
 */
void Highs_destroy(void* highs);

/*
 * @brief
 */
HighsInt Highs_readModel(void* highs,
                         const char* filename  //!< filename
);

/*
 * @brief
 */
HighsInt Highs_writeModel(void* highs,
                          const char* filename  //!< filename
);

/*
 * @brief
 */
HighsInt Highs_clearModel(void* highs);

/*
 * @brief Runs HiGHS
 */
HighsInt Highs_run(void* highs);

/*
 * @brief Reports the solution and basis status
 */
HighsInt Highs_writeSolution(const void* highs,
                             const char* filename  //!< filename
);

/*
 * @brief Reports the solution and basis status in a human-readable fashion
 */
HighsInt Highs_writeSolutionPretty(const void* highs,
                                   const char* filename  //!< filename
);

/*
 * @brief pass an LP to HiGHS
 */
HighsInt Highs_passLp(
    void* highs,
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt* astart,  //!< array of length [numcol or numrow if(rowwise)]
                             //!< with start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with indices of matrix entries
    const double*
        avalue  //!< array of length [numnz] with value of matrix entries
);

/*
 * @brief pass a MIP to HiGHS
 */
HighsInt Highs_passMip(
    void* highs,
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt* astart,  //!< array of length [numcol or numrow if(rowwise)]
                             //!< with start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with indices of matrix entries
    const double*
        avalue,  //!< array of length [numnz] with value of matrix entries
    const HighsInt*
        integrality  //!< array of length [numcol] indicating whether
                     //!< variables are continuous (0) or integer (1)
);

/*
 * @brief pass a general model to HiGHS
 */
HighsInt Highs_passModel(
    void* highs,
    const HighsInt numcol,    //!< number of columns
    const HighsInt numrow,    //!< number of rows
    const HighsInt numnz,     //!< number of entries in the constraint matrix
    const HighsInt q_num_nz,  //!< number of nonzeros in Hessian
    const HighsInt a_format,  //!< format of the constraint matrix (2 =>
                              //!< Row-wise; Otherise column-wise)
    const HighsInt q_format,  //!< format of the Hessian matrix (if numnz>0 it
                              //!< must be 1; Otherwise N/A)
    const HighsInt sense,   //!< sense of the optimization (-1 => Maximization;
                            //!< Otherwise minimization)
    const double offset,    //!< objective constant
    const double* colcost,  //!< array of length [numcol] with column costs
    const double*
        collower,  //!< array of length [numcol] with lower column bounds
    const double*
        colupper,  //!< array of length [numcol] with upper column bounds
    const double* rowlower,  //!< array of length [numrow] with lower row bounds
    const double* rowupper,  //!< array of length [numrow] with upper row bounds
    const HighsInt* astart,  //!< array of length [numcol or numrow if(rowwise)]
                             //!< with start indices
    const HighsInt*
        aindex,  //!< array of length [numnz] with indices of matrix entries
    const double*
        avalue,  //!< array of length [numnz] with value of matrix entries
    const HighsInt* qstart,  //!< array of length [numcol] with Hessian start
                             //!< indices - or NULL if model is linear
    const HighsInt* qindex,  //!< array of length [q_num_nz] with indices of
                             //!< Hessian entries - or NULL if model is linear
    const double* qvalue,    //!< array of length [q_num_nz] with values of
                             //!< Hessian entries - or NULL if model is linear
    const HighsInt*
        integrality  //!< array of length [numcol] indicating whether
                     //!< variables are continuous (0) or integer (1) - or NULL
                     //!< if model is continuous
);

/*
 * @brief pass a hessian to HiGHS
 */
HighsInt Highs_passHessian(void* highs, const HighsInt dim,
                           const HighsInt num_nz, const HighsInt format,
                           const HighsInt* start, const HighsInt* index,
                           const double* value);

HighsInt Highs_setBoolOptionValue(void* highs,
                                  const char* option,   //!< name of the option
                                  const HighsInt value  //!< new value of option
);

HighsInt Highs_setIntOptionValue(void* highs,
                                 const char* option,   //!< name of the option
                                 const HighsInt value  //!< new value of option
);

HighsInt Highs_setDoubleOptionValue(void* highs,
                                    const char* option,  //!< name of the option
                                    const double value  //!< new value of option
);

HighsInt Highs_setStringOptionValue(void* highs,
                                    const char* option,  //!< name of the option
                                    const char* value  //!< new value of option
);

HighsInt Highs_getBoolOptionValue(const void* highs,
                                  const char* option,  //!< name of the option
                                  HighsInt* value      //!< value of option
);

HighsInt Highs_getIntOptionValue(const void* highs,
                                 const char* option,  //!< name of the option
                                 HighsInt* value      //!< value of option
);

HighsInt Highs_getDoubleOptionValue(const void* highs,
                                    const char* option,  //!< name of the option
                                    double* value        //!< value of option
);

HighsInt Highs_getStringOptionValue(
    const void* highs,
    const char* option,  //!< name of the option
    char* value  //!< pointer to allocated memory to store value of option
);

/*
 * @brief Get the type expected by an option
 */
HighsInt Highs_getOptionType(const void* highs,
                             const char* option,  //!< The name of the option
                             HighsInt* type       //!< The type of the option.
);

/*
 * @brief
 */
HighsInt Highs_resetOptions(void* highs);

HighsInt Highs_writeOptions(const void* highs, const char* filename);

HighsInt Highs_writeOptionsDeviations(const void* highs, const char* filename);

/*
 * @brief
 */
HighsInt Highs_getIntInfoValue(const void* highs,
                               const char* info,  //!< The info name
                               HighsInt* value    //!< The info value
);

/*
 * @brief
 */
HighsInt Highs_getDoubleInfoValue(const void* highs,
                                  const char* info,  //!< The info name
                                  double* value      //!< The info value
);

/*
 * Get an int64-valued info value.
 *
 * @param highs a pointer to the HiGHS model object
 * @param info  the name of the info
 * @param value a reference to a int64 that the result will be stored in
 *
 * @returns A non-zero return value indicates that a problem occured
 */
HighsInt Highs_getInt64InfoValue(const void* highs, const char* info,
                                 int64_t* value);

/*
 * @brief
 */
HighsInt Highs_getSolution(
    const void* highs,
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* coldual,   //!< array of length [numcol], filled with column duals
    double* rowvalue,  //!< array of length [numrow], filled with row values
    double* rowdual    //!< array of length [numrow], filled with row duals
);

/*
 * @brief
 */
HighsInt Highs_getBasis(
    const void* highs,
    HighsInt* colstatus,  //!< array of length [numcol], filled
                          //!< with column basis stati
    HighsInt* rowstatus   //!< array of length [numrow], filled
                          //!< with row basis stati
);

/**
 * @brief Returns the status of the model
 */
HighsInt Highs_getModelStatus(const void* highs);

/**
 * @brief Returns the status of the scaled model
 */
HighsInt Highs_getScaledModelStatus(const void* highs);

/**
 * @brief Returns an unbounded dual ray that is a certificate of primal
 * infeasibility.
 */
HighsInt Highs_getDualRay(
    const void* highs,
    HighsInt* has_dual_ray,  //!< TRUE if the dual ray exists
    double* dual_ray_value   //!< array of length [numrow],
                             //!< filled with an unbounded ray
);

/**
 * @brief Returns an unbounded primal ray that is a certificate of dual
 * infeasibility.
 */
HighsInt Highs_getPrimalRay(
    const void* highs,
    HighsInt* has_primal_ray,  //!< TRUE if the primal ray exists
    double* primal_ray_value   //!< array of length [numcol], filled with an
                               //!< unbounded ray
);

/**
 * @brief Returns the objective function value (if known)
 */
double Highs_getObjectiveValue(const void* highs);

/**
 * @brief Gets the basic variables in the order corresponding to
 * calls to getBasisInverseRow, getBasisInverseCol, getBasisSolve,
 * getBasisTransposeSolve, getReducedRow and getReducedColumn. As
 * required by SCIP, non-negative entries are indices of columns,
 * and negative entries are -(row_index+1).
 */
HighsInt Highs_getBasicVariables(const void* highs,
                                 HighsInt* basic_variables  //!< Basic variables
);

/**
 * @brief Gets a row of \f$B^{-1}\f$ for basis matrix \f$B\f$
 */
HighsInt Highs_getBasisInverseRow(
    const void* highs,
    const HighsInt row,    //!< Index of row required
    double* row_vector,    //!< Row required
    HighsInt* row_num_nz,  //!< Number of nonzeros
    HighsInt* row_indices  //!< Indices of nonzeros
);

/**
 * @brief Gets a column of \f$B^{-1}\f$ for basis matrix \f$B\f$
 */
HighsInt Highs_getBasisInverseCol(
    const void* highs,
    const HighsInt col,    //!< Index of column required
    double* col_vector,    //!< Column required
    HighsInt* col_num_nz,  //!< Number of nonzeros
    HighsInt* col_indices  //!< Indices of nonzeros
);

/**
 * @brief Forms \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$
 */
HighsInt Highs_getBasisSolve(
    const void* highs,
    const double* rhs,          //!< RHS \f$\mathbf{b}\f$
    double* solution_vector,    //!< Solution \f$\mathbf{x}\f$
    HighsInt* solution_num_nz,  //!< Number of nonzeros
    HighsInt* solution_indices  //!< Indices of nonzeros
);

/**
 * @brief Forms \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$
 */
HighsInt Highs_getBasisTransposeSolve(
    const void* highs,
    const double* rhs,          //!< RHS \f$\mathbf{b}\f$
    double* solution_vector,    //!< Solution  \f$\mathbf{x}\f$
    HighsInt* solution_nz,      //!< Number of nonzeros
    HighsInt* solution_indices  //!< Indices of nonzeros
);

/**
 * @brief Forms a row of \f$B^{-1}A\f$
 */
HighsInt Highs_getReducedRow(const void* highs,
                             const HighsInt row,    //!< Index of row required
                             double* row_vector,    //!< Row required
                             HighsInt* row_num_nz,  //!< Number of nonzeros
                             HighsInt* row_indices  //!< Indices of nonzeros
);

/**
 * @brief Forms a column of \f$B^{-1}A\f$
 */
HighsInt Highs_getReducedColumn(
    const void* highs,
    const HighsInt col,    //!< Index of column required
    double* col_vector,    //!< Column required
    HighsInt* col_num_nz,  //!< Number of nonzeros
    HighsInt* col_indices  //!< Indices of nonzeros
);

/**
 * @brief Passes a basis to HiGHS
 */
HighsInt Highs_setBasis(void* highs,
                        const HighsInt* colstatus,  //!< Column status
                        const HighsInt* rowstatus   //!< Row status
);

/**
 * @brief Sets up a logical basis in HiGHS
 */
HighsInt Highs_setLogicalBasis(void* highs);

/**
 * @brief Returns the cumulative wall-clock time spent in Highs_run();
 */
double Highs_getRunTime(const void* highs);

/**
 * @brief Adds a row to the model
 */
HighsInt Highs_addRow(
    void* highs,
    const double lower,         //!< Lower bound of the row
    const double upper,         //!< Upper bound of the row
    const HighsInt num_new_nz,  //!< Number of nonzeros in the row
    const HighsInt* indices,  //!< Array of size num_new_nz with column indices
    const double* values      //!< Array of size num_new_nz with column values
);

/**
 * @brief Adds multiple rows to the model
 */
HighsInt Highs_addRows(
    void* highs,
    const HighsInt num_new_row,  //!< Number of new rows
    const double* lower,        //!< Array of size num_new_row with lower bounds
    const double* upper,        //!< Array of size num_new_row with upper bounds
    const HighsInt num_new_nz,  //!< Number of new nonzeros
    const HighsInt*
        starts,  //!< Array of size num_new_row with start indices of the rows
    const HighsInt*
        indices,  //!< Array of size num_new_nz with column indices for all rows
    const double*
        values  //!< Array of size num_new_nz with column values for all rows
);

/**
 * @brief Adds a column to the model
 */
HighsInt Highs_addCol(
    void* highs,
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
HighsInt Highs_addCols(
    void* highs,
    const HighsInt num_new_col,  //!< Number of new columns
    const double* costs,         //!< Array of size num_new_col with costs
    const double* lower,        //!< Array of size num_new_col with lower bounds
    const double* upper,        //!< Array of size num_new_col with upper bounds
    const HighsInt num_new_nz,  //!< Number of new nonzeros
    const HighsInt* starts,   //!< Array of size num_new_row with start indices
                              //!< of the columns
    const HighsInt* indices,  //!< Array of size num_new_nz with row indices for
                              //!< all columns
    const double*
        values  //!< Array of size num_new_nz with row values for all columns
);

/**
 * @brief Change the objective sense of the model
 */
HighsInt Highs_changeObjectiveSense(
    void* highs,
    const HighsInt sense  //!< New objective sense
);

/**
 * @brief Change the objective offset of the model
 */
HighsInt Highs_changeObjectiveOffset(
    void* highs,
    const double offset  //!< New objective offset
);

/**
 * @brief Change the integrality of a column
 */
HighsInt Highs_changeColIntegrality(
    void* highs,
    const HighsInt
        col,  //!< The index of the column whose integrality is to change
    const HighsInt integrality  //!< The new integrality
);

/**
 * @brief Change the integrality of multiple columns given by an interval
 */
HighsInt Highs_changeColsIntegralityByRange(
    void* highs,
    const HighsInt
        from_col,  //!< The index of the first column whose integrality changes
    const HighsInt to_col,  //!< One more than the index of the last column
                            //!< whose integrality changes
    const HighsInt*
        integrality  //!< Array of size num_set_entries with new integralitys
);

/**
 * @brief Change the integrality of multiple columns given by a set of indices
 */
HighsInt Highs_changeColsIntegralityBySet(
    void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set,  //!< Array of size num_set_entries with indices of
                          //!< columns whose integralitys change
    const HighsInt*
        integrality  //!< Array of size num_set_entries with new integralitys
);

/**
 * @brief Change the integrality of multiple columns given by a mask
 */
HighsInt Highs_changeColsIntegralityByMask(
    void* highs,
    const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
    const HighsInt* integrality  //!< Full length array of new integralitys
);

/**
 * @brief Change the cost of a column
 */
HighsInt Highs_changeColCost(
    void* highs,
    const HighsInt col,  //!< The index of the column whose cost is to change
    const double cost    //!< The new cost
);

/**
 * @brief Change the cost of multiple columns given by an interval
 */
HighsInt Highs_changeColsCostByRange(
    void* highs,
    const HighsInt
        from_col,  //!< The index of the first column whose cost changes
    const HighsInt to_col,  //!< One more than the index of the last column
                            //!< whose cost changes
    const double* cost      //!< Array of size num_set_entries with new costs
);

/**
 * @brief Change the cost of multiple columns given by a set of indices
 */
HighsInt Highs_changeColsCostBySet(
    void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set,  //!< Array of size num_set_entries with indices of
                          //!< columns whose costs change
    const double* cost    //!< Array of size num_set_entries with new costs
);

/**
 * @brief Change the cost of multiple columns given by a mask
 */
HighsInt Highs_changeColsCostByMask(
    void* highs,
    const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
    const double* cost     //!< Full length array of new costs
);

/**
 * @brief Change the bounds of a column
 */
HighsInt Highs_changeColBounds(
    void* highs,
    const HighsInt col,  //!< The index of the column whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
);

/**
 * @brief Change the bounds of multiple columns given by an interval
 */
HighsInt Highs_changeColsBoundsByRange(
    void* highs,
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
HighsInt Highs_changeColsBoundsBySet(
    void* highs,
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
HighsInt Highs_changeColsBoundsByMask(
    void* highs,
    const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
    const double* lower,   //!< Full length array of new lower bounds
    const double* upper    //!< Full length array of new upper bounds
);

/**
 * @brief Change the bounds of a row
 */
HighsInt Highs_changeRowBounds(
    void* highs,
    const HighsInt row,  //!< The index of the row whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
);

/**
 * @brief Change the bounds of multiple rows given by a set of indices
 */
HighsInt Highs_changeRowsBoundsBySet(
    void* highs,
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
HighsInt Highs_changeRowsBoundsByMask(
    void* highs,
    const HighsInt* mask,  //!< Full length array with 1 => change; 0 => not
    const double* lower,   //!< Full length array of new lower bounds
    const double* upper    //!< Full length array of new upper bounds
);

/**
 * @brief Change a coefficient in the constraint matrix.
 */
HighsInt Highs_changeCoeff(
    void* highs,
    const HighsInt row,  //!< The index of the row to change
    const HighsInt col,  //!< The index of the column to change
    const double value   //!< The new coefficient
);

/**
 * @brief Get the objective sense
 */
HighsInt Highs_getObjectiveSense(const void* highs, HighsInt* sense);

/**
 * @brief Get the objective offset
 */
HighsInt Highs_getObjectiveOffset(const void* highs, double* offset);

/**
 * @brief Get multiple columns from the model given by an interval
 */
HighsInt Highs_getColsByRange(
    const void* highs,
    const HighsInt from_col,  //!< The index of the first column to
                              //!< get from the model
    const HighsInt to_col,    //!< One more than the last column to get
                              //!< from the model
    HighsInt* num_col,        //!< Number of columns got from the model
    double* costs,            //!< Array of size num_col with costs
    double* lower,            //!< Array of size num_col with lower bounds
    double* upper,            //!< Array of size num_col with upper bounds
    HighsInt* num_nz,         //!< Number of nonzeros got from the model
    HighsInt* matrix_start,   //!< Array of size num_col with start
                              //!< indices of the columns
    HighsInt* matrix_index,   //!< Array of size num_nz with row
                              //!< indices for the columns
    double* matrix_value      //!< Array of size num_nz with row
                              //!< values for the columns
);

/**
 * @brief Get multiple columns from the model given by a set
 */
HighsInt Highs_getColsBySet(
    const void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set,     //!< Array of size num_set_entries with indices
                             //!< of columns to get
    HighsInt* num_col,       //!< Number of columns got from the model
    double* costs,           //!< Array of size num_col with costs
    double* lower,           //!< Array of size num_col with lower bounds
    double* upper,           //!< Array of size num_col with upper bounds
    HighsInt* num_nz,        //!< Number of nonzeros got from the model
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
HighsInt Highs_getColsByMask(
    const void* highs,
    const HighsInt* mask,    //!< Full length array with 1 => get; 0 => not
    HighsInt* num_col,       //!< Number of columns got from the model
    double* costs,           //!< Array of size num_col with costs
    double* lower,           //!< Array of size num_col with lower bounds
    double* upper,           //!< Array of size num_col with upper bounds
    HighsInt* num_nz,        //!< Number of nonzeros got from the model
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
HighsInt Highs_getRowsByRange(
    const void* highs,
    const HighsInt
        from_row,  //!< The index of the first row to get from the model
    const HighsInt to_row,   //!< One more than the last row get from the model
    HighsInt* num_row,       //!< Number of rows got from the model
    double* lower,           //!< Array of size num_row with lower bounds
    double* upper,           //!< Array of size num_row with upper bounds
    HighsInt* num_nz,        //!< Number of nonzeros got from the model
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
HighsInt Highs_getRowsBySet(
    const void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set,     //!< Array of size num_set_entries with indices
                             //!< of rows to get
    HighsInt* num_row,       //!< Number of rows got from the model
    double* lower,           //!< Array of size num_row with lower bounds
    double* upper,           //!< Array of size num_row with upper bounds
    HighsInt* num_nz,        //!< Number of nonzeros got from the model
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
HighsInt Highs_getRowsByMask(
    const void* highs,
    const HighsInt* mask,    //!< Full length array with 1 => get; 0 => not
    HighsInt* num_row,       //!< Number of rows got from the model
    double* lower,           //!< Array of size num_row with lower bounds
    double* upper,           //!< Array of size num_row with upper bounds
    HighsInt* num_nz,        //!< Number of nonzeros got from the model
    HighsInt* matrix_start,  //!< Array of size num_row with start indices
                             //!< of the rows
    HighsInt* matrix_index,  //!< Array of size num_nz with column indices
                             //!< for the rows
    double* matrix_value     //!< Array of size num_nz with column
                             //!< values for the rows
);

/**
 * @brief Delete multiple columns from the model given by an interval
 */
HighsInt Highs_deleteColsByRange(
    void* highs,
    const HighsInt from_col,  //!< The index of the first column
                              //!< to delete from the model
    const HighsInt to_col     //!< One more than the last column to
                              //!< delete from the model
);

/**
 * @brief Delete multiple columns from the model given by a set
 */
HighsInt Highs_deleteColsBySet(
    void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set  //!< Array of size num_set_entries with indices of
                         //!< columns to delete
);

/**
 * @brief Delete multiple columns from the model given by a mask
 */
HighsInt Highs_deleteColsByMask(
    void* highs,
    HighsInt* mask  //!< Full length array with 1 => delete; 0 => not
);

/**
 * @brief Delete multiple rows from the model given by an interval
 */
HighsInt Highs_deleteRowsByRange(
    void* highs,
    const int
        from_row,  //!< The index of the first row to delete from the model
    const HighsInt to_row  //!< One more than the last row delete from the model
);

/**
 * @brief Delete multiple rows from the model given by a set
 */
HighsInt Highs_deleteRowsBySet(
    void* highs,
    const HighsInt num_set_entries,  //!< The number of indides in the set
    const HighsInt* set  //!< Array of size num_set_entries with indices of
                         //!< columns to delete
);

/**
 * @brief Delete multiple rows from the model given by a mask
 */
HighsInt Highs_deleteRowsByMask(
    void* highs,
    HighsInt* mask  //!< Full length array with 1 => delete; 0 => not
);

/**
 * @brief Scale a matrix column (and cost) by a constant - flipping
 * bounds if the constant is negative
 */
HighsInt Highs_scaleCol(void* highs,
                        const HighsInt col,    //!< Column to scale
                        const double scaleval  //!< Value to scale by
);

/**
 * @brief Scale a matrix row by a constant - flipping bounds if the
 * constant is negative
 */
HighsInt Highs_scaleRow(void* highs,
                        const HighsInt row,    //!< Row to scale
                        const double scaleval  //!< Value to scale by
);

/**
 * @brief Returns the value of infinity used by HiGHS
 */
double Highs_getInfinity(const void* highs);

/**
 * @brief Returns the number of columns of the current model
 */
HighsInt Highs_getNumCol(const void* highs);

/**
 * @brief Returns the number of rows of the current model
 */
HighsInt Highs_getNumRow(const void* highs);

/**
 * @brief Returns the number of nonzeros of the current model
 */
HighsInt Highs_getNumNz(const void* highs);

/**
 * @brief Returns the number of nonzeroes of the current Hessian
 */
HighsInt Highs_getHessianNumNz(const void* highs);

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

HighsInt Highs_crossover(void* highs  //!< HiGHS object reference
);

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

HighsInt Highs_setHighsLogfile(
    void* highs,
    const void* logfile  //!< File handle of the logfile
);

HighsInt Highs_setHighsOutput(
    void* highs,
    const void* outputfile  //!< File handle of the output file
);

HighsInt Highs_getIterationCount(const void* highs);

HighsInt Highs_getSimplexIterationCount(const void* highs);

HighsInt Highs_setHighsBoolOptionValue(
    void* highs,
    const char* option,   //!< name of the option
    const HighsInt value  //!< new value of option
);

// /**
//  * @brief Returns the current model
//  */
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
HighsInt Highs_setHighsIntOptionValue(
    void* highs,
    const char* option,   //!< name of the option
    const HighsInt value  //!< new value of option
);

HighsInt Highs_setHighsDoubleOptionValue(
    void* highs,
    const char* option,  //!< name of the option
    const double value   //!< new value of option
);

HighsInt Highs_setHighsStringOptionValue(
    void* highs,
    const char* option,  //!< name of the option
    const char* value    //!< new value of option
);

HighsInt Highs_setHighsOptionValue(void* highs,
                                   const char* option,  //!< name of the option
                                   const char* value    //!< new value of option
);

HighsInt Highs_getHighsBoolOptionValue(
    const void* highs,
    const char* option,  //!< name of the option
    HighsInt* value      //!< value of option
);

HighsInt Highs_getHighsIntOptionValue(
    const void* highs,
    const char* option,  //!< name of the option
    HighsInt* value      //!< value of option
);

HighsInt Highs_getHighsDoubleOptionValue(
    const void* highs,
    const char* option,  //!< name of the option
    double* value        //!< value of option
);

HighsInt Highs_getHighsStringOptionValue(
    const void* highs,
    const char* option,  //!< name of the option
    char* value  //!< pointer to allocated memory to store value of option
);

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
