# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE, tmpfile

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, make_unique

cimport numpy as np
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import OptimizeResult
from warnings import warn

from HConst cimport (
    # Verbosity levels
    ML_DETAILED,
    ML_NONE,
    ML_VERBOSE,
    ML_MINIMAL,

    # HighsBasisStatus
    HighsBasisStatus,
    LOWER,
    BASIC,
    UPPER,
    ZERO,
    NONBASIC,
    SUPER)
from Highs cimport Highs
from HighsLp cimport (
    HighsSolution,
    HighsBasis,
    HighsModelStatus)
from HighsInfo cimport HighsInfo
from highs_c_api cimport Highs_passLp

cdef int Highs_call(int numcol, int numrow, int numnz, double* colcost,
                    double* collower, double* colupper, double* rowlower,
                    double* rowupper, int* astart, int* aindex, double* avalue,
                    double* colvalue, double* coldual, double* rowvalue,
                    double* rowdual, int* colbasisstatus, int* rowbasisstatus,
                    int* modelstatus, Highs & highs):
    # cdef Highs highs
    cdef int status = Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                                   rowlower, rowupper, astart, aindex, avalue)

    # Customize sense : MIN or MAX
    # This API is not currently working, do it manually in caller
    # highs.changeObjectiveSense(sense)

    if (status != 0):
        return status
    status = <int>highs.run()

    cdef unique_ptr[HighsSolution] solution
    cdef HighsBasis basis
    if (status == 0):
        solution = make_unique[HighsSolution](highs.getSolution())
        basis = highs.getBasis()
        modelstatus[0] = <int>highs.getModelStatus()

        for ii in range(numcol):
            colvalue[ii] = solution.get().col_value[ii]
            coldual[ii] = solution.get().col_dual[ii]
            colbasisstatus[ii] = <int>basis.col_status[ii]

        for ii in range(numrow):
            rowvalue[ii] = solution.get().row_value[ii]
            rowdual[ii] = solution.get().row_dual[ii]
            rowbasisstatus[ii] = <int>basis.row_status[ii]

    return status

cdef apply_options(dict options, Highs & highs):
    '''Take options from dictionary and apply to HiGHS object.'''

    # Make sure we have a good verbosity level
    disp = options.get('disp', ML_NONE)
    if disp not in [ML_NONE, ML_VERBOSE, ML_DETAILED, ML_MINIMAL]:
        warn('disp level not one of {0, 1, 2, 4}! Choosing disp=0.')
        disp = ML_NONE

    # Set verbosity level for logging
    highs.setHighsOptionValue(b'message_level', disp)

    # Send logging to dummy file to get rid of output from stdout
    cdef FILE * f
    if disp == ML_NONE:
        f = tmpfile()
        highs.setHighsLogfile(f)

    # Set the presolve option
    presolve = options.get('presolve', None)
    if presolve is None:
        highs.setHighsOptionValueStr(b'presolve', b'choose')
    elif not presolve:
        highs.setHighsOptionValueStr(b'presolve', b'off')

    # Set the solver to use
    solver = options.get('solver', None)
    if solver is None:
        highs.setHighsOptionValueStr(b'solver', b'choose')
    elif solver == 'simplex':
        highs.setHighsOptionValueStr(b'solver', b'simplex')
    elif solver == 'ipm':
        highs.setHighsOptionValueStr(b'solver', b'ipm')
    else:
        warn('%s is not a recognized solver. Using default.' % solver)

    # Choose parallel or serial
    parallel = options.get('parallel', None)
    if parallel is None:
        highs.setHighsOptionValueStr(b'parallel', b'choose')
    elif parallel:
        highs.setHighsOptionValueStr(b'parallel', b'on')
    else:
        highs.setHighsOptionValueStr(b'parallel', b'off')

    # Set a time limit
    time_limit = options.get('time_limit', None)
    if time_limit is not None:
        highs.setHighsOptionValueDbl(b'time_limit', time_limit)

def highs_wrapper(
        double[::1] c,
        A,
        double[::1] rhs,
        double[::1] lhs=None,
        double[::1] lb=None,
        double[::1] ub=None,
        dict options=None):
    '''Solve linear programs using HiGHS.

    Assume problems of the form:

        MIN/MAX c.T @ x
        s.t. lhs <= A @ x <= rhs
             lb <= x <= ub

    Default is MIN (for MAX set `sense=-1`).

    Parameters
    ----------
    c : 1-D array, (n,)
        Array of objective value coefficients.
    A : 2-D array, (m, n)
        Sparse (or dense) matrix of constraint coefficients.
    rhs : 1-D array, (m,)
        Array of right hand side values of the inequality constraints.
    lhs : 1-D array (or None), (m,)
        Array of left hand side values of the inequality constraints.
        If `lhs=None`, then an array of `-inf` is assumed.
    lb : 1-D array (or None), (n,)
        Lower bounds on solution variables x.  If `lb=None`, then an
        array of all `0` is assumed.
    ub : 1-D array (or None), (n,)
        Upper bounds on solution variables x.  If `ub=None`, then an
        array of `inf` is assumed.
    options : dict
        A dictionary of solver options with the following fields:

            - disp : int {0, 1, 2, 4}
                Verbosity level, corresponds to:

                    - `0`: ML_NONE
                    - `1`: ML_VERBOSE
                    - `2`: ML_DETAILED
                    - `4`: ML_MINIMAL

            - dual_feasibility_tolerance : double
                Dual feasibility tolerance
            - dual_objective_value_upper_bound : double
                Upper bound on objective value for dual simplex:
                algorithm terminates if reached
            - infinite_bound : double
                Limit on abs(constraint bound): values larger than
                this will be treated as infinite
            - infinite_cost : double
                Limit on cost coefficient: values larger than this
                will be treated as infinite.
            - large_matrix_value : double
                Upper limit on abs(matrix entries): values larger than
                this will be treated as infinite
            - max_threads : int
                Maximum number of threads in parallel execution.
            - min_threads : int
                Minimum number of threads in parallel execution.
            - parallel : bool
                Run the solver in serial (False) or parallel (True).
            - presolve : bool
                Run the presolve or not (or if `None`, then choose).
            - primal_feasibility_tolerance : double
                Primal feasibility tolerance.
            - sense : int {1, -1}
                `sense=1` corresponds to the MIN problem, `sense=-1`
                corresponds to the MAX problem.
            - simplex_crash_strategy : int {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
                Strategy for simplex crash: off / LTSSF / Bixby (0/1/2).
                Default is `0`.  Corresponds to the following:

                    - `0`: `SIMPLEX_CRASH_STRATEGY_OFF`
                    - `1`: `SIMPLEX_CRASH_STRATEGY_LTSSF_K`
                    - `2`: `SIMPLEX_CRASH_STRATEGY_BIXBY`
                    - `3`: `SIMPLEX_CRASH_STRATEGY_LTSSF_PRI`
                    - `4`: `SIMPLEX_CRASH_STRATEGY_LTSF_K`
                    - `5`: `SIMPLEX_CRASH_STRATEGY_LTSF_PRI`
                    - `6`: `SIMPLEX_CRASH_STRATEGY_LTSF`
                    - `7`: `SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS`
                    - `8`: `SIMPLEX_CRASH_STRATEGY_BASIC`
                    - `9`: `SIMPLE_CRASH_STRATEGY_TEST_SING`

            - simplex_dual_edge_weight_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex dual edge weights:
                Dantzig / Devex / Steepest Edge. Corresponds
                to the following:

                    - `0`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX`
                    - `2`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH`
                    - `3`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE`
                    - `4`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL`

            - simplex_iteration_limit : int
                Iteration limit for simplex solver.

            - simplex_primal_edge_weight_strategy : int {0, 1}
                Strategy for simplex primal edge weights:
                Dantzig / Devex.  Corresponds to the following:

                    - `0`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX`

            - simplex_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex solver. Default: 1. Corresponds
                to the following:

                    - `0`: `SIMPLEX_STRATEGY_MIN`
                    - `1`: `SIMPLEX_STRATEGY_DUAL`
                    - `2`: `SIMPLEX_STRATEGY_DUAL_TASKS`
                    - `3`: `SIMPLEX_STRATEGY_DUAL_MULTI`
                    - `4`: `SIMPLEX_STRATEGY_PRIMAL`

            - simplex_update_limit : int
                Limit on the number of simplex UPDATE operations.
            - small_matrix_value : double
                Lower limit on abs(matrix entries): values smaller
                than this will be treated as zero.
            - solver : str {'simplex', 'ipm'}
                Choose which solver to use.  If `solver='simplex'`
                and `parallel=True` then PAMI will be used.
            - time_limit : double
                Max number of seconds to run the solver for.

    Returns
    -------
    res : OptimizeResult

        - col_basis_status : dict
            Key: `'statuses'` contains `n` status codes corresponding
                 to the `n` columns.
            Key: `'messages'` contains the `n` messages corresponding
                 to each status.
        - col_dual : 1-D array, (n,)
            The dual solution.
        - col_value : 1-D array, (n,)
            Solution variables.
        - crossover_nit : int
            Number of iterations taken to transform the interior
            solution produced by barrier into a basic solution
        - dual_status : dict
            Key: `'status'` contains the dual solution status code.
            Key: `'message'` contains the corresponding message.
        - fun : double
            The final objective value.
        - ipm_nit : int
            Number of iterations taken by IPM (interior-point solver).
        - max_dual_infeasibility : double
        - max_primal_infeasibility : double
        - model_status : dict
            Key: `'status'` contains the status code of the LP model.
            Key: `'message'` contains the corresponding message.
        - num_dual_infeasibilities : int
        - num_primal_infeasibilities : int
        - primal_status : dict
            Key: `'status'` contains the primal solution status code.
            Key: `'message'` contains the corresponding message.
        - row_basis_status : dict
            Key: `'statuses'` contains `m` status codes corresponding
                 to the `m` rows.
            Key: `'messages'` contains the `m` messages corresponding
                 to each status.
        - row_dual : 1-D array, (m,)
        - simplex_nit : int
            Number of iterations taken by the simplex solver.
        - sum_dual_infeasibilities : double
        - sum_primal_infeasibilities : double
    '''

    # Try to cast, it'll raise a type error if it don't work
    if not isinstance(A, csc_matrix):
        A = csc_matrix(A)

    # Get dimensions of problem
    cdef int numrow = A.shape[0]
    cdef int numcol = A.shape[1]
    cdef int numnz = A.nnz

    # Objective function coefficients
    # Do MIN/MAX conversion here because API not working for HiGHS
    cdef double[::1] cc = c.copy()
    if options.get('sense', 1) == -1:
        for ii in range(numcol):
            cc[ii] *= -1
    cdef double * colcost = &cc[0]

    # Bounds on variables
    cdef double * collower
    cdef double * colupper
    if lb is None:
        # Default is lower bound of 0
        lb = np.zeros(numcol, dtype='double')
    if ub is None:
        # Default is upper bound of inf
        ub = 1e20*np.ones(numcol, dtype='double')
    collower = &lb[0]
    colupper = &ub[0]

    # LHS/RHS constraints
    cdef double * rowlower
    if lhs is None:
        # Default to no LHS (all -Inf)
        lhs = -1e20*np.ones(numrow, dtype='double')
    rowlower = &lhs[0]
    cdef double * rowupper = &rhs[0]

    # Contents of constraint matrices as memoryviews
    cdef int[::1] astart = A.indptr
    cdef int[::1] aindex = A.indices
    cdef double[::1] avalue = A.data

    # Allocate memoryviews to hold results
    cdef double[::1] colvalue = np.empty(numcol, dtype='double')
    cdef double[::1] coldual = np.empty(numcol, dtype='double')
    cdef double[::1] rowvalue = np.empty(numcol, dtype='double')
    cdef double[::1] rowdual = np.empty(numcol, dtype='double')

    # Result status flags
    cdef int[::1] colbasisstatus = np.empty(numcol, dtype=np.int32)
    cdef int[::1] rowbasisstatus = np.empty(numrow, dtype=np.int32)
    cdef int modelstatus = 0

    cdef Highs highs

    # Apply any options
    apply_options(options, highs)

    # Call the solver
    cdef int ret = Highs_call(
        numcol, numrow, numnz,
        colcost, collower, colupper,
        &rowlower[0], rowupper,
        &astart[0], &aindex[0], &avalue[0],
        &colvalue[0], &coldual[0], &rowvalue[0], &rowdual[0],
        &colbasisstatus[0], &rowbasisstatus[0], &modelstatus,
        highs)

    # Decode HighsBasisStatus:
    HighsBasisStatusToStr = {
        <int>LOWER: 'LOWER: (slack) variable is at its lower bound [including fixed variables]',
        <int>BASIC: 'BASIC: (slack) variable is basic',
        <int>UPPER: 'UPPER: (slack) variable is at its upper bound',
        <int>ZERO: 'ZERO: free variable is non-basic and set to zero',
        <int>NONBASIC: 'NONBASIC: nonbasic with no specific bound information - useful for users and postsolve',
        <int>SUPER: 'SUPER: Super-basic variable: non-basic and either free and nonzero or not at a bound. No SCIP equivalent',
    }

    # Pull info out of out of highs
    cdef HighsInfo info = highs.getHighsInfo()
    return OptimizeResult({
        # From HighsInfo
        'fun': info.objective_function_value,
        'simplex_nit': info.simplex_iteration_count,
        'ipm_nit': info.ipm_iteration_count,
        'crossover_nit': info.crossover_iteration_count,
        'primal_status': {
            'status': info.primal_status,
            'message': highs.highsPrimalDualStatusToString(info.primal_status).decode(),
        },
        'dual_status': {
            'status': info.dual_status,
            'message': highs.highsPrimalDualStatusToString(info.dual_status).decode(),
        },
        'num_primal_infeasibilities': info.num_primal_infeasibilities,
        'max_primal_infeasibility': info.max_primal_infeasibility,
        'sum_primal_infeasibilities': info.sum_primal_infeasibilities,
        'num_dual_infeasibilities': info.num_dual_infeasibilities,
        'max_dual_infeasibility': info.max_dual_infeasibility,
        'sum_dual_infeasibilities': info.sum_dual_infeasibilities,

        # From C API
        'col_value': np.array(colvalue),
        'col_dual': np.array(coldual),
        'row_value': np.array(rowvalue),
        'row_dual': np.array(rowdual),
        'col_basis_status': {
            'statuses': [colbasisstatus[ii] for ii in range(numcol)],
            'messages': [HighsBasisStatusToStr[colbasisstatus[ii]] for ii in range(numcol)],
        },
        'row_basis_status': {
            'statuses': [rowbasisstatus[ii] for ii in range(numrow)],
            'messages': [HighsBasisStatusToStr[rowbasisstatus[ii]] for ii in range(numrow)],
        },
        'model_status': {
            'status': modelstatus,
            'message': highs.highsModelStatusToString(<HighsModelStatus>modelstatus).decode(),
        },
    })
