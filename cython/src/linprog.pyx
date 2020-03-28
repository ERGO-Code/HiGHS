# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE, tmpfile

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, make_unique

cimport numpy as np
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import OptimizeResult

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
                    int* modelstatus, int sense, Highs & highs):
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


def linprog(double[::1] c, A, double[::1] b, double[::1] lhs=None, int sense=1, bool disp=False):
    '''Solve linear programs.

    Assume the form:

        min c.T @ x
        s.t. lhs <= A @ x  <= b

    Still working on bounds on x, currently assumed to be (0, inf).

    Parameters
    ----------
    c : 1-D array, (n,)
        Array of objective value coefficients.
    A : 2-D array, (m, n)
        Sparse (or dense) matrix of constraint coefficients.
    b : 1-D array, (m,)
        Array of right hand side values of the inequality constraints.
    lhs : 1-D array (or None), (m,)
        Array of left hand side values of the inequality constraints.
        If `lhs=None`, then an array of `-inf` is assumed.
    sense : int {1, -1}, optional
        `sense=1` corresponds to the MIN problem, `sense=-1`
        corresponds to the MAX problem.
    disp : bool, optional
        Write information about solver to stdio.

    Returns
    -------
    res : OptimizeResult

        - x : 1-D array, (n,)
        - fun : double
        - slack : 1-D array, (n,)
        - rowvalue : 1-D array, (m,)
        - rowdual : 1-D array, (m,)
        - colbasisstatus : 1-D array, (n,)
        - rowbasisstatus : 1-D array, (m,)
        - modelstatus : int
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
    if sense == -1:
        for ii in range(numcol):
            cc[ii] *= -1
    cdef double * colcost = &cc[0]

    # Bounds on variables
    cdef double[::1] collower_memview = np.zeros(numcol, dtype='double')
    cdef double[::1] colupper_memview = 1e20*np.ones(numcol, dtype='double')
    cdef double * collower = &collower_memview[0]
    cdef double * colupper = &colupper_memview[0]

    # LHS/RHS constraints
    cdef double * rowlower
    if lhs is None:
        # Default to no LHS (all -Inf)
        lhs = -1e20*np.ones(numrow, dtype='double')
    rowlower = &lhs[0]
    cdef double * rowupper = &b[0]

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
    cdef FILE * f
    if not disp:
        # set output to tmp to remove stdio spam
        f = tmpfile()
        highs.setHighsOutput(f)

    # Call the solver
    cdef int ret = Highs_call(
        numcol, numrow, numnz,
        colcost, collower, colupper,
        &rowlower[0], rowupper,
        &astart[0], &aindex[0], &avalue[0],
        &colvalue[0], &coldual[0], &rowvalue[0], &rowdual[0],
        &colbasisstatus[0], &rowbasisstatus[0], &modelstatus,
        sense, highs)

    # Pull info out of out of highs
    cdef HighsInfo info = highs.getHighsInfo()
    return OptimizeResult({
        # From HighsInfo
        'fun': info.objective_function_value,
        'simplex_nit': info.simplex_iteration_count,
        'ipm_nit': info.ipm_iteration_count,
        'crossover_nit': info.crossover_iteration_count,
        'primal_status': info.primal_status,
        'dual_status': info.dual_status,
        'num_primal_infeasibilities': info.num_primal_infeasibilities,
        'max_primal_infeasibility': info.max_primal_infeasibility,
        'sum_primal_infeasibilities': info.sum_primal_infeasibilities,
        'num_dual_infeasibilities': info.num_dual_infeasibilities,
        'max_dual_infeasibility': info.max_dual_infeasibility,
        'sum_dual_infeasibilities': info.sum_dual_infeasibilities,

        # From C API
        'x': np.array(colvalue),
        'slack': np.array(coldual),
        'row_value': np.array(rowvalue),
        'row_dual': np.array(rowdual),
        'col_basis_status': np.array(colbasisstatus),
        'row_basis_status': np.array(rowbasisstatus),
        'model_status': highs.highsModelStatusToString(<HighsModelStatus>modelstatus).decode(),
    })
