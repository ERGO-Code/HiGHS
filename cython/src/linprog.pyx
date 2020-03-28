# distutils: language=c++
# cython: language_level=3

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, allocator

from cython.operator cimport dereference

cimport numpy as np
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import OptimizeResult

cdef extern from "highs_c_api.h" nogil:
    int Highs_call(
        int numcol,
        int numrow,
        int numnz,
        double* colcost,
        double* collower,
        double* colupper,
        double* rowlower,
        double* rowupper,
        int* astart,
        int* aindex,
        double* avalue,
        double* colvalue,
        double* coldual,
        double* rowvalue,
        double* rowdual,
        int* colbasisstatus,
        int* rowbasisstatus,
        int* modelstatus)

def linprog(double[::1] c, A, double[::1] b, double[::1] lhs=None):
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
    cdef double * colcost = &c[0]

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

    cdef int ret = Highs_call(
        numcol, numrow, numnz,
        colcost, collower, colupper,
        &rowlower[0], rowupper,
        &astart[0], &aindex[0], &avalue[0],
        &colvalue[0], &coldual[0], &rowvalue[0], &rowdual[0],
        &colbasisstatus[0], &rowbasisstatus[0], &modelstatus)

    return OptimizeResult({
        'fun': np.sum(c*np.array(colvalue)), # There's a way to get this, just haven't found it yet
        'x': np.array(colvalue),
        'slack': np.array(coldual),
        'rowvalue': np.array(rowvalue),
        'rowdual': np.array(rowdual),
        'colbasisstatus': np.array(colbasisstatus),
        'rowbasisstatus': np.array(rowbasisstatus),
        'modelstatus': modelstatus,
    })
