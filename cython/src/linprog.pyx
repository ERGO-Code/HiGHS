# distutils: language=c++
# cython: language_level=3

from libcpp.memory cimport unique_ptr, allocator

from cython.operator cimport dereference

cimport numpy as np
import numpy as np

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

def linprog(double[::1] c, A, double[::1] b):
    '''Solve linear programs using dense matrices.'''

    cdef int numrow = A.shape[0]
    cdef int numcol = A.shape[1]
    cdef int numnz = A.nnz

    # Allocate pointers to hold results;
    # smart pointers will manage their own memory.
    cdef allocator[double] double_al
    cdef unique_ptr[double] colvalue
    cdef unique_ptr[double] coldual
    cdef unique_ptr[double] rowvalue
    cdef unique_ptr[double] rowdual
    colvalue.reset(double_al.allocate(numcol))
    coldual.reset(double_al.allocate(numcol))
    rowvalue.reset(double_al.allocate(numrow))
    rowdual.reset(double_al.allocate(numrow))

    # Objective function coefficients
    cdef double * colcost = &c[0]

    # Bounds on variables
    cdef double[::1] collower_memview = np.zeros(numcol, dtype='double')
    cdef double[::1] colupper_memview = 1e20*np.ones(numcol, dtype='double')
    cdef double * collower = &collower_memview[0]
    cdef double * colupper = &colupper_memview[0]

    # LHS/RHS constraints
    cdef double[::1] rowlower_memview = -1e20*np.ones(numrow, dtype='double')
    cdef double * rowlower = &rowlower_memview[0]
    cdef double * rowupper = &b[0]

    # Contents of constraint matrices
    cdef int[::1] astart_memview = A.indptr
    cdef int[::1] aindex_memview = A.indices
    cdef double[::1] avalue_memview = A.data
    cdef int * astart = &astart_memview[0]
    cdef int * aindex = &aindex_memview[0]
    cdef double * avalue = &avalue_memview[0]

    # Result status flags
    cdef allocator[int] int_al
    cdef unique_ptr[int] colbasisstatus
    cdef unique_ptr[int] rowbasisstatus
    colbasisstatus.reset(int_al.allocate(numcol))
    rowbasisstatus.reset(int_al.allocate(numrow))
    cdef int modelstatus = 0

    cdef int ret = Highs_call(
        numcol, numrow, numnz,
        colcost, collower, colupper,
        rowlower, rowupper,
        astart, aindex, avalue,
        colvalue.get(), coldual.get(), rowvalue.get(), rowdual.get(),
        colbasisstatus.get(), rowbasisstatus.get(), &modelstatus)

    print([colvalue.get()[ii] for ii in range(numcol)])
    print([coldual.get()[ii] for ii in range(numcol)])
    print([rowvalue.get()[ii] for ii in range(numrow)])
    print([rowdual.get()[ii] for ii in range(numrow)])
