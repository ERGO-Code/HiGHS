#ifndef FACTOR_HIGHS_C_API_H
#define FACTOR_HIGHS_C_API_H

#include "util/HighsInt.h"

/* C API to HiPO linear solver

Consider a sparse symmetric matrix M in CSC format.
Only its lower triangular part is used; entries in the upper triangle are
ignored.
The matrix has n rows/cols and nz nonzero entries.
It is stored using three arrays:
- ptr, column pointers, of length n+1;
- rows, row indices, of length nz;
- vals, values, of length nz.


The direct solver uses the following objects:
- Symbolic, to store the symbolic factorization;
- FHsolver, to perform analyse and factorise phases.
They are both accessed through void* pointers.

Define a vector signs that contains the expected sign of each pivot (1 or -1).
Define a right-hand side rhs, which will be overwritten with the solution of
M^{-1} * rhs. The pre-computed fill-reducing ordering to use is stored in the
vector perm.

Then, the factorization is performed as follows.

    void* S = FactorHighs_symbolic_create();
    void* FH = FactorHighs_create();
    FactorHighs_analyse(FH, S, n, nz, rows, ptr, signs, perm);
    FactorHighs_factorise(FH, S, n, nz, rows, ptr, val);
    FactorHighs_solve(FH, x);
    FactorHighs_symbolic_destroy(S);
    FactorHighs_destroy(FH);

Printing to screen is not available through the C API for now.

To add static regularisation when the pivots are selected, use
FactorHighs_setRegularisation(reg_p, reg_d) to choose values of primal and dual
regularisation. If regularisation is already added to the matrix, ignore.

The default block size is 128. Changing this options is not available through
the C API for now.

*/

#ifdef __cplusplus
extern "C" {
#endif

void* FactorHighs_symbolic_create(void);
void FactorHighs_symbolic_destroy(void* S);
void* FactorHighs_create(void);
void FactorHighs_destroy(void* FH);
HighsInt FactorHighs_analyse(void* FH, void* S, HighsInt n, HighsInt nz,
                             const HighsInt* rows, const HighsInt* ptr,
                             const HighsInt* signs, const HighsInt* perm);
HighsInt FactorHighs_factorise(void* FH, const void* S, HighsInt n, HighsInt nz,
                               const HighsInt* rows, const HighsInt* ptr,
                               const double* vals);
HighsInt FactorHighs_solve(void* FH, double* x);
void FactorHighs_setRegularisation(void* FH, double reg_p, double reg_d);
void FactorHighs_getRegularisation(void* FH, double* reg);
void FactorHighs_newIter(void* FH);
void FactorHighs_setBlockSize(void* FH, HighsInt nb);

#ifdef __cplusplus
}
#endif

#endif