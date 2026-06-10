#ifndef FACTOR_HIGHS_C_API_H
#define FACTOR_HIGHS_C_API_H

#include "util/HighsInt.h"

/* C API to HiPO linear solver

Refer to FactorHighs.h for a description of the interface.

    void* S = FactorHighs_symbolic_create();
    void* FH = FactorHighs_create();
    FactorHighs_analyse(FH, S, n, nz, rows, ptr, signs, perm);
    FactorHighs_factorise(FH, S, n, nz, rows, ptr, val);
    FactorHighs_solve(FH, x);
    FactorHighs_symbolic_destroy(S);
    FactorHighs_destroy(FH);

Printing to screen is not available through the C API for now.
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
void FactorHighs_setPivoting(void* FH, HighsInt pivoting);
void FactorHighs_inertia(void* FH, HighsInt& pos, HighsInt& neg,
                         HighsInt& zero, double tol);

#ifdef __cplusplus
}
#endif

#endif