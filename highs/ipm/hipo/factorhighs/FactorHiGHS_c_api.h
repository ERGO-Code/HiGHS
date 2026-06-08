#ifndef FACTOR_HIGHS_C_API_H
#define FACTOR_HIGHS_C_API_H

#include "util/HighsInt.h"

// C API to HiPO linear solver

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif