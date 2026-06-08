#include "FactorHighs_c_api.h"

#include "FactorHighs.h"

void* FactorHighs_create(void) { return new hipo::FHsolver(); }

void FactorHighs_destroy(void* FH) { delete (hipo::FHsolver*)FH; }

HighsInt FactorHighs_analyse(void* FH, void* S, HighsInt n, HighsInt nz,
                             const HighsInt* rows, const HighsInt* ptr,
                             const HighsInt* signs, const HighsInt* perm) {
  return ((hipo::FHsolver*)FH)
      ->analyse(*(hipo::Symbolic*)S, n, nz, rows, ptr, signs, perm);
}

HighsInt FactorHighs_factorise(void* FH, const void* S, HighsInt n, HighsInt nz,
                               const HighsInt* rows, const HighsInt* ptr,
                               const double* vals) {
  return ((hipo::FHsolver*)FH)
      ->factorise(*(hipo::Symbolic*)S, n, nz, rows, ptr, vals);
}

HighsInt FactorHighs_solve(void* FH, double* x) {
  return ((hipo::FHsolver*)FH)->solve(x);
}

void FactorHighs_setRegularisation(void* FH, double reg_p, double reg_d) {
  ((hipo::FHsolver*)FH)->setRegularisation(reg_p, reg_d);
}

void FactorHighs_getRegularisation(void* FH, double* reg) {
  ((hipo::FHsolver*)FH)->getRegularisation(reg);
}

void FactorHighs_newIter(void* FH) { ((hipo::FHsolver*)FH)->newIter(); }
