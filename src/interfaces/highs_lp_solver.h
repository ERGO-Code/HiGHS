#ifndef HIGHS_LP_SOLVER
#define HIGHS_LP_SOLVER

#ifdef __cplusplus
extern "C" {
#endif
void callhighs(int numcol, int numrow, int numnz, double* colcost,
                          double* collower, double* colupper, double* rowlower,
                          double* rowupper, int* astart, int* aindex,
                          double* avalue, double* col_value, double* col_dual,
                          double* row_value, double* row_dual,
                          int* col_basis_status, int* row_basis_status);

#ifdef __cplusplus
}
#endif

#endif
