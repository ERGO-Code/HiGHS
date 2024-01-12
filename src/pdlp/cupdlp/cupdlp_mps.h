#ifndef lp_mps_h
#define lp_mps_h

#include "cupdlp_defs.h"

/* Implement an LP mps file reader */
cupdlp_int cupdlpMpsRead(char *fname, char *name, int *pnRow, int *pnEqRow,
                         int *pnInEqRow, int *pnCol, int *pnElem,
                         int **pfullMatBeg, int **pfullMatIdx,
                         double **pfullMatElem, int **peqMatBeg,
                         int **peqMatIdx, double **peqMatElem,
                         int **pIneqMatBeg, int **pIneqMatIdx,
                         double **pIneqMatElem, double **prowRHS,
                         double **pcolObj, int *pnColUb, int **pcolUbIdx,
                         double **pcolUbElem);

#endif /* lp_mps_h */
