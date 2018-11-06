#ifndef SIMPLEX_HAPI_H_
#define SIMPLEX_HAPI_H_

const int HiGHS_probStatusUnset = -1;
const int HiGHS_probStatusOptimal = 0;
const int HiGHS_probStatusInfeasible = 1;
const int HiGHS_probStatusUnbounded = 2;
const int HiGHS_probStatusSingular = 3;
const int HiGHS_probStatusFailed = 4;
const int HiGHS_probStatusObjUB = 5;
const int HiGHS_probStatusOutOfTime = 6;

const int HiGHS_basisStatus_no = 0;
const int HiGHS_basisStatus_yes = 1;

void solve_fromArrays_dense(int *probStatus, int *basisStatus,
                            const int XnumCol, const int XnumRow,
                            const int XobjSense, const int XobjOffset,
                            const double *XcolCost, const double *XcolLower,
                            const double *XcolUpper, const double *XrowLower,
                            const double *XrowUpper, const double *XAmatrix,
                            double *colPrimalValues, double *colDualValues,
                            double *rowPrimalValues, double *rowDualValues,
                            int *basicVariables);

void solve_fromArrays(int *probStatus, int *basisStatus, const int XnumCol,
                      const int XnumRow, const int XnumNz, const int XobjSense,
                      const int XobjOffset, const double *XcolCost,
                      const double *XcolLower, const double *XcolUpper,
                      const double *XrowLower, const double *XrowUpper,
                      const int *XAstart, const int *XAindex,
                      const double *XAvalue, double *colPrimalValues,
                      double *colDualValues, double *rowPrimalValues,
                      double *rowDualValues, int *basicVariables);
#endif
