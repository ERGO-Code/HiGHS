#ifndef HCONST_H_
#define HCONST_H_

const int HIGHS_CONST_I_INF = 32767;
const double HIGHS_CONST_INF = 1e200;
const double HIGHS_CONST_TINY = 1e-14;
const double HIGHS_CONST_ZERO = 1e-50;

// Until everyone is weaned off HSOL_CONST*
const int HSOL_CONST_I_INF = HIGHS_CONST_I_INF;
const double HSOL_CONST_INF = HIGHS_CONST_INF;
const double HSOL_CONST_TINY = HIGHS_CONST_TINY;
const double HSOL_CONST_ZERO = HIGHS_CONST_ZERO;

#endif /* HCONST_H_ */
