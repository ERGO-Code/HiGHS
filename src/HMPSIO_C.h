#ifndef HMPSIO_C_H_
#define HMPSIO_C_H_

int readMPS_LP_dense_c(const char *filename, int mxNumRow, int mxNumCol,
                       int *numRow_p, int *numCol_p, int *objSense_p, double *objOffset_p,
                       double **A_rw_p,
                       double **rhs_p, double **cost_p, double **lb_p, double **ub_p);

#endif
