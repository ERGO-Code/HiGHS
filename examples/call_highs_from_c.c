#include "interfaces/highs_lp_solver.h"

#include <stdlib.h>
#include <stdio.h>

int main() {
  int numcol = 2;
  int numrow = 2;
  int nnz = 4;
  int i;

  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  int astart[3] = {0, 2, 4};
  int aindex[4] = {0, 1, 0, 1};
  double avalue[4] = {1.0, 2.0, 1.0, 3.0};

  double* cv = (double*)malloc(sizeof(double) * numcol);
  double* cd = (double*)malloc(sizeof(double) * numcol);
  double* rv = (double*)malloc(sizeof(double) * numrow);
  double* rd = (double*)malloc(sizeof(double) * numrow);

  int* cbs = (int*)malloc(sizeof(int) * numcol);
  int* rbs = (int*)malloc(sizeof(int) * numrow);

  callhighs(numcol, numrow, nnz, cc, cl, cu, rl, ru, astart, aindex, avalue, cv,
            cd, rv, rd, cbs, rbs);

  for (i=0; i<numcol; i++) {
     printf("x%d = %lf\n", i, cv[i]);
  }

  free(cv);
  free(cd);
  free(rv);
  free(rd);
  free(cbs);
  free(rbs);

  return 0;
}