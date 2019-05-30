#include "interfaces/highs_lp_solver.h"
#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>

// gcc call_highs_from_c.c -o highstest -I ../build/install_folder/include/ -L ../build/install_folder/lib/ -lhighs

void minimal_api() {
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

  for (i = 0; i < numcol; i++) {
    printf("x%d = %lf\n", i, cv[i]);
  }

  free(cv);
  free(cd);
  free(rv);
  free(rd);
  free(cbs);
  free(rbs);
}

void full_api() {
  void* highs;

  highs = Highs_create();
  Highs_loadFromFile(highs, "/home/s1613957/lpInstances/qap04.mps");
  Highs_run(highs);
  Highs_destroy(highs);
}

int main() { 
  minimal_api();
  full_api();
  return 0; 
}