#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
// Force asserts to be checked always.
#undef NDEBUG
#include <assert.h>

void minimal_api() {
  HighsInt numcol = 2;
  HighsInt numrow = 2;
  HighsInt nnz = 4;
  HighsInt i;

  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  HighsInt astart[3] = {0, 2, 4};
  HighsInt aindex[4] = {0, 1, 0, 1};
  double avalue[4] = {1.0, 2.0, 1.0, 3.0};

  double* cv = (double*)malloc(sizeof(double) * numcol);
  double* cd = (double*)malloc(sizeof(double) * numcol);
  double* rv = (double*)malloc(sizeof(double) * numrow);
  double* rd = (double*)malloc(sizeof(double) * numrow);

  HighsInt* cbs = (HighsInt*)malloc(sizeof(HighsInt) * numcol);
  HighsInt* rbs = (HighsInt*)malloc(sizeof(HighsInt) * numrow);

  int modelstatus;

  int status = Highs_call(numcol, numrow, nnz, cc, cl, cu, rl, ru, astart, aindex, avalue, cv,
            cd, rv, rd, cbs, rbs, &modelstatus);
  assert(status == 0);

  for (i = 0; i < numcol; i++) {
    printf("x%"HIGHSINT_FORMAT" = %lf\n", i, cv[i]);
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

  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  HighsInt astart[3] = {0, 2, 4};
  HighsInt aindex[4] = {0, 1, 0, 1};
  double avalue[4] = {1.0, 2.0, 1.0, 3.0};

  assert( Highs_addCols(highs, 2, cc, cl, cu, 0, NULL, NULL, NULL) );
  assert( Highs_addRows(highs, 2, rl, ru,  4, astart, aindex, avalue) );

  Highs_run(highs);
  Highs_destroy(highs);
}

void options() {
  void* highs = Highs_create();

<<<<<<< HEAD
  HighsInt simplex_scale_strategy;
  Highs_setHighsIntOptionValue(highs, "simplex_scale_strategy", 0);
  Highs_getHighsIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
=======
  int simplex_scale_strategy;
  Highs_setIntOptionValue(highs, "simplex_scale_strategy", 0);
  Highs_getIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
>>>>>>> f90-interface
  assert( simplex_scale_strategy == 0 );

  double primal_feasibility_tolerance;
  Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", 2.0);
  Highs_getDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  assert( primal_feasibility_tolerance == 2.0 );

  Highs_destroy(highs);
}

void test_getColsByRange() {
    void* highs = Highs_create();
    Highs_addCol(highs, -1.0, 0.0, 1.0, 0, NULL, NULL);
    Highs_addCol(highs, -1.0, 0.0, 1.0, 0, NULL, NULL);
    HighsInt aindex[2] = {0, 1};
    double avalue[2] = {1.0, -1.0};
    Highs_addRow(highs, 0.0, 0.0, 2, aindex, avalue);
    HighsInt num_cols;
    HighsInt num_nz;
    HighsInt matrix_start[2] = {-1, -1};
    Highs_getColsByRange(highs, 0, 1, &num_cols, NULL, NULL, NULL, &num_nz,
                         matrix_start, NULL, NULL);
    assert( num_cols == 2 );
    assert( num_nz == 2 );
    assert( matrix_start[0] == 0 );
    assert( matrix_start[1] == 1 );
    HighsInt matrix_indices[2] = {-1, -1};
    double matrix_values[2] = {0.0, 0.0};
    Highs_getColsByRange(highs, 0, 1, &num_cols, NULL, NULL, NULL, &num_nz,
                         matrix_start, matrix_indices, matrix_values);
    assert( matrix_indices[0] == 0 );
    assert( matrix_indices[1] == 0 );
    assert( matrix_values[0] == 1.0 );
    assert( matrix_values[1] == -1.0 );
    Highs_destroy(highs);
}

int main() {
  minimal_api();
  full_api();
  options();
  test_getColsByRange();
  return 0;
}
