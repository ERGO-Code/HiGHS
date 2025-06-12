#ifndef FACTOR_HIGHS_H
#define FACTOR_HIGHS_H

#include "Analyse.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "Factorise.h"
#include "Numeric.h"
#include "ReturnValues.h"
#include "Symbolic.h"

/*

Direct solver for IPM matrices.
It requires Metis and BLAS.

Consider a sparse symmetric matrix M in CSC format.
Only its lower triangular part is used.
The matrix has n rows/cols and nz nonzero entries in the lower triangle.
It is stored using three arrays:
- ptr, column pointers, of length n+1;
- rows, row indices, of length nz;
- val, values, of length nz.

The direct solver uses the following objects:
- Symbolic, to store the symbolic factorization;
- Numeric, to store the numeric factorization;
- Analyse, to perform the analyse phase;
- Factorise, to perform the factorise phase.

Define the integer neg, which contains the size of the (1,1)-block for the
augmented system, or 0 for the normal equations.

Define a right-hand side rhs, which will be overwritten with the solution of
M^{-1} * rhs.

Then, the factorization is performed as follows.

    Symbolic S;
    Numeric N(S);

    Analyse analyse(S, rows, ptr, neg);
    analyse.run();

    Factorise factorise(S, rows, ptr, val);
    factorise.run(N);

    N.solve(rhs);

*/

#endif