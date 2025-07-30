#ifndef FACTOR_HIGHS_H
#define FACTOR_HIGHS_H

#include "FactorHiGHSSettings.h"
#include "Numeric.h"
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
- vals, values, of length nz.

The direct solver uses the following objects:
- Symbolic, to store the symbolic factorization;
- Numeric, to store the numeric factorization;
- FHsolver, to perform analyse and factorise phases.

Define the integer neg, which contains the size of the (1,1)-block for the
augmented system, or 0 for the normal equations.

Define a right-hand side rhs, which will be overwritten with the solution of
M^{-1} * rhs.

Then, the factorization is performed as follows.

    Symbolic S;
    Numeric N(S);

    FHsolver();
    FHsolver.analyse(S, rows, ptr, neg);
    FHsolver.factorise(N, S, rows, ptr, val);

    N.solve(rhs);

Printing to screen is achieved using the interface in auxiliary/Log.h, which
uses Highs logging. To use the linear solver outside of Highs, redefine the
functions in auxiliary/Log.h to print in other ways.

*/

namespace hipo {

class FHsolver {
 public:
  // Create object and initialise DataCollector
  FHsolver();

  // Print collected data (if any) and terminate DataCollector
  ~FHsolver();

  // Perform analyse phase of matrix with sparsity pattern given by rows and
  // ptr, and store symbolic factorisation in object S.
  // See ReturnValues.h for errors.
  Int analyse(Symbolic& S, const std::vector<Int>& rows,
              const std::vector<Int>& ptr, Int negative_pivots = 0) const;

  // Perform factorise phase of matrix given by rows, ptr, vals, and store
  // numerical factorisation in object N. Matrix is moved into the object, so
  // rows, ptr, vals are invalid afterwards.
  // See ReturnValues.h for errors.
  Int factorise(Numeric& N, const Symbolic& S, const std::vector<Int>& rows,
                const std::vector<Int>& ptr,
                const std::vector<double>& vals) const;

  // If multiple factorisation are performed, call newIter() before each
  // factorisation. This is used only to collect data for debugging, if
  // expensive data collection is turned on at compile time.
  void newIter() const;
};

}  // namespace hipo

#endif