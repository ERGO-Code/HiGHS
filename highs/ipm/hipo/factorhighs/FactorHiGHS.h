#ifndef FACTOR_HIGHS_H
#define FACTOR_HIGHS_H

#include "DataCollector.h"
#include "Numeric.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/Log.h"

/*

Direct solver for IPM matrices.
It requires Metis and BLAS.

Consider a sparse symmetric matrix M in CSC format.
Only its lower triangular part is used; entries in the upper triangle are
ignored.
The matrix has n rows/cols and nz nonzero entries in the lower triangle.
It is stored using three arrays:
- ptr, column pointers, of length n+1;
- rows, row indices, of length nz;
- vals, values, of length nz.

The direct solver uses the following objects:
- Symbolic, to store the symbolic factorization;
- Numeric, to store the numeric factorization;
- FHsolver, to perform analyse and factorise phases.

Define a vector signs that contains the expected sign of each pivot (1 or -1).
Define a right-hand side rhs, which will be overwritten with the solution of
M^{-1} * rhs.

Then, the factorization is performed as follows.

    Symbolic S;
    Numeric N(S);

    FHsolver FH;
    FH.analyse(S, rows, ptr, signs);
    FH.factorise(N, S, rows, ptr, val);

    N.solve(rhs);

Printing to screen is achieved using the interface in auxiliary/Log.h. Pass an
object of type Log for normal printing:
    ...
    Log log;
    FHsolver FH(&log);
    ...
Pass an object of type LogHighs for Highs logging.
Pass nothing to suppress all logging.

To add static regularisation when the pivots are selected, use
setRegularisation(reg_p,reg_d) to choose values of primal and dual
regularisation. If regularisation is already added to the matrix, ignore.

The default block size is 128. To set a different block size, pass it as second
input to the constructor.

*/

namespace hipo {

class FHsolver {
  const Log* log_;
  DataCollector data_;
  Regul regul_;

  const Int nb_;  // block size
  static const Int default_nb_ = 128;

  // columns of factorisation, stored by supernode
  std::vector<std::vector<double>> sn_columns_;

 public:
  // Create object and initialise DataCollector
  FHsolver(const Log* log = nullptr, Int block_size = default_nb_);

  // Print collected data (if any) and terminate DataCollector
  ~FHsolver();

  // Perform analyse phase of matrix with sparsity pattern given by rows and
  // ptr, and store symbolic factorisation in object S.
  // See ReturnValues.h for errors.
  Int analyse(Symbolic& S, const std::vector<Int>& rows,
              const std::vector<Int>& ptr, const std::vector<Int>& signs);

  // Perform factorise phase of matrix given by rows, ptr, vals, and store
  // numerical factorisation in object N. Matrix is moved into the object, so
  // rows, ptr, vals are invalid afterwards.
  // See ReturnValues.h for errors.
  Int factorise(Numeric& N, const Symbolic& S, const std::vector<Int>& rows,
                const std::vector<Int>& ptr, const std::vector<double>& vals);

  // If multiple factorisation are performed, call newIter() before each
  // factorisation. This is used only to collect data for debugging, if
  // expensive data collection is turned on at compile time.
  void newIter();

  // Set values for static regularisation to be added when a pivot is selected.
  // If regularisation is already added to the matrix, ignore.
  void setRegularisation(double reg_p, double reg_d);
};

/* To do

- At the moment, the format handler allocates new space for frontal each time,
  and then moves the values into sn_columns_, with the previous space being
  deallocated. This is very inefficient, because it causes many memory
  allocation/deallocation. It may also cause fragmentation.
- Have sn_columns_ outside of factorise, into FHsolver. Then, pass it into
  Factorise and assemble the contributions directly in place. Numeric then
  receives a pointer to the data, so that it can access the data in the same way
  as now. Format handler needs to assign zeros during initialization each time.
- In this way, the space for the columns of the factor is persistent throughout
  the ipm iterations, no extra allocation/deallocation is needed, no copies are
  needed. Potentially, there may be more contention of the cache lines, but this
  is to be seen.

- For the moment, the same can be done with cliques. No need to
  allocate/deallocate them. They can stay allocated and be reused. This would
  probably considerably increase the memory usage though and the stack approach
  is to be preferred. However, the latter requires a different parallelisation.

*/

}  // namespace hipo

#endif