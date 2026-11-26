#ifndef FACTORHIGHS_PERMUTED_MATRIX_H
#define FACTORHIGHS_PERMUTED_MATRIX_H

#include <vector>

#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

class ScratchArena;

// Matrix to factorise, permuted as in the symbolic factorisation, lower
// triangular, with the rows of each column sorted. The pattern is fixed by the
// symbolic factorisation, so it is built once for each Symbolic, and each
// factorisation only writes the values again.
class PermutedMatrix {
  const Symbolic* symbolic_ = nullptr;
  Int64 symbolic_id_ = 0;
  Int nz_input_ = 0;
  bool one_indexing_ = false;

  // column pointers of the permuted matrix
  std::vector<Int> ptr_{};

  // map_[k] is the position in val_ of entry k of the input, or -1 if the entry
  // is in the upper triangle and is ignored
  std::vector<Int> map_{};

  std::vector<double> val_{};

 public:
  bool matches(const Symbolic& S, Int nz, bool one_indexing) const;

  // build the pattern of the input matrix rows, ptr, with nz entries
  void build(const Symbolic& S, Int nz, const Int* rows, const Int* ptr,
             bool one_indexing);

  // write the values of the input matrix, which must have the pattern given to
  // build, and return the 1-norm of the (symmetric) matrix
  double fill(const Int* rows, const Int* ptr, const double* vals,
              ScratchArena& scratch);

  Int nz() const { return val_.size(); }
  const std::vector<Int>& ptr() const { return ptr_; }
  const std::vector<double>& val() const { return val_; }
};

}  // namespace hipo

#endif
