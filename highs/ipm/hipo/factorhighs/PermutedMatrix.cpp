#include "PermutedMatrix.h"

#include <algorithm>
#include <cmath>

#include "ScratchArena.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

bool PermutedMatrix::matches(const Symbolic& S, Int nz,
                             bool one_indexing) const {
  return symbolic_ == &S && symbolic_id_ == S.id() && nz_input_ == nz &&
         one_indexing_ == one_indexing;
}

void PermutedMatrix::build(const Symbolic& S, Int nz, const Int* rows,
                           const Int* ptr, bool one_indexing) {
  const Int n = S.size();
  const Int offset = one_indexing ? 1 : 0;
  const std::vector<Int>& iperm = S.iperm();

  // 1. Position of each input entry in the permuted lower triangle; -1 for the
  // entries in the upper triangle of the input, which are ignored.
  std::vector<Int> perm_row(nz, -1);
  std::vector<Int> perm_col(nz, -1);
  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j] - offset; el < ptr[j + 1] - offset; ++el) {
      const Int i = rows[el] - offset;
      if (i < j) continue;
      perm_row[el] = std::max(iperm[i], iperm[j]);
      perm_col[el] = std::min(iperm[i], iperm[j]);
    }
  }

  // 2. Counting sort of the entries by row.
  std::vector<Int> next(n, 0);
  for (Int entry = 0; entry < nz; ++entry)
    if (perm_row[entry] >= 0) ++next[perm_row[entry]];
  std::vector<Int> row_start(n + 1);
  counts2Ptr(n, row_start.data(), next.data());

  std::vector<Int> sorted_by_row(row_start[n]);
  for (Int entry = 0; entry < nz; ++entry)
    if (perm_row[entry] >= 0) sorted_by_row[next[perm_row[entry]]++] = entry;

  // 3. Counting sort by column: the rows of each column end up sorted, as with
  // a double transpose.
  next.assign(n, 0);
  for (Int entry : sorted_by_row) ++next[perm_col[entry]];
  ptr_.resize(n + 1);
  counts2Ptr(n, ptr_.data(), next.data());

  map_.assign(nz, -1);
  for (Int entry : sorted_by_row) map_[entry] = next[perm_col[entry]]++;

  val_.resize(ptr_[n]);

  symbolic_ = &S;
  symbolic_id_ = S.id();
  nz_input_ = nz;
  one_indexing_ = one_indexing;
}

double PermutedMatrix::fill(const Int* rows, const Int* ptr, const double* vals,
                            ScratchArena& scratch) {
  const Int n = ptr_.empty() ? 0 : ptr_.size() - 1;
  const Int offset = one_indexing_ ? 1 : 0;

  ScratchBuffer<double> col_norm(scratch, n);
  std::fill_n(col_norm.data(), n, 0.0);

  for (Int j = 0; j < n; ++j) {
    for (Int el = ptr[j] - offset; el < ptr[j + 1] - offset; ++el) {
      if (map_[el] < 0) continue;
      const Int i = rows[el] - offset;
      const double val = vals[el];
      val_[map_[el]] = val;
      col_norm[j] += std::abs(val);
      if (i != j) col_norm[i] += std::abs(val);
    }
  }

  double norm = 0.0;
  for (Int j = 0; j < n; ++j) norm = std::max(norm, col_norm[j]);
  return norm;
}

}  // namespace hipo
