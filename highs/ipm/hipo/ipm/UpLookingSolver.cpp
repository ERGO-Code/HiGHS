#include "UpLookingSolver.h"

#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

UpLookingSolver::UpLookingSolver(KktMatrix& kkt)
    : kkt_{kkt},
      ptr_{kkt_.ptrAS.empty() ? kkt_.ptrNE : kkt_.ptrAS},
      rows_{kkt_.rowsAS.empty() ? kkt_.rowsNE : kkt_.rowsAS},
      val_{kkt_.valAS.empty() ? kkt_.valNE : kkt_.valAS} {}

void UpLookingSolver::etreeAndCounts() {
  // compute elimination tree and column counts
  const Int n = ptr_.size() - 1;
  parent_.assign(n, -1);
  colcount_.assign(n, 1);
  std::vector<Int> ancestor(n, 0);

  for (Int col = 0; col < n; ++col) {
    ancestor[col] = col;
    for (Int el = ptr_[col]; el < ptr_[col + 1]; ++el) {
      Int row = rows_[el];
      while (ancestor[row] != col) {
        if (parent_[row] < 0) parent_[row] = col;
        ++colcount_[row];
        ancestor[row] = col;
        row = parent_[row];
      }
    }
  }
}

Int UpLookingSolver::setup() {
  // KktMatrix stores lower triangle, unpermuted.
  // Transform it into upper triangle and permute it.
  std::vector<Int> ptrT(ptr_.size()), rowsT(rows_.size());

  transpose(ptr_, rows_, ptrT, rowsT);
  ptr_ = std::move(ptrT);
  rows_ = std::move(rowsT);
  std::vector<double> empty;
  permuteSym(kkt_.iperm, ptr_, rows_, empty, false);

  etreeAndCounts();

  return 0;
}

}  // namespace hipo