#include "Numeric.h"

#include <algorithm>

#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "HybridSolveHandler.h"
#include "ReturnValues.h"
#include "ScratchArena.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace hipo {

void Numeric::computeAnySwaps() {
  // compute which blocks of columns require swaps
  // one entry per column, of which only the first n_blocks of each supernode
  // are used
  any_swaps_.assign(swaps_.size(), 0);
  const Int nb = options_->nb;
  for (Int sn = 0; sn < S_->sn(); ++sn) {
    const Int sn_start = S_->snStart(sn);
    const Int sn_size = S_->snStart(sn + 1) - sn_start;
    const Int n_blocks = (sn_size - 1) / nb + 1;

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      for (Int i = 0; i < jb; ++i) {
        if (swaps_[sn_start + nb * j + i] != i) {
          any_swaps_[sn_start + j] = 1;
          break;
        }
      }
    }
  }
}

Int Numeric::prepare() {
  if (!factor_ || !S_ || !data_ || !options_ || !scratch_)
    return kRetInvalidPointer;

  // Rebuild the solve handler only if needed
  const bool reuse = SH_ && SH_symbolic_ == S_ && SH_symbolic_id_ == S_->id() &&
                     SH_parallel_forward_ == options_->parallel_forward &&
                     SH_parallel_backward_ == options_->parallel_backward;
  if (!reuse) {
    SH_.reset(new HybridSolveHandler(*S_, *factor_, swaps_, any_swaps_,
                                     pivot_2x2_, *data_, *options_, *scratch_));
    SH_symbolic_ = S_;
    SH_symbolic_id_ = S_->id();
    SH_parallel_forward_ = options_->parallel_forward;
    SH_parallel_backward_ = options_->parallel_backward;
  }

  if (options_->pivoting) computeAnySwaps();
  return kRetOk;
}

// The solves work on a copy of the rhs in scratch memory, permuted, instead of
// on x: BLAS kernels may round differently depending on the alignment of their
// operands, and the scratch memory is aligned, while x may not be.

Int Numeric::solve(double* x) const {
  // Return the number of solves performed

  // scratch_ is set whenever SH_ is; it must have a stack for each thread
  if (!SH_ || !scratch_->matchesThreads()) return kRetGeneric;

  HIPO_CLOCK_CREATE;

  const std::vector<Int>& iperm = S_->iperm();
  const Int n = iperm.size();
  ScratchBuffer<double> y(*scratch_, n);

  // permute rhs
  HIPO_CLOCK_START(2);
  for (Int i = 0; i < n; ++i) y[iperm[i]] = x[i];
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  // solve
  HIPO_CLOCK_START(2);
  SH_->forwardSolve(y.data());
  SH_->diagSolve(y.data());
  SH_->backwardSolve(y.data());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolveSolve);

  // unpermute solution
  HIPO_CLOCK_START(2);
  for (Int i = 0; i < n; ++i) x[i] = y[iperm[i]];
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  HIPO_CLOCK_STOP(1, *data_, kTimeSolve);
  return kRetOk;
}

Int Numeric::forwardSolve(double* x) const {
  if (!SH_ || !scratch_->matchesThreads()) return kRetGeneric;
  const std::vector<Int>& iperm = S_->iperm();
  const Int n = iperm.size();
  ScratchBuffer<double> y(*scratch_, n);
  for (Int i = 0; i < n; ++i) y[iperm[i]] = x[i];
  SH_->forwardSolve(y.data());
  std::copy_n(y.data(), n, x);
  return kRetOk;
}
Int Numeric::diagSolve(double* x) const {
  if (!SH_) return kRetGeneric;
  SH_->diagSolve(x);
  return kRetOk;
}
Int Numeric::backwardSolve(double* x) const {
  if (!SH_ || !scratch_->matchesThreads()) return kRetGeneric;
  const std::vector<Int>& iperm = S_->iperm();
  const Int n = iperm.size();
  ScratchBuffer<double> y(*scratch_, n);
  std::copy_n(x, n, y.data());
  SH_->backwardSolve(y.data());
  for (Int i = 0; i < n; ++i) x[i] = y[iperm[i]];
  return kRetOk;
}

// Serial: the solve handler workspace is shared, so rhs cannot run in parallel
#define SOLVE_MULTIPLE(f)       \
  if (k == 1) return f(x);      \
  const Int n = S_->size();     \
  for (Int i = 0; i < k; ++i) { \
    Int status = f(&x[i * n]);  \
    if (status) return status;  \
  }                             \
  return kRetOk;

Int Numeric::solve(double* x, Int k) const { SOLVE_MULTIPLE(solve); }
Int Numeric::forwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(forwardSolve);
}
Int Numeric::diagSolve(double* x, Int k) const { SOLVE_MULTIPLE(diagSolve); }
Int Numeric::backwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(backwardSolve);
}

void Numeric::getReg(double* reg) {
  std::memcpy(reg, total_reg_.data(), total_reg_.size() * sizeof(double));
}

void Numeric::getUnregPivots(std::vector<double>& unreg_pivots) const {
  unreg_pivots = unreg_pivots_;
}

Int Numeric::num2x2() const {
  // only the first column of a 2x2 pivot stores a nonzero entry
  Int num = 0;
  for (double v : pivot_2x2_)
    if (v != 0.0) ++num;
  return num;
}

void Numeric::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  if (!SH_) return;
  SH_->inertia(pos, neg, zero, tol);
}

}  // namespace hipo