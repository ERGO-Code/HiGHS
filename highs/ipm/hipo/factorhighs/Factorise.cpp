#include "Factorise.h"

#include <algorithm>
#include <fstream>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "FormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "ReturnValues.h"
#include "SymScaling.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Factorise::Factorise(const Symbolic& S, const std::vector<Int>& rowsA,
                     const std::vector<Int>& ptrA,
                     const std::vector<double>& valA, const Log* log)
    : S_{S}, log_{log} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyse.
  // Only the lower triangular part of the matrix is used.
  // The Factorise object takes ownership of the matrix; rowsA, ptrA and valA
  // are not valid anymore.

  n_ = ptrA.size() - 1;

  if (n_ != S_.size()) {
    if (log_)
      log_->printDevInfo(
          "Matrix provided to Factorise has size incompatible with symbolic "
          "object.\n");
    return;
  }

  // take ownership of the matrix
  rowsA_ = std::move(rowsA);
  valA_ = std::move(valA);
  ptrA_ = std::move(ptrA);

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  permute(S_.iperm());

  nzA_ = ptrA_.back();

  // Double transpose to sort columns
  std::vector<Int> temp_ptr(n_ + 1);
  std::vector<Int> temp_rows(nzA_);
  std::vector<double> temp_val(nzA_);
  transpose(ptrA_, rowsA_, valA_, temp_ptr, temp_rows, temp_val);
  transpose(temp_ptr, temp_rows, temp_val, ptrA_, rowsA_, valA_);

  // create linked lists of children in supernodal elimination tree
  childrenLinkedList(S_.snParent(), first_child_, next_child_);

  if (S_.parTree()) {
    // create reverse linked lists of children
    first_child_reverse_ = first_child_;
    next_child_reverse_ = next_child_;
    reverseLinkedList(first_child_reverse_, next_child_reverse_);
  }

  // compute largest diagonal entry in absolute value
  max_diag_ = 0.0;
  min_diag_ = std::numeric_limits<double>::max();
  for (Int col = 0; col < n_; ++col) {
    double val = std::abs(valA_[ptrA_[col]]);
    max_diag_ = std::max(max_diag_, val);
    min_diag_ = std::min(min_diag_, val);
  }

  // infinity norm of columns of A
  inf_norm_cols_.assign(n_, 0.0);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      Int row = rowsA_[el];
      double val = valA_[el];
      inf_norm_cols_[col] = std::max(inf_norm_cols_[col], std::abs(val));
      if (row != col)
        inf_norm_cols_[row] = std::max(inf_norm_cols_[row], std::abs(val));
    }
  }

  // one norm of columns of A
  one_norm_cols_.assign(n_, 0.0);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      Int row = rowsA_[el];
      double val = valA_[el];
      one_norm_cols_[col] += std::abs(val);
      if (row != col) one_norm_cols_[row] += std::abs(val);
    }
  }
  A_norm1_ = *std::max_element(one_norm_cols_.begin(), one_norm_cols_.end());

  DataCollector::get()->setNorms(A_norm1_, max_diag_);
}

void Factorise::permute(const std::vector<Int>& iperm) {
  // Symmetric permutation of the lower triangular matrix A based on inverse
  // permutation iperm.
  // The resulting matrix is lower triangular, regardless of the input matrix.

  std::vector<Int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (Int j = 0; j < n_; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const Int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const Int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      Int actual_col = std::min(row, col);
      ++work[actual_col];
    }
  }

  std::vector<Int> new_ptr(n_ + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(new_ptr, work);

  std::vector<Int> new_rows(new_ptr.back());
  std::vector<double> new_val(new_ptr.back());

  // go through the columns to assign row indices
  for (Int j = 0; j < n_; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const Int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const Int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      const Int actual_col = std::min(row, col);
      const Int actual_row = std::max(row, col);

      Int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      new_val[pos] = valA_[el];
    }
  }

  ptrA_ = std::move(new_ptr);
  rowsA_ = std::move(new_rows);
  valA_ = std::move(new_val);
}

std::unique_ptr<FormatHandler> getFormatHandler(const Symbolic& S, Int sn) {
  std::unique_ptr<FormatHandler> ptr;
  ptr.reset(new HybridHybridFormatHandler(S, sn));
  return ptr;
}

void Factorise::processSupernode(Int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.

  if (flag_stop_) return;

  if (S_.parTree()) {
    // spawn children of this supernode in reverse order
    Int child_to_spawn = first_child_reverse_[sn];
    while (child_to_spawn != -1) {
      highs::parallel::spawn([=]() { processSupernode(child_to_spawn); });
      child_to_spawn = next_child_reverse_[child_to_spawn];
    }

    // wait for first child to finish, before starting the parent (if there is a
    // first child)
    if (first_child_reverse_[sn] != -1) highs::parallel::sync();
  }

#if HIPO_TIMING_LEVEL >= 2
  Clock clock;
#endif
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  const Int sn_begin = S_.snStart(sn);
  const Int sn_end = S_.snStart(sn + 1);
  const Int sn_size = sn_end - sn_begin;

  // initialise the format handler
  // this also allocates space for the frontal matrix and schur complement
  std::unique_ptr<FormatHandler> FH = getFormatHandler(S_, sn);

#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeFactorisePrepare, clock.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock.start();
#endif
  // ===================================================
  // Assemble original matrix A into frontal
  // ===================================================
  // j is relative column index in the frontal matrix
  for (Int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    const Int col = sn_begin + j;

    // go through the column
    for (Int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      // relative row index in the frontal matrix
      const Int i = S_.relindCols(el);

      FH->assembleFrontal(i, j, valA_[el]);
    }
  }
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeFactoriseAssembleOriginal, clock.stop());
#endif

  // ===================================================
  // Assemble frontal matrices of children
  // ===================================================
  Int child_sn = first_child_[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    std::vector<double>& child_clique = schur_contribution_[child_sn];

    if (S_.parTree()) {
      // sync with spawned child, apart from the first one
      if (child_sn != first_child_[sn]) highs::parallel::sync();

      if (flag_stop_) return;

      if (child_clique.size() == 0) {
        if (log_) log_->printDevInfo("Missing child supernode contribution\n");
        flag_stop_ = true;
        return;
      }
    }

    // determine size of clique of child
    const Int child_begin = S_.snStart(child_sn);
    const Int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const Int child_size = child_end - child_begin;

    // size of clique of child sn
    const Int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

// ASSEMBLE INTO FRONTAL
#if HIPO_TIMING_LEVEL >= 2
    clock.start();
#endif
    // go through the columns of the contribution of the child
    for (Int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      Int j = S_.relindClique(child_sn, col);

      if (j < sn_size) {
        // assemble into frontal

        // go through the rows of the contribution of the child
        Int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix frontal
          const Int i = S_.relindClique(child_sn, row);

          // how many entries to sum
          const Int consecutive = S_.consecutiveSums(child_sn, row);

          FH->assembleFrontalMultiple(consecutive, child_clique, nc, child_sn,
                                      row, col, i, j);

          row += consecutive;
        }
      }
    }
#if HIPO_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeFactoriseAssembleChildrenFrontal,
                                  clock.stop());
#endif

// ASSEMBLE INTO CLIQUE
#if HIPO_TIMING_LEVEL >= 2
    clock.start();
#endif
    FH->assembleClique(child_clique, nc, child_sn);
#if HIPO_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeFactoriseAssembleChildrenClique,
                                  clock.stop());
#endif

    // Schur contribution of the child is no longer needed
    // Swap with temporary empty vector to deallocate memory
    std::vector<double> temp_empty;
    schur_contribution_[child_sn].swap(temp_empty);

    // move on to the next child
    child_sn = next_child_[child_sn];
  }

  if (flag_stop_) return;

    // ===================================================
    // Partial factorisation
    // ===================================================
#if HIPO_TIMING_LEVEL >= 2
  clock.start();
#endif
  // threshold for regularisation
  // const double reg_thresh = max_diag_ * kDynamicDiagCoeff;
  const double reg_thresh = A_norm1_ * kDynamicDiagCoeff;

  if (Int flag = FH->denseFactorise(reg_thresh)) {
    flag_stop_ = true;

    if (log_ && flag == kRetInvalidInput)
      log_->printDevInfo("DenseFact: invalid input\n");
    else if (log_ && flag == kRetInvalidPivot)
      log_->printDevInfo("DenseFact: invalid pivot\n");

    return;
  }
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeFactoriseDenseFact, clock.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock.start();
#endif
  // compute largest elements in factorisation
  FH->extremeEntries();

  // terminate the format handler
  FH->terminate(sn_columns_[sn], schur_contribution_[sn], total_reg_,
                swaps_[sn], pivot_2x2_[sn]);
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeFactoriseTerminate, clock.stop());
#endif
}

bool Factorise::run(Numeric& num) {
#if HIPO_TIMING_LEVEL >= 1
  Clock clock;
#endif

  total_reg_.assign(n_, 0.0);

  // allocate space for list of generated elements and columns of L
  schur_contribution_.resize(S_.sn());
  sn_columns_.resize(S_.sn());
  swaps_.resize(S_.sn());
  pivot_2x2_.resize(S_.sn());

  if (S_.parTree()) {
    Int spawned_roots{};
    // spawn tasks for root supernodes
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      if (S_.snParent(sn) == -1) {
        highs::parallel::spawn([=]() { processSupernode(sn); });
        ++spawned_roots;
      }
    }

    // sync tasks for root supernodes
    for (Int root = 0; root < spawned_roots; ++root) {
      highs::parallel::sync();
    }
  } else {
    // go through each supernode serially
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      processSupernode(sn);
    }
  }

  if (flag_stop_) return true;

  // move factorisation to numerical object
  num.sn_columns_ = std::move(sn_columns_);
  num.total_reg_ = std::move(total_reg_);
  num.swaps_ = std::move(swaps_);
  num.pivot_2x2_ = std::move(pivot_2x2_);
  num.ptrA_ = std::move(ptrA_);
  num.rowsA_ = std::move(rowsA_);
  num.valA_ = std::move(valA_);
  num.one_norm_cols_ = std::move(one_norm_cols_);
  num.inf_norm_cols_ = std::move(inf_norm_cols_);

#if HIPO_TIMING_LEVEL >= 1
  DataCollector::get()->sumTime(kTimeFactorise, clock.stop());
#endif

  return false;
}

}  // namespace hipo
