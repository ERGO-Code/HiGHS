#include "Factorise.h"

#include <algorithm>
#include <fstream>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "FormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "ReturnValues.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Factorise::Factorise(const Symbolic& S, const std::vector<Int>& rowsM,
                     const std::vector<Int>& ptrM,
                     const std::vector<double>& valM, const Regul& regul,
                     const Log* log, DataCollector& data,
                     std::vector<std::vector<double>>& sn_columns,
                     CliqueStack* stack)
    : S_{S},
      sn_columns_{sn_columns},
      regul_{regul},
      log_{log},
      data_{data},
      stack_{stack} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyse.
  // Only the lower triangular part of the matrix is used.

  n_ = ptrM.size() - 1;

  if (n_ != S_.size()) {
    if (log_)
      log_->printDevInfo(
          "Matrix provided to Factorise has size incompatible with symbolic "
          "object.\n");
    return;
  }

  // Make a copy of the matrix to be factorised
  rowsM_ = rowsM;
  valM_ = valM;
  ptrM_ = ptrM;

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  permute(S_.iperm());

  nzM_ = ptrM_.back();

  // Double transpose to sort columns
  std::vector<Int> temp_ptr(n_ + 1);
  std::vector<Int> temp_rows(nzM_);
  std::vector<double> temp_val(nzM_);
  transpose(ptrM_, rowsM_, valM_, temp_ptr, temp_rows, temp_val);
  transpose(temp_ptr, temp_rows, temp_val, ptrM_, rowsM_, valM_);

  // create linked lists of children in supernodal elimination tree
  childrenLinkedList(S_.snParent(), first_child_, next_child_);

  // create reverse linked lists of children
  first_child_reverse_ = first_child_;
  next_child_reverse_ = next_child_;
  reverseLinkedList(first_child_reverse_, next_child_reverse_);

  // compute largest diagonal entry in absolute value
  max_diag_ = 0.0;
  min_diag_ = kHighsInf;
  for (Int col = 0; col < n_; ++col) {
    double val = std::abs(valM_[ptrM_[col]]);
    max_diag_ = std::max(max_diag_, val);
    min_diag_ = std::min(min_diag_, val);
  }

  // one norm of columns of M
  std::vector<double> one_norm_cols(n_, 0.0);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = ptrM_[col]; el < ptrM_[col + 1]; ++el) {
      Int row = rowsM_[el];
      double val = valM_[el];
      one_norm_cols[col] += std::abs(val);
      if (row != col) one_norm_cols[row] += std::abs(val);
    }
  }
  M_norm1_ = *std::max_element(one_norm_cols.begin(), one_norm_cols.end());

  data_.setNorms(M_norm1_, max_diag_);
}

void Factorise::permute(const std::vector<Int>& iperm) {
  // Symmetric permutation of the lower triangular matrix M based on inverse
  // permutation iperm.
  // The resulting matrix is lower triangular, regardless of the input matrix.

  std::vector<Int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (Int j = 0; j < n_; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptrM_[j]; el < ptrM_[j + 1]; ++el) {
      const Int i = rowsM_[el];

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
    for (Int el = ptrM_[j]; el < ptrM_[j + 1]; ++el) {
      const Int i = rowsM_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const Int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      const Int actual_col = std::min(row, col);
      const Int actual_row = std::max(row, col);

      Int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      new_val[pos] = valM_[el];
    }
  }

  ptrM_ = std::move(new_ptr);
  rowsM_ = std::move(new_rows);
  valM_ = std::move(new_val);
}

class TaskGroupSpecial : public highs::parallel::TaskGroup {
  // Using TaskGroup may throw an exception when tasks are cancelled. Not sure
  // exactly why this happens, but for now this fix seems to work.

 public:
  ~TaskGroupSpecial() {
    // No virtual destructor in TaskGroup. Do not call this class via pointer to
    // the base!

    cancel();

    // re-call taskWait if it throws, until it succeeds
    while (true) {
      try {
        taskWait();
        break;
      } catch (HighsTask::Interrupt) {
        continue;
      }
    }
  }
};

void Factorise::processSupernode(Int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.

  TaskGroupSpecial tg;
  HIPO_CLOCK_CREATE;

  const bool parallel = S_.parTree();
  const bool serial = !parallel;

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  if (parallel) {
    // spawn children of this supernode in forward order
    Int child_to_spawn = first_child_[sn];
    while (child_to_spawn != -1) {
      tg.spawn([=]() { processSupernode(child_to_spawn); });
      child_to_spawn = next_child_[child_to_spawn];
    }

    // wait for first child to finish, before starting the parent (if there is a
    // first child)
    if (first_child_[sn] != -1) tg.sync();
  }

  // ===================================================
  // Supernode information
  // ===================================================
  HIPO_CLOCK_START(2);
  // first and last+1 column of the supernodes
  const Int sn_begin = S_.snStart(sn);
  const Int sn_end = S_.snStart(sn + 1);
  const Int sn_size = sn_end - sn_begin;

  // When the tree is processed in serial, use CliqueStack to store the cliques.
  // Otherwise, use local storage in FormatHandler.
  double* clique_ptr = nullptr;
  if (serial) {
    bool reallocation = false;
    clique_ptr = stack_->setup(S_.cliqueSize(sn), reallocation);
    if (reallocation && log_)
      log_->printDevInfo("Reallocation of CliqueStack\n");
  }

  // initialise the format handler
  // this also allocates space for the frontal matrix and schur complement
  std::unique_ptr<FormatHandler> FH(new HybridHybridFormatHandler(
      S_, sn, regul_, data_, sn_columns_[sn], clique_ptr));

  HIPO_CLOCK_STOP(2, data_, kTimeFactorisePrepare);

  // ===================================================
  // Assemble original matrix M into frontal
  // ===================================================
  HIPO_CLOCK_START(2);
  // j is relative column index in the frontal matrix
  for (Int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    const Int col = sn_begin + j;

    // go through the column
    for (Int el = ptrM_[col]; el < ptrM_[col + 1]; ++el) {
      // relative row index in the frontal matrix
      const Int i = S_.relindCols(el);

      FH->assembleFrontal(i, j, valM_[el]);
    }
  }
  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleOriginal);

  // ===================================================
  // Assemble frontal matrices of children
  // ===================================================
  Int child_sn = first_child_reverse_[sn];
  while (child_sn != -1) {
    // Child contribution is found:
    // - in cliquestack, if we are processing the tree in serial.
    // - in schur_contribution_ if we are processing the tree in parallel.
    // Children are always summed from last to first.

    const double* child_clique;

    if (parallel) {
      // sync with spawned child, apart from the first one
      if (child_sn != first_child_reverse_[sn]) tg.sync();

      if (flag_stop_.load(std::memory_order_relaxed)) return;

      child_clique = schur_contribution_[child_sn].data();

      if (!child_clique) {
        if (log_) log_->printDevInfo("Missing child supernode contribution\n");
        flag_stop_.store(true, std::memory_order_relaxed);
        return;
      }
    } else {
      Int child;
      child_clique = stack_->getChild(child);
      assert(child == child_sn);
    }

    // determine size of clique of child
    const Int child_begin = S_.snStart(child_sn);
    const Int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const Int child_size = child_end - child_begin;

    // size of clique of child sn
    const Int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

    // ASSEMBLE INTO FRONTAL
    HIPO_CLOCK_START(2);
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
    HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenFrontal);

    // ASSEMBLE INTO CLIQUE
    HIPO_CLOCK_START(2);
    FH->assembleClique(child_clique, nc, child_sn);
    HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenClique);

    // Schur contribution of the child is no longer needed
    if (parallel) {
      freeVector(schur_contribution_[child_sn]);
    } else {
      stack_->popChild();
    }

    // move on to the next child
    child_sn = next_child_reverse_[child_sn];
  }

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  // ===================================================
  // Partial factorisation
  // ===================================================
  HIPO_CLOCK_START(2);
  // threshold for regularisation
  // const double reg_thresh = max_diag_ * kDynamicDiagCoeff;
  const double reg_thresh = M_norm1_ * kDynamicDiagCoeff;

  if (Int flag = FH->denseFactorise(reg_thresh)) {
    flag_stop_.store(true, std::memory_order_relaxed);

    if (log_ && flag == kRetInvalidInput)
      log_->printDevInfo("DenseFact: invalid input\n");
    else if (log_ && flag == kRetInvalidPivot)
      log_->printDevInfo("DenseFact: invalid pivot\n");

    return;
  }
  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseDenseFact);

  HIPO_CLOCK_START(2);
  // compute largest elements in factorisation
  FH->extremeEntries();

  // terminate the format handler
  FH->terminate(schur_contribution_[sn], total_reg_, swaps_[sn],
                pivot_2x2_[sn]);

  if (serial) stack_->pushWork(sn);

  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseTerminate);
}

bool Factorise::run(Numeric& num) {
  HIPO_CLOCK_CREATE;

  TaskGroupSpecial tg;

  total_reg_.assign(n_, 0.0);

  // allocate space
  schur_contribution_.resize(S_.sn());
  swaps_.resize(S_.sn());
  pivot_2x2_.resize(S_.sn());

  // This should actually allocate only the first time, then sn_columns_ reuses
  // the memory of previous factorisations.
  sn_columns_.resize(S_.sn());

  if (S_.parTree()) {
    Int spawned_roots{};
    // spawn tasks for root supernodes
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      if (S_.snParent(sn) == -1) {
        tg.spawn([=]() { processSupernode(sn); });
        ++spawned_roots;
      }
    }

    // sync tasks for root supernodes
    tg.taskWait();
  } else {
    // processing the tree in serial requires a CliqueStack
    if (!stack_) return true;
    if (stack_->empty()) stack_->init(S_.maxStackSize());

    // go through each supernode serially
    for (Int sn = 0; sn < S_.sn(); ++sn) {
      processSupernode(sn);
    }
  }

  if (flag_stop_.load(std::memory_order_relaxed)) return true;

  // move factorisation to numerical object
  num.S_ = &S_;
  num.sn_columns_ = &sn_columns_;
  num.total_reg_ = std::move(total_reg_);
  num.swaps_ = std::move(swaps_);
  num.pivot_2x2_ = std::move(pivot_2x2_);
  num.data_ = &data_;

  HIPO_CLOCK_STOP(1, data_, kTimeFactorise);

  return false;
}

}  // namespace hipo
