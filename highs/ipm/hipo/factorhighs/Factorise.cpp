#include "Factorise.h"

#include <algorithm>
#include <fstream>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "HybridHybridFormatHandler.h"
#include "ReturnValues.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Factorise::Factorise(const Symbolic& S, const std::vector<Int>& rowsM,
                     const std::vector<Int>& ptrM,
                     const std::vector<double>& valM, const Regul& regul,
                     const Logger* logger, DataCollector& data,
                     std::vector<std::vector<double>>& sn_columns,
                     CliqueStack* stack)
    : S_{S},
      sn_columns_{sn_columns},
      regul_{regul},
      logger_{logger},
      data_{data},
      stack_{stack} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyse.
  // Only the lower triangular part of the matrix is used.

  n_ = ptrM.size() - 1;

  if (n_ != S_.size()) {
    if (logger_)
      logger_->printInfo(
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
  permuteSym(iperm, ptrM_, rowsM_, valM_, true);
}

void Factorise::assembleInitial(Int sn) {
  // This can execute while the children are still running.

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  HIPO_CLOCK_CREATE;

  const bool parallel = S_.parTree();
  const bool serial = !parallel;

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
    if (reallocation && logger_)
      logger_->printInfo("Reallocation of CliqueStack\n");
  }

  // initialise the format handler
  // this also allocates space for the frontal matrix and schur complement
  std::unique_ptr<FormatHandler>& FH = sn_FH_[sn];
  FH.reset(new HybridHybridFormatHandler(S_, sn, regul_, data_, sn_columns_[sn],
                                         clique_ptr));

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
}

void Factorise::assembleChild(Int sn, Int child_sn) {
  // This must execute after syncing with the given child.

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  HIPO_CLOCK_CREATE;

  std::unique_ptr<FormatHandler>& FH = sn_FH_[sn];
  const bool parallel = S_.parTree();
  const bool serial = !parallel;
  const Int sn_begin = S_.snStart(sn);
  const Int sn_end = S_.snStart(sn + 1);
  const Int sn_size = sn_end - sn_begin;

  // Child contribution is found:
  // - in cliquestack, if we are processing the tree in serial.
  // - in schur_contribution_ if we are processing the tree in parallel.
  // Children are always summed from last to first.

  const double* child_clique;

  if (parallel) {
    child_clique = schur_contribution_[child_sn].data();
    if (!child_clique) {
      if (logger_)
        logger_->printInfo("%d: Missing child supernode contribution\n",
                           child_sn);
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
        Int consecutive = S_.consecutiveSums(child_sn, row);

        FH->assembleFrontalMultiple(consecutive, child_clique, nc, child_sn,
                                    row, col, i, j);

        row += consecutive;
      }
    } else
      break;
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
}

void Factorise::denseFactorise(Int sn) {
  // This must execute after syncing with all children.

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  HIPO_CLOCK_CREATE;

  std::unique_ptr<FormatHandler>& FH = sn_FH_[sn];
  const bool parallel = S_.parTree();
  const bool serial = !parallel;

  // ===================================================
  // Partial factorisation
  // ===================================================
  HIPO_CLOCK_START(2);
  // threshold for regularisation
  // const double reg_thresh = max_diag_ * kDynamicDiagCoeff;
  const double reg_thresh = M_norm1_ * kDynamicDiagCoeff;

  if (Int flag = FH->denseFactorise(reg_thresh)) {
    flag_stop_.store(true, std::memory_order_relaxed);

    if (logger_ && flag == kRetInvalidInput)
      logger_->printInfo("DenseFact: invalid input\n");
    else if (logger_ && flag == kRetInvalidPivot)
      logger_->printInfo("DenseFact: invalid pivot\n");

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

  FH.reset();

  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseTerminate);
}

void Factorise::processSupernode(Int sn) {
  // Process a supernode with no parallelism

  assembleInitial(sn);

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  Int child_sn = first_child_reverse_[sn];
  while (child_sn != -1) {
    assembleChild(sn, child_sn);
    child_sn = next_child_reverse_[child_sn];
  }

  if (flag_stop_.load(std::memory_order_relaxed)) return;

  denseFactorise(sn);
}

void Factorise::processTask(Int task) {
  const std::vector<Int>& sns = S_.schedule().sn_per_task[task];

  // Spawn all children tasks
  highs::parallel::TaskGroup tg;
  Int child_task = first_task_child_[task];
  while (child_task != -1) {
    tg.spawn([=]() { processTask(child_task); });
    child_task = next_task_child_[child_task];
  }

  for (Int sn : sns) {
    if (first_child_[sn] == -1)
      processSupernode(sn);  // leaf sn
    else
      assembleInitial(sn);  // sn dependent on other sn
  }

  tg.taskWait();

  for (Int sn : sns) {
    if (first_child_[sn] == -1) continue;

    Int child_sn = first_child_[sn];
    while (child_sn != -1) {
      assembleChild(sn, child_sn);
      child_sn = next_child_[child_sn];
    }

    denseFactorise(sn);
  }
}

void Factorise::processParallelTree() {
  const TreeSchedule& sched = S_.schedule();
  childrenLinkedList(sched.task_parent, first_task_child_, next_task_child_);

  first_task_child_reverse_ = first_task_child_;
  next_task_child_reverse_ = next_task_child_;
  reverseLinkedList(first_task_child_reverse_, next_task_child_reverse_);

  highs::parallel::TaskGroup tg;
  for (Int task = 0; task < sched.count(); ++task) {
    if (sched.task_parent[task] == -1) {
      tg.spawn([=]() { processTask(task); });
    }
  }
  tg.taskWait();
}

void Factorise::processSerialTree() {
  // processing the tree in serial requires a CliqueStack
  if (!stack_) {
    flag_stop_.store(true, std::memory_order_relaxed);
    return;
  }

  if (stack_->empty()) stack_->init(S_.maxStackSize());

  for (Int sn = 0; sn < S_.sn(); ++sn) {
    processSupernode(sn);
  }
}

bool Factorise::run(Numeric& num) {
  HIPO_CLOCK_CREATE;

  highs::parallel::TaskGroup tg;

  total_reg_.assign(n_, 0.0);

  // allocate space
  schur_contribution_.resize(S_.sn());
  swaps_.resize(S_.sn());
  pivot_2x2_.resize(S_.sn());
  sn_FH_.resize(S_.sn());

  // This should actually allocate only the first time, then sn_columns_ reuses
  // the memory of previous factorisations.
  sn_columns_.resize(S_.sn());

  const bool use_parallel_tree = S_.parTree() && S_.schedule().valid;

  if (use_parallel_tree)
    processParallelTree();
  else
    processSerialTree();

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
