#include "Factorise.h"

#include <algorithm>
#include <fstream>

#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "FormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "ReturnValues.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Factorise::Factorise(const Symbolic& S, Numeric& num, Int n, Int nz,
                     const Int* rowsM, const Int* ptrM, const double* valM,
                     const std::vector<double>& static_reg,
                     const FHoptions& FH_opt, const Logger& logger,
                     DataCollector& data, FactorStorage& factor,
                     CliqueStack* stack, ScratchArena& scratch)
    : S_{S},
      num_{num},
      factor_{factor},
      swaps_{num.swaps_},
      pivot_2x2_{num.pivot_2x2_},
      total_reg_{num.total_reg_},
      unreg_pivots_{num.unreg_pivots_},
      logger_{logger},
      data_{data},
      FH_opt_{FH_opt},
      stack_{stack},
      scratch_{scratch} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyse.
  // Only the lower triangular part of the matrix is used.
  // Diagonal entries must be stored for each column, even if zero.

  n_ = n;

  if (n_ != S_.size()) {
    logger_.printInfo(
        "Matrix provided to Factorise has size incompatible with symbolic "
        "object.\n");
    return;
  }

  // Make a copy of the matrix to be factorised
  rowsM_ = std::vector<Int>(rowsM, rowsM + nz);
  valM_ = std::vector<double>(valM, valM + nz);
  ptrM_ = std::vector<Int>(ptrM, ptrM + n + 1);

  // adjust data if one-based indexing is used
  if (FH_opt_.one_indexing) {
    for (Int& i : rowsM_) --i;
    for (Int& i : ptrM_) --i;
  }

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  permute(S_.iperm());

  // static regularisation in the same order as the permuted matrix: the given
  // values, or reg_p for negative pivots and reg_d for the others
  if (!static_reg.empty()) {
    static_reg_ = static_reg;
    permuteVectorInverse(static_reg_, S_.iperm());
  } else {
    const std::vector<Int>& sign = S_.pivotSign();
    static_reg_.resize(n_);
    for (Int i = 0; i < n_; ++i)
      static_reg_[i] = sign[i] < 0 ? FH_opt_.reg_p : FH_opt_.reg_d;
  }

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
  M_norm1_ = one_norm_cols.empty() ? 0.0
                                   : *std::max_element(one_norm_cols.begin(),
                                                       one_norm_cols.end());
}

void Factorise::permute(const std::vector<Int>& iperm) {
  permuteSym(iperm, ptrM_, rowsM_, valM_, true);
}

void Factorise::processSupernode(Int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.

  highs::parallel::TaskGroup tg;
  HIPO_CLOCK_CREATE;

  const bool parallel = FH_opt_.parallel_tree;
  const bool serial = !parallel;

  if (flag_stop_.load(std::memory_order_relaxed)) return;
  // Called from parallel tasks, so the stop callback must be thread safe
  if (FH_opt_.stopRequested()) {
    flag_stop_.store(true, std::memory_order_relaxed);
    return;
  }

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
    if (reallocation) logger_.printInfo("Reallocation of CliqueStack\n");
  }

  // initialise the format handler
  // this also initialises the frontal matrix, and allocates space for the schur
  // complement if CliqueStack is not used
  HybridHybridFormatHandler FH_obj(
      S_, sn, data_, factor_[sn], clique_ptr, FH_opt_, total_reg_.data(),
      unreg_pivots_.data(), swaps_.data(), pivot_2x2_.data());
  FormatHandler* FH = &FH_obj;

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

      TSAN_ANNOTATE_HAPPENS_AFTER(&schur_contribution_[child_sn]);

      child_clique = schur_contribution_[child_sn].data();

      if (!child_clique) {
        logger_.printInfo("Missing child supernode contribution\n");
        flag_stop_.store(true, std::memory_order_relaxed);
        return;
      }
    } else {
      Int child;
      child_clique = stack_->getChild(child);
      assert(child == child_sn);
    }

    FH->assembleChild(child_sn, child_clique);

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
  const double reg_thresh = M_norm1_ * kDynamicDiagCoeff;

  if (Int flag = FH->denseFactorise(reg_thresh, &static_reg_[S_.snStart(sn)],
                                    scratch_)) {
    flag_stop_.store(true, std::memory_order_relaxed);

    if (flag == kRetInvalidInput)
      logger_.printInfo("DenseFact: invalid input\n");
    else if (flag == kRetInvalidPivot)
      logger_.printInfo("DenseFact: invalid pivot\n");

    return;
  }
  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseDenseFact);

  HIPO_CLOCK_START(2);
  // terminate the format handler
  FH->terminate(schur_contribution_[sn]);

  if (serial) stack_->pushWork(sn);

  TSAN_ANNOTATE_HAPPENS_BEFORE(&schur_contribution_[sn]);

  HIPO_CLOCK_STOP(2, data_, kTimeFactoriseTerminate);
}

void Factorise::processTreeSerial() {
  if (!stack_) {
    // processing the tree in serial requires a CliqueStack
    flag_stop_.store(true, std::memory_order_relaxed);
    return;
  }
  if (stack_->empty()) stack_->init(S_.maxStackSize());
  for (Int sn = 0; sn < S_.sn(); ++sn) {
    processSupernode(sn);
  }
}

void Factorise::processTreeParallel() {
  highs::parallel::TaskGroup tg;
  // spawn roots
  for (Int sn = 0; sn < S_.sn(); ++sn) {
    if (S_.snParent(sn) == -1) {
      tg.spawn([=]() { processSupernode(sn); });
    }
  }
  tg.taskWait();
}

bool Factorise::run() {
  HIPO_CLOCK_CREATE;

  // The outputs reuse the memory of the previous factorisation. The values are
  // reset here, as with fresh vectors.
  total_reg_.assign(n_, 0.0);
  unreg_pivots_.assign(n_, 0.0);
  swaps_.assign(n_, 0);
  pivot_2x2_.assign(n_, 0.0);

  // allocate space
  schur_contribution_.resize(S_.sn());

  // Layout of the factor: the frontal matrix of each supernode, followed by
  // extra_space_frontal entries, one after the other. The memory is allocated
  // only the first time, then factor_ reuses the memory of previous
  // factorisations.
  factor_.setLayout(S_.sn(), [&](Int sn) {
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
    const Int ldf = (Int)(S_.ptr(sn + 1) - S_.ptr(sn));
    return HybridHybridFormatHandler::frontalSize(ldf, sn_size, FH_opt_.nb) +
           extra_space_frontal;
  });

  if (FH_opt_.parallel_tree) {
    processTreeParallel();
  } else {
    processTreeSerial();
  }

  if (flag_stop_.load(std::memory_order_relaxed)) return true;

  permuteVector(total_reg_, S_.iperm());
  permuteVector(unreg_pivots_, S_.iperm());

  // complete the numerical object; the outputs are already in it
  num_.S_ = &S_;
  num_.factor_ = &factor_;
  num_.data_ = &data_;
  num_.options_ = &FH_opt_;
  num_.scratch_ = &scratch_;

  if (num_.prepare()) return true;

  HIPO_CLOCK_STOP(1, data_, kTimeFactorise);

  return false;
}

}  // namespace hipo
