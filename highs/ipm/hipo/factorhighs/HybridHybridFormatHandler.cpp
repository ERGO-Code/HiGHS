#include "HybridHybridFormatHandler.h"

#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

HybridHybridFormatHandler::HybridHybridFormatHandler(const Symbolic& S, Int sn,
                                                     const Regul& regul,
                                                     DataCollector& data)
    : FormatHandler(S, sn, regul), data_{data} {
  // initialise frontal and clique
  initFrontal();
  initClique();
}

void HybridHybridFormatHandler::initFrontal() {
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;
  diag_start_.resize(n_blocks);
  Int frontal_size = getDiagStart(ldf_, sn_size_, nb_, n_blocks, diag_start_);
  frontal_.resize(frontal_size + extra_space);
  // NB: the plus 10 is not needed, but it avoids weird problems later on.
}

void HybridHybridFormatHandler::initClique() {
  clique_.resize(S_->cliqueSize(sn_));
}

void HybridHybridFormatHandler::assembleFrontal(Int i, Int j, double val) {
  Int block = j / nb_;
  Int ldb = ldf_ - block * nb_;
  Int ii = i - block * nb_;
  Int jj = j - block * nb_;
  frontal_[diag_start_[block] + ii + ldb * jj] = val;
}

void HybridHybridFormatHandler::assembleFrontalMultiple(
    Int num, const std::vector<double>& child, Int nc, Int child_sn, Int row,
    Int col, Int i, Int j) {
  const Int jblock = col / nb_;
  const Int jb = std::min(nb_, nc - nb_ * jblock);
  const Int row_ = row - jblock * nb_;
  const Int col_ = col - jblock * nb_;
  const Int start_block = S_->cliqueBlockStart(child_sn, jblock);

  Int block = j / nb_;
  Int ldb = ldf_ - block * nb_;
  Int ii = i - block * nb_;
  Int jj = j - block * nb_;

  callAndTime_daxpy(num, 1.0, &child[start_block + col_ + jb * row_], jb,
                    &frontal_[diag_start_[block] + ii + ldb * jj], 1, data_);
}

Int HybridHybridFormatHandler::denseFactorise(double reg_thresh) {
  Int status;

  status = denseFactFP2FH(frontal_.data(), ldf_, sn_size_, nb_, data_);
  if (status) return status;

  // find the position within pivot_sign corresponding to this supernode
  Int sn_start = S_->snStart(sn_);
  const Int* pivot_sign = &S_->pivotSign().data()[sn_start];

  status = denseFactFH('H', ldf_, sn_size_, S_->blockSize(), frontal_.data(),
                       clique_.data(), pivot_sign, reg_thresh, regul_,
                       local_reg_.data(), swaps_.data(), pivot_2x2_.data(),
                       S_->parNode(), data_);

  return status;
}

void HybridHybridFormatHandler::assembleClique(const std::vector<double>& child,
                                               Int nc, Int child_sn) {
  // assemble the child clique into the current clique by blocks of columns.
  // within a block, assemble by rows.

  const Int n_blocks = (nc - 1) / nb_ + 1;

  Int row_start{};

  // go through the blocks of columns of the child sn
  for (Int b = 0; b < n_blocks; ++b) {
    const Int b_start = S_->cliqueBlockStart(child_sn, b);

    const Int col_start = row_start;
    const Int col_end = std::min(col_start + nb_, nc);

    // go through the rows within this block
    for (Int row = row_start; row < nc; ++row) {
      const Int i = S_->relindClique(child_sn, row) - sn_size_;

      // already assembled into frontal
      if (i < 0) continue;

      // go through the columns of the block
      Int col = col_start;
      while (col < col_end) {
        Int j = S_->relindClique(child_sn, col);
        if (j < sn_size_) {
          ++col;
          continue;
        }
        j -= sn_size_;

        // information and sizes of child sn
        const Int jblock_c = b;
        const Int jb_c = std::min(nb_, nc - nb_ * jblock_c);
        const Int row_ = row - jblock_c * nb_;
        const Int col_ = col - jblock_c * nb_;
        const Int start_block_c = b_start;

        // sun consecutive entries in a row.
        // consecutive need to be reduced, to account for edge of the block
        const Int zeros_stored_row =
            std::max((Int)0, jb_c - (row - row_start) - 1);
        Int consecutive = S_->consecutiveSums(child_sn, col);
        const Int left_in_child = col_end - col - zeros_stored_row;
        consecutive = std::min(consecutive, left_in_child);

        // consecutive need to account also for edge of block in parent
        const Int block_in_parent = j / nb_;
        const Int col_end_parent = std::min((block_in_parent + 1) * nb_, ldc_);
        const Int left_in_parent = col_end_parent - j;
        consecutive = std::min(consecutive, left_in_parent);

        // needed to deal with zeros stored in upper right part of block
        if (consecutive == 0) break;

        // information and sizes of current sn
        const Int jblock = block_in_parent;
        const Int jb = std::min(nb_, ldc_ - nb_ * jblock);
        const Int i_ = i - jblock * nb_;
        const Int j_ = j - jblock * nb_;
        const Int start_block = S_->cliqueBlockStart(sn_, jblock);

        const double d_one = 1.0;
        const Int i_one = 1;
        callAndTime_daxpy(consecutive, 1.0,
                          &child[start_block_c + col_ + jb_c * row_], 1,
                          &clique_[start_block + j_ + jb * i_], 1, data_);

        col += consecutive;
      }
    }

    row_start += nb_;
  }
}

void HybridHybridFormatHandler::extremeEntries() {
  double minD = 1e100;
  double maxD = 0.0;
  double minoffD = 1e100;
  double maxoffD = 0.0;

  // number of blocks of columns
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;

  // index to access frontal
  Int index{};

  // go through blocks of columns for this supernode
  for (Int j = 0; j < n_blocks; ++j) {
    // number of columns in the block
    const Int jb = std::min(nb_, sn_size_ - nb_ * j);

    for (Int k = 0; k < jb; ++k) {
      // off diagonal entries
      for (Int i = 0; i < k; ++i) {
        if (frontal_[index] != 0.0) {
          minoffD = std::min(minoffD, std::abs(frontal_[index]));
          maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
        }
        index++;
      }

      // diagonal entry
      minD = std::min(minD, std::abs(1.0 / frontal_[index]));
      maxD = std::max(maxD, std::abs(1.0 / frontal_[index]));

      index += jb - k;
    }

    const Int entries_left = (ldf_ - nb_ * j - jb) * jb;

    for (Int i = 0; i < entries_left; ++i) {
      if (frontal_[index] != 0.0) {
        minoffD = std::min(minoffD, std::abs(frontal_[index]));
        maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
      }
      index++;
    }
  }

  data_.setExtremeEntries(minD, maxD, minoffD, maxoffD);
}

}  // namespace hipo