#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"
#include "ReturnValues.h"
#include "Swaps.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

// Factorisation with "hybrid formats".

Int denseFactFH(Int n, Int k, double* A, double* B, const Int* pivot_sign,
                double thresh, double* totalreg, Int* swaps, double* pivot_2x2,
                DataCollector& data, const FHoptions& options) {
  // ===========================================================================
  // Partial blocked factorisation
  // Matrix A is in format FH
  // Matrix B is in format FH
  // BLAS calls: dcopy, dscal, daxpy, dgemm, dtrsm
  // ===========================================================================
  // While processing a block of columns:
  // - D is the diagonal block
  // - R is the collection of blocks below D
  // - while using R to perform a right-looking update to a later block:
  //    * Rj is the portion of R being used
  //    * Pj is the portion of R being used, before being solved with the pivots
  //    * Qj is the block of columns being updated
  //   so that the update is Qj -= Rj * Pj^T
  // ===========================================================================

  HIPO_CLOCK_CREATE;

  if (n < 0 || k < 0 || !A || (k < n && !B)) return kRetInvalidInput;
  if (n == 0) return kRetOk;

  const Int nb = options.nb;
  const Int blocks_in_A = (k - 1) / nb + 1;

  std::vector<Int64> blocks_diag_start(blocks_in_A);
  getDiagStart(n, k, nb, blocks_in_A, blocks_diag_start);

  std::vector<double> buffer(n * nb);

  const Int B_size = n - k;
  const Int blocks_in_B = (B_size - 1) / nb + 1;

  // ===========================================================================
  // LOOP OVER BLOCKS
  // ===========================================================================
  for (Int block_id = 0; block_id < blocks_in_A; ++block_id) {
    HIPO_CLOCK_START(2);

    const Int this_block_col = std::min(nb, k - nb * block_id);
    const Int jb = this_block_col;

    const Int diag_block_entries = jb * jb;
    const Int rows_below = n - nb * block_id - jb;

    double* D = &A[blocks_diag_start[block_id]];
    double* R = &A[blocks_diag_start[block_id] + diag_block_entries];

    // ===========================================================================
    // FACTORISE DIAGONAL BLOCK
    // ===========================================================================
    double* this_block_regularisation = &totalreg[block_id * nb];
    std::vector<Int> this_block_pivot_sign(&pivot_sign[block_id * nb],
                                           &pivot_sign[block_id * nb] + jb);
    Int* this_block_swaps = &swaps[block_id * nb];
    double* this_block_pivot_2x2 = &pivot_2x2[block_id * nb];
    Int status = denseFactK('U', jb, D, jb, this_block_pivot_sign.data(),
                            thresh, this_block_regularisation, this_block_swaps,
                            this_block_pivot_2x2, data, options);
    if (status != 0) return status;

    if (options.pivoting) {
      applySwaps(this_block_swaps, rows_below, jb, R, data);
      permuteWithSwaps(this_block_regularisation, this_block_swaps, jb, true);
    }

    if (rows_below > 0) {
      // ===========================================================================
      // SOLVE COLUMNS
      // ===========================================================================
      // solve block R with D
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, rows_below, 1.0, D, jb, R, jb,
                        data);

      // make copy of partially solved columns
      callAndTime_dcopy(jb * rows_below, R, 1, buffer.data(), 1, data);

      // solve block R with pivots
      Int step = 1;
      for (Int col = 0; col < jb; col += step) {
        if (this_block_pivot_2x2[col] == 0.0) {
          // 1x1 pivots
          step = 1;
          callAndTime_dscal(rows_below, D[col + jb * col], &R[col], jb, data);
        } else {
          // 2x2 pivots
          step = 2;

          // columns affected
          double* c1 = &R[col];
          double* c2 = &R[col + 1];

          // inverse of 2x2 pivot
          double i_d1 = D[col + jb * col];
          double i_d2 = D[col + 1 + jb * (col + 1)];
          double i_off = this_block_pivot_2x2[col];

          // copy of original col1
          std::vector<double> c1_temp(rows_below);
          callAndTime_dcopy(rows_below, c1, jb, c1_temp.data(), 1, data);

          // solve col and col+1
          callAndTime_dscal(rows_below, i_d1, c1, jb, data);
          callAndTime_daxpy(rows_below, i_off, c2, jb, c1, jb, data);
          callAndTime_dscal(rows_below, i_d2, c2, jb, data);
          callAndTime_daxpy(rows_below, i_off, c1_temp.data(), 1, c2, jb, data);
        }
      }

      // ===========================================================================
      // UPDATE
      // ===========================================================================

      auto split_gemm = [&](Int num_row, Int num_col, const double* Rj,
                            const double* Pj, double* Qj) {
        // Qj -= Rj * Pj^T

        const bool do_split = options.parallel_node && num_col > nb / 2 &&
                              jb > nb / 2 &&
                              num_row >= kBlockParallelThreshold * nb;

        if (do_split) {
          auto do_gemm_block_rows = [=, &data](Int start_row, Int end_row) {
            const Int Rj_offset = jb * start_row;
            const Int Qj_offset = num_col * start_row;
            const Int row_count = end_row - start_row;

            const double* Rj_block = &Rj[Rj_offset];
            double* Qj_block = &Qj[Qj_offset];

            callAndTime_dgemm('T', 'N', num_col, row_count, jb, -1.0, Pj, jb,
                              Rj_block, jb, 1.0, Qj_block, num_col, data);
          };

          highs::parallel::for_each(0, num_row, do_gemm_block_rows, nb);
        } else {
          auto do_gemm_full = [=, &data]() {
            callAndTime_dgemm('T', 'N', num_col, num_row, jb, -1.0, Pj, jb, Rj,
                              jb, 1.0, Qj, num_col, data);
          };
          do_gemm_full();
        }
      };

      // ===========================================================================
      // UPDATE FRONTAL
      // ===========================================================================
      Int64 R_offset{};
      for (Int j = block_id + 1; j < blocks_in_A; ++j) {
        const Int col_block_j = std::min(nb, k - nb * j);
        const Int row_block_j = n - nb * j;

        const double* Pj = &buffer[R_offset];
        double* Qj = &A[blocks_diag_start[j]];
        const double* Rj = &R[R_offset];

        split_gemm(row_block_j, col_block_j, Rj, Pj, Qj);

        R_offset += jb * col_block_j;
      }

      // ===========================================================================
      // UPDATE SCHUR COMPLEMENT
      // ===========================================================================
      HIPO_CLOCK_START(2);
      if (k < n) {
        Int64 B_offset{};

        for (Int j = 0; j < blocks_in_B; ++j) {
          const Int row_block_j = B_size - nb * j;
          const Int col_block_j = std::min(nb, row_block_j);

          const double* Pj = &buffer[R_offset];
          double* Qj = &B[B_offset];
          const double* Rj = &R[R_offset];

          split_gemm(row_block_j, col_block_j, Rj, Pj, Qj);

          B_offset += row_block_j * col_block_j;
          R_offset += jb * col_block_j;
        }
      }

      HIPO_CLOCK_STOP(2, data, kTimeDenseFact_update);
    }
  }

  return kRetOk;
}

Int denseFactFP2FH(double* A, Int nrow, Int ncol, Int nb, DataCollector& data) {
  // ===========================================================================
  // Packed to Hybrid conversion
  // Matrix A on  input is in format FP
  // Matrix A on output is in format FH
  // BLAS calls: dcopy
  // ===========================================================================

  HIPO_CLOCK_CREATE;

  std::vector<double> buffer(nrow * nb);

  Int64 offset_for_read = 0;
  Int64 offset_for_write = 0;

  for (Int k = 0; k <= (ncol - 1) / nb; ++k) {
    const Int this_block_col = std::min(nb, ncol - k * nb);
    const Int this_block_row = nrow - k * nb;

    // Copy block into buffer
    callAndTime_dcopy(this_block_row * this_block_col, &A[offset_for_read], 1,
                      buffer.data(), 1, data);
    offset_for_read += this_block_row * this_block_col;

    // Copy columns back into A, row by row.
    // One call of dcopy_ for each row of the block of columns.
    for (Int i = 0; i < this_block_row; ++i) {
      callAndTime_dcopy(this_block_col, &buffer[i], this_block_row,
                        &A[offset_for_write], 1, data);
      offset_for_write += this_block_col;
    }
  }

  HIPO_CLOCK_STOP(2, data, kTimeDenseFact_convert);

  return kRetOk;
}

}  // namespace hipo
