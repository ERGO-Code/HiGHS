#ifndef FACTORHIGHS_DENSE_FACT_H
#define FACTORHIGHS_DENSE_FACT_H

#include "ipm/hpm/auxiliary/IntConfig.h"

namespace hipo {

/*
  Names:
  denseFact:
  - K : factorisation kernel for diagonal blocks
  - F : blocked factorisation in full format
  - FP: blocked factorisation in format FP
  - FH : blocked factorisation in "hybrid formats"

  Formats used:
  - F : Full format
  - P : lower Packed format
  - FP: lower Packed format with Full diagonal blocks
  - H : lower-blocked-Hybrid format
  - FH: lower-blocked-Hybrid format with Full diagonal blocks

  F, P do not use blocks. FP, H, FH use blocks.
  Blocks are always blocks of columns.
  F, P store by columns.
  FP stores by columns within the blocks. H, FH store by rows within the blocks.
  See report for more details.
*/

// dense factorisation kernel
Int denseFactK(char uplo, Int n, double* A, Int lda, Int* pivot_sign,
               double thresh, double* regul, Int* swaps, double* pivot_2x2,
               Int sn, Int bl, double max_in_R = -1);

// dense partial factorisation, in full format
Int denseFactF(Int n, Int k, Int nb, double* A, Int lda, double* B, Int ldb,
               const Int* pivot_sign, double thresh, double* regul, Int sn);

// dense partial factorisation, in packed format with full diagonal blocks
Int denseFactFP(Int n, Int k, Int nb, double* A, double* B,
                const Int* pivot_sign, double thresh, double* regul, Int sn);

// dense partial factorisation, in "hybrid formats"
Int denseFactFH(char format, Int n, Int k, Int nb, double* A, double* B,
                const Int* pivot_sign, double thresh, double* regul, Int* swaps,
                double* pivot_2x2, Int sn, bool parnode);

// function to convert A from lower packed, to lower-blocked-hybrid format
Int denseFactFP2FH(double* A, Int nrow, Int ncol, Int nb);

}  // namespace hipo

#endif