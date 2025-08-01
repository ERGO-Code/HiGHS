#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "ReturnValues.h"
#include "Swaps.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "util/HighsRandom.h"

namespace hipo {

// Dense Factorisation kernel

std::pair<Int, double> maxInCol(Int j, Int n, Int m, double* A, Int lda) {
  // Given the symemtric matrix A, of size nxn, accessed with leading dimension
  // lda, in upper triangular format, ignoring rows 0:j-1, find the maximum in
  // row/col m.

  double maxval = -1.0;
  Int r = -1;
  // scan column m, rows j:m-1
  for (Int row = j; row < m; ++row) {
    double val = std::abs(A[row + lda * m]);
    if (val > maxval) {
      maxval = val;
      r = row;
    }
  }
  // scan row m, columns m+1:n
  for (Int col = m + 1; col < n; ++col) {
    double val = std::abs(A[m + lda * col]);
    if (val > maxval) {
      maxval = val;
      r = col;
    }
  }
  return {r, maxval};
}

void staticReg(double& pivot, Int sign, double& regul) {
  // apply static regularisation

  double old_pivot = pivot;
  if (sign > 0)
    pivot += kDualStaticRegularisation;
  else
    pivot -= kPrimalStaticRegularisation;
  regul = pivot - old_pivot;
}

bool blockBunchKaufman(Int j, Int n, double* A, Int lda, Int* swaps, Int* sign,
                       double thresh, double* regul, DataCollector& data) {
  // Perform Bunch-Kaufman pivoting within a block of the supernode (see Schenk,
  // Gartner, ETNA 2006).
  // It works only for upper triangular A.
  // Return true if 2x2 pivot should be used, false otherwise.
  // Swap of columns may be performed.
  // Regularisation of pivot may be performed.

#if HIPO_TIMING_LEVEL >= 2
  Clock clock;
#endif
  bool flag_2x2 = false;

  // Find largest diagonal entry in the residual part of the block
  Int ind_max_diag = -1;
  double max_diag = -1.0;
  for (Int i = j; i < n; ++i) {
    double val = std::abs(A[i + lda * i]);
    if (val > max_diag) {
      max_diag = val;
      ind_max_diag = i;
    }
  }
  assert(ind_max_diag >= j);

  // put column with max pivot as first column
  swapCols('U', n, A, lda, j, ind_max_diag, swaps, sign, data);

  // Max in column j of diagonal block
  auto res = maxInCol(j, n, j, A, lda);
  double gamma_j = res.second;
  Int r = res.first;
  double Ajj = sign[j] > 0 ? A[j + lda * j] + kDualStaticRegularisation
                           : A[j + lda * j] - kPrimalStaticRegularisation;

  if (std::max(std::abs(Ajj), gamma_j) <= thresh || sign[j] * Ajj < 0 ||
      j == n - 1) {
    // Must accept current pivot
    double old_pivot = A[j + lda * j];
    staticReg(A[j + lda * j], sign[j], regul[j]);

    if (sign[j] * A[j + lda * j] < 0) {
      data.setWrongSign(A[j + lda * j]);
      // A[j + lda * j] = sign[j] * thresh;
    }

    if (std::max(std::abs(Ajj), gamma_j) < thresh) {
      // perturbe pivot
      A[j + lda * j] = sign[j] * thresh;
      data.countRegPiv();
    }
    regul[j] = A[j + lda * j] - old_pivot;
    data.setMaxReg(std::abs(regul[j]));

  } else {
    // Max in column r of diagonal block
    assert(r >= 0);
    res = maxInCol(j, n, r, A, lda);
    double gamma_r = res.second;
    double Arr = sign[r] > 0 ? A[r + lda * r] + kDualStaticRegularisation
                             : A[r + lda * r] - kPrimalStaticRegularisation;

    if ((std::abs(Ajj) >= kAlphaBK * gamma_j ||
         std::abs(Ajj) * gamma_r >= kAlphaBK * gamma_j * gamma_j)) {
      // Accept current pivot
      staticReg(A[j + lda * j], sign[j], regul[j]);

      if (sign[j] * A[j + lda * j] < 0) data.setWrongSign(A[j + lda * j]);

    } else if (std::abs(Arr) >= kAlphaBK * gamma_r) {
      // Use pivot r

      swapCols('U', n, A, lda, j, r, swaps, sign, data);
      staticReg(A[j + lda * j], sign[j], regul[j]);

      if (sign[j] * A[j + lda * j] < 0) data.setWrongSign(A[j + lda * j]);

    } else {
      // Use 2x2 pivot (j,r)

      swapCols('U', n, A, lda, j + 1, r, swaps, sign, data);
      flag_2x2 = true;
      staticReg(A[j + lda * j], sign[j], regul[j]);
      staticReg(A[j + 1 + lda * (j + 1)], sign[j + 1], regul[j + 1]);

      if (sign[j] * A[j + lda * j] < 0) data.setWrongSign(A[j + lda * j]);
      if (sign[j + 1] * A[j + 1 + lda * (j + 1)] < 0)
        data.setWrongSign(A[j + 1 + lda * (j + 1)]);
    }
  }

#if HIPO_TIMING_LEVEL >= 2
  data.sumTime(kTimeDenseFact_pivoting, clock.stop());
#endif
  return flag_2x2;
}

Int denseFactK(char uplo, Int n, double* A, Int lda, Int* pivot_sign,
               double thresh, double* regul, Int* swaps, double* pivot_2x2,
               DataCollector& data) {
  // ===========================================================================
  // Factorisation kernel
  // Matrix A is in format F
  // BLAS calls: dscal, dcopy, dger
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) return kRetInvalidInput;

  // quick return
  if (n == 0) return kRetOk;

#if HIPO_TIMING_LEVEL >= 2
  Clock clock;
#endif

  if (uplo == 'L') {
    assert(1 == 0);
  } else {
    if (!swaps || !pivot_2x2) return kRetInvalidInput;

    // initialise order of pivots
    for (Int i = 0; i < n; ++i) swaps[i] = i;

    // allocate space for copy of col(s)
    std::vector<double> temp(n - 1);
    std::vector<double> temp2(n - 1);

    Int step = 1;

    for (Int j = 0; j < n; j += step) {
      bool flag_2x2 = false;

#ifdef HIPO_PIVOTING
      flag_2x2 = blockBunchKaufman(j, n, A, lda, swaps, pivot_sign, thresh,
                                   regul, data);
#endif

      // cannot do 2x2 pivoting on last column
      assert(j < n - 1 || flag_2x2 == false);

      if (!flag_2x2) {
        // 1x1 pivots
        step = 1;

        // diagonal element
        double Ajj = A[j + lda * j];

        if (std::isnan(Ajj)) return kRetInvalidPivot;

        // add regularisation
        staticReg(Ajj, pivot_sign[j], regul[j]);

#ifndef HIPO_PIVOTING
        // add static regularisation
        staticReg(Ajj, pivot_sign[j], regul[j]);
        data.setMaxReg(std::abs(regul[j]));
#endif

        // save reciprocal of pivot
        A[j + lda * j] = 1.0 / Ajj;

        const Int M = n - j - 1;
        if (M > 0) {
          // make copy of row
          callAndTime_dcopy(M, &A[j + (j + 1) * lda], lda, temp.data(), 1,
                            data);

          // scale row j
          callAndTime_dscal(M, 1.0 / Ajj, &A[j + (j + 1) * lda], lda, data);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + (j + 1) * lda],
                           lda, &A[j + 1 + (j + 1) * lda], lda, data);
        }
      } else {
        // 2x2 pivots
        step = 2;

        // diagonal pivot elements
        const double d1 = A[j + lda * j];
        const double d2 = A[j + 1 + lda * (j + 1)];

        if (std::isnan(d1) || std::isnan(d2)) return kRetInvalidPivot;

        // off-diagonal pivot element
        const double offd = A[j + lda * (j + 1)];
        A[j + lda * (j + 1)] = 0.0;

        // compute coefficients of 2x2 inverse
        const double denom = d1 * d2 - offd * offd;
        const double i_d1 = d2 / denom;
        const double i_d2 = d1 / denom;
        const double i_off = -offd / denom;

        // save them in place of pivots
        A[j + lda * j] = i_d1;
        A[j + 1 + lda * (j + 1)] = i_d2;
        pivot_2x2[j] = i_off;

        const Int M = n - j - 2;
        if (M > 0) {
          double* r1 = &A[j + (j + 2) * lda];
          double* r2 = &A[j + 1 + (j + 2) * lda];

          // make a copy of first row
          callAndTime_dcopy(M, r1, lda, temp.data(), 1, data);

          // make a copy of second row
          callAndTime_dcopy(M, r2, lda, temp2.data(), 1, data);

          // solve rows j,j+1
          callAndTime_dscal(M, i_d1, r1, lda, data);
          callAndTime_daxpy(M, i_off, temp2.data(), 1, r1, lda, data);
          callAndTime_dscal(M, i_d2, r2, lda, data);
          callAndTime_daxpy(M, i_off, temp.data(), 1, r2, lda, data);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, r1, lda,
                           &A[j + 2 + (j + 2) * lda], lda, data);
          callAndTime_dger(M, M, -1.0, temp2.data(), 1, r2, lda,
                           &A[j + 2 + (j + 2) * lda], lda, data);
        }

        data.count2x2();
      }
    }
  }

#if HIPO_TIMING_LEVEL >= 2
  data.sumTime(kTimeDenseFact_kernel, clock.stop());
#endif

  return kRetOk;
}

}  // namespace hipo
