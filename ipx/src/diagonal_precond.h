// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_DIAGONAL_PRECOND_H_
#define IPX_DIAGONAL_PRECOND_H_

#include "linear_operator.h"
#include "model.h"
#include "sparse_matrix.h"

namespace ipx {

// DiagonalPrecond provides inverse operations with the diagonal matrix
//
//   diag(AI*W*AI')                                 (1)
//
// or, if AI = [As Ad I] is split into sparse and dense part, with
//
//   diag(As*Ws*As') + W[n+1:n+m] + Ad*Wd*Ad'.      (2)
//
// Here AI is the m-by-(n+m) matrix defined by the model, and W is a diagonal
// (weight) matrix that is provided by the user. If the form (2) is used, then
// the splitting into As and Ad is defined by the model, and the part of W
// corresponding to Ad (denoted Wd) must be invertible (not checked).

class DiagonalPrecond : public LinearOperator {
public:
    // Constructor stores a reference to the model. No data is copied. The model
    // must be valid as long as the preconditioner is used.
    explicit DiagonalPrecond(const Model& model);

    // Factorizes the preconditioner. W must either hold n+m entries, or be
    // NULL, in which case the first n entries are assumed 1.0 and the last
    // m entries are assumed 0.0. If precond_dense_cols is true, then (2)
    // becomes the preconditioner, otherwise (1).
    // On return info->errflag is 0 on succcess and IPX_ERROR_lapack_chol
    // if the LAPACK Cholesky factorization failed (used for (2) only).
    void Factorize(const double* W, bool precond_dense_cols, Info* info);

    // Returns computation time for calls to Apply() since last reset_time().
    double time() const;
    void reset_time();

private:
    void _Apply(const Vector& rhs, Vector& lhs, double* rhs_dot_lhs) override;

    const Model& model_;
    bool factorized_{false};    // preconditioner factorized?
    Vector diagonal_;           // diagonal of (sparse part of) normal matrix
    SparseMatrix Atdense_;      // dense columns of A, stored in CSR format
    Vector chol_factor_;        // Cholesky factor of pivotal matrix from SMW
    Vector work_;               // size ndense workspace for Apply()
    double time_{0.0};
};

}  // namespace ipx

#endif  // IPX_DIAGONAL_PRECOND_H_
