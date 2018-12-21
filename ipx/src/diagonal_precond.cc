// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "diagonal_precond.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "lapack.h"
#include "timer.h"

namespace ipx {

DiagonalPrecond::DiagonalPrecond(const Model& model) : model_(model) {
    const Int m = model_.rows();
    diagonal_.resize(m);
}

void DiagonalPrecond::Factorize(const double* W, bool precond_dense_cols,
                               Info* info) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Int num_dense = model_.num_dense_cols();
    const SparseMatrix& AI = model_.AI();

    factorized_ = false;

    // Build diagonal of normal matrix, excluding dense columns if
    // precond_dense_cols is true.
    if (W) {
        for (Int i = 0; i < m; i++)
            diagonal_[i] = W[n+i];
        for (Int j = 0; j < n; j++) {
            if (precond_dense_cols && model_.IsDenseColumn(j))
                continue;
            double w = W[j];
            for (Int p = AI.begin(j); p < AI.end(j); p++)
                diagonal_[AI.index(p)] += AI.value(p) * w * AI.value(p);
        }
    } else {
        diagonal_ = 0.0;        // rightmost m columns have weight zero
        for (Int j = 0; j < n; j++) {
            if (precond_dense_cols && model_.IsDenseColumn(j))
                continue;
            for (Int p = AI.begin(j); p < AI.end(j); p++)
                diagonal_[AI.index(p)] += AI.value(p) * AI.value(p);
        }
    }

    if (precond_dense_cols && num_dense > 0) {
        // Compute a representation of the inverse of the preconditioner (2)
        // from the Sherman-Morrison-Woodbury (SMW) formula. Let P denote the
        // preconditioner and E = P - Ad*Wd*Ad'. Then E is diagonal matrix and
        //
        //   inv(P) = inv(E) - inv(E)*Ad * inv(S) * Ad'*inv(E),
        //
        // where S = inv(Wd)+Ad'*inv(E)*Ad is a Schur complement of dimension
        // ndense-by-ndense. This representation requires Wd to be invertible.

        // Build dense part of A rowwise.
        std::vector<Int> dense_cols;
        for (Int j = 0; j < n; j++)
            if (model_.IsDenseColumn(j))
                dense_cols.push_back(j);
        assert(dense_cols.size() == num_dense);
        Atdense_ = Transpose(CopyColumns(AI, dense_cols));
        assert(Atdense_.rows() == num_dense);

        // Build Schur complement for SMW formula.
        chol_factor_.resize(num_dense*num_dense);
        chol_factor_ = 0.0;
        Int k = 0;
        for (Int j : dense_cols) {
            // Compute column k of S:
            // S[:,k] = inv(Wd)[:,k] + Ad'*inv(E)*Ad[:,k]
            double* col = &chol_factor_[k*num_dense];
            for (Int p = AI.begin(j); p < AI.end(j); p++) {
                Int i = AI.index(p);
                double alpha = AI.value(p) / diagonal_[i];
                assert(std::isfinite(alpha));
                for (Int pp = Atdense_.begin(i); pp < Atdense_.end(i); pp++)
                    col[Atdense_.index(pp)] += alpha * Atdense_.value(pp);
            }
            double w = W ? W[j] : 1.0;
            col[k] += 1.0 / w;  // add to diagonal
            k++;
        }

        // Cholesky factorization of Schur complement.
        Int err = Lapack_dpotrf('L', num_dense, &chol_factor_[0], num_dense);
        if (err) {
            info->errflag = IPX_ERROR_lapack_chol;
            return;
        }
        // double rcondest = Lapack_dtrcon('1', 'L', 'N', num_dense,
        //                                 &chol_factor_[0], num_dense);
        // if (rcondest < 1e-6) {
        //     std::cout << " Cholesky factor condest = "
        //                    << sci2(1./rcondest) << '\n';
        // }

        // Allocate workspace required to apply preconditioner.
        work_.resize(num_dense);
    } else {
        // At least we have to resize Atdense to zero rows, so that Apply()
        // knows that no dense column preconditioning is used. We also free
        // memory that might have been allocated in a previous factorization.
        Atdense_.clear();
        chol_factor_.resize(0);
        work_.resize(0);
    }
    factorized_ = true;
}

double DiagonalPrecond::time() const {
    return time_;
}

void DiagonalPrecond::reset_time() {
    time_ = 0.0;
}

void DiagonalPrecond::_Apply(const Vector& rhs, Vector& lhs,
                              double* rhs_dot_lhs) {
    const Int m = model_.rows();
    const Int ndense = Atdense_.rows();
    double rldot = 0.0;
    Timer timer;

    assert(factorized_);
    assert(lhs.size() == m);
    assert(rhs.size() == m);
    assert(work_.size() >= ndense);

    if (ndense > 0) {
        // Compute Ad'*inv(E)*rhs in workspace.
        work_ = 0.0;
        for (Int i = 0; i < m; i++)
            ScatterColumn(Atdense_, i, rhs[i]/diagonal_[i], work_);

        // Solve with Cholesky factorization of Schur complement.
        Int err = Lapack_dpotrs('L', ndense, 1, &chol_factor_[0], ndense,
                                &work_[0], ndense);
        assert(err == 0);

        // Compute lhs = inv(E)*rhs - inv(E)*Ad*work.
        for (Int i = 0; i < m; i++) {
            double d = rhs[i] - DotColumn(Atdense_, i, work_);
            lhs[i] = d / diagonal_[i];
            rldot += lhs[i] * rhs[i];
        }
    } else {
        for (Int i = 0; i < m; i++) {
            lhs[i] = rhs[i] / diagonal_[i];
            rldot += lhs[i] * rhs[i];
        }
    }
    if (rhs_dot_lhs)
        *rhs_dot_lhs = rldot;
    time_ += timer.Elapsed();
}

}  // namespace ipx
