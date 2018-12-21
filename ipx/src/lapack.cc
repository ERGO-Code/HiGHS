// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "lapack.h"
#include <stdexcept>
#include <vector>

namespace ipx {

#ifdef BLAS64
#include <cstdint>
#define BLASINT int64_t
#else
#define BLASINT int
#endif

extern"C" void dpotrf_(const char *, const BLASINT *, double *, const BLASINT *,
                       BLASINT *);
extern"C" void dpotrs_(const char *, const BLASINT *, const BLASINT *, const
                       double *, const BLASINT *, double *, const BLASINT *,
                       BLASINT *);
extern"C" void dtrcon_(const char *, const char *, const char *, const
                       BLASINT *, const double *, const BLASINT *, double *,
                       double *, BLASINT *, BLASINT *);

Int Lapack_dpotrf(char uplo, Int n, double* a, Int lda) {
    if (n == 0)
        return 0;
    BLASINT N = n;
    BLASINT LDA = lda;
    BLASINT INFO = 0;
    if (N != n || LDA != lda)
        throw std::overflow_error("BLAS int overflow");
    dpotrf_(&uplo, &N, a, &LDA, &INFO);
    return static_cast<Int>(INFO);
}

Int Lapack_dpotrs(char uplo, Int n, Int nrhs, const double* a, Int lda,
                  double* b, Int ldb) {
    if (n == 0)
        return 0;
    BLASINT N = n;
    BLASINT NRHS = nrhs;
    BLASINT LDA = lda;
    BLASINT LDB = ldb;
    if (N != n || NRHS != nrhs || LDA != lda || LDB != ldb)
        throw std::overflow_error("BLAS int overflow");
    BLASINT INFO = 0;
    dpotrs_(&uplo, &N, &NRHS, a, &LDA, b, &LDB, &INFO);
    return static_cast<Int>(INFO);
}

double Lapack_dtrcon(char norm, char uplo, char diag, Int n, const double* a,
                     Int lda) {
    if (n == 0)
        return 0.0;
    BLASINT N = n;
    BLASINT LDA = lda;
    if (N != n || LDA != lda)
        throw std::overflow_error("BLAS int overflow");
    BLASINT INFO = 0;
    double rcond = 0.0;
    std::vector<double> work(3*n);
    std::vector<BLASINT> iwork(n);
    dtrcon_(&norm, &uplo, &diag, &N, a, &LDA, &rcond, work.data(), iwork.data(),
            &INFO);
    if (INFO != 0)
        throw std::logic_error("invalid input to dtrcon");
    return rcond;
}

}  // namespace ipx
