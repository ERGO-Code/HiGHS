// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_LAPACK_H_
#define IPX_LAPACK_H_

#include "ipx_internal.h"

namespace ipx {

// Wrapper functions around LAPACK that convert between integer types and throw
// std::overflow_error if an integer overflow occurs.

Int Lapack_dpotrf(char uplo, Int n, double* a, Int lda);

Int Lapack_dpotrs(char uplo, Int n, Int nrhs, const double* a, Int lda,
                  double* b, Int ldb);

double Lapack_dtrcon(char norm, char uplo, char diag, Int n, const double* a,
                     Int lda);

}  // namespace ipx

#endif  // IPX_LAPACK_H_
