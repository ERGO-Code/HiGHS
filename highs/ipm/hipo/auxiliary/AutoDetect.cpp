#include "AutoDetect.h"

#include <stdint.h>

// declare Metis with 32-bit integers for the test
#define IDXTYPEWIDTH 32
#include "metis.h"

// Weird tricks to detect the integer type width used by BLAS and Metis.
// These are technically undefined behaviour, because they rely on using a
// function declaration that involves a certain integer type, while the actual
// implementation may use a different one. Their behaviour may depend on the
// endianness of the CPU (?).

namespace hipo {

extern "C" {
int64_t cblas_isamax(const int64_t N, const float* X, const int64_t incX);
}
IntegerModel getBlasIntegerModel() {
  // Test if BLAS uses 32-bit integers (LP64) or 64-bit integers (ILP64) at
  // runtime. Inspired by libblastrampoline's autodetection.c

  // Even though isamax is declared to use 64-bit integers, it may actually
  // use 32-bit integers. If a negative number is passed as first argument,
  // isamax returns 0. If the correct value of 3 is passed, it returns 2
  // instead.

  static IntegerModel blas_model = IntegerModel::not_set;

  if (blas_model == IntegerModel::not_set) {
    // This is a very negative 64-bit number, but it is just 3 if only the lower
    // 32 bits are used.
    const int64_t n = 0xffffffff00000003;

    const float X[3] = {1.0f, 2.0f, 3.0f};

    const int64_t incx = 1;
    int64_t max_idx = cblas_isamax(n, X, incx);

    // Ignore potential upper 32 bits of the result
    max_idx = max_idx & 0xffffffff;

    if (max_idx == 0) {
      // isamax read negative n and returned 0, so it's using ilp64
      blas_model = IntegerModel::ilp64;

    } else if (max_idx == 2) {
      // isamax read correct n and returned 2, so it's using lp64
      blas_model = IntegerModel::lp64;

    } else {
      // something went wrong
      blas_model = IntegerModel::unknown;
    }
  }

  return blas_model;
}

std::string getIntegerModelString(IntegerModel i) {
  switch (i) {
    case IntegerModel::not_set:
      return "Not set";

    case IntegerModel::unknown:
      return "Unknown";

    case IntegerModel::lp64:
      return "LP64";

    case IntegerModel::ilp64:
      return "ILP64";
  }
}

IntegerModel getMetisIntegerModel() {
  static IntegerModel metis_model = IntegerModel::not_set;

  if (metis_model == IntegerModel::not_set) {
    idx_t options[METIS_NOPTIONS];
    for (int i = 0; i < METIS_NOPTIONS; ++i) options[i] = -1;

    // if 32 bits are used, this sets iptype to 2, otherwise it sets objtype to
    // a wrong value
    options[METIS_OPTION_IPTYPE] = 2;

    idx_t n = 3;
    idx_t ptr[4] = {0, 2, 4, 6};
    idx_t rows[6] = {1, 2, 0, 2, 0, 1};
    idx_t perm[3], iperm[3];

    idx_t status = METIS_NodeND(&n, ptr, rows, NULL, options, perm, iperm);

    if (status == METIS_OK) {
      if (perm[0] != 0 || perm[1] != 1 || perm[2] != 2)
        metis_model = IntegerModel::unknown;
      else
        metis_model = IntegerModel::lp64;
    } else if (status == METIS_ERROR_INPUT) {
      metis_model = IntegerModel::ilp64;
    } else {
      metis_model = IntegerModel::unknown;
    }
  }

  return metis_model;
}
}  // namespace hipo