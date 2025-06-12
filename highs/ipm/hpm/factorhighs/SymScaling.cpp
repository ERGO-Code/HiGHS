
#include "SymScaling.h"

#include <limits>

#include "ipm/hpm/auxiliary/VectorOperations.h"

namespace hipo {

void product(const std::vector<double>& x, std::vector<double>& y,
             const std::vector<Int>& ptr, const std::vector<Int>& rows,
             const std::vector<double>& N) {
  // Multiply by matrix E, i.e. matrix A with all entries equal to one, in lower
  // triangular form, and sum component-wise product of N with x.
  // E * x + N .* x= y

  Int n = x.size();

  y.assign(n, 0.0);

  // multiply by E
  for (Int col = 0; col < n; ++col) {
    for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
      Int row = rows[el];
      y[row] += x[col];
      if (row != col) y[col] += x[row];
    }
  }

  // multiply by N
  for (Int i = 0; i < n; ++i) y[i] += N[i] * x[i];
}
void CG_for_CR_scaling(const std::vector<double>& b, std::vector<double>& x,
                       const std::vector<double>& N,
                       const std::vector<Int>& ptr,
                       const std::vector<Int>& rows) {
  Int n = N.size();

  // initial residual
  std::vector<double> r = b;

  // initial approximation
  x.assign(n, 0.0);

  // direction
  std::vector<double> p = r;

  Int iter{};
  std::vector<double> Ap(n);

  while (iter < 100) {
    product(p, Ap, ptr, rows, N);

    double norm_r = dotProd(r, r);
    double alpha = norm_r / dotProd(p, Ap);

    // x = x + alpha * p
    vectorAdd(x, p, alpha);

    // r = r - alpha * Ap;
    vectorAdd(r, Ap, -alpha);

    // exit test
    if (norm2(r) / norm2(b) < 1e-6) break;

    double beta = dotProd(r, r) / norm_r;

    // p = r + beta * p
    vectorAdd(p, r, 1.0, beta);

    ++iter;
  }
}
void CurtisReidScalingSym(const std::vector<Int>& ptr,
                          const std::vector<Int>& rows,
                          const std::vector<double>& val,
                          std::vector<double>& colscale) {
  // Takes as input the CSC matrix A.
  // Computes symmetric Curtis-Reid scaling of the matrix, using powers of 2.

  Int n = ptr.size() - 1;

  colscale.assign(n, 0.0);

  // rhs for CG
  std::vector<double> logsumcol(n, 0.0);

  // number of entries in each column
  std::vector<double> col_entries(n, 0.0);

  // log A_ij
  for (Int col = 0; col < n; ++col) {
    for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
      Int row = rows[el];
      if (val[el] != 0.0) {
        double temp = log2(std::abs(val[el]));
        logsumcol[col] += temp;
        col_entries[col] += 1.0;

        // only lower triangle is used, so add components corresponding to the
        // upper triangle
        if (col != row) {
          logsumcol[row] += temp;
          col_entries[row] += 1.0;
        }
      }
    }
  }

  // solve linear system with CG
  std::vector<double> exponents(n);
  CG_for_CR_scaling(logsumcol, exponents, col_entries, ptr, rows);

  // compute scaling using exact powers of 2
  for (Int i = 0; i < n; ++i)
    colscale[i] = std::ldexp(1.0, std::round(-exponents[i]));
}

void RuizScalingSym(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                    const std::vector<double>& val,
                    std::vector<double>& colscale) {
  // Compute symmetric Ruiz scaling, i.e. simultaneous row and column iterative
  // scaling, using powers of 2.

  const Int n = ptr.size() - 1;

  colscale.assign(n, 1.0);

  for (Int iter = 0; iter < 10; ++iter) {
    // inf norm of columns
    std::vector<double> col_norms(n);
    for (Int col = 0; col < n; ++col) {
      for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
        Int row = rows[el];
        double v = std::abs(val[el]);
        v *= colscale[col] * colscale[row];

        if (v != 0.0) {
          col_norms[col] = std::max(col_norms[col], v);
          if (row != col) col_norms[row] = std::max(col_norms[row], v);
        }
      }
    }

    double max_error{};
    for (double d : col_norms) max_error = std::max(max_error, std::abs(1 - d));

    if (max_error < 1e-3) break;

    for (Int i = 0; i < n; ++i) colscale[i] *= 1.0 / sqrt(col_norms[i]);
  }

  // round scaling to a power of 2
  for (Int i = 0; i < n; ++i) {
    int exp;  // needs to remain int
    std::frexp(colscale[i], &exp);
    colscale[i] = std::ldexp(1.0, exp);
  }
}

void JacekScalingSym(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                     const std::vector<double>& val,
                     std::vector<double>& colscale) {
  // Compute scaling symilar to Ruiz, but with sqrt(max*min) instead of
  // sqrt(max), rounded to a power of 2.

  const Int n = ptr.size() - 1;

  colscale.assign(n, 1.0);

  for (Int iter = 0; iter < 10; ++iter) {
    // inf norm of columns
    std::vector<double> col_norms(n);
    std::vector<double> col_min(n, std::numeric_limits<double>::max());
    for (Int col = 0; col < n; ++col) {
      for (Int el = ptr[col]; el < ptr[col + 1]; ++el) {
        Int row = rows[el];
        double v = std::abs(val[el]);
        v *= colscale[col] * colscale[row];

        if (v != 0.0) {
          col_norms[col] = std::max(col_norms[col], v);
          col_min[col] = std::min(col_min[col], v);
          if (row != col) {
            col_norms[row] = std::max(col_norms[row], v);
            col_min[row] = std::min(col_min[row], v);
          }
        }
      }
    }

    double max_error{};
    for (double d : col_norms) max_error = std::max(max_error, std::abs(1 - d));

    if (max_error < 1e-3) break;

    for (Int i = 0; i < n; ++i)
      colscale[i] *= 1.0 / sqrt(col_norms[i] * col_min[i]);
  }

  // round scaling to a power of 2
  for (Int i = 0; i < n; ++i) {
    int exp;  // needs to remain int
    std::frexp(colscale[i], &exp);
    colscale[i] = std::ldexp(1.0, exp);
  }
}

}  // namespace hipo