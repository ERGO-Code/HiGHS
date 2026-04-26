#include "perf_spmv.h"

#include <vector>

#include "perf_fixtures.h"
#include "perf_harness.h"
#include "util/HVector.h"
#include "util/HighsRandom.h"
#include "util/HighsSparseMatrix.h"

namespace highs_perf {

namespace {

void make_dense_random(std::vector<double>& v, HighsRandom& rng) {
  for (auto& x : v) x = 2.0 * rng.fraction() - 1.0;
}

void make_sparse_column_hvector(HVector& col, HighsInt num_row, double density,
                                HighsRandom& rng) {
  col.clear();
  HighsInt count = 0;
  for (HighsInt i = 0; i < num_row; ++i) {
    if (rng.fraction() < density) {
      double v = 2.0 * rng.fraction() - 1.0;
      if (v == 0.0) v = 1e-3;
      col.array[i] = v;
      col.index[count++] = i;
    }
  }
  if (count == 0) {
    col.array[0] = 1.0;
    col.index[0] = 0;
    count = 1;
  }
  col.count = count;
}

}  // namespace

void bench_spmv(Fixture& fx, std::vector<ResultRow>& out) {
  HighsRandom rng(0xBEEF);

  // ---- HighsSparseMatrix::product ---------------------------------------
  // Colwise product: result of size num_row, input of size num_col.
  {
    std::vector<double> in(fx.num_col), result(fx.num_row);
    make_dense_random(in, rng);
    auto stats = time_kernel(
        [&]() { fx.matrix_colwise.product(result, in); },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL, /*max_reps=*/200000);
    out.push_back({fx.name, "product", "colwise", stats});
  }

  // ---- HighsSparseMatrix::productTranspose ------------------------------
  // Colwise transpose-product: result of size num_col, input of size num_row.
  {
    std::vector<double> in(fx.num_row), result(fx.num_col);
    make_dense_random(in, rng);
    auto stats = time_kernel(
        [&]() { fx.matrix_colwise.productTranspose(result, in); },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL, /*max_reps=*/200000);
    out.push_back({fx.name, "productTranspose", "colwise", stats});
  }

  // ---- HighsSparseMatrix::priceByColumn ---------------------------------
  // Two variants of input column density.
  struct Variant {
    const char* name;
    double density;
  };
  const Variant variants[] = {{"sparse_5pct", 0.05}, {"dense", 1.0}};

  for (const auto& v : variants) {
    HVector column;
    column.setup(fx.num_row);
    make_sparse_column_hvector(column, fx.num_row, v.density, rng);

    HVector result;
    result.setup(fx.num_col);

    auto stats = time_kernel(
        [&]() {
          // priceByColumn writes into result.array[iCol] for nonzero entries
          // and accumulates result.index[result.count++] = iCol. Reset between
          // calls so each rep does the same work.
          for (HighsInt k = 0; k < result.count; ++k)
            result.array[result.index[k]] = 0.0;
          result.count = 0;
          fx.matrix_colwise.priceByColumn(/*quad=*/false, result, column,
                                          kDebugReportOff);
        },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL, /*max_reps=*/200000);
    out.push_back({fx.name, "priceByColumn", v.name, stats});
  }

  // ---- HighsSparseMatrix::priceByRow ------------------------------------
  for (const auto& v : variants) {
    HVector column;
    column.setup(fx.num_row);
    make_sparse_column_hvector(column, fx.num_row, v.density, rng);

    HVector result;
    result.setup(fx.num_col);

    auto stats = time_kernel(
        [&]() {
          for (HighsInt k = 0; k < result.count; ++k)
            result.array[result.index[k]] = 0.0;
          result.count = 0;
          fx.matrix_rowwise.priceByRow(/*quad=*/false, result, column,
                                       kDebugReportOff);
        },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL, /*max_reps=*/200000);
    out.push_back({fx.name, "priceByRow", v.name, stats});
  }
}

}  // namespace highs_perf
