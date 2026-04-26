#include "perf_saxpy.h"

#include <vector>

#include "perf_fixtures.h"
#include "perf_harness.h"
#include "util/HVector.h"
#include "util/HighsRandom.h"

namespace highs_perf {

namespace {

void fill_hvector(HVector& v, HighsInt num_row, double density,
                  HighsRandom& rng) {
  v.clear();
  HighsInt count = 0;
  for (HighsInt i = 0; i < num_row; ++i) {
    if (rng.fraction() < density) {
      double x = 2.0 * rng.fraction() - 1.0;
      if (x == 0.0) x = 1e-3;
      v.array[i] = x;
      v.index[count++] = i;
    }
  }
  if (count == 0) {
    v.array[0] = 1.0;
    v.index[0] = 0;
    count = 1;
  }
  v.count = count;
}

}  // namespace

void bench_saxpy(Fixture& fx, std::vector<ResultRow>& out) {
  HighsRandom rng(0xDEADBEEF);

  struct Variant {
    const char* name;
    double a_density;
    double b_density;
  };
  const Variant variants[] = {
      {"a25_b1",  0.25, 0.01},
      {"a25_b5",  0.25, 0.05},
      {"a50_b25", 0.50, 0.25},
  };

  for (const auto& v : variants) {
    HVector a, b;
    a.setup(fx.num_row);
    b.setup(fx.num_row);
    fill_hvector(a, fx.num_row, v.a_density, rng);
    fill_hvector(b, fx.num_row, v.b_density, rng);

    // Snapshot of `a` so each repetition starts from the same input.
    std::vector<HighsInt> a_idx0 = a.index;
    std::vector<double> a_arr0 = a.array;
    HighsInt a_count0 = a.count;

    auto stats = time_kernel(
        [&]() {
          a.saxpy(0.7, &b);
          // Restore (excluded from timing).
          a.index = a_idx0;
          a.array = a_arr0;
          a.count = a_count0;
        },
        /*min_reps=*/100, /*min_wall_ns=*/100000000LL, /*max_reps=*/500000);
    out.push_back({fx.name, "saxpy", v.name, stats});
  }
}

}  // namespace highs_perf
