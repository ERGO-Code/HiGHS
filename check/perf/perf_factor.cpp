#include "perf_factor.h"

#include <cmath>
#include <cstdio>

#include "perf_fixtures.h"
#include "perf_harness.h"
#include "util/HFactor.h"
#include "util/HVector.h"
#include "util/HighsRandom.h"

namespace highs_perf {

namespace {

// Populate rhs as a hyper-sparse unit vector e_k for a deterministic k.
void make_rhs_unit(HVector& rhs, HighsInt num_row, HighsRandom& rng) {
  rhs.clear();
  HighsInt k = rng.integer(num_row);
  rhs.array[k] = 1.0;
  rhs.index[0] = k;
  rhs.count = 1;
}

// Populate rhs with approximately `density` fraction of nonzeros, random
// values in (-1, 1).
void make_rhs_random(HVector& rhs, HighsInt num_row, double density,
                     HighsRandom& rng) {
  rhs.clear();
  HighsInt count = 0;
  for (HighsInt i = 0; i < num_row; ++i) {
    if (rng.fraction() < density) {
      double v = 2.0 * rng.fraction() - 1.0;
      if (v == 0.0) v = 1e-3;
      rhs.array[i] = v;
      rhs.index[count++] = i;
    }
  }
  if (count == 0) {
    rhs.array[0] = 1.0;
    rhs.index[0] = 0;
    count = 1;
  }
  rhs.count = count;
}

// Save / restore HVector contents around the timed call so each iteration
// sees the same input. We snapshot only what `clear()` resets.
struct HVectorSnapshot {
  std::vector<HighsInt> index;
  std::vector<double> array;
  HighsInt count;
};

HVectorSnapshot snapshot(const HVector& v) {
  return {v.index, v.array, v.count};
}

void restore(HVector& v, const HVectorSnapshot& s) {
  v.index = s.index;
  v.array = s.array;
  v.count = s.count;
}

}  // namespace

void bench_factor(Fixture& fx, std::vector<ResultRow>& out) {
  if (!fx.factor_ok) return;

  // ---- HFactor::build ----------------------------------------------------
  // Each rep needs a fresh HFactor because build() mutates internal L/U state
  // and a second build on the same factor isn't equivalent to a from-scratch
  // factorisation.
  {
    auto stats = time_kernel(
        [&]() {
          HFactor f;
          // Local copy of basic_set because HFactor::setup keeps a non-const
          // reference and may reorder it during build.
          std::vector<HighsInt> bs = fx.basic_set;
          f.setup(fx.lp.a_matrix_, bs);
          f.build();
        },
        /*min_reps=*/10, /*min_wall_ns=*/200000000LL,
        /*max_reps=*/2000);
    out.push_back({fx.name, "build", "-", stats});
  }

  // ---- HFactor::ftranCall / btranCall at 3 densities --------------------
  HighsRandom rng(0xC0FFEE);
  HVector rhs;
  rhs.setup(fx.num_row);

  struct Variant {
    const char* name;
    double density;       // approximate input density
    double expected;      // expected_density passed to ftran/btran
  };
  const Variant variants[] = {
      {"hyper_sparse", 0.0, 0.001},  // unit vector
      {"sparse",       0.05, 0.05},
      {"dense",        1.0, 0.5},
  };

  // FTRAN
  for (const auto& v : variants) {
    if (v.density == 0.0)
      make_rhs_unit(rhs, fx.num_row, rng);
    else
      make_rhs_random(rhs, fx.num_row, v.density, rng);
    auto snap = snapshot(rhs);
    auto stats = time_kernel(
        [&]() {
          fx.factor.ftranCall(rhs, v.expected);
          // Reset rhs so the next iteration starts from the same input.
          // Restore is excluded from timing because steady_clock::now() bracket
          // captures only the call itself.
          restore(rhs, snap);
        },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL,
        /*max_reps=*/200000);
    out.push_back({fx.name, "ftranCall", v.name, stats});
  }

  // BTRAN
  for (const auto& v : variants) {
    if (v.density == 0.0)
      make_rhs_unit(rhs, fx.num_row, rng);
    else
      make_rhs_random(rhs, fx.num_row, v.density, rng);
    auto snap = snapshot(rhs);
    auto stats = time_kernel(
        [&]() {
          fx.factor.btranCall(rhs, v.expected);
          restore(rhs, snap);
        },
        /*min_reps=*/50, /*min_wall_ns=*/100000000LL,
        /*max_reps=*/200000);
    out.push_back({fx.name, "btranCall", v.name, stats});
  }
}

}  // namespace highs_perf
