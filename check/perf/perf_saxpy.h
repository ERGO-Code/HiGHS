#ifndef HIGHS_PERF_SAXPY_H_
#define HIGHS_PERF_SAXPY_H_

#include <vector>

#include "perf_harness.h"

namespace highs_perf {
struct Fixture;
void bench_saxpy(Fixture& fx, std::vector<ResultRow>& out);
}  // namespace highs_perf

#endif
