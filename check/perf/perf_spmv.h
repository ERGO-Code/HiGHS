#ifndef HIGHS_PERF_SPMV_H_
#define HIGHS_PERF_SPMV_H_

#include <vector>

#include "perf_harness.h"

namespace highs_perf {
struct Fixture;
void bench_spmv(Fixture& fx, std::vector<ResultRow>& out);
}  // namespace highs_perf

#endif
