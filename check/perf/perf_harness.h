#ifndef HIGHS_PERF_HARNESS_H_
#define HIGHS_PERF_HARNESS_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace highs_perf {

struct PerfStats {
  std::int64_t reps = 0;
  std::int64_t min_ns = 0;
  std::int64_t median_ns = 0;
  double mean_ns = 0.0;
  double stddev_ns = 0.0;
  double total_ns = 0.0;
};

template <typename Fn>
inline PerfStats time_kernel(Fn&& fn,
                             std::int64_t min_reps = 50,
                             std::int64_t min_wall_ns = 100000000LL,
                             std::int64_t max_reps = 1000000LL) {
  using clock = std::chrono::steady_clock;
  std::vector<std::int64_t> samples;
  samples.reserve(static_cast<std::size_t>(min_reps));
  std::int64_t total_ns = 0;
  std::int64_t reps = 0;
  while (reps < max_reps && (reps < min_reps || total_ns < min_wall_ns)) {
    auto t0 = clock::now();
    fn();
    auto t1 = clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)
                  .count();
    samples.push_back(dt);
    total_ns += dt;
    ++reps;
  }
  PerfStats s;
  s.reps = reps;
  s.total_ns = static_cast<double>(total_ns);
  if (samples.empty()) return s;
  std::sort(samples.begin(), samples.end());
  s.min_ns = samples.front();
  s.median_ns = samples[samples.size() / 2];
  double mean = 0.0;
  for (auto v : samples) mean += static_cast<double>(v);
  mean /= static_cast<double>(samples.size());
  double sq = 0.0;
  for (auto v : samples) {
    double d = static_cast<double>(v) - mean;
    sq += d * d;
  }
  s.mean_ns = mean;
  s.stddev_ns = std::sqrt(sq / static_cast<double>(samples.size()));
  return s;
}

struct ResultRow {
  std::string instance;
  std::string kernel;
  std::string variant;
  PerfStats stats;
};

inline void print_header() {
  std::printf(
      "%-14s %-18s %-14s %8s %12s %12s %12s %12s\n", "instance", "kernel",
      "variant", "reps", "min_ns", "median_ns", "mean_ns", "stddev_ns");
}

inline void print_row(const ResultRow& r) {
  std::printf("%-14s %-18s %-14s %8lld %12lld %12lld %12.0f %12.0f\n",
              r.instance.c_str(), r.kernel.c_str(), r.variant.c_str(),
              static_cast<long long>(r.stats.reps),
              static_cast<long long>(r.stats.min_ns),
              static_cast<long long>(r.stats.median_ns), r.stats.mean_ns,
              r.stats.stddev_ns);
}

inline void write_csv(const std::string& path,
                      const std::vector<ResultRow>& rows) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) return;
  std::fprintf(f,
               "instance,kernel,variant,reps,min_ns,median_ns,mean_ns,stddev_"
               "ns\n");
  for (const auto& r : rows) {
    std::fprintf(f, "%s,%s,%s,%lld,%lld,%lld,%.1f,%.1f\n", r.instance.c_str(),
                 r.kernel.c_str(), r.variant.c_str(),
                 static_cast<long long>(r.stats.reps),
                 static_cast<long long>(r.stats.min_ns),
                 static_cast<long long>(r.stats.median_ns), r.stats.mean_ns,
                 r.stats.stddev_ns);
  }
  std::fclose(f);
}

}  // namespace highs_perf

#endif
