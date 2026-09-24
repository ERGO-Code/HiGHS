#include <features.h>

#ifndef __GLIBC__
#error "The solve heap audit requires Linux/glibc"
#endif

#include <execinfo.h>
#include <gnu/libc-version.h>
#include <malloc.h>
#include <pthread.h>

#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <thread>

#include "HConfig.h"
#include "Highs.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/api/options.hpp"

#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif

namespace solve_heap_audit {

enum HeapField {
  kMallocCalls,
  kCallocCalls,
  kReallocCalls,
  kAlignedCalls,
  kFreeCalls,
  kFailedCalls,
  kRequestedBytes,
  kMaxRequestBytes,
  kFreeUsableBytes,
  kHeapFieldCount
};

const std::size_t kPhaseCount =
    static_cast<std::size_t>(hipdlp::algorithm::SolvePhase::kCount);
const std::size_t kPreamble = kPhaseCount;
const std::size_t kSelfTest = kPhaseCount + 1;
const std::size_t kTeardown = kPhaseCount + 2;
const std::size_t kBucketCount = kPhaseCount + 3;
static_assert(kBucketCount < 32, "Heap phase mask is too small");
static_assert(ATOMIC_LLONG_LOCK_FREE == 2,
              "Heap auditing requires lock-free 64-bit atomics");

struct HeapSample {
  std::uint64_t values[kHeapFieldCount]{};

  std::uint64_t allocations() const {
    return values[kMallocCalls] + values[kCallocCalls] + values[kReallocCalls] +
           values[kAlignedCalls];
  }

  void add(const HeapSample& other) {
    for (std::size_t field = 0; field < kHeapFieldCount; ++field) {
      if (field == kMaxRequestBytes) {
        if (other.values[field] > values[field])
          values[field] = other.values[field];
      } else {
        values[field] += other.values[field];
      }
    }
  }
};

struct AtomicHeapSample {
  std::atomic<unsigned long long> values[kHeapFieldCount]{};
};

AtomicHeapSample samples[kBucketCount];
std::atomic<unsigned int> active_mask{0};

struct HeapTrace {
  unsigned int phase_mask = 0;
  HeapField kind = kMallocCalls;
  std::size_t bytes = 0;
  long long iteration = 0;
  long long elapsed_ns = 0;
  char thread_name[16]{};
  void* frames[32]{};
  int depth = 0;
  std::atomic<bool> ready{false};
};

const unsigned int kTraceCapacity = 256;
HeapTrace traces[kTraceCapacity];
std::atomic<unsigned int> trace_count{0};
std::atomic<bool> trace_enabled{false};
std::atomic<long long> trace_iteration{0};
std::atomic<long long> trace_epoch_ns{0};
thread_local bool tracing = false;

void captureHeapCall(const unsigned int mask, const HeapField kind,
                     const std::size_t bytes) noexcept {
  if ((mask & (1u << static_cast<unsigned int>(
                   hipdlp::algorithm::SolvePhase::kIterations))) == 0 ||
      !trace_enabled.load(std::memory_order_relaxed) || tracing)
    return;
  const int saved_errno = errno;
  tracing = true;
  const unsigned int index =
      trace_count.fetch_add(1, std::memory_order_relaxed);
  if (index < kTraceCapacity) {
    HeapTrace& trace = traces[index];
    trace.phase_mask = mask;
    trace.kind = kind;
    trace.bytes = bytes;
    trace.iteration = trace_iteration.load(std::memory_order_relaxed);
    trace.elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::steady_clock::now().time_since_epoch())
                           .count() -
                       trace_epoch_ns.load(std::memory_order_relaxed);
    pthread_getname_np(pthread_self(), trace.thread_name,
                       sizeof(trace.thread_name));
    trace.depth = backtrace(trace.frames, 32);
    trace.ready.store(true, std::memory_order_release);
  }
  tracing = false;
  errno = saved_errno;
}

void enableTraces() {
  void* frames[32];
  backtrace(frames, 32);
  char thread_name[16];
  pthread_getname_np(pthread_self(), thread_name, sizeof(thread_name));
  trace_epoch_ns.store(std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::steady_clock::now().time_since_epoch())
                           .count(),
                       std::memory_order_relaxed);
  trace_enabled.store(true, std::memory_order_relaxed);
}

void writeTraces() {
  const unsigned int count = trace_count.load(std::memory_order_relaxed);
  std::fprintf(stderr, "heap_traces total=%u capacity=%u\n", count,
               kTraceCapacity);
  const char* const names[] = {"malloc", "calloc", "realloc", "aligned",
                               "free"};
  for (unsigned int index = 0; index < count && index < kTraceCapacity;
       ++index) {
    const HeapTrace& trace = traces[index];
    if (!trace.ready.load(std::memory_order_acquire)) continue;
    std::fprintf(stderr,
                 "heap_trace index=%u phase_mask=%u api=%s bytes=%zu "
                 "iteration=%lld elapsed_ns=%lld thread=%s\n",
                 index, trace.phase_mask, names[trace.kind], trace.bytes,
                 trace.iteration, trace.elapsed_ns, trace.thread_name);
    backtrace_symbols_fd(trace.frames, trace.depth, 2);
  }
}

void recordAllocation(unsigned int mask, const HeapField kind,
                      const std::size_t bytes, const bool failed) noexcept {
  captureHeapCall(mask, kind, bytes);
  for (std::size_t bucket = 0; mask != 0; ++bucket, mask >>= 1) {
    if ((mask & 1u) == 0) continue;
    AtomicHeapSample& sample = samples[bucket];
    sample.values[kind].fetch_add(1, std::memory_order_relaxed);
    sample.values[kRequestedBytes].fetch_add(bytes, std::memory_order_relaxed);
    unsigned long long maximum =
        sample.values[kMaxRequestBytes].load(std::memory_order_relaxed);
    while (maximum < bytes &&
           !sample.values[kMaxRequestBytes].compare_exchange_weak(
               maximum, bytes, std::memory_order_relaxed)) {
    }
    if (failed)
      sample.values[kFailedCalls].fetch_add(1, std::memory_order_relaxed);
  }
}

void recordFree(unsigned int mask, const std::size_t usable_bytes) noexcept {
  captureHeapCall(mask, kFreeCalls, usable_bytes);
  for (std::size_t bucket = 0; mask != 0; ++bucket, mask >>= 1) {
    if ((mask & 1u) == 0) continue;
    samples[bucket].values[kFreeCalls].fetch_add(1, std::memory_order_relaxed);
    samples[bucket].values[kFreeUsableBytes].fetch_add(
        usable_bytes, std::memory_order_relaxed);
  }
}

void begin(const std::size_t bucket) {
  for (std::size_t field = 0; field < kHeapFieldCount; ++field)
    samples[bucket].values[field].store(0, std::memory_order_relaxed);
  active_mask.fetch_or(1u << bucket, std::memory_order_relaxed);
}

HeapSample snapshot(const std::size_t bucket) {
  HeapSample result;
  for (std::size_t field = 0; field < kHeapFieldCount; ++field)
    result.values[field] =
        samples[bucket].values[field].load(std::memory_order_relaxed);
  return result;
}

HeapSample end(const std::size_t bucket) {
  active_mask.fetch_and(~(1u << bucket), std::memory_order_relaxed);
  return snapshot(bucket);
}

bool selfTest() {
  void* (*volatile allocate)(std::size_t) = &std::malloc;
  void* (*volatile allocate_zero)(std::size_t, std::size_t) = &std::calloc;
  void* (*volatile resize)(void*, std::size_t) = &std::realloc;
  void (*volatile release)(void*) = &std::free;
  void* (*volatile allocate_new)(std::size_t) = &::operator new;
  void (*volatile release_new)(void*) = &::operator delete;
  begin(kSelfTest);
  void* first = allocate(17);
  void* second = allocate_zero(3, 19);
  first = resize(first, 31);
  void* third = allocate_new(43);
  release(first);
  release(second);
  release_new(third);
  const HeapSample direct = end(kSelfTest);
  const bool direct_ok =
      direct.values[kMallocCalls] == 2 && direct.values[kCallocCalls] == 1 &&
      direct.values[kReallocCalls] == 1 && direct.values[kFreeCalls] == 3 &&
      direct.values[kRequestedBytes] == 148 &&
      direct.values[kMaxRequestBytes] == 57;

  void* (*volatile allocate_aligned)(std::size_t, std::size_t) =
      &::aligned_alloc;
  void* (*volatile allocate_memalign)(std::size_t, std::size_t) = &::memalign;
  int (*volatile allocate_posix)(void**, std::size_t, std::size_t) =
      &::posix_memalign;
  begin(kSelfTest);
  void* aligned = allocate_aligned(64, 128);
  void* legacy_aligned = allocate_memalign(64, 64);
  void* posix_aligned = nullptr;
  const int posix_status = allocate_posix(&posix_aligned, 64, 192);
  void* unchanged = &posix_aligned;
  errno = EDOM;
  const int invalid_status = allocate_posix(&unchanged, 3, 17);
  const bool alignment_ok =
      aligned != nullptr && legacy_aligned != nullptr && posix_status == 0 &&
      invalid_status == EINVAL && errno == EDOM &&
      unchanged == &posix_aligned &&
      reinterpret_cast<std::uintptr_t>(aligned) % 64 == 0 &&
      reinterpret_cast<std::uintptr_t>(legacy_aligned) % 64 == 0 &&
      reinterpret_cast<std::uintptr_t>(posix_aligned) % 64 == 0;
  release(aligned);
  release(legacy_aligned);
  release(posix_aligned);
  const HeapSample aligned_sample = end(kSelfTest);
  const bool aligned_ok = alignment_ok &&
                          aligned_sample.values[kAlignedCalls] == 4 &&
                          aligned_sample.values[kFailedCalls] == 1 &&
                          aligned_sample.values[kRequestedBytes] == 401 &&
                          aligned_sample.values[kMaxRequestBytes] == 192;

  bool child_seen = false;
  begin(kSelfTest);
  std::thread child([&]() {
    const HeapSample before = snapshot(kSelfTest);
    void* allocation = allocate(7919);
    const HeapSample after = snapshot(kSelfTest);
    child_seen =
        after.values[kMallocCalls] > before.values[kMallocCalls] &&
        after.values[kRequestedBytes] >= before.values[kRequestedBytes] + 7919;
    release(allocation);
  });
  child.join();
  end(kSelfTest);

  begin(kSelfTest);
  {
    Highs highs;
    highs.setOptionValue("output_flag", false);
  }
  const HeapSample library = end(kSelfTest);
  const bool library_ok = library.values[kMallocCalls] > 0;
  std::printf(
      "{\"allocator_self_test\":%s,\"new_via_malloc\":%s,"
      "\"aligned_self_test\":%s,\"child_thread_observed\":%s,"
      "\"libhighs_interposition\":%s,\"direct_calls\":%llu,"
      "\"direct_requested_bytes\":%llu,\"library_malloc_calls\":%llu}\n",
      direct_ok ? "true" : "false",
      direct.values[kMallocCalls] == 2 ? "true" : "false",
      aligned_ok ? "true" : "false", child_seen ? "true" : "false",
      library_ok ? "true" : "false",
      static_cast<unsigned long long>(direct.allocations()),
      static_cast<unsigned long long>(direct.values[kRequestedBytes]),
      static_cast<unsigned long long>(library.values[kMallocCalls]));
  return direct_ok && aligned_ok && child_seen && library_ok;
}

}  // namespace solve_heap_audit

extern "C" {
void* __libc_malloc(std::size_t) noexcept;
void* __libc_calloc(std::size_t, std::size_t) noexcept;
void* __libc_realloc(void*, std::size_t) noexcept;
void __libc_free(void*) noexcept;
void* __libc_memalign(std::size_t, std::size_t) noexcept;

void* malloc(const std::size_t bytes) noexcept {
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  void* allocation = __libc_malloc(bytes);
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kMallocCalls,
                                     bytes, allocation == nullptr);
  return allocation;
}

void* calloc(const std::size_t count, const std::size_t size) noexcept {
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  void* allocation = __libc_calloc(count, size);
  const std::size_t bytes =
      size != 0 && count > std::numeric_limits<std::size_t>::max() / size
          ? std::numeric_limits<std::size_t>::max()
          : count * size;
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kCallocCalls,
                                     bytes, allocation == nullptr);
  return allocation;
}

void* realloc(void* previous, const std::size_t bytes) noexcept {
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  void* allocation = __libc_realloc(previous, bytes);
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kReallocCalls,
                                     bytes,
                                     allocation == nullptr && bytes != 0);
  return allocation;
}

void free(void* allocation) noexcept {
  const int saved_errno = errno;
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  const std::size_t bytes =
      mask != 0 && allocation != nullptr ? malloc_usable_size(allocation) : 0;
  __libc_free(allocation);
  solve_heap_audit::recordFree(mask, bytes);
  errno = saved_errno;
}

void* memalign(const std::size_t alignment, const std::size_t bytes) noexcept {
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  void* allocation = __libc_memalign(alignment, bytes);
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kAlignedCalls,
                                     bytes, allocation == nullptr);
  return allocation;
}

void* aligned_alloc(const std::size_t alignment,
                    const std::size_t bytes) noexcept {
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  void* allocation = nullptr;
  if (alignment == 0 || (alignment & (alignment - 1)) != 0)
    errno = EINVAL;
  else
    allocation = __libc_memalign(alignment, bytes);
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kAlignedCalls,
                                     bytes, allocation == nullptr);
  return allocation;
}

int posix_memalign(void** output, const std::size_t alignment,
                   const std::size_t bytes) noexcept {
  const int saved_errno = errno;
  const unsigned int mask =
      solve_heap_audit::active_mask.load(std::memory_order_relaxed);
  int result = EINVAL;
  if (alignment >= sizeof(void*) && (alignment & (alignment - 1)) == 0) {
    void* allocation = __libc_memalign(alignment, bytes);
    result = allocation != nullptr ? 0 : ENOMEM;
    if (result == 0) *output = allocation;
  }
  solve_heap_audit::recordAllocation(mask, solve_heap_audit::kAlignedCalls,
                                     bytes, result != 0);
  errno = saved_errno;
  return result;
}
}

namespace {

using hipdlp::algorithm::SolveDiagnostics;
using hipdlp::algorithm::SolveObservation;
using hipdlp::algorithm::SolvePhase;
using solve_heap_audit::HeapSample;
using Clock = std::chrono::steady_clock;

enum CudaField {
  kDeviceAllocations,
  kDeviceAllocationBytes,
  kDeviceFrees,
  kDeviceFreeBytes,
  kPinnedAllocations,
  kPinnedAllocationBytes,
  kPinnedFrees,
  kPinnedFreeBytes,
  kH2dCopies,
  kH2dBytes,
  kD2hCopies,
  kD2hBytes,
  kD2dCopies,
  kD2dBytes,
  kSynchronizations,
  kStreamCreations,
  kStreamDestructions,
  kBlasCreations,
  kBlasDestructions,
  kSparseCreations,
  kSparseDestructions,
  kDescriptorCreations,
  kDescriptorDestructions,
  kPrimitiveEnqueues,
  kGraphLaunches,
  kCudaFieldCount
};

const char* const kCudaNames[] = {"device_allocations",
                                  "device_allocation_bytes",
                                  "device_frees",
                                  "device_free_bytes",
                                  "pinned_allocations",
                                  "pinned_allocation_bytes",
                                  "pinned_frees",
                                  "pinned_free_bytes",
                                  "h2d_copies",
                                  "h2d_bytes",
                                  "d2h_copies",
                                  "d2h_bytes",
                                  "d2d_copies",
                                  "d2d_bytes",
                                  "synchronizations",
                                  "stream_creations",
                                  "stream_destructions",
                                  "cublas_creations",
                                  "cublas_destructions",
                                  "cusparse_creations",
                                  "cusparse_destructions",
                                  "descriptor_creations",
                                  "descriptor_destructions",
                                  "primitive_enqueues",
                                  "graph_launches"};
static_assert(sizeof(kCudaNames) / sizeof(kCudaNames[0]) == kCudaFieldCount,
              "CUDA counter labels must match the counters");

struct CudaSample {
  std::uint64_t values[kCudaFieldCount]{};

  void add(const CudaSample& other) {
    for (std::size_t field = 0; field < kCudaFieldCount; ++field)
      values[field] += other.values[field];
  }

  CudaSample since(const CudaSample& before) const {
    CudaSample delta;
    for (std::size_t field = 0; field < kCudaFieldCount; ++field)
      delta.values[field] = values[field] - before.values[field];
    return delta;
  }

  bool noResourcesOrUploads() const {
    if (values[kDeviceAllocations] || values[kDeviceFrees] ||
        values[kPinnedAllocations] || values[kPinnedFrees] ||
        values[kH2dCopies] || values[kH2dBytes])
      return false;
    for (std::size_t field = kStreamCreations; field <= kDescriptorDestructions;
         ++field)
      if (values[field]) return false;
    return true;
  }
};

CudaSample cudaSnapshot(const bool enabled) {
  CudaSample result;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  if (enabled) {
    using hipdlp::backend::cuda::testing::ResourceAudit;
    static std::size_t ResourceAudit::* const fields[] = {
        &ResourceAudit::allocation_calls,
        &ResourceAudit::allocation_bytes,
        &ResourceAudit::free_calls,
        &ResourceAudit::free_bytes,
        &ResourceAudit::pinned_allocation_calls,
        &ResourceAudit::pinned_allocation_bytes,
        &ResourceAudit::pinned_free_calls,
        &ResourceAudit::pinned_free_bytes,
        &ResourceAudit::host_to_device_transfers,
        &ResourceAudit::host_to_device_bytes,
        &ResourceAudit::device_to_host_transfers,
        &ResourceAudit::device_to_host_bytes,
        &ResourceAudit::device_to_device_transfers,
        &ResourceAudit::device_to_device_bytes,
        &ResourceAudit::synchronization_calls,
        &ResourceAudit::stream_creations,
        &ResourceAudit::stream_destructions,
        &ResourceAudit::cublas_creations,
        &ResourceAudit::cublas_destructions,
        &ResourceAudit::cusparse_creations,
        &ResourceAudit::cusparse_destructions,
        &ResourceAudit::descriptor_creations,
        &ResourceAudit::descriptor_destructions,
        &ResourceAudit::direct_primitive_enqueue_calls,
        &ResourceAudit::graph_launch_calls};
    static_assert(sizeof(fields) / sizeof(fields[0]) == kCudaFieldCount,
                  "CUDA snapshot fields must match the counters");
    const ResourceAudit audit = hipdlp::backend::cuda::testing::resourceAudit();
    for (std::size_t field = 0; field < kCudaFieldCount; ++field)
      result.values[field] = audit.*fields[field];
  }
#else
  (void)enabled;
#endif
  return result;
}

struct AuditRow {
  std::uint64_t calls = 0;
  std::uint64_t violations = 0;
  double seconds = 0.0;
  HeapSample heap;
  CudaSample cuda;

  void add(const AuditRow& other) {
    calls += other.calls;
    violations += other.violations;
    seconds += other.seconds;
    heap.add(other.heap);
    cuda.add(other.cuda);
  }

  bool noAllocationOrFree() const {
    return heap.allocations() == 0 &&
           heap.values[solve_heap_audit::kFreeCalls] == 0 &&
           cuda.noResourcesOrUploads();
  }
};

struct PhaseFrame {
  SolvePhase phase;
  HighsInt iteration;
  Clock::time_point start;
  CudaSample cuda;
};

class AuditObserver : public hipdlp::algorithm::SolveObserver {
 public:
  explicit AuditObserver(const bool cuda_enabled)
      : cuda_enabled_(cuda_enabled) {}

  AuditRow phases[solve_heap_audit::kPhaseCount];
  AuditRow before_iterations;
  AuditRow after_iterations;
  AuditRow ordinary_boundaries;
  AuditRow special_boundaries;
  AuditRow terminal_boundaries;
  AuditRow candidate_downloads;
  AuditRow first_candidate_download;
  AuditRow complete_blocks;
  std::uint64_t violations = 0;
  std::size_t depth = 0;

  void onEvent(const SolveObservation& event) noexcept override {
    const std::size_t phase = static_cast<std::size_t>(event.phase);
    if (event.begin) {
      if (solve_heap_audit::trace_enabled.load(std::memory_order_relaxed))
        solve_heap_audit::trace_iteration.store(event.iteration,
                                                std::memory_order_relaxed);
      if (depth == 32 ||
          (solve_heap_audit::active_mask.load(std::memory_order_relaxed) &
           (1u << phase))) {
        ++violations;
        return;
      }
      if (event.phase == SolvePhase::kIterations) {
        before_iterations.calls = 1;
        before_iterations.heap = solve_heap_audit::snapshot(
            static_cast<std::size_t>(SolvePhase::kSolve));
        before_iterations.cuda =
            cudaSnapshot(cuda_enabled_).since(frames_[0].cuda);
        before_iterations.seconds =
            std::chrono::duration<double>(Clock::now() - frames_[0].start)
                .count();
        iterating_ = true;
        had_iterations_ = true;
      }
      if (event.phase == SolvePhase::kBoundary) {
        boundary_active_ = true;
        boundary_special_ = false;
      }
      if (boundary_active_ && (event.phase == SolvePhase::kDownload ||
                               event.phase == SolvePhase::kHostVerification ||
                               event.phase == SolvePhase::kFixedPointFallback))
        boundary_special_ = true;
      PhaseFrame& frame = frames_[depth++];
      frame.phase = event.phase;
      frame.iteration = event.iteration;
      frame.cuda = cudaSnapshot(cuda_enabled_);
      solve_heap_audit::begin(phase);
      frame.start = Clock::now();
      return;
    }

    const Clock::time_point stop = Clock::now();
    if (depth == 0 || frames_[depth - 1].phase != event.phase) {
      ++violations;
      return;
    }
    const PhaseFrame& frame = frames_[--depth];
    AuditRow delta;
    delta.calls = 1;
    delta.seconds = std::chrono::duration<double>(stop - frame.start).count();
    delta.heap = solve_heap_audit::end(phase);
    delta.cuda = cudaSnapshot(cuda_enabled_).since(frame.cuda);

    if (event.phase == SolvePhase::kStep &&
        (!delta.noAllocationOrFree() || delta.cuda.values[kD2hCopies] ||
         delta.cuda.values[kSynchronizations]))
      ++delta.violations;
    if (iterating_ &&
        (event.phase == SolvePhase::kDownload ||
         event.phase == SolvePhase::kHostPreparation ||
         event.phase == SolvePhase::kReferenceFixedPoint) &&
        !delta.noAllocationOrFree())
      ++delta.violations;
    if (event.phase == SolvePhase::kDownload && iterating_ &&
        event.iteration > 0) {
      candidate_downloads.add(delta);
      if (first_candidate_download.calls == 0) first_candidate_download = delta;
    }
    if (event.phase == SolvePhase::kBoundary) {
      if (iterating_) {
        if (!event.continued) {
          terminal_boundaries.add(delta);
        } else if (boundary_special_) {
          special_boundaries.add(delta);
        } else {
          if (!delta.noAllocationOrFree() ||
              (cuda_enabled_ && (delta.cuda.values[kD2hCopies] != 1 ||
                                 delta.cuda.values[kD2hBytes] != 112 ||
                                 delta.cuda.values[kSynchronizations] != 1)))
            ++delta.violations;
          ordinary_boundaries.add(delta);
        }
      }
      boundary_active_ = false;
    }
    if (event.phase == SolvePhase::kIterationBlock &&
        event.iteration - frame.iteration == 40)
      complete_blocks.add(delta);
    violations += delta.violations;
    phases[phase].add(delta);
    if (event.phase == SolvePhase::kIterations) {
      iterating_ = false;
      teardown_start_ = Clock::now();
      teardown_cuda_ = cudaSnapshot(cuda_enabled_);
      solve_heap_audit::begin(solve_heap_audit::kTeardown);
    }
    if (event.phase == SolvePhase::kSolve) {
      if (had_iterations_) {
        after_iterations.calls = 1;
        after_iterations.seconds =
            std::chrono::duration<double>(stop - teardown_start_).count();
        after_iterations.heap =
            solve_heap_audit::end(solve_heap_audit::kTeardown);
        after_iterations.cuda =
            cudaSnapshot(cuda_enabled_).since(teardown_cuda_);
      } else {
        before_iterations = delta;
      }
    }
  }

 private:
  bool cuda_enabled_;
  bool iterating_ = false;
  bool had_iterations_ = false;
  bool boundary_active_ = false;
  bool boundary_special_ = false;
  PhaseFrame frames_[32];
  Clock::time_point teardown_start_;
  CudaSample teardown_cuda_;
};

struct PhaseDescription {
  const char* name;
  double SolveDiagnostics::* seconds;
};

const PhaseDescription kPhases[] = {
    {"solve", &SolveDiagnostics::total_wall_seconds},
    {"preprocessing", &SolveDiagnostics::preprocessing_seconds},
    {"backend_setup", &SolveDiagnostics::backend_setup_seconds},
    {"operator_norm", &SolveDiagnostics::operator_norm_seconds},
    {"iterations", &SolveDiagnostics::iterations_seconds},
    {"iteration_block", &SolveDiagnostics::iteration_block_seconds},
    {"step", &SolveDiagnostics::step_seconds},
    {"boundary", &SolveDiagnostics::boundary_seconds},
    {"download", &SolveDiagnostics::download_seconds},
    {"host_preparation", &SolveDiagnostics::host_preparation_seconds},
    {"host_verification", &SolveDiagnostics::host_verification_seconds},
    {"reference_fixed_point", &SolveDiagnostics::reference_fixed_point_seconds},
    {"fixed_point_fallback", &SolveDiagnostics::fixed_point_fallback_seconds}};
static_assert(sizeof(kPhases) / sizeof(kPhases[0]) ==
                  solve_heap_audit::kPhaseCount,
              "Phase labels must match the solve observer");

void jsonString(const char* value) {
  std::putchar('"');
  for (const unsigned char* character =
           reinterpret_cast<const unsigned char*>(value);
       *character; ++character) {
    if (*character == '"' || *character == '\\') {
      std::putchar('\\');
      std::putchar(*character);
    } else if (*character < 32) {
      std::printf("\\u%04x", static_cast<unsigned int>(*character));
    } else {
      std::putchar(*character);
    }
  }
  std::putchar('"');
}

void writeRow(const AuditRow& row, const bool cuda_enabled) {
  static const char* const heap_names[] = {
      "malloc_calls",    "calloc_calls",      "realloc_calls",
      "aligned_calls",   "free_calls",        "failed_calls",
      "requested_bytes", "max_request_bytes", "free_usable_bytes"};
  std::printf(
      "{\"calls\":%llu,\"violations\":%llu,\"wall_seconds\":%.17g,\"heap\":{",
      static_cast<unsigned long long>(row.calls),
      static_cast<unsigned long long>(row.violations), row.seconds);
  for (std::size_t field = 0; field < solve_heap_audit::kHeapFieldCount;
       ++field)
    std::printf("%s\"%s\":%llu", field == 0 ? "" : ",", heap_names[field],
                static_cast<unsigned long long>(row.heap.values[field]));
  std::printf("},\"cuda\":");
  if (cuda_enabled) {
    std::putchar('{');
    for (std::size_t field = 0; field < kCudaFieldCount; ++field)
      std::printf("%s\"%s\":%llu", field == 0 ? "" : ",", kCudaNames[field],
                  static_cast<unsigned long long>(row.cuda.values[field]));
    std::putchar('}');
  } else {
    std::printf("null");
  }
  std::putchar('}');
}

struct Arguments {
  const char* instance = nullptr;
  const char* backend = "cpu";
  const char* step_size = "fixed";
  const char* scaling = "auto";
  HighsInt iteration_limit = 400;
  double time_limit = kHighsInf;
  double kkt_tolerance = 1e-4;
  bool require_ordinary = false;
  bool heap_stacks = false;
};

bool parseArguments(const int argc, char** argv, Arguments& arguments) {
  for (int index = 1; index < argc; ++index) {
    const char* flag = argv[index];
    if (std::strcmp(flag, "--require-ordinary") == 0) {
      arguments.require_ordinary = true;
      continue;
    }
    if (std::strcmp(flag, "--heap-stacks") == 0) {
      arguments.heap_stacks = true;
      continue;
    }
    if (flag[0] != '-' && arguments.instance == nullptr) {
      arguments.instance = flag;
      continue;
    }
    if (++index == argc) return false;
    const char* value = argv[index];
    if (std::strcmp(flag, "--instance") == 0)
      arguments.instance = value;
    else if (std::strcmp(flag, "--backend") == 0)
      arguments.backend = value;
    else if (std::strcmp(flag, "--step-size") == 0)
      arguments.step_size = value;
    else if (std::strcmp(flag, "--scaling") == 0)
      arguments.scaling = value;
    else {
      errno = 0;
      char* end = nullptr;
      if (std::strcmp(flag, "--iteration-limit") == 0) {
        const long long limit = std::strtoll(value, &end, 10);
        if (errno || end == value || *end || limit < 0 ||
            limit > std::numeric_limits<HighsInt>::max())
          return false;
        arguments.iteration_limit = static_cast<HighsInt>(limit);
      } else if (std::strcmp(flag, "--time-limit") == 0 ||
                 std::strcmp(flag, "--kkt-tolerance") == 0) {
        const double number = std::strtod(value, &end);
        if (errno || end == value || *end || !std::isfinite(number) ||
            number < 0)
          return false;
        if (std::strcmp(flag, "--time-limit") == 0)
          arguments.time_limit = number;
        else if (number > 0)
          arguments.kkt_tolerance = number;
        else
          return false;
      } else
        return false;
    }
  }
  return arguments.instance != nullptr &&
         (std::strcmp(arguments.backend, "cpu") == 0 ||
          std::strcmp(arguments.backend, "cuda") == 0) &&
         (std::strcmp(arguments.step_size, "fixed") == 0 ||
          std::strcmp(arguments.step_size, "pid") == 0);
}

const char* statusName(const hipdlp::api::SolverStatus status) {
  switch (status) {
    case hipdlp::api::SolverStatus::kNotSet:
      return "not_set";
    case hipdlp::api::SolverStatus::kOptimal:
      return "optimal";
    case hipdlp::api::SolverStatus::kInfeasible:
      return "infeasible";
    case hipdlp::api::SolverStatus::kUnbounded:
      return "unbounded";
    case hipdlp::api::SolverStatus::kUnboundedOrInfeasible:
      return "unbounded_or_infeasible";
    case hipdlp::api::SolverStatus::kIterationLimit:
      return "iteration_limit";
    case hipdlp::api::SolverStatus::kTimeLimit:
      return "time_limit";
    case hipdlp::api::SolverStatus::kError:
      return "error";
  }
  return "unknown";
}

void writeReport(const Arguments& arguments, const HighsLp& lp,
                 const hipdlp::api::SolveResult& result,
                 const SolveDiagnostics& diagnostics,
                 const AuditObserver& observer, const AuditRow& preamble,
                 const bool cuda_enabled) {
  std::printf("{\"schema_version\":1,\"instance\":");
  jsonString(arguments.instance);
  std::printf(",\"backend\":");
  jsonString(arguments.backend);
  std::printf(",\"step_size\":");
  jsonString(arguments.step_size);
  std::printf(",\"scaling\":");
  jsonString(arguments.scaling);
  std::printf(",\"glibc\":");
  jsonString(gnu_get_libc_version());
  std::printf(",\"compiler\":");
  jsonString(__VERSION__);
  std::printf(
      ",\"rows\":%lld,\"columns\":%lld,\"nonzeros\":%llu,"
      "\"iteration_limit\":%lld,\"iterations\":%lld,\"kkt_tolerance\":%.17g,"
      "\"status\":\"%s\",\"solve_ok\":%s,\"error\":",
      static_cast<long long>(lp.num_row_), static_cast<long long>(lp.num_col_),
      static_cast<unsigned long long>(lp.a_matrix_.value_.size()),
      static_cast<long long>(arguments.iteration_limit),
      static_cast<long long>(result.iteration_count), arguments.kkt_tolerance,
      statusName(result.status), result.ok() ? "true" : "false");
  jsonString(result.error.message.c_str());
  std::printf(
      ",\"audit_metadata\":{\"observer_storage_bytes\":%llu,"
      "\"heap_counter_storage_bytes\":%llu,\"cuda_counters_available\":%s,"
      "\"heap_scope\":\"process-wide interposable glibc entry points\","
      "\"heap_apis\":\"malloc,calloc,realloc,free,memalign,aligned_alloc,posix_"
      "memalign\","
      "\"new_included_via_malloc\":true,\"loader_and_reporting_excluded\":true,"
      "\"vendor_hidden_allocators_and_direct_driver_calls\":\"unmeasured\","
      "\"phase_times\":\"inclusive host wall; step is enqueue time; blocks "
      "include check waits\","
      "\"baseline\":\"not measured\"},\n\"preamble\":",
      static_cast<unsigned long long>(sizeof(observer)),
      static_cast<unsigned long long>(sizeof(solve_heap_audit::samples)),
      cuda_enabled ? "true" : "false");
  writeRow(preamble, false);
  std::printf(
      ",\n\"diagnostics\":{\"restarts\":%llu,\"primal_weight_updates\":%llu,"
      "\"check_boundaries\":%llu,\"cadence_check_boundaries\":%llu,"
      "\"cadence_packet_d2h_bytes\":%llu,\"cadence_packet_copies\":%llu,"
      "\"cadence_packet_synchronizations\":%llu,\"candidate_downloads\":%llu,"
      "\"host_oracle_checks\":%llu,\"certificate_checks\":%llu,"
      "\"restart_reference_checks\":%llu,\"fixed_point_evaluations\":%llu},"
      "\n\"contracts\":{\"violations\":%llu,\"balanced_events\":%s,"
      "\"step_heap_observed\":%s,\"ordinary_cuda_observed\":%s},\n\"phases\":{",
      static_cast<unsigned long long>(diagnostics.restarts),
      static_cast<unsigned long long>(diagnostics.primal_weight_updates),
      static_cast<unsigned long long>(diagnostics.check_boundaries),
      static_cast<unsigned long long>(diagnostics.ordinary_check_boundaries),
      static_cast<unsigned long long>(diagnostics.ordinary_check_d2h_bytes),
      static_cast<unsigned long long>(diagnostics.ordinary_check_transfers),
      static_cast<unsigned long long>(
          diagnostics.ordinary_check_synchronizations),
      static_cast<unsigned long long>(diagnostics.candidate_downloads),
      static_cast<unsigned long long>(diagnostics.host_oracle_checks),
      static_cast<unsigned long long>(diagnostics.certificate_checks),
      static_cast<unsigned long long>(diagnostics.restart_reference_checks),
      static_cast<unsigned long long>(diagnostics.fixed_point_evaluations),
      static_cast<unsigned long long>(observer.violations),
      observer.depth == 0 ? "true" : "false",
      observer.phases[static_cast<std::size_t>(SolvePhase::kStep)].calls
          ? "true"
          : "false",
      cuda_enabled && observer.ordinary_boundaries.calls ? "true" : "false");
  for (std::size_t phase = 0; phase < solve_heap_audit::kPhaseCount; ++phase) {
    if (phase) std::putchar(',');
    jsonString(kPhases[phase].name);
    std::putchar(':');
    AuditRow row = observer.phases[phase];
    row.seconds = diagnostics.*kPhases[phase].seconds;
    writeRow(row, cuda_enabled);
  }
  std::printf("},\n\"before_iterations\":");
  writeRow(observer.before_iterations, cuda_enabled);
  std::printf(",\n\"after_iterations\":");
  writeRow(observer.after_iterations, cuda_enabled);
  std::printf(",\n\"ordinary_boundaries\":");
  writeRow(observer.ordinary_boundaries, cuda_enabled);
  std::printf(",\n\"special_continuing_boundaries\":");
  writeRow(observer.special_boundaries, cuda_enabled);
  std::printf(",\n\"terminal_boundaries\":");
  writeRow(observer.terminal_boundaries, cuda_enabled);
  std::printf(",\n\"candidate_downloads\":");
  writeRow(observer.candidate_downloads, cuda_enabled);
  std::printf(",\n\"first_candidate_download\":");
  writeRow(observer.first_candidate_download, cuda_enabled);
  std::printf(",\n\"complete_40_step_blocks\":");
  writeRow(observer.complete_blocks, cuda_enabled);
  std::printf(",\n\"cuda_resource_maxima\":");
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  if (cuda_enabled) {
    const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
    std::printf(
        "{\"device_allocation_bytes\":%llu,\"pinned_allocation_bytes\":%llu}",
        static_cast<unsigned long long>(audit.max_allocation_bytes),
        static_cast<unsigned long long>(audit.max_pinned_allocation_bytes));
  } else
#endif
    std::printf("null");
  std::printf("}\n");
}

}  // namespace

int main(const int argc, char** argv) {
  solve_heap_audit::begin(solve_heap_audit::kPreamble);
  const Clock::time_point preamble_start = Clock::now();
  if (argc == 2 && std::strcmp(argv[1], "--self-test") == 0) {
    solve_heap_audit::end(solve_heap_audit::kPreamble);
    return solve_heap_audit::selfTest() ? 0 : 1;
  }
  Arguments arguments;
  if (!parseArguments(argc, argv, arguments)) {
    std::fprintf(
        stderr,
        "usage: hipdlp_solve_audit --instance MODEL --backend cpu|cuda\n"
        "       [--step-size fixed|pid] [--iteration-limit 400]\n"
        "       [--scaling auto|off|curtis_reid|legacy] [--time-limit "
        "SECONDS]\n"
        "       [--kkt-tolerance TOL] [--require-ordinary]\n"
        "       [--heap-stacks] (bounded allocation backtraces on stderr)\n"
        "       hipdlp_solve_audit --self-test\n");
    return 2;
  }
  if (arguments.heap_stacks) solve_heap_audit::enableTraces();
  const bool use_cuda = std::strcmp(arguments.backend, "cuda") == 0;
  bool cuda_counters = false;
#if defined(HIGHS_HIPDLP_CUDA_DETECTED) && defined(HIGHS_HIPDLP_TESTING)
  cuda_counters = use_cuda;
#endif
  if (use_cuda && !cuda_counters) {
    std::fprintf(stderr,
                 "CUDA solve auditing is not compiled in this binary\n");
    return 2;
  }
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(arguments.instance) == HighsStatus::kError) {
    std::fprintf(stderr, "failed to read model: %s\n", arguments.instance);
    return 3;
  }
  const HighsInt step_strategy = std::strcmp(arguments.step_size, "pid") == 0
                                     ? kPdlpStepSizeStrategyPid
                                     : kPdlpStepSizeStrategyFixed;
  if (highs.setOptionValue("solver", kHiPdlpString) != HighsStatus::kOk ||
      highs.setOptionValue("presolve", "off") != HighsStatus::kOk ||
      highs.setOptionValue("hipdlp_backend", arguments.backend) !=
          HighsStatus::kOk ||
      highs.setOptionValue("hipdlp_scaling", arguments.scaling) !=
          HighsStatus::kOk ||
      highs.setOptionValue("pdlp_step_size_strategy", step_strategy) !=
          HighsStatus::kOk ||
      highs.setOptionValue("pdlp_iteration_limit", arguments.iteration_limit) !=
          HighsStatus::kOk ||
      highs.setOptionValue("time_limit", arguments.time_limit) !=
          HighsStatus::kOk ||
      highs.setOptionValue("kkt_tolerance", arguments.kkt_tolerance) !=
          HighsStatus::kOk) {
    std::fprintf(stderr, "invalid solve option\n");
    return 2;
  }
  const hipdlp::api::MigrationOptions migration(
      hipdlp::api::SolverImplementation::kRewrite,
      use_cuda ? hipdlp::backend::Kind::kCuda : hipdlp::backend::Kind::kCpu);
  const hipdlp::api::OptionTranslationResult translated =
      hipdlp::api::translateOptions(highs.getOptions(), migration);
  if (!translated.ok()) {
    std::fprintf(stderr, "%s\n", translated.error.message.c_str());
    return 2;
  }
  SolveDiagnostics diagnostics;
  AuditObserver observer(cuda_counters);
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  context.observer = &observer;
  AuditRow preamble;
  preamble.calls = 1;
  preamble.seconds =
      std::chrono::duration<double>(Clock::now() - preamble_start).count();
  preamble.heap = solve_heap_audit::end(solve_heap_audit::kPreamble);
  const HighsLp& lp = highs.getLp();
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      translated.options, hipdlp::model::OriginalModel(lp), context);
  if (observer.depth != 0 || solve_heap_audit::active_mask.load() != 0)
    ++observer.violations;
  solve_heap_audit::active_mask.store(0);
  if (observer.phases[static_cast<std::size_t>(SolvePhase::kSolve)]
          .heap.allocations() == 0)
    ++observer.violations;
  if (arguments.require_ordinary && observer.ordinary_boundaries.calls == 0)
    ++observer.violations;
  writeReport(arguments, lp, result, diagnostics, observer, preamble,
              cuda_counters);
  if (arguments.heap_stacks) solve_heap_audit::writeTraces();
  return !result.ok() ? 3 : observer.violations != 0 ? 4 : 0;
}