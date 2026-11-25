#ifndef HIPO_ALIGNED_MEMORY_H
#define HIPO_ALIGNED_MEMORY_H

#include <cstddef>
#include <memory>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Cache-line alignment of the memory of the operands of BLAS kernels. BLAS
// kernels may round differently depending on the alignment of their operands,
// so a fixed alignment keeps the results independent of where the memory is.
// The alignment also affects the speed of BLAS: 128 bytes is the largest cache
// line of common machines, and a smaller alignment made the dense kernels
// measurably slower.
const std::size_t kCacheAlignBytes = 128;

// Allocate into alloc uninitialised memory for size entries of T, and return a
// pointer to the first of them, which is on an align_bytes boundary.
// align_bytes must be a power of two and a multiple of sizeof(T).
template <typename T>
T* allocateAligned(std::unique_ptr<T[]>& alloc, Int64 size,
                   std::size_t align_bytes) {
  // new T[] is aligned at least to a T, so the first boundary is at most
  // align_bytes / sizeof(T) - 1 entries after its start
  const Int64 alloc_size = size + align_bytes / sizeof(T) - 1;
  alloc.reset(new T[alloc_size]);

  // std::align moves ptr to the first boundary with room for size entries
  // after it
  void* ptr = alloc.get();
  std::size_t space = alloc_size * sizeof(T);
  std::align(align_bytes, size * sizeof(T), ptr, space);
  return static_cast<T*>(ptr);
}

}  // namespace hipo

#endif
