#ifndef FACTORHIGHS_SCRATCH_ARENA_H
#define FACTORHIGHS_SCRATCH_ARENA_H

#include <memory>
#include <new>
#include <type_traits>
#include <vector>

#include "ipm/hipo/auxiliary/AlignedMemory.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

// Scratch memory kept between factorisations and solves, with one stack per
// thread used in LIFO order. A thread that waits inside a task may run other
// tasks, which finish first, so LIFO order holds also with nested parallelism.
// Memory in use is never moved, so pointers stay valid until popped.
//
// Use ScratchBuffer to get memory from the arena: it is released when the
// buffer goes out of scope.

namespace hipo {

class ScratchArena {
  // Memory as allocated, and base, its first kAlignment boundary; capacity and
  // used are in bytes from base.
  struct Chunk {
    std::unique_ptr<unsigned char[]> data{};
    unsigned char* base = nullptr;
    Int64 capacity = 0;
    Int64 used = 0;
  };
  struct Stack {
    std::vector<Chunk> chunks{};
    Int top = 0;
  };
  std::vector<Stack> stacks_{};

 public:
  // Every buffer starts on a kAlignment boundary (see AlignedMemory.h).
  static constexpr Int64 kAlignment = kCacheAlignBytes;

  // Make room for one stack per thread of the task executor. Call it from the
  // main thread, before any task that uses the arena is spawned.
  void init();

  // True if init made a stack for each thread of the current task executor. It
  // may not have, if the executor changed since the last call of init.
  bool matchesThreads() const;

  // Index of the stack of the calling thread.
  Int thread() const;

  // Uninitialised memory for bytes bytes (a multiple of kAlignment) on the
  // stack of thread.
  void* push(Int thread, Int64 bytes);

  // Release ptr, of bytes bytes, which must be the last buffer pushed on the
  // stack of thread.
  void pop(Int thread, const void* ptr, Int64 bytes);
};

template <typename T>
class ScratchBuffer {
  static_assert(std::is_same<T, double>::value || std::is_same<T, Int>::value,
                "ScratchBuffer holds double or Int");
  static_assert(alignof(T) <= ScratchArena::kAlignment,
                "ScratchBuffer type needs more alignment than the arena gives");

  ScratchArena& arena_;
  Int thread_;
  Int64 bytes_;
  T* data_;

  // Bytes of size entries of T, rounded up to a multiple of the arena
  // alignment, so that the next buffer is aligned too.
  static Int64 numBytes(Int64 size) {
    const Int64 bytes = size * (Int64)sizeof(T);
    const Int64 align = ScratchArena::kAlignment;
    return (bytes + align - 1) / align * align;
  }

 public:
  // Uninitialised memory for size entries of type T.
  ScratchBuffer(ScratchArena& arena, Int64 size)
      : arena_{arena},
        thread_{arena.thread()},
        bytes_{numBytes(size)},
        data_{static_cast<T*>(arena.push(thread_, bytes_))} {
    // Start the lifetime of the entries, without initialising them; it
    // generates no code. T() instead of T would zero them.
    for (Int64 i = 0; i < size; ++i) ::new (static_cast<void*>(data_ + i)) T;
  }
  ~ScratchBuffer() { arena_.pop(thread_, data_, bytes_); }

  ScratchBuffer(const ScratchBuffer&) = delete;
  ScratchBuffer& operator=(const ScratchBuffer&) = delete;

  T* data() const { return data_; }
  T& operator[](Int64 i) const { return data_[i]; }
};

}  // namespace hipo

#endif
