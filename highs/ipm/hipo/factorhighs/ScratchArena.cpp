#include "ScratchArena.h"

#include <cassert>

#include "parallel/HighsParallel.h"

namespace hipo {

void ScratchArena::init() {
  // Without a task executor, everything runs on the calling thread.
  HighsSplitDeque* deque = HighsTaskExecutor::getThisWorkerDeque();
  const Int threads = deque ? deque->getNumWorkers() : 1;
  if ((Int)stacks_.size() < threads) stacks_.resize(threads);
}

Int ScratchArena::thread() const {
  HighsSplitDeque* deque = HighsTaskExecutor::getThisWorkerDeque();
  const Int thread = deque ? deque->getOwnerId() : 0;
  assert(thread < (Int)stacks_.size());
  return thread;
}

void* ScratchArena::push(Int thread, Int64 bytes) {
  assert(bytes >= 0 && bytes % kAlignment == 0);
  if (bytes == 0) return nullptr;

  Stack& stack = stacks_[thread];
  if (stack.chunks.empty()) stack.chunks.resize(1);

  // Chunks below top are in use, chunks above top are empty.
  Chunk* chunk = &stack.chunks[stack.top];
  if (chunk->capacity - chunk->used < bytes) {
    // No room: move up if this chunk is in use, then grow the (empty) chunk.
    if (chunk->used > 0) {
      ++stack.top;
      if ((Int)stack.chunks.size() == stack.top)
        stack.chunks.resize(stack.top + 1);
      chunk = &stack.chunks[stack.top];
    }
    if (chunk->capacity < bytes) {
      // An array of unsigned char from new is aligned for any fundamental type
      // that fits in it, and offsets are multiples of kAlignment, so every
      // buffer is aligned.
      // capacity is reset first, so that it is correct if new throws
      chunk->data.reset();
      chunk->capacity = 0;
      chunk->data.reset(new unsigned char[bytes]);
      chunk->capacity = bytes;
    }
  }

  void* ptr = chunk->data.get() + chunk->used;
  chunk->used += bytes;
  return ptr;
}

void ScratchArena::pop(Int thread, const void* ptr, Int64 bytes) {
  if (bytes == 0) return;

  Stack& stack = stacks_[thread];
  Chunk& chunk = stack.chunks[stack.top];
  // ptr must be the last buffer pushed (LIFO order)
  assert(chunk.used >= bytes);
  assert(ptr == chunk.data.get() + chunk.used - bytes);
  (void)ptr;
  chunk.used -= bytes;
  if (chunk.used == 0 && stack.top > 0) --stack.top;
}

}  // namespace hipo
