#include "CliqueStack.h"

#include <cassert>
#include <cstring>

namespace hipo {

void CliqueStack::reserve(Int64 size) {
  // the old stack is kept until its content is copied
  std::unique_ptr<double[]> alloc;
  double* stack = allocateAligned(alloc, size, kCacheAlignBytes);

  if (top_ > 0) std::memcpy(stack, stack_, top_ * sizeof(double));
  alloc_ = std::move(alloc);
  stack_ = stack;
  capacity_ = size;
}

void CliqueStack::init(Int64 stack_size) {
  top_ = 0;
  workspace_ = nullptr;
  worksize_ = 0;
  sn_pushed_ = std::stack<std::pair<Int, Int64>>();

  if (capacity_ < stack_size) reserve(stack_size);
}

double* CliqueStack::setup(Int64 clique_size, bool& reallocation) {
  // Clear workspace

  assert(!workspace_ && !worksize_);
  reallocation = false;

  if (clique_size > 0) {
    // This should not trigger reallocation, because init reserves the maximum
    // possible size of the stack.
    if (top_ + clique_size > capacity_) {
      reallocation = true;
      reserve(top_ + clique_size);
    }

    workspace_ = stack_ + top_;
    worksize_ = clique_size;

    // initialize workspace to zero
    std::memset(workspace_, 0, worksize_ * sizeof(double));
  }

  return workspace_;
}

const double* CliqueStack::getChild(Int& child_sn) const {
  // Get the top of the stack, in terms of supernode ID of the child and pointer
  // to its data.

  child_sn = sn_pushed_.top().first;
  Int64 child_size = sn_pushed_.top().second;
  const double* child = stack_ + top_ - child_size;

  return child;
}

void CliqueStack::popChild() {
  // Remove top child from the stack

  Int64 child_size = sn_pushed_.top().second;
  sn_pushed_.pop();

  top_ -= child_size;
}

void CliqueStack::pushWork(Int sn) {
  // Put the content of the workspace at the top of the stack

  if (worksize_ > 0) {
    // the top of the stack has lower address than workspace
    std::memmove(stack_ + top_, workspace_, worksize_ * sizeof(double));

    top_ += worksize_;

    // keep track of supernodes pushed
    sn_pushed_.push({sn, worksize_});
  }

  worksize_ = 0;
  workspace_ = nullptr;
}

}  // namespace hipo