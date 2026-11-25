#ifndef FACTORHIGHS_SN_COLUMNS_H
#define FACTORHIGHS_SN_COLUMNS_H

#include <cstddef>
#include <memory>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Columns of the factorisation, stored by supernode in a single buffer:
// supernode sn occupies size(sn) entries from start_[sn]. The buffer is reused
// by later factorisations and grows only when a layout needs more space.
//
// Every supernode starts on a 64-byte boundary: BLAS kernels may round
// differently depending on the alignment of their operands, so this keeps the
// results independent of where a supernode lies in the buffer.
class SnColumns {
  // alignment of each supernode, in bytes and in entries
  static constexpr std::size_t kAlignBytes = 64;
  static constexpr Int64 kAlignEntries = kAlignBytes / sizeof(double);

  // smallest multiple of kAlignEntries that is at least n
  static Int64 roundUpToAlignment(Int64 n) {
    return (n + kAlignEntries - 1) / kAlignEntries * kAlignEntries;
  }

  // memory as allocated, and data_, its first 64-byte boundary
  std::unique_ptr<double[]> alloc_{};
  double* data_ = nullptr;
  Int64 capacity_ = 0;
  std::vector<Int64> start_{};
  std::vector<Int64> size_{};

 public:
  // Set the layout of n_sn supernodes, where supernode sn has size(sn) entries.
  // The entries are left uninitialised.
  template <typename F>
  void setLayout(Int n_sn, F size) {
    // each supernode starts where the previous one ends, rounded up to the
    // next 64-byte boundary
    start_.resize(n_sn + 1);
    size_.resize(n_sn);
    start_[0] = 0;
    for (Int sn = 0; sn < n_sn; ++sn) {
      size_[sn] = size(sn);
      start_[sn + 1] = start_[sn] + roundUpToAlignment(size_[sn]);
    }

    const Int64 total = start_[n_sn];
    if (total > capacity_) {
      // capacity is reset first, so that it is correct if new throws
      alloc_.reset();
      data_ = nullptr;
      capacity_ = 0;

      // new double[] is aligned at least to a double, so the first 64-byte
      // boundary is at most kAlignEntries - 1 entries after its start
      const Int64 alloc_entries = total + kAlignEntries - 1;
      alloc_.reset(new double[alloc_entries]);

      // std::align moves ptr to the first 64-byte boundary with room for total
      // entries after it
      void* ptr = alloc_.get();
      std::size_t space = alloc_entries * sizeof(double);
      std::align(kAlignBytes, total * sizeof(double), ptr, space);
      data_ = static_cast<double*>(ptr);
      capacity_ = total;
    }
  }

  double* operator[](Int sn) { return data_ + start_[sn]; }
  const double* operator[](Int sn) const { return data_ + start_[sn]; }

  // number of entries of supernode sn
  Int64 size(Int sn) const { return size_[sn]; }
};

}  // namespace hipo

#endif
