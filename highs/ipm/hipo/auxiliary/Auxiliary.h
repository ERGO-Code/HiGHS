#ifndef HIPO_AUXILIARY_H
#define HIPO_AUXILIARY_H

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

void inversePerm(const std::vector<Int>& perm, std::vector<Int>& iperm);
void subtreeSize(const std::vector<Int>& parent, std::vector<Int>& sizes);
void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               std::vector<Int>& ptrT, std::vector<Int>& rowsT);
void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               const std::vector<double>& val, std::vector<Int>& ptrT,
               std::vector<Int>& rowsT, std::vector<double>& valT);
void childrenLinkedList(const std::vector<Int>& parent, std::vector<Int>& head,
                        std::vector<Int>& next);
void reverseLinkedList(std::vector<Int>& head, std::vector<Int>& next);
void dfsPostorder(Int node, Int& start, std::vector<Int>& head,
                  const std::vector<Int>& next, std::vector<Int>& order);
void processEdge(Int j, Int i, const std::vector<Int>& first,
                 std::vector<Int>& maxfirst, std::vector<Int>& delta,
                 std::vector<Int>& prevleaf, std::vector<Int>& ancestor);
Int64 getDiagStart(Int n, Int k, Int nb, Int n_blocks,
                   std::vector<Int64>& start, bool triang = false);

template <typename T>
void counts2Ptr(std::vector<T>& ptr, std::vector<T>& w) {
  // Given the column counts in the vector w (of size n),
  // compute the column pointers in the vector ptr (of size n+1),
  // and copy the first n pointers back into w.

  T temp_nz{};
  T n = w.size();
  for (T j = 0; j < n; ++j) {
    ptr[j] = temp_nz;
    temp_nz += w[j];
    w[j] = ptr[j];
  }
  ptr[n] = temp_nz;
}

template <typename T>
void permuteVector(std::vector<T>& v, const std::vector<Int>& perm) {
  // Permute vector v according to permutation perm.
  std::vector<T> temp_v(v);
  for (Int i = 0; i < v.size(); ++i) v[i] = temp_v[perm[i]];
}

template <typename T>
void permuteVectorInverse(std::vector<T>& v, const std::vector<Int>& iperm) {
  // Permute vector v according to inverse permutation iperm.
  std::vector<T> temp_v(v);
  for (Int i = 0; i < v.size(); ++i) v[iperm[i]] = temp_v[i];
}

template <typename T>
void printTest(const std::vector<T>& v, const std::string s) {
  std::ofstream out_file;
  char name[80];
  snprintf(name, 80, "%s.txt", s.c_str());
  out_file.open(name);
  for (T i : v) {
    out_file << std::setprecision(16) << i << '\n';
  }
  out_file.close();
}

template <typename T>
void freeVector(std::vector<T>& v) {
  // Give up memory allocated to v.
  // (technically shrink_to_fit does not guarantee to deallocate)
  v.clear();
  v.shrink_to_fit();
}

class Clock {
  std::chrono::high_resolution_clock::time_point t0;

 public:
  Clock();
  void start();
  double stop() const;
};

}  // namespace hipo

#endif
