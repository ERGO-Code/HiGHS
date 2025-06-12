#ifndef HIPO_AUXILIARY_H
#define HIPO_AUXILIARY_H

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "ipm/hpm/auxiliary/IntConfig.h"
#include "util/HighsCDouble.h"

namespace hipo {

void counts2Ptr(std::vector<Int>& ptr, std::vector<Int>& w);
void inversePerm(const std::vector<Int>& perm, std::vector<Int>& iperm);
void subtreeSize(const std::vector<Int>& parent, std::vector<Int>& sizes);
void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               std::vector<Int>& ptrT, std::vector<Int>& rowsT);
void transpose(const std::vector<Int>& ptr, const std::vector<Int>& rows,
               const std::vector<double>& val, std::vector<Int>& ptrT,
               std::vector<Int>& rowsT, std::vector<double>& valT);
void symProduct(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                const std::vector<double>& vals, const std::vector<double>& x,
                std::vector<double>& y, double alpha = 1.0);
void symProductQuad(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                    const std::vector<double>& vals,
                    const std::vector<double>& x, std::vector<HighsCDouble>& y,
                    double alpha);
void childrenLinkedList(const std::vector<Int>& parent, std::vector<Int>& head,
                        std::vector<Int>& next);
void reverseLinkedList(std::vector<Int>& head, std::vector<Int>& next);
void dfsPostorder(Int node, Int& start, std::vector<Int>& head,
                  const std::vector<Int>& next, std::vector<Int>& order);
void processEdge(Int j, Int i, const std::vector<Int>& first,
                 std::vector<Int>& maxfirst, std::vector<Int>& delta,
                 std::vector<Int>& prevleaf, std::vector<Int>& ancestor);
double getDiagStart(Int n, Int k, Int nb, Int n_blocks, std::vector<Int>& start,
                    bool triang = false);

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

class Clock {
  std::chrono::high_resolution_clock::time_point t0;

 public:
  Clock();
  void start();
  double stop() const;
};

}  // namespace hipo

#endif
