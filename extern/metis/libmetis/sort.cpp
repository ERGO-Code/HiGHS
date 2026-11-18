#include <algorithm>
#include <functional>

#include "metislib.h"

extern "C" void isortd(size_t n, idx_t* base) {
  std::sort(base, base + n, std::greater<idx_t>());
}

extern "C" void ikvsorti(size_t n, ikv_t* base) {
  std::sort(base, base + n,
            [](const ikv_t& a, const ikv_t& b) { return a.key < b.key; });
}