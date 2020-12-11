#ifndef HIGHS_UTIL_HASH_H_
#define HIGHS_UTIL_HASH_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>

inline constexpr size_t fibonacci_muliplier() {
  return sizeof(size_t) == 4 ? size_t(0x9e3779b9) : size_t(0x9e3779b97f4a7c15);
}

inline constexpr size_t rotate_left(size_t x, int n) {
  return (x << n) | (x >> (sizeof(size_t) * 8 - n));
}

inline void hash_combine(size_t& hash, size_t val) {
  hash = (rotate_left(hash, 5) ^ val) * fibonacci_muliplier();
}

struct HighsPairHasher {
  template <typename U, typename V>
  size_t operator()(const std::pair<U, V>& pair) const {
    size_t hash = std::hash<U>()(pair.first);
    hash_combine(hash, std::hash<V>()(pair.second));
    return hash;
  }
};

#endif