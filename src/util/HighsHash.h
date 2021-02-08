#ifndef HIGHS_UTIL_HASH_H_
#define HIGHS_UTIL_HASH_H_

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <type_traits>
#include <utility>

struct HighsHashHelpers {
  using u32 = std::uint32_t;
  using u64 = std::uint64_t;
  /// mersenne prime 2^31 - 1
  static constexpr u32 M31() { return u32{0x7fffffff}; };

  /// mersenne prime 2^61 - 1
  static constexpr u64 M61() { return u64{0x1fffffffffffffff}; };

  /// compute a * b mod 2^31-1
  static u32 multiply_modM31(u32 a, u32 b) {
    u64 result = a;
    result *= b;
    result = (result >> 31) + (result & (u64)M31());
    if (result >= M31()) result -= M31();
    return result;
  }

  /// compute a * b mod 2^31-1
  static u64 multiply_modM61(u64 a, u64 b) {
    u64 ahi = a >> 32;
    u64 bhi = b >> 32;
    u64 alo = a & 0xffffffffu;
    u64 blo = b & 0xffffffffu;

    // compute the different order terms with adicities 2^64, 2^32, 2^0
    u64 term_64 = ahi * bhi;
    u64 term_32 = ahi * blo + bhi * alo;
    u64 term_0 = blo * blo;

    // now extract the upper 61 and the lower 61 bits of the result a * b
    u64 ab61 = term_64 << 3 | (term_32 + (term_0 >> 32)) >> 61;
    u64 ab0 = (term_0 + term_32) & M61();

    // finally take the result modulo M61 which is computed by exploiting
    // that M61 is a mersenne prime
    u64 result = ab0 + ab61;
    if (result >= M61()) result -= M61();
    return result;
  }

  // take the given 32bit number a, which is assumed to be smaller than the
  // mersenne prime M31=2^31-1, to the power of e modulo the M31
  static u32 modexp_M31(u32 a, u64 e) {
    u32 result = a;

    while (e != 1) {
      // square
      result = multiply_modM31(result, result);

      // multiply with a if exponent is odd
      if (e & 1) result = multiply_modM31(result, a);

      // shift to next bit
      e = e >> 1;
    }

    return result;
  }

  static u64 modexp_M61(u64 a, u64 e) {
    u64 result = a;

    while (e != 1) {
      // square
      result = multiply_modM61(result, result);

      // multiply with a if exponent is odd
      if (e & 1) result = multiply_modM61(result, a);

      // shift to next bit
      e = e >> 1;
    }

    return result;
  }

  static constexpr u32 fibonacci_muliplier32() { return 0x9e3779b9u; }

  static constexpr u64 fibonacci_muliplier64() {
    return uint64_t{0x9e3779b97f4a7c15};
  }

  static uint64_t pair_hash_1(u32 a, u32 b) {
    constexpr static u64 c1 = uint64_t{0xc8497d2a400d9551};
    constexpr static u64 c2 = uint64_t{0x80c8963be3e4c2f3};
    return (a + c1) * (b + c2);
  }

  static uint64_t pair_hash_2(uint32_t a, uint32_t b) {
    constexpr static u64 c1 = uint64_t{0x042d8680e260ae5b};
    constexpr static u64 c2 = uint64_t{0x8a183895eeac1536};
    return (a + c1) * (b + c2);
  }

  static uint64_t pair_hash_3(uint32_t a, uint32_t b) {
    constexpr static u64 c1 = uint64_t{0xa94e9c75f80ad6de};
    constexpr static u64 c2 = uint64_t{0x7e92251dec62835e};
    return (a + c1) * (b + c2);
  }

  static uint64_t pair_hash_4(uint32_t a, uint32_t b) {
    constexpr static u64 c1 = uint64_t{0x07294165cb671455};
    constexpr static u64 c2 = uint64_t{0x89b0f6212b0a4292};
    return (a + c1) * (b + c2);
  }

  static constexpr size_t fibonacci_muliplier() {
    return sizeof(size_t) == 4 ? fibonacci_muliplier32()
                               : fibonacci_muliplier64();
  }

  static constexpr size_t rotate_left(size_t x, int n) {
    return (x << n) | (x >> (sizeof(size_t) * 8 - n));
  }

  static void sparse_combine(u32& hash, int index, u32 value) {
    // we take each value of the sparse hash as coefficient for a polynomial
    // of the finite field modulo the mersenne prime 2^31-1 where the monomial
    // for a sparse entry has the degree of its index. We evaluate at the
    // polynomial at a random constnat. This allows to compute the hashes of
    // sparse vectors independently of each others nonzero contribution and
    // therefore allows to use the order of best access patterns for cache
    // performance. E.g. we can compute a strong hash value for parallel row and
    // column detection and only need to loop over the nonzeros once in
    // arbitrary order. This comes at the expense of more expensive hash
    // calculations as it would be more efficient to evaluate the polynomial
    // with horners scheme, but allows for parallelization and arbitrary order.
    constexpr u32 a = 0x5f4dd81u;
    hash += HighsHashHelpers::multiply_modM31(
        value, HighsHashHelpers::modexp_M31(a, index + 1));
    hash = (hash >> 31) + (hash & HighsHashHelpers::M31());
    if (hash >= HighsHashHelpers::M31()) hash -= HighsHashHelpers::M31();
  }

  static void sparse_combine(u64& hash, int index, u64 value) {
    // 64bit version of sparse_combine using the mersenne prime 2^61 - 1. Less
    // efficient as it needs to compute a 122bit product requiring 3
    // multiplications instead of one and combine it using more integer
    // operations.
    constexpr u64 a = u64{0x9da3e9244a42ca9};
    hash += multiply_modM61(value, modexp_M61(a, index + 1));
    hash = (hash >> 61) + (hash & HighsHashHelpers::M61());
    if (hash >= HighsHashHelpers::M61()) hash -= HighsHashHelpers::M61();
  }

  static void combine(size_t& hash, size_t val) {
    // very fast but decent order dependent hash combine function
    hash = (rotate_left(hash, 5) ^ val) * fibonacci_muliplier();
  }

  template <typename T>
  static size_t vector_hash(const T* vals, size_t numvals) {
    std::array<uint32_t, 2> pair;
    uint64_t hash = 0;
    constexpr uint64_t a = uint64_t{0x45694844fe407};

    const char* dataptr = (const char*)vals;
    const char* dataend = (const char*)(vals + numvals);

    while (dataptr != dataend) {
      size_t chunksize = std::min(int(dataend - dataptr), 32);
      size_t chunkhash = 0;
      switch (chunksize) {
        case 32:
          if (hash != 0) hash = multiply_modM61(hash, a);
          // fall through
        case 31:
        case 30:
        case 29:
        case 28:
        case 27:
        case 26:
        case 25:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash_4(pair[0], pair[1]);
          dataptr += 8;
          chunksize -= 8;
          // fall through
        case 24:
        case 23:
        case 22:
        case 21:
        case 20:
        case 19:
        case 18:
        case 17:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash_3(pair[0], pair[1]);
          dataptr += 8;
          chunksize -= 8;
          // fall through
        case 16:
        case 15:
        case 14:
        case 13:
        case 12:
        case 11:
        case 10:
        case 9:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash_2(pair[0], pair[1]);
          dataptr += 8;
          chunksize -= 8;
          // fall through
        case 8:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 1:

          std::memcpy(&pair[0], dataptr, chunksize);
          chunkhash += pair_hash_1(pair[0], pair[1]);
          dataptr += chunksize;
      }

      hash += chunkhash >> 32;
    }

    return hash;
  }

  template <typename T,
            typename std::enable_if<(sizeof(T) <= 8) && (sizeof(T) >= 4),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<uint32_t, 2> bytes;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return pair_hash_1(bytes[0], bytes[1]) >> 32;
  }

  template <typename T,
            typename std::enable_if<(sizeof(T) >= 12) && (sizeof(T) <= 16),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<uint32_t, 4> bytes;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash_1(bytes[0], bytes[1]) +
            pair_hash_2(bytes[2], bytes[3])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<(sizeof(T) >= 20) && (sizeof(T) <= 24),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<uint32_t, 6> bytes;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash_1(bytes[0], bytes[1]) + pair_hash_2(bytes[2], bytes[3]) +
            pair_hash_3(bytes[4], bytes[5])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<(sizeof(T) >= 28) && (sizeof(T) <= 32),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<uint32_t, 8> bytes;
    std::memcpy(&bytes[0], &val);
    return (pair_hash_1(bytes[0], bytes[1]) + pair_hash_2(bytes[2], bytes[3]) +
            pair_hash_3(bytes[4], bytes[5]) +
            pair_hash_4(bytes[6], bytes[7])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<(sizeof(T) > 32), int>::type = 0>
  static size_t hash(const T& val) {
    return vector_hash(&val, 1);
  }

  static constexpr double golden_ratio_reciprocal() {
    return 0.61803398874989484;
  }

  static uint32_t double_hash_code(double val) {
    // we multiply by some irrational number, so that the buckets in which we
    // put the real numbers do not break on a power of two pattern. E.g.
    // consider the use case for detecting parallel rows when we have two
    // parallel rows scaled to have their largest coefficient 1.0 and another
    // coefficient which is 0.5
    // +- epsilon. Clearly we want to detect those rows as parallel and give
    // them the same hash value for small enough epsilon. The exponent,
    // however will switch to -2 for the value just below 0.5 and the hashcodes
    // will differ. when multiplying with the reciprocal of the golden ratio the
    // exact 0.5 will yield 0.30901699437494742 and 0.5 - 1e-9 will yield
    // 0.3090169937569134 which has the same exponent and matches in the most
    // significant bits. Hence it yields the same hashcode. Obviously there will
    // now be different values which exhibit the same pattern as the 0.5 case,
    // but they do not have a small denominator like 1/2 in their rational
    // representation but are power of two multiples of the golden ratio and
    // therefore irrational, which we do not expect in non-artifical input data.
    int exponent;
    double hashbits = std::frexp(val * golden_ratio_reciprocal(), &exponent);

    // some extra casts to be more verbose about what is happening.
    // We want the exponent to use only 16bits so that the remaining 16 bits
    // are used for the most significant bits of the mantissa and the sign bit.
    // casting to unsigned 16bits first ensures that the value after the cast is
    // defined to be UINT16_MAX - |exponent| when the exponent is negative.
    // casting the exponent to a uint32_t directly would give wrong promotion
    // of negative exponents as UINT32_MAX - |exponent| and take up to many bits
    // or possibly loose information after the 16 bit shift. For the mantissa we
    // take the 15 most significant bits, even though we could squeeze out a few
    // more of the exponent. We don't need more bits as this would make the
    // buckets very small and might miss more values that are equal within
    // epsilon. Therefore the most significant 15 bits of the mantissa and the
    // sign is encoded in the 16 lower bits of the hashcode and the upper 16bits
    // encode the sign and value of the exponent.
    uint32_t hashvalue = (uint32_t)(uint16_t)(int16_t)exponent;
    hashvalue = (hashvalue << 16) |
                (uint32_t)(uint16_t)(int16_t)std::ldexp(hashbits, 15);

    return hashvalue;
  }
};

struct HighsHasher {
  template <typename T>
  size_t operator()(const T& x) const {
    return HighsHashHelpers::hash(x);
  }
};

struct HighsVectorHasher {
  template <typename T>
  size_t operator()(const std::vector<T>& vec) const {
    return HighsHashHelpers::vector_hash(vec.data(), vec.size());
  }
};

struct HighsVectorEqual {
  template <typename T>
  bool operator()(const std::vector<T>& vec1,
                  const std::vector<T>& vec2) const {
    if (vec1.size() != vec2.size()) return false;
    return std::equal(vec1.begin(), vec1.end(), vec2.begin());
  }
};

#endif