/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_HASH_H_
#define HIGHS_UTIL_HASH_H_

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
struct HighsHashable : public std::is_trivially_copyable<T> {};

template <typename U, typename V>
struct HighsHashable<std::pair<U, V>>
    : public std::integral_constant<bool, HighsHashable<U>::value &&
                                              HighsHashable<V>::value> {};

template <typename U, typename V>
struct HighsHashable<std::tuple<U, V>> : public HighsHashable<std::pair<U, V>> {
};
template <typename U, typename V, typename W, typename... Args>
struct HighsHashable<std::tuple<U, V, W, Args...>>
    : public std::integral_constant<
          bool, HighsHashable<U>::value &&
                    HighsHashable<std::tuple<V, W, Args...>>::value> {};

struct HighsHashHelpers {
  using u8 = std::uint8_t;
  using i8 = std::int8_t;

  using u16 = std::uint16_t;
  using i16 = std::int16_t;

  using u32 = std::uint32_t;
  using i32 = std::int32_t;

  using u64 = std::uint64_t;
  using i64 = std::uint64_t;

  static constexpr u64 c[] = {
      u64{0xc8497d2a400d9551}, u64{0x80c8963be3e4c2f3}, u64{0x042d8680e260ae5b},
      u64{0x8a183895eeac1536}, u64{0xa94e9c75f80ad6de}, u64{0x7e92251dec62835e},
      u64{0x07294165cb671455}, u64{0x89b0f6212b0a4292}, u64{0x31900011b96bf554},
      u64{0xa44540f8eee2094f}, u64{0xce7ffd372e4c64fc}, u64{0x51c9d471bfe6a10f},
      u64{0x758c2a674483826f}, u64{0xf91a20abe63f8b02}, u64{0xc2a069024a1fcc6f},
      u64{0xd5bb18b70c5dbd59}};

  /// mersenne prime 2^61 - 1
  static constexpr u64 M61() { return u64{0x1fffffffffffffff}; };

  /// compute a * b mod 2^61-1
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
    // that M61 is a mersenne prime, particularly, if a * b = q * 2^61 + r
    // then a * b = (q + r) (mod 2^61 - 1)
    u64 result = ab0 + ab61;
    if (result >= M61()) result -= M61();
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

  template <int k>
  static u64 pair_hash(u32 a, u32 b) {
    return (a + c[2 * k]) * (b + c[2 * k + 1]);
  }

  static void sparse_combine(u64& hash, int index, u64 value) {
    // we take each value of the sparse hash as coefficient for a polynomial
    // of the finite field modulo the mersenne prime 2^61-1 where the monomial
    // for a sparse entry has the degree of its index. We evaluate the
    // polynomial at a random constant. This allows to compute the hashes of
    // sparse vectors independently of each others nonzero contribution and
    // therefore allows to use the order of best access patterns for cache
    // performance. E.g. we can compute a strong hash value for parallel row and
    // column detection and only need to loop over the nonzeros once in
    // arbitrary order. This comes at the expense of more expensive hash
    // calculations as it would be more efficient to evaluate the polynomial
    // with horners scheme, but allows for parallelization and arbitrary order.
    // Since we have 16 random constants available, we slightly improve
    // the scheme by using a lower degree polynomial with 16 variables
    // which we evaluate at the random vector of 16.

    // make sure that the constant has at most 61 bits, as otherwise the modulo
    // algorithm for multiplication mod M61 might not work properly due to
    // overflow
    u64 a = c[index & 15] & M61();
    int degree = (index / 16) + 1;

    hash += multiply_modM61(value, modexp_M61(a, degree));
    hash = (hash >> 61) + (hash & M61());
    if (hash >= M61()) hash -= M61();
  }

  static constexpr u64 fibonacci_muliplier() { return u64{0x9e3779b97f4a7c15}; }

  static constexpr size_t rotate_left(size_t x, int n) {
    return (x << n) | (x >> (sizeof(size_t) * 8 - n));
  }

  static void combine(size_t& hash, size_t val) {
    // very fast but decent order dependent hash combine function
    hash = (rotate_left(hash, 5) ^ val) * fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static u64 vector_hash(const T* vals, size_t numvals) {
    std::array<u32, 2> pair{};
    u64 hash = 0;
    int k = 0;

    const char* dataptr = (const char*)vals;
    const char* dataend = (const char*)(vals + numvals);

    while (dataptr != dataend) {
      size_t numBytes = std::min(size_t(dataend - dataptr), size_t{64});
      size_t numPairs = (numBytes + 7) / 8;
      size_t lastPairBytes = numBytes - (numPairs - 1) * 8;
      size_t chunkhash = 0;
      switch (numPairs) {
        case 8:
          if (hash != 0) hash = multiply_modM61(hash, c[k++ & 15] & M61());
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<0>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 7:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<1>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 6:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<2>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 5:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<3>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 4:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<4>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 3:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<5>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 2:
          std::memcpy(&pair[0], dataptr, 8);
          chunkhash += pair_hash<6>(pair[0], pair[1]);
          dataptr += 8;
          // fall through
        case 1:
          std::memcpy(&pair[0], dataptr, lastPairBytes);
          chunkhash += pair_hash<7>(pair[0], pair[1]);
          dataptr += lastPairBytes;
      }

      hash += chunkhash >> 32;
    }

    return hash;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) <= 8) && (sizeof(T) >= 1),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 2> bytes;
    if (sizeof(T) < 4) bytes[0] = 0;
    if (sizeof(T) < 8) bytes[1] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return pair_hash<0>(bytes[0], bytes[1]) >> 32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 9) && (sizeof(T) <= 16),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 4> bytes;
    if (sizeof(T) < 12) bytes[4] = 0;
    if (sizeof(T) < 16) bytes[5] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 17) && (sizeof(T) <= 24),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 6> bytes;
    if (sizeof(T) < 20) bytes[4] = 0;
    if (sizeof(T) < 24) bytes[5] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 25) && (sizeof(T) <= 32),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 8> bytes;
    if (sizeof(T) < 28) bytes[6] = 0;
    if (sizeof(T) < 32) bytes[7] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5]) +
            pair_hash<3>(bytes[6], bytes[7])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 33) && (sizeof(T) <= 40),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 10> bytes;
    if (sizeof(T) < 36) bytes[8] = 0;
    if (sizeof(T) < 40) bytes[9] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5]) +
            pair_hash<3>(bytes[6], bytes[7]) +
            pair_hash<4>(bytes[8], bytes[9])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 41) && (sizeof(T) <= 48),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 12> bytes;
    if (sizeof(T) < 44) bytes[10] = 0;
    if (sizeof(T) < 48) bytes[11] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5]) +
            pair_hash<3>(bytes[6], bytes[7]) +
            pair_hash<4>(bytes[8], bytes[9]) +
            pair_hash<5>(bytes[10], bytes[11])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 49) && (sizeof(T) <= 56),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 14> bytes;
    if (sizeof(T) < 52) bytes[12] = 0;
    if (sizeof(T) < 56) bytes[13] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5]) +
            pair_hash<3>(bytes[6], bytes[7]) +
            pair_hash<4>(bytes[8], bytes[9]) +
            pair_hash<5>(bytes[10], bytes[11]) +
            pair_hash<6>(bytes[12], bytes[13])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 57) && (sizeof(T) <= 64),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    std::array<u32, 16> bytes;
    if (sizeof(T) < 60) bytes[14] = 0;
    if (sizeof(T) < 64) bytes[15] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) +
            pair_hash<1>(bytes[2], bytes[3]) +
            pair_hash<2>(bytes[4], bytes[5]) +
            pair_hash<3>(bytes[6], bytes[7]) +
            pair_hash<4>(bytes[8], bytes[9]) +
            pair_hash<5>(bytes[10], bytes[11]) +
            pair_hash<6>(bytes[12], bytes[13]) +
            pair_hash<7>(bytes[14], bytes[15])) >>
           32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value && (sizeof(T) > 64),
                                    int>::type = 0>
  static size_t hash(const T& val) {
    return vector_hash(&val, 1);
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static size_t hash(const std::vector<T>& val) {
    return vector_hash(val.data(), val.size());
  }

  template <typename T, typename std::enable_if<
                            std::is_same<decltype(*reinterpret_cast<T*>(0) ==
                                                  *reinterpret_cast<T*>(0)),
                                         bool>::value,
                            int>::type = 0>
  static bool equal(const T& a, const T& b) {
    return a == b;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static bool equal(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) return false;
    return std::memcmp(a.data(), b.data(), sizeof(T) * a.size()) == 0;
  }

  static constexpr double golden_ratio_reciprocal() {
    return 0.61803398874989484;
  }

  static u32 double_hash_code(double val) {
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
    u32 hashvalue = (u32)(u16)(i16)exponent;
    hashvalue = (hashvalue << 16) | (u32)(u16)(i16)std::ldexp(hashbits, 15);

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

template <typename K, typename V>
struct HighsHashTableEntry {
 private:
  K key_;
  V value_;

 public:
  template <typename K_>
  HighsHashTableEntry(K_&& k) : key_(k), value_() {}
  template <typename K_, typename V_>
  HighsHashTableEntry(K_&& k, V_&& v) : key_(k), value_(v) {}

  const K& key() const { return key_; }
  const V& value() const { return value_; }
  V& value() { return value_; }
};

template <typename T>
struct HighsHashTableEntry<T, void> {
 private:
  T value_;

 public:
  template <typename... Args>
  HighsHashTableEntry(Args&&... args) : value_(std::forward<Args>(args)...) {}

  const T& key() const { return value_; }
  const T& value() const { return value_; }
};

template <typename K, typename V = void>
class HighsHashTable {
  struct OpNewDeleter {
    void operator()(void* ptr) { ::operator delete(ptr); }
  };

 public:
  using u8 = std::uint8_t;
  using i8 = std::int8_t;

  using u16 = std::uint16_t;
  using i16 = std::int16_t;

  using u32 = std::uint32_t;
  using i32 = std::int32_t;

  using u64 = std::uint64_t;
  using i64 = std::uint64_t;

  using Entry = HighsHashTableEntry<K, V>;
  using KeyType = K;
  using ValueType = typename std::remove_reference<decltype(
      reinterpret_cast<Entry*>(0)->value())>::type;

  std::unique_ptr<Entry, OpNewDeleter> entries;
  std::unique_ptr<u8[]> metadata;
  u32 tableSizeMask;
  u32 numElements = 0;

  template <typename IterType>
  class HashTableIterator {
    u8* pos;
    u8* end;
    Entry* entryEnd;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = IterType;
    using pointer = IterType*;
    using reference = IterType&;
    using iterator_category = std::forward_iterator_tag;
    HashTableIterator(u8* pos, u8* end, Entry* entryEnd)
        : pos(pos), end(end), entryEnd(entryEnd) {}
    HashTableIterator() = default;

    HashTableIterator<IterType> operator++(int) {
      // postfix
      HashTableIterator<IterType> oldpos = *this;
      for (++pos; pos != end; ++pos)
        if ((*pos) & 0x80u) break;

      return oldpos;
    }

    HashTableIterator<IterType>& operator++() {
      // prefix
      for (++pos; pos != end; ++pos)
        if ((*pos) & 0x80u) break;

      return *this;
    }

    reference operator*() const { return *(entryEnd - (end - pos)); }
    pointer operator->() const { return (entryEnd - (end - pos)); }
    HashTableIterator<IterType> operator+(difference_type v) const {
      for (difference_type k = 0; k != v; ++k) ++(*this);
    }

    bool operator==(const HashTableIterator<IterType>& rhs) const {
      return pos == rhs.pos;
    }
    bool operator!=(const HashTableIterator<IterType>& rhs) const {
      return pos != rhs.pos;
    }
  };

  using const_iterator = HashTableIterator<const Entry>;
  using iterator = HashTableIterator<Entry>;

  HighsHashTable() { makeEmptyTable(128); }
  HighsHashTable(u32 minCapacity) {
    u32 initCapacity =
        1u << (u32)std::ceil(std::log2(std::max(128.0, 8 * minCapacity / 7.0)));
    makeEmptyTable(initCapacity);
  }

  iterator end() {
    u32 capacity = tableSizeMask + 1;
    return iterator{metadata.get() + capacity, metadata.get() + capacity,
                    entries.get() + capacity};
  };

  const_iterator end() const {
    u32 capacity = tableSizeMask + 1;
    return const_iterator{metadata.get() + capacity, metadata.get() + capacity,
                          entries.get() + capacity};
  };

  const_iterator begin() const {
    if (numElements == 0) return end();
    u32 capacity = tableSizeMask + 1;
    const_iterator iter{metadata.get(), metadata.get() + capacity,
                        entries.get() + capacity};
    if (!occupied(metadata[0])) ++iter;

    return iter;
  };

  iterator begin() {
    if (numElements == 0) return end();
    u32 capacity = tableSizeMask + 1;
    iterator iter{metadata.get(), metadata.get() + capacity,
                  entries.get() + capacity};
    if (!occupied(metadata[0])) ++iter;

    return iter;
  };

 private:
  static constexpr u8 toMetadata(u32 hash) { return (hash & 0x7Fu) | 0x80u; }

  static constexpr u32 maxDistance() { return 127; }

  void makeEmptyTable(u32 capacity) {
    tableSizeMask = capacity - 1;
    numElements = 0;

    metadata = decltype(metadata)(new u8[capacity]{});
    entries =
        decltype(entries)((Entry*)::operator new(sizeof(Entry) * capacity));
  }

  bool occupied(u8 meta) const { return meta & 0x80; }

  u32 distanceFromIdealSlot(u32 pos) const {
    // we store 7 bits of the hash in the metadata. Assuming a decent
    // hashfunction it is practically never happening that an item travels more
    // then 127 slots from its ideal position, therefore, we can compute the
    // distance from the ideal position just as it would normally be done
    // assuming there is at most one overflow. Consider using 3 bits which gives
    // values from 0 to 7. When an item is at a position with lower bits 7 and
    // is placed 3 positions after its ideal position, the lower bits of the
    // hash value will overflow and yield the value 2. With the assumption that
    // an item never cycles through one full cycle of the range 0 to 7, its
    // position would never be placed in a position with lower bits 7 other than
    // its ideal position. This allows us to compute the distance from its ideal
    // position by simply ignoring an overflow. In our case the correct answer
    // would be 3, but we get (2 - 7)=-5. This, however, is the correct result 3
    // when promoting to an unsigned value and looking at the lower 3 bits.

    return ((pos - metadata[pos])) & 0x7f;
  }

  void growTable() {
    decltype(entries) oldEntries = std::move(entries);
    decltype(metadata) oldMetadata = std::move(metadata);
    u32 oldCapactiy = tableSizeMask + 1;

    makeEmptyTable(2 * oldCapactiy);

    for (u32 i = 0; i != oldCapactiy; ++i)
      if (occupied(oldMetadata[i])) insert(std::move(oldEntries.get()[i]));
  }

  void shrinkTable() {
    decltype(entries) oldEntries = std::move(entries);
    decltype(metadata) oldMetadata = std::move(metadata);
    u32 oldCapactiy = tableSizeMask + 1;

    makeEmptyTable(oldCapactiy / 2);

    for (u32 i = 0; i != oldCapactiy; ++i)
      if (occupied(oldMetadata[i])) insert(std::move(oldEntries.get()[i]));
  }

  bool findPosition(const KeyType& key, u8& meta, u32& startPos, u32& maxPos,
                    u32& pos) const {
    u32 hash = HighsHashHelpers::hash(key);
    startPos = hash & tableSizeMask;
    maxPos = (startPos + maxDistance()) & tableSizeMask;
    meta = toMetadata(hash);

    const Entry* entryArray = entries.get();
    pos = startPos;
    do {
      if (!occupied(metadata[pos])) return false;
      if (metadata[pos] == meta &&
          HighsHashHelpers::equal(key, entryArray[pos].key()))
        return true;

      u32 currentDistance = (pos - startPos) & tableSizeMask;
      if (currentDistance > distanceFromIdealSlot(pos)) return false;

      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    return false;
  }

 public:
  void clear() { makeEmptyTable(128); }

  const ValueType* find(const KeyType& key) const {
    u32 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return &(entries.get()[pos].value());

    return nullptr;
  }

  ValueType* find(const KeyType& key) {
    u32 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return &(entries.get()[pos].value());

    return nullptr;
  }

  ValueType& operator[](const KeyType& key) {
    Entry* entryArray = entries.get();
    u32 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return entryArray[pos].value();

    if (numElements == ((tableSizeMask + 1) * 7) / 8 || pos == maxPos) {
      growTable();
      return (*this)[key];
    }

    using std::swap;
    ValueType& insertLocation = entryArray[pos].value();
    Entry entry(key, ValueType());
    ++numElements;

    do {
      if (!occupied(metadata[pos])) {
        metadata[pos] = meta;
        new (&entryArray[pos]) Entry{std::move(entry)};
        return insertLocation;
      }

      u32 currentDistance = (pos - startPos) & tableSizeMask;
      u32 distanceOfCurrentOccupant = distanceFromIdealSlot(pos);
      if (currentDistance > distanceOfCurrentOccupant) {
        // steal the position
        swap(entry, entryArray[pos]);
        swap(meta, metadata[pos]);

        startPos = (pos - distanceOfCurrentOccupant) & tableSizeMask;
        maxPos = (startPos + maxDistance()) & tableSizeMask;
      }
      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    growTable();
    insert(std::move(entry));
    return (*this)[key];
  }

  template <typename... Args>
  bool insert(Args&&... args) {
    Entry entry(std::forward<Args>(args)...);

    u32 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(entry.key(), meta, startPos, maxPos, pos)) return false;

    if (numElements == ((tableSizeMask + 1) * 7) / 8 || pos == maxPos) {
      growTable();
      return insert(std::move(entry));
    }

    using std::swap;
    Entry* entryArray = entries.get();
    ++numElements;

    do {
      if (!occupied(metadata[pos])) {
        metadata[pos] = meta;
        new (&entryArray[pos]) Entry{std::move(entry)};
        return true;
      }

      u32 currentDistance = (pos - startPos) & tableSizeMask;
      u32 distanceOfCurrentOccupant = distanceFromIdealSlot(pos);
      if (currentDistance > distanceOfCurrentOccupant) {
        // steal the position
        swap(entry, entryArray[pos]);
        swap(meta, metadata[pos]);

        startPos = (pos - distanceOfCurrentOccupant) & tableSizeMask;
        maxPos = (startPos + maxDistance()) & tableSizeMask;
      }
      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    growTable();
    insert(std::move(entry));
    return true;
  }

  bool erase(const KeyType& key) {
    u32 pos, startPos, maxPos;
    u8 meta;
    if (!findPosition(key, meta, startPos, maxPos, pos)) return false;
    // delete element at position pos
    Entry* entryArray = entries.get();
    entryArray[pos].~Entry();
    metadata[pos] = 0;

    // retain at least a quarter of slots occupied, otherwise shrink the table
    // if its not at its minimum size already
    --numElements;
    u32 capacity = tableSizeMask + 1;
    if (capacity != 128 && numElements < capacity / 4) {
      shrinkTable();
      return true;
    }

    // shift elements after pos backwards
    while (true) {
      u32 shift = (pos + 1) & tableSizeMask;
      if (!occupied(metadata[shift])) return true;

      u32 dist = distanceFromIdealSlot(shift);
      if (dist == 0) return true;

      entryArray[pos] = std::move(entryArray[shift]);
      metadata[pos] = metadata[shift];
      metadata[shift] = 0;
      pos = shift;
    }
  }

  size_t size() const { return numElements; }
};

#endif
