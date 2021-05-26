/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsRandom.h
 * @brief Random number generators for HiGHS
 */
#ifndef UTIL_HIGHSRANDOM_H_
#define UTIL_HIGHSRANDOM_H_

#include <cassert>

#include "util/HighsHash.h"


/**
 * @brief Class for HiGHS random number generators
 */
class HighsRandom {
 private:
  // integer log2 algorithm without floating point arithmetic. It uses an
  // unrolled loop and requires few instructions that can be well optimized.

  int log2i(uint64_t n) {
    int x = 0;

    auto log2Iteration = [&](int p) {
      if (n >= uint64_t{1} << p) {
        x += p;
        n >>= p;
      }
    };

    log2Iteration(32);
    log2Iteration(16);
    log2Iteration(8);
    log2Iteration(4);
    log2Iteration(2);
    log2Iteration(1);

    return x;
  }

  int log2i(uint32_t n) {
    int x = 0;

    auto log2Iteration = [&](int p) {
      if (n >= 1u << p) {
        x += p;
        n >>= p;
      }
    };

    log2Iteration(16);
    log2Iteration(8);
    log2Iteration(4);
    log2Iteration(2);
    log2Iteration(1);

    return x;
  }

  uint32_t drawUniform(uint32_t sup, int nbits) {
    // draw uniformly in interval [0,sup) where nbits is the maximal number of
    // bits the results can have we draw random numbers with nbits many bits
    // until we get one that is in the desired range we first use all available
    // output functions for the same state before we advance the state again. We
    // expect nbits to be at most 32 for this 32 bit version.
    assert(sup <= uint32_t{1} << nbits);
    while (true) {
      advance();
      uint32_t lo = state;
      uint32_t hi = state >> 32;

      uint64_t val = HighsHashHelpers::pair_hash<0>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<1>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<2>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<3>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<4>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<5>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<6>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<7>(lo, hi) >> (64 - nbits);
      if (val < sup) return val;
    }
  }

  uint64_t drawUniform(uint64_t sup, int nbits) {
    // 64 bit version for drawUniform. If result fits in 32 bits use 32 bit
    // version above so that it is only called when we actually need more than
    // 32 bits
    if (nbits <= 32) return drawUniform(uint32_t(sup), nbits);

    assert(sup <= uint64_t{1} << nbits);
    while (true) {
      advance();
      uint32_t lo = state;
      uint32_t hi = state >> 32;

      uint64_t val = HighsHashHelpers::pair_hash<0>(lo, hi) >> (64 - nbits) ^
                     HighsHashHelpers::pair_hash<1>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<1>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<2>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<2>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<3>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<3>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<4>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<4>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<5>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<5>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<6>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<6>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<7>(lo, hi) >> 32;
      if (val < sup) return val;

      val = HighsHashHelpers::pair_hash<7>(lo, hi) >> (64 - nbits) ^
            HighsHashHelpers::pair_hash<0>(lo, hi) >> 32;
      if (val < sup) return val;
    }
  }

 public:
  /**
   * @brief Initialisations
   */
  HighsRandom(HighsUInt seed = 0) { initialise(seed); }

  /**
   * @brief (Re-)initialise the random number generator
   */
  void initialise(HighsUInt seed = 0) {
    state = seed;
    do {
      state = HighsHashHelpers::pair_hash<0>(state, state >> 32);
      state ^= (HighsHashHelpers::pair_hash<1>(state >> 32, seed) >> 32);
    } while (state == 0);
  }

  void advance() {
    // advance state with simple xorshift the outputs are produced by applying
    // strongly universal hash functions to the state so that the lower
    // order bits are as strong as the upper order bits
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
  }

  /**
   * @brief Return a random integer between 0 and 2147483647
   */
  HighsInt integer() {
    advance();
#ifdef HIGHSINT64
    // use 63 bits of the first hash result and use second hash for lower half
    // of bits
    return (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> 1) ^
           (HighsHashHelpers::pair_hash<1>(state, state >> 32) >> 32);
#else
    // use 31 bits of the 64 bit result
    return HighsHashHelpers::pair_hash<0>(state, state >> 32) >> 33;
#endif
  }

  /**
   * @brief Return a random integer between [0,sup)
   */
  HighsInt integer(HighsInt sup) {  // let overload resolution select the 32bit
                                    // or the 64bit version
    if (sup <= 1) return 0;
    int nbits = log2i(HighsUInt(sup - 1)) + 1;
    return drawUniform(HighsUInt(sup), nbits);
  }

  /**
   * @brief Return a random integer between [min,sup)
   */
  HighsInt integer(HighsInt min, HighsInt sup) {
    return min + integer(sup - min);
  }

  /**
   * @brief Return a random fraction - real in (0, 1)
   */
  double fraction() {
    advance();
    // 52 bit output is in interval [0,2^52-1]
    uint64_t output =
        (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> (64 - 52)) ^
        (HighsHashHelpers::pair_hash<1>(state, state >> 32) >> (64 - 26));
    // compute (1+output) / (2^52+1) which is strictly between 0 and 1
    return (1 + output) * 2.2204460492503125e-16;
  }

  /**
   * @brief Return a random fraction - real in [0, 1]
   */
  double closedFraction() {
    advance();
    // 53 bit result is in interval [0,2^53-1]
    uint64_t output =
        (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> (64 - 53)) ^
        (HighsHashHelpers::pair_hash<1>(state, state >> 32) >> 32);
    // compute output / (2^53-1) in double precision which is in the closed
    // interval [0,1]
    return output * 1.1102230246251566e-16;
  }

  /**
   * @brief Return a random real value in the interval [a,b]
   */
  double real(double a, double b) { return a + (b - a) * closedFraction(); }

  /**
   * @brief Return a random bit
   */
  bool bit() {
    advance();
    return HighsHashHelpers::pair_hash<0>(state, state >> 32) >> 63;
  }

  /**
   * @brief shuffle the given data array
   */
  template <typename T>
  void shuffle(T* data, HighsInt N) {
    for (HighsInt i = N; i > 1; --i) {
      HighsInt pos = integer(i);
      std::swap(data[pos], data[i - 1]);
    }
  }

 private:
  uint64_t state;
};

#endif
