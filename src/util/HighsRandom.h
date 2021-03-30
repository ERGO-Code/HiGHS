/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsRandom.h
 * @brief Random number generators for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSRANDOM_H_
#define UTIL_HIGHSRANDOM_H_

#include "util/HighsHash.h"

/**
 * @brief Class for HiGHS random number generators
 */
class HighsRandom {
 public:
  /**
   * @brief Initialisations
   */
  HighsRandom(unsigned seed = 0) { initialise(seed); }

  /**
   * @brief (Re-)initialise the random number generator
   */
  void initialise(unsigned seed = 0) {
    state = seed;
    do {
      state = (HighsHashHelpers::pair_hash<0>(state, 0) >> 32);
      state =
          (state << 32) ^ (HighsHashHelpers::pair_hash<1>(state, seed) >> 32);
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
  int integer() {
    // use 31 bits of the 64 bit result
    advance();
    return HighsHashHelpers::pair_hash<0>(state, state >> 32) >> 33;
  }

  /**
   * @brief Return a random integer between [0,sup)
   */
  int integer(int sup) { return fractionOrZero() * sup; }

  /**
   * @brief Return a random integer between [min,sup)
   */
  int integer(int min, int sup) { return min + integer(sup - min); }

  /**
   * @brief Return a random fraction - real in (0, 1)
   */
  double fraction() {
    advance();
    // 52 bit output is in interval [0,2^52-1]
    uint64_t output =
        (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> (64 - 26))
            << 26 |
        HighsHashHelpers::pair_hash<1>(state, state >> 32) >> (64 - 26);
    // compute (1+output) / (2^52+1) which is strictly between 0 and 1
    return (1 + output) * 2.2204460492503125e-16;
  }

  /**
   * @brief Return a random fraction - real in [0, 1)
   */
  double fractionOrZero() {
    advance();
    // 52 bit result is in interval [0,2^52-1]
    uint64_t output =
        (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> (64 - 26))
            << 26 |
        HighsHashHelpers::pair_hash<1>(state, state >> 32) >> (64 - 26);
    // compute output / (2^52) which is in the half-open interval [0,1)
    return output * 2.22044604925031308e-16;
  }

  /**
   * @brief Return a random real value in the half-open interval [a,b)
   */
  double real(double a, double b) { return a + (b - a) * fractionOrZero(); }

  /**
   * @brief Return a random bit
   */
  bool bit() {
    advance();
    return state >> 63;
  }

  /**
   * @brief shuffle the given data array
   */
  template <typename T>
  void shuffle(T* data, int N) {
    for (int i = N; i > 1; --i) {
      int pos = integer(i);
      std::swap(data[pos], data[i - 1]);
    }
  }

 private:
  uint64_t state;
};

#endif
