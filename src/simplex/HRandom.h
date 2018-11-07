/**@file HRandom.h
 * @brief Random number generator for HiGHS
 * @author Qi Huangfu
 */
#ifndef SIMPLEX_HRANDOM_H_
#define SIMPLEX_HRANDOM_H_

/**
 * @brief Class for the random number generator for HiGHS
 */
class HRandom {
 public:
  /**
   * @brief Initialise the two seeds
   */
  HRandom() {
    random_mw = 1985;
    random_mz = 2012;
  }

  /**
   * @brief Return a random integer between 0 and 2147483647
   */
  int intRandom() {
    random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
    random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
    unsigned result = (random_mz << 16) + random_mw;
    return result >> 1;
  }

  /**
   * @brief Return a random real in (0, 1)
   */
  double dblRandom() {
    random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
    random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
    unsigned result = (random_mz << 16) + random_mw;
    double returnValue = (result + 1.0) * 2.328306435454494e-10;
    return returnValue;
  }

 private:
  unsigned random_mw;
  unsigned random_mz;
};

#endif /* SIMPLEX_HRANDOM_H_ */
