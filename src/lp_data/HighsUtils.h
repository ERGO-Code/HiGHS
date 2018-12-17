/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSUTILS_H_
#define LP_DATA_HIGHSUTILS_H_

#include "HConfig.h"
#include "HighsLp.h"
#include "HConst.h"
//#include "HModel.h"
//#include "HighsModelObject.h"

/**
 * @brief Logical check of double being +Infinity
 */
bool highs_isInfinity(
		      double val //!< Value being tested against +Infinity
		      );

const int initial_random_mw = 1985;
const int initial_random_mz = 2012;

/**
 * @brief Class for HiGHS utilities
 */
class HighsUtils {
 public:

  /**
   * @brief Initialisations
   */
  HighsUtils() {
  /**
   * @brief Initialise the two seeds to default values
   */
    random_mw = initial_random_mw;
    random_mz = initial_random_mz;
  }

  /**
   * @brief (Re-)initialise the random number generator
   */
  void initialiseRandom() {
    random_mw = initial_random_mw;
    random_mz = initial_random_mz;
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

#ifdef HiGHSDEV
  // Analyse the values of a vector, assessing how many are in each power of ten
  void util_anVecV(const char* message, int vecDim, std::vector<double>& vec, bool anVLs);
#endif

  // Methods taking HighsLp as an argument
  /**
   * @brief Report the data of an LP
   */
  void reportLp(
		HighsLp lp //!< LP whose data are to be reported
		);
  /**
   * @brief Report the brief data of an LP 
   */
  void reportLpBrief(
		     HighsLp lp //!< LP whose data are to be reported
		     );
  /**
   * @brief Report the data of an LP
   */
  void reportLpDimensions(
			  HighsLp lp //!< LP whose data are to be reported
			  );
  /**
   * @brief Report the data of an LP
   */
  void reportLpObjSense(
			HighsLp lp //!< LP whose data are to be reported
			);
  /**
   * @brief Report the data of an LP
   */
  void reportLpColVec(
		      HighsLp lp //!< LP whose data are to be reported
		      );
  /**
   * @brief Report the data of an LP
   */
  void reportLpRowVec(
		      HighsLp lp //!< LP whose data are to be reported
		      );
  /**
   * @brief Report the data of an LP
   */
  void reportLpColMtx(
		      HighsLp lp //!< LP whose data are to be reported
		      );


  /*
void reportLpSolution(
			HighsModelObject highs_model //!< Model object whose LP solution is to be reported
			);
  */

#ifdef HiGHSDEV
  
  void util_anMl(HighsLp lp, const char* message);
  void util_anMlBd(const char* message, int numBd, std::vector<double>& lower, std::vector<double>& upper);
#endif
 private:
  unsigned random_mw;
  unsigned random_mz;

};
#endif // LP_DATA_HIGHSUTILS_H_
