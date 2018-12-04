/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSensitivity.h
 * @brief Sensitivity and ranging information for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSENSITIVITY_H_
#define SIMPLEX_HSENSITIVITY_H_

#include <vector>

#include "HModel.h"
//#using namespace std;

/**
 * @brief Sensitivity and ranging information for HiGHS
 */
class HSensitivity {
 public:
  /**
   * @brief Compute sensitivity and ranging information
   */
  int getSensitivityData(
      HModel *model  //!< Instance of HModel class for which sensitivity and
                     //!< ranging data are to be generated
  );

  /**
   * @brief Check sensitivity and ranging information
   */
  int checkSensitivityData(
      HModel *model  //!< Instance of HModel class for which sensitivity and
                     //!< ranging data are to be checked
  );

 private:
  void checkSensitivityZeroMlFg(
      HModel *model   //!< Instance of HModel class for which model
		      //!< flags are to be zeroed
				 );

  void checkSensitivityDataSolve(
      HModel *model,  //!< Instance of HModel class for which test
                      //!< problem is to be solved when checking
                      //!< sensitivity and ranging data
      bool rp //!< Report when model is tested
				 );
};

#endif /* SIMPLEX_HSENSITIVITY_H_ */
