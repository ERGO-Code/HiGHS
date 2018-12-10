/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HRanging.h
 * @brief Compute and check ranging information for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HRANGING_H_
#define SIMPLEX_HRANGING_H_

#include <vector>

#include "HModel.h"
//#using namespace std;

/**
 * @brief Compute and check ranging information for HiGHS
 */
class HRanging {
 public:
  /**
   * @brief Compute ranging information
   */
  int computeData(
		  HModel *model  //!< Instance of HModel class for which ranging
		                 //!< data are to be generated
  );

  /**
   * @brief Check ranging information
   */
  int checkData(
      HModel *model  //!< Instance of HModel class for which ranging
		     //!< data are to be checked
  );

 private:
  void checkDataZeroMlFg(
			 HModel *model   //!< Instance of HModel class for which model flags are to be zeroed
			 );

  void checkDataSolve(
		      HModel *model,  //!< Instance of HModel class for which test
		                      //!< problem is to be solved when checking ranging
		                      //!< data
		      bool rp         //!< Report when model is tested
		      );
};

#endif /* SIMPLEX_HRANGING_H_ */
