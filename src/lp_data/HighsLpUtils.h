/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLpUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSLPUTILS_H_
#define LP_DATA_HIGHSLPUTILS_H_

#include "HConfig.h"
#include "HighsLp.h"

class HighsLp;

// Methods taking HighsLp as an argument
/**
 * @brief Report the data of an LP
 */
void reportLp(
	      HighsLp &lp //!< LP whose data are to be reported
	      );
/**
 * @brief Report the brief data of an LP 
 */
void reportLpBrief(
		   HighsLp &lp //!< LP whose data are to be reported
		   );
/**
 * @brief Report the data of an LP
 */
void reportLpDimensions(
			HighsLp &lp //!< LP whose data are to be reported
			);
/**
 * @brief Report the data of an LP
 */
void reportLpObjSense(
		      HighsLp &lp //!< LP whose data are to be reported
		      );
/**
 * @brief Report the data of an LP
 */
void reportLpColVec(
		    HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpRowVec(
		    HighsLp &lp //!< LP whose data are to be reported
		    );
/**
 * @brief Report the data of an LP
 */
void reportLpColMtx(
		    HighsLp &lp //!< LP whose data are to be reported
		    );


/*
  void reportLpSolution(
  HighsModelObject &highs_model //!< Model object whose LP solution is to be reported
  );
*/
#ifdef HiGHSDEV
// Analyse the values of a vector, assessing how many are in each power of ten
void util_anMl(HighsLp &lp, const char* message);
void util_anMlBd(const char* message, int numBd, std::vector<double>& lower, std::vector<double>& upper);
#endif
#endif // LP_DATA_HIGHSLPUTILS_H_
