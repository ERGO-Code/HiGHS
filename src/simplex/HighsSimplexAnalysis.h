/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.h
 * @brief Analyse simplex iterations, both for run-time control and data gathering
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSIMPLEXANALYSIS_H_
#define SIMPLEX_HIGHSSIMPLEXANALYSIS_H_

//#include <vector>

#include "HConfig.h"
//#include "simplex/HVector.h"

/**
 * @brief Analyse simplex iterations, both for run-time control and data gathering
 */
class HighsSimplexAnalysis {
 public:
  void setup(int numCol_,            //!< Number of columns
             int numRow_            //!< Number of rows
	     );

  void updateOperationResultDensity(const double local_density,
				    double& density
				    );

  void equalDensity(const double density0, const double density1);

  int numRow;
  int numCol;
  const double running_average_multiplier = 0.05;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;

 private:
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
