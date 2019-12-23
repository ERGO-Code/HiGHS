/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include <cmath>
#include <cstdio>
#include "simplex/HighsSimplexAnalysis.h"

void HighsSimplexAnalysis::setup(const int numCol_, const int numRow_) {
  
  // Copy Problem size
  numRow = numRow_;
  numCol = numCol_;
  col_aq_density = 0;
  row_ep_density = 0;
  row_ap_density = 0;
  row_DSE_density = 0;
}

void HighsSimplexAnalysis::updateOperationResultDensity(const double local_density, double& density) {
  density = (1 - running_average_multiplier) * density +
    running_average_multiplier * local_density;
}

void HighsSimplexAnalysis::equalDensity(const double density0, const double density1) {
  const double delta_density = std::fabs(density1-density0);
  if (delta_density>1e-15) {
    printf("ERROR: Difference %g in density0 - %g and density1 = %g\n", delta_density, density0, density1);
  }
}
