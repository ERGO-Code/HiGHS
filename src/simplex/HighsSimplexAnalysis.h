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

#include <vector>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HVector.h"

/**
 * @brief Analyse simplex iterations, both for run-time control and data gathering
 */
class HighsSimplexAnalysis {
 public:
  HighsSimplexAnalysis(HighsModelObject& highs_model_object)
      : highs_model_object(highs_model_object) {}

  HighsModelObject& highs_model_object;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_dse_density;

};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
