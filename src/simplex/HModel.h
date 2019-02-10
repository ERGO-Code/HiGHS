/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.h
 * @brief LP model representation and management for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HMODEL_H_
#define SIMPLEX_HMODEL_H_

#include "HFactor.h"
#include "HMatrix.h"
#include "HighsLp.h"
#include "HighsTimer.h" //For timer_
#include "HighsRandom.h"

//class HVector;

#include <sstream>
#include <string>
#include <vector>

// After removing HTimer.h add the following

class HModel {
 public:
  HModel();
  // Methods which load whole models, initialise the basis then
  // allocate and populate (where possible) work* arrays and
  // allocate basis* arrays
  int load_fromToy(const char* filename);

  // Methods which initialise the basis then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  void rp_basis();
  void replaceFromNonbasic();

  // ???? Housekeeping done from here down ????
#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void changeUpdate(int updateMethod);
#endif

  int writeToMPS(const char* filename);
  void util_getBasicIndexNonbasicFlag(
				      vector<int> &XbasicIndex,
				      vector<int> &XnonbasicFlag
				      );
#ifdef HiGHSDEV
  void util_anMlLargeCo(HighsLp lp, const char* message);
#endif

 public:
  
  // The scaled model
  HighsLp *solver_lp_;
  HMatrix *matrix_;
  HFactor *factor_;
  HighsSimplexInfo *simplex_info_;
  HighsBasis *basis_;
  HighsScale *scale_;
  HighsRanging *ranging_;
  HighsRandom *random_;
  HighsTimer *timer_;

};

#endif /* SIMPLEX_HMODEL_H_ */
