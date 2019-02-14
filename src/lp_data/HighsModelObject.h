/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef LP_DATA_HIGHS_MODEL_OBJECT_H_
#define LP_DATA_HIGHS_MODEL_OBJECT_H_

#include "lp_data/HighsLp.h"
#include "util/HighsTimer.h"
#include "util/HighsRandom.h"
#include "simplex/HMatrix.h"
#include "simplex/HFactor.h"
// include Sensitivity(or Ranging) header

struct BasisInfo {
  std::vector<int> basis_index;
  std::vector<int> nonbasic_flag;
  std::vector<int> nonbasic_move;
};

// Class to communicate data between the simplex solver and the class
// Highs below. Sensitivity data structure would be added here. Only
// include essential data.
class HighsModelObject {
public:
 HighsModelObject(HighsLp& lp, HighsTimer& timer) : lp_(lp),
    timer_(timer) {}

  HighsLp& lp_;
  HighsBasis basis_;
  HighsSolution solution_;
  HighsTimer& timer_;

  HighsLp solver_lp_;
  HighsSimplexInfo simplex_info_;
  //  HighsRanging ranging_;
  HighsScale scale_;
  HMatrix matrix_;
  HFactor factor_;
  HighsRandom random_;

  bool reportModelOperationsClock = false;

  BasisInfo basis_info_;

};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_
