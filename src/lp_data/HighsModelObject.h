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

#include "HighsLp.h"
#include "HighsTimer.h"
#include "HighsRandom.h"
#include "HModel.h"
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
  HighsModelObject(HighsLp& lp) : lp_(lp) {}

  HighsLp& lp_;

  HighsLp solver_lp_;
  HighsSimplexInfo simplex_info_;
  HighsSolution solution_;
  HighsRanging ranging_;
  HighsBasis basis_;
  HighsScale scale_;
  HMatrix matrix_;
  HFactor factor_;
  HighsTimer timer_;
  HighsRandom random_;


  bool reportModelOperationsClock = false;

  //
  // Basis consists of basicIndex, nonbasicFlag and nonbasicMove. To
  // have them means that they correspond to a consistent basis
  // logically, but B is not necessarily nonsingular.
  int has_basis;
  int has_some_solver_lp_data;

  BasisInfo basis_info_;

  // the vector below either contains one vector or zero.
  std::vector<HModel> hmodel_;

};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_
