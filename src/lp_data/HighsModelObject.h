#ifndef LP_DATA_HIGHS_MODEL_OBJECT_H_
#define LP_DATA_HIGHS_MODEL_OBJECT_H_

#include "HighsLp.h"
#include "HighsTimer.h"
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

  int problemStatus;
  string modelName;
  
  HighsLp lp_scaled_;
  HighsSimplexInfo simplex_;
  HighsSolution solution_;
  HighsRanging ranging_;
  HighsBasis basis_;
  HighsScale scale_;
  HMatrix matrix_;
  HFactor factor_;
  HighsTimer timer_;

  BasisInfo basis_info_;

  // the vector below either contains one vector or zero.
  std::vector<HModel> hmodel_;

};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_
