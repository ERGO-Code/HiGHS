#ifndef LP_DATA_HIGHS_MODEL_OBJECT_H_
#define LP_DATA_HIGHS_MODEL_OBJECT_H_

#include "HighsLp.h"
#include "HModel.h"
// include Sensitivity(or Ranging) header

struct BasisInfo {
  std::vector<int> basis_index;
  std::vector<int> nonbasic_flag;
};

// Class to communicate data between the simplex solver and the class
// Highs below. Sensitivity data structure would be added here. Only
// include essential data.
class HighsModelObject {
public:
  HighsModelObject(HighsLp& lp) : lp_(lp) {}

  HighsLp& lp_;
  HighsSolution solution_;
  BasisInfo basis_info_;
  HighsBasis basis_;

  // the vector below either contains one vector or zero.
  std::vector<HModel> hmodel_;

  int* getBaseIndex() { return &basis_.basicIndex_[0]; }
  int* getNonbasicFlag() { return &basis_.nonbasicFlag_[0]; }
  int* getNonbasicMove() { return &basis_.nonbasicMove_[0]; }
};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_
