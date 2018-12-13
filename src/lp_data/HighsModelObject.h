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
  HighsSimplexInfo simplex_;
  HighsSolution solution_;
  HighsRanging ranging_;
  HighsBasis basis_;
  HighsScale scale_;

  BasisInfo basis_info_;

  // the vector below either contains one vector or zero.
  std::vector<HModel> hmodel_;

  int* getBaseIndex() { return &basis_.basicIndex_[0]; }
  int* getNonbasicFlag() { return &basis_.nonbasicFlag_[0]; }
  int* getNonbasicMove() { return &basis_.nonbasicMove_[0]; }
  double* getWorkCost() { return &simplex_.workCost_[0]; }
  double* getWorkDual() { return &simplex_.workDual_[0]; }
  double* getWorkShift() { return &simplex_.workShift_[0]; }
  double* getWorkLower() { return &simplex_.workLower_[0]; }
  double* getWorkUpper() { return &simplex_.workUpper_[0]; }
  double* getWorkRange() { return &simplex_.workRange_[0]; }
  double* getWorkValue() { return &simplex_.workValue_[0]; }
  double* getBaseLower() { return &simplex_.baseLower_[0]; }
  double* getBaseUpper() { return &simplex_.baseUpper_[0]; }
  double* getBaseValue() { return &simplex_.baseValue_[0]; }
};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_
