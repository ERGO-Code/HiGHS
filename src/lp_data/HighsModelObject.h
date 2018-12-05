#ifndef LP_DATA_HIGHS_MODEL_OBJECT_H_
#define LP_DATA_HIGHS_MODEL_OBJECT_H_

#include "HighsLp.h"
// include Sensitivity(or Ranging) header

// Class to communicate data between the simplex solver and the class
// Highs below. Sensitivity data structure would be added here. Only
// include essential data.
class HighsModelObject {
public:
  HighsModelObject(const HighsLp& lp) : lp_(lp) {}

  const HighsLp& lp_;
  HighsSolution solution_;
};

#endif // LP_DATA_HIGHS_MODEL_OBJECT_H_