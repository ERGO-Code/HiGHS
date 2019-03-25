#ifndef HIGHS_QP_H
#define HIGHS_QP_H

#include "HighsLp.h"

class HighsQp : public HighsLp {
 public:
  std::vector<int> Hstart_;
  std::vector<int> Hindex_;
  std::vector<double> Hvalue_;
};

#endif