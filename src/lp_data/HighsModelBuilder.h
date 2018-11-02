#ifndef HIGHS_MODEL_BUILDER_H
#define HIGHS_MODEL_BUILDER_H

#include "HighsLp.h"

#include <list>

class HighsVar {
  char* name;
  double lowerBound;
  double upperBound;
};

class HighsLinearConsCoef {
  HighsVar* var;
  double coef;
};

class QuadraticConsCoef {
  HighsVar var1;
  HighsVar var2;
  double coef;
};

class HighsCons {
  char* name;
  double lowerBound;
  double upperBound;
};

class HighsLinearCons : public HighsCons {
  std::list<HighsLinearConsCoef> linearCoefs;
};

class HighsQuadraticCons : public HighsLinearCons {
  std::list<HighsQuadraticCons> quadraticCoefs;
};

class HighsModelBuilder {
 public:
  void HighsAddVar(HighsVar& var);
  void HighsAddCons(HighsCons& cons);
  void HighsCreateLp(HighsLp& lp);

 private:
  std::list<HighsCons> constraints;
  std::list<HighsVar> variables;
  HighsCons objective;
};

#endif