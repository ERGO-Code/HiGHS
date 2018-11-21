#ifndef HIGHS_MODEL_H
#define HIGHS_MODEL_H

#include <string.h>
#include <list>
#include <map>

#include "HighsLp.h"

enum class HighsVarType { CONT, SEMI, BIN, GEN };

class HighsVar {
 public:
  char* name;
  double lowerBound;
  double upperBound;
  double obj;
  HighsVarType type;

  HighsVar(const char* name, double lo, double hi, double obj,
           HighsVarType type);
  ~HighsVar();
};

class HighsCons {
 public:
  char* name;
  double lowerBound;
  double upperBound;

  HighsCons(const char* name, double lo, double hi);
  ~HighsCons();
};

class HighsLinearConsCoef {
 public:
  HighsVar* var;
  double coef;

  HighsLinearConsCoef(HighsVar* var, double coef);
  ~HighsLinearConsCoef();  
};

struct char_cmp {
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
};

typedef std::map<const char*, HighsVar*, char_cmp> VarMap;
typedef std::map<HighsVar*, HighsLinearConsCoef*> VarConsCoefMap;
typedef std::map<HighsVar*, std::list<HighsLinearConsCoef*>*> VarConsCoefsMap;

class HighsLinearCons : public HighsCons {
 public:
  VarConsCoefMap linearCoefs;

  HighsLinearCons(const char* name, double lo, double hi);
  ~HighsLinearCons();
};

typedef std::map<const char*, HighsLinearCons*, char_cmp> ConsMap;
typedef std::map<HighsVar*, std::list<HighsLinearCons*>*> VarConsMap;

class HighsModel {
 public:
  double objOffset;

  void HighsCreateVar(const char* name, double lo, double hi, double obj,
                      HighsVarType type, HighsVar** var);
  void HighsCreateVar(const char* name, HighsVar** var);
  void HighsCreateVar(HighsVar** var);
  void HighsGetVarByName(const char* name, HighsVar** var);
  void HighsRemoveVar(HighsVar* var);

  void HighsCreateLinearConsCoef(HighsVar* var, double coef, HighsLinearConsCoef** consCoef);
  void HighsAddLinearConsCoefToCons(HighsLinearCons* cons, HighsLinearConsCoef* coef);
  //void HighsDestroyLinearConsCoef();

  void HighsCreateLinearCons(const char* name, double lo, double hi,
                             HighsLinearCons** cons);
  void HighsCreateLinearCons(const char* name, HighsLinearCons** cons);
  void HighsCreateLinearCons(HighsLinearCons** cons);
  void HighsGetLinearConsByName(const char* name, HighsLinearCons** cons);
  void HighsDestroyLinearCons();

  // conversion from/to technical Lp representation
  HighsModel(HighsLp* lp);
  void HighsBuildModel(HighsLp* lp);

 private:
  // major data structures
  std::list<HighsLinearCons*> linearConstraints;
  std::list<HighsVar*> variables;

  // helper data structures
  VarMap variableMap;
  ConsMap constraintMap;
  VarConsMap variableConstraintMap;
  VarConsCoefsMap variableConstraintCoefficientMap;
};

#endif