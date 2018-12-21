#include "HighsModel.h"

HighsModel::~HighsModel() {
  while (this->variables.size() > 0) {
    HighsVar* variable;
    variable = this->variables.front();
    this->variables.pop_front();

    // find all coefficients corresponding to the variable
    VarConsCoefsMap::iterator it =
        this->variableConstraintCoefficientMap.find(variable);
    if (it != this->variableConstraintCoefficientMap.end()) {
      std::list<HighsLinearConsCoef*>* coefficients = it->second;

      while(coefficients->size() > 0) {
        HighsLinearConsCoef* coef = coefficients->front();
        coefficients->pop_front();
        // remove coefficient from constraint
        
        CoefConsMap::iterator iter = this->coefficientConstraintMap.find(coef);
        assert(iter != this->coefficientConstraintMap.end());
        HighsLinearCons* constraint = iter->second;
        VarConsCoefMap::iterator iterator = constraint->linearCoefs.find(variable);
        assert(iterator != constraint->linearCoefs.end());
        constraint->linearCoefs.erase(iterator);
        this->coefficientConstraintMap.erase(iter);

       

        delete coef;
      }
      VarConsMap::iterator iter = this->variableConstraintMap.find(variable);
      if(iter != variableConstraintMap.end()) {
        std::list<HighsLinearCons*>* conslist = iter->second;
        assert(conslist->empty());
        this->variableConstraintMap.erase(iter);
        delete conslist;
      }

      this->variableConstraintCoefficientMap.erase(it);
      delete coefficients;
    }

    delete variable;
  }

  while (this->linearConstraints.size() > 0) {
    HighsLinearCons* constraint;
    constraint = this->linearConstraints.front();
    this->linearConstraints.pop_front();

    delete constraint;
  }
}

#pragma region HighsVar
HighsVar::HighsVar(const char* name, double lo, double hi, double obj,
                   HighsVarType type) {
  // create a copy of the name
  if (name != NULL) {
    int namelen = strlen(name);
    this->name = new char[namelen + 1];
    strcpy(this->name, name);
  } else {
    this->name = NULL;
  }

  // copy all remaining data
  this->lowerBound = lo;
  this->upperBound = hi;
  this->obj = obj;
  this->type = type;
}

HighsVar::~HighsVar() {
  if (this->name != NULL) {
    delete[] this->name;
  }
}

#pragma endregion

#pragma region HighsCons

HighsCons::HighsCons(const char* name, double lo, double hi) {
  // create a copy of the name
  if (name != NULL) {
    int namelen = strlen(name);
    this->name = new char[namelen + 1];
    strcpy(this->name, name);
  } else {
    this->name = NULL;
  }

  // copy all remaining data
  this->lowerBound = lo;
  this->upperBound = hi;
}

HighsCons::~HighsCons() {
  if (this->name != NULL) {
    delete[] this->name;
  }
}

#pragma endregion

#pragma region HighsLinearCons

HighsLinearCons::HighsLinearCons(const char* name, double lo, double hi)
    : HighsCons(name, lo, hi) {}

HighsLinearCons::~HighsLinearCons() {}

#pragma endregion

#pragma region HighsLinearConsCoef

HighsLinearConsCoef::HighsLinearConsCoef(HighsVar* var, double coef) {
  this->var = var;
  this->coef = coef;
}

HighsLinearConsCoef::~HighsLinearConsCoef() {}

#pragma endregion

#pragma region HighsModel
#pragma region HighsModel Variables

void HighsModel::HighsCreateVar(const char* name, double lo, double hi,
                                double obj, HighsVarType type, HighsVar** var) {
  if (name != NULL) {
    // make sure name is available
    VarMap::iterator it = this->variableMap.find(name);
    if (it != this->variableMap.end()) {
      // name already in use
      // TODO: Error Message
      return;
    }
  }

  // create the new variable and add it to the model
  *var = new HighsVar(name, lo, hi, obj, type);
  this->variables.push_back(*var);
  if (name != NULL) {
    this->variableMap.insert(VarMap::value_type((*var)->name, *var));
  }
}

void HighsModel::HighsCreateVar(const char* name, HighsVar** var) {
  this->HighsCreateVar(name, 0.0, __DBL_MAX__, 0.0, HighsVarType::CONT, var);
}

void HighsModel::HighsGetOrCreateVarByName(const char* name, HighsVar** var) {
  this->HighsGetVarByName(name, var);
  if (*var == NULL) {
    this->HighsCreateVar(name, var);
  }
}

void HighsModel::HighsCreateVar(HighsVar** var) {
  this->HighsCreateVar(NULL, var);
}

void HighsModel::HighsGetVarByName(const char* name, HighsVar** var) {
  VarMap::iterator it = this->variableMap.find(name);
  if (it != this->variableMap.end()) {
    *var = it->second;
  } else {
    // variable not found
    // TODO: Error Message
    *var = NULL;
  }
}

void HighsModel::HighsRemoveVar(HighsVar* var) {
  // check that variable is no longer used in any constraints
  // TODO

  // remove variable from map
  VarMap::iterator it = this->variableMap.find(var->name);
  if (it == this->variableMap.end()) {
    // variable no longer in Model?
    // TODO: Error Message
    return;
  }
  this->variableMap.erase(var->name);

  // remove variable from list
  // TODO
  return;
}

#pragma endregion

#pragma region HighsModel Constraints

void HighsModel::HighsCreateLinearCons(const char* name, double lo, double hi,
                                       HighsLinearCons** cons) {
  if (name != NULL) {
    // make sure name is available
    ConsMap::iterator it = this->constraintMap.find(name);
    if (it != this->constraintMap.end()) {
      // name already in use
      // TODO: Error Message
      return;
    }
  }

  // create the new constraint and add it to the model
  *cons = new HighsLinearCons(name, lo, hi);
  this->linearConstraints.push_back(*cons);
  if (name != NULL) {
    this->constraintMap.insert(ConsMap::value_type((*cons)->name, *cons));
  }
}

void HighsModel::HighsCreateLinearCons(const char* name,
                                       HighsLinearCons** cons) {
  this->HighsCreateLinearCons(name, -__DBL_MAX__, __DBL_MAX__, cons);
}

void HighsModel::HighsCreateLinearCons(HighsLinearCons** cons) {
  this->HighsCreateLinearCons(NULL, cons);
}

void HighsModel::HighsGetLinearConsByName(const char* name,
                                          HighsLinearCons** cons) {}

void HighsModel::HighsDestroyLinearCons() {}

#pragma endregion

#pragma region HighsModel Coefficients

void HighsModel::HighsCreateLinearConsCoef(HighsVar* var, double coef,
                                           HighsLinearConsCoef** consCoef) {
  *consCoef = new HighsLinearConsCoef(var, coef);
  VarConsCoefsMap::iterator it =
      this->variableConstraintCoefficientMap.find(var);
  if (it != this->variableConstraintCoefficientMap.end()) {
    it->second->push_back(*consCoef);
    ;
  } else {
    std::list<HighsLinearConsCoef*>* coefList =
        new std::list<HighsLinearConsCoef*>;
    coefList->push_back(*consCoef);
    this->variableConstraintCoefficientMap.insert(
        VarConsCoefsMap::value_type(var, coefList));
  }
}

void HighsModel::HighsAddLinearConsCoefToCons(HighsLinearCons* cons,
                                              HighsLinearConsCoef* coef) {
  VarConsCoefMap::iterator it = cons->linearCoefs.find(coef->var);
  if (it != cons->linearCoefs.end()) {
    // constraint already has a coefficient for this variable
  } else {
    coefficientConstraintMap.insert(CoefConsMap::value_type(coef, cons));
    cons->linearCoefs.insert(VarConsCoefMap::value_type(coef->var, coef));
    VarConsMap::iterator it = this->variableConstraintMap.find(coef->var);
    if (it != this->variableConstraintMap.end()) {
      it->second->push_back(cons);
    } else {
      std::list<HighsLinearCons*>* consList = new std::list<HighsLinearCons*>;
      consList->push_back(cons);
      this->variableConstraintMap.insert(
          VarConsMap::value_type(coef->var, consList));
    }
  }
}

#pragma endregion

void HighsModel::HighsBuildTechnicalModel(HighsLp* lp) {
  lp->numCol_ = this->variables.size();
  lp->numRow_ = this->linearConstraints.size();

  // determine order of variables
  HighsVar** variables = new HighsVar*[lp->numCol_];
  for (int i = 0; i < lp->numCol_; i++) {
    HighsVar* front = this->variables.front();
    this->variables.pop_front();
    this->variables.push_back(front);
    variables[i] = front;
    lp->colCost_.push_back(this->objSense * front->obj);
    lp->colLower_.push_back(front->lowerBound);
    lp->colUpper_.push_back(front->upperBound);
  }

  // determine order of constraints
  HighsLinearCons** constraints = new HighsLinearCons*[lp->numRow_];
  for (int i = 0; i < lp->numRow_; i++) {
    HighsLinearCons* front = this->linearConstraints.front();
    this->linearConstraints.pop_front();
    this->linearConstraints.push_back(front);
    constraints[i] = front;
    lp->rowLower_.push_back(front->lowerBound);
    lp->rowUpper_.push_back(front->upperBound);
  }

  // handle constraints
  lp->Astart_.clear();
  lp->Astart_.push_back(0);
  for (int var = 0; var < lp->numCol_; var++) {
    VarConsCoefsMap::iterator iter =
        this->variableConstraintCoefficientMap.find(variables[var]);
    if (iter != this->variableConstraintCoefficientMap.end()) {
      std::list<HighsLinearConsCoef*>* coefs = iter->second;
      int numberOfCoefficients = coefs->size();

      lp->Astart_.push_back(lp->Astart_[var] + numberOfCoefficients);

      for (int coef = 0; coef < numberOfCoefficients; coef++) {
        HighsLinearConsCoef* front = coefs->front();
        coefs->pop_front();
        coefs->push_back(front);
        lp->Avalue_.push_back(front->coef);
        CoefConsMap::iterator it = this->coefficientConstraintMap.find(front);
        if (it != this->coefficientConstraintMap.end()) {
          // find index of constraint
          HighsCons* currentCons = it->second;
          for (int cons = 0; cons < lp->numRow_; cons++) {
            if (constraints[cons] == currentCons) {
              lp->Aindex_.push_back(cons);
              break;
            }
          }
        } else {
          // ERROR
        }
      }
    }
  }

  delete[] variables;
  delete[] constraints;
}

#pragma endregion
