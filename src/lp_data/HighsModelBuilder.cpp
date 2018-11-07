/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModelBuilder.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HighsModelBuilder.h"

void HighsModelBuilder::HighsAddVar(HighsVar& var) {
  this->variables.push_back(var);
}

void HighsModelBuilder::HighsAddCons(HighsCons& cons) {
  this->constraints.push_back(cons);
}

void HighsModelBuilder::HighsCreateLp(HighsLp& lp) {}
