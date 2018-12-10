/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "FilereaderLp.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const char* filename, HighsLp& model) {
  return FilereaderRetcode::PARSERERROR;
}

FilereaderRetcode FilereaderLp::writeModelToFile(const char* filename, HighsLp& model) {
  return FilereaderRetcode::PARSERERROR;
}