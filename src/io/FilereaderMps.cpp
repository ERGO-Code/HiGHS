/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderMps.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "FilereaderMps.h"

FilereaderRetcode FilereaderMps::readModelFromFile(const char filename,
                                                   HighsLp& model) {
  // todo

  // Which parser
  // will be (options.getValue("parser") == MpsParser::new)

  // Initialize arrays

  // call MPSParser::loadProblem(arrays of HighsLp object)

  return FilereaderRetcode::OKAY;
}
