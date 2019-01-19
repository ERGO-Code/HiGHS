/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderEms.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "FilereaderEms.h"
#include "HConst.h"

FilereaderRetcode FilereaderMps::readModelFromFile(const char* filename,
                                                   HighsLp& model) {
  return FilereaderRetcode::OKAY;
}


FilereaderRetcode FilereaderMps::writeModelToFile(const char* filename,
                                                  HighsLp& model) {
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderMps::readModelFromFile(const char* filename,
                                                   HighsModel& model) {
  return FilereaderRetcode::OKAY;
}