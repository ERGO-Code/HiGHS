/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderEms.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#ifndef IO_FILEREADER_EMS_H_
#define IO_FILEREADER_EMS_H_

#include <list>

#include "Filereader.h"
#include "HighsIO.h"  // For messages.

class FilereaderEms : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const char* filename, HighsLp& model);
  FilereaderRetcode writeModelToFile(const char* filename, HighsLp& model);
  FilereaderRetcode readModelFromFile(const char* filename, HighsModel& model);
};

#endif

