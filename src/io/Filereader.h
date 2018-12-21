/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/Filereader.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_FILEREADER_H_
#define IO_FILEREADER_H_

#include "HighsLp.h"
#include "HighsModel.h"

enum class FilereaderRetcode { OKAY = 0, FILENOTFOUND = 1, PARSERERROR = 2 };

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsLp& model) = 0;
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsModel& model) = 0;
  virtual FilereaderRetcode writeModelToFile(const char* filename,
                                             HighsLp& model) = 0;
  static Filereader* getFilereader(const char* filename);

  virtual ~Filereader() {};

 private:
  static void readLineFromFile(FILE* file, char* buffer, int buffersize);
};

#endif
