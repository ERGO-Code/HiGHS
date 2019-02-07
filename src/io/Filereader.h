/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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

enum class HighsInputStatus {
  OK,
  FileNotFound,
  ErrorMatrixDimensions,
  ErrorMatrixIndices,
  ErrorMatrixStart,
  ErrorMatrixValue,
  ErrorColBounds,
  ErrorRowBounds,
  ErrorObjective
};

enum class FilereaderRetcode {
  OKAY = 0,
  FILENOTFOUND = 1,
  PARSERERROR = 2,
  NOT_IMPLEMENTED = 3
};

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsLp& model) = 0;
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsModel& model) = 0;
  virtual FilereaderRetcode writeModelToFile(const char* filename,
                                             HighsLp& model) = 0;
  static Filereader* getFilereader(const char* filename);

  virtual ~Filereader(){};
};

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status);

HighsInputStatus checkLp(const HighsLp& lp);

#endif
