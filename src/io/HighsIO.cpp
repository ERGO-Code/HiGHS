/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.cpp
 * @brief IO methods for HiGHS - currently just print/log messages
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HighsIO.h"

#include <stdarg.h>
#include <stdio.h>
#include <time.h>

void HighsPrintMessage(unsigned int level, const char* format, ...) {
  FILE* output = stdout; // TODO: read from options
  int messageLevel = 1+2+4; // TODO: read from options

  if (messageLevel & level) {
    va_list argptr;
    va_start(argptr, format);
    vfprintf(output, format, argptr);
    va_end(argptr);
  }
}

void HighsLogMessage(HighsMessageType type, const char* format, ...) {
  FILE* logfile = stdout; // TODO: read from options
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  va_list argptr;
  va_start(argptr, format);

  fprintf(logfile, "%02d:%02d:%02d [%s] ", timeinfo->tm_hour, timeinfo->tm_min,
          timeinfo->tm_sec, HighsMessageTypeTag[type]);
  vfprintf(logfile, format, argptr);
  fprintf(logfile, "\n");

  va_end(argptr);
}
