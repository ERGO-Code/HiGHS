/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.h
 * @brief IO methods for HiGHS - currently just print/log messages
 */
#ifndef HIGHS_IO_H
#define HIGHS_IO_H

#include <iostream>

#include "util/HighsInt.h"

class HighsOptions;

/**
 * @brief IO methods for HiGHS - currently just print/log messages
 */
enum class HighsLogType { kInfo = 1, kDetailed, kVerbose, kWarning, kError };
const char* const HighsLogTypeTag[] = {"", "",          "",
                                       "", "WARNING: ", "ERROR:   "};
enum LogDevLevel {
  kHighsLogDevLevelMin = 0,
  kHighsLogDevLevelNone = kHighsLogDevLevelMin,  // 0
  kHighsLogDevLevelInfo,                         // 1
  kHighsLogDevLevelDetailed,                     // 2
  kHighsLogDevLevelVerbose,                      // 3
  kHighsLogDevLevelMax = kHighsLogDevLevelVerbose
};

struct HighsLogOptions {
  FILE* log_file_stream;
  bool* output_flag;
  bool* log_to_console;
  HighsInt* log_dev_level;
};

/**
 * @brief For _single-line_ user logging with message type notification
 */
// Printing format: must contain exactly one "\n" at end of format
void highsLogUser(const HighsLogOptions& log_options_, const HighsLogType type,
                  const char* format, ...);

/**
 * @brief For development logging
 */
void highsLogDev(const HighsLogOptions& log_options_, const HighsLogType type,
                 const char* format, ...);

/*
 * @brief sets the callbacks used to print output and and log
 *
 * Set to NULL to reset to default, which is to print to logfile and output file
 */
void highsSetLogCallback(void (*printmsgcb_)(HighsInt level, const char* msg,
                                             void* msgcb_data),
                         void (*logmsgcb_)(HighsLogType type, const char* msg,
                                           void* msgcb_data),
                         void* msgcb_data_);

/*
 * @brief sets callbacks from options
 */
void highsSetLogCallback(HighsOptions& options  //!< the options
);

void highsReportLogOptions(const HighsLogOptions& log_options_);

std::string highsFormatToString(const char* format, ...);

const std::string highsBoolToString(const bool b);

#endif
