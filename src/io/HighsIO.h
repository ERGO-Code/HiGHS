/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.h
 * @brief IO methods for HiGHS - currently just print/log messages
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef HIGHS_IO_H
#define HIGHS_IO_H

#include <iostream>

#include "util/HighsInt.h"

class HighsOptions;

/**
 * @brief IO methods for HiGHS - currently just print/log messages
 */
enum class HighsLogType { INFO = 1, DETAILED, VERBOSE, WARNING, ERROR };
const char* const HighsLogTypeTag[] = {"", "",          "",
                                       "", "WARNING: ", "ERROR:   "};
enum LogDevLevel {
  LOG_DEV_LEVEL_MIN = 0,
  LOG_DEV_LEVEL_NONE = LOG_DEV_LEVEL_MIN,  // 0
  LOG_DEV_LEVEL_INFO,                      // 1
  LOG_DEV_LEVEL_DETAILED,                  // 2
  LOG_DEV_LEVEL_VERBOSE,                   // 3
  LOG_DEV_LEVEL_MAX = LOG_DEV_LEVEL_VERBOSE
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
void highsSetLogCallback(void (*printmsgcb_)(int level, const char* msg,
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
