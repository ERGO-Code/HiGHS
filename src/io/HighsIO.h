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

class HighsOptions;

/**
 * @brief IO methods for HiGHS - currently just print/log messages
 */
enum class HighsMessageType { INFO = 1, DETAILED, VERBOSE, WARNING, ERROR };
const char* const HighsMessageTypeTag[] = {"", "", "", "WARNING: ", "ERROR: "};
enum OutputDevLevel {
		     OUTPUT_DEV_MIN = 0,
		     OUTPUT_DEV_NONE = OUTPUT_DEV_MIN, // 0
		     OUTPUT_DEV_INFO,                  // 1
		     OUTPUT_DEV_DETAILED,              // 2
		     OUTPUT_DEV_VERBOSE,               // 3
		     OUTPUT_DEV_MAX = OUTPUT_DEV_VERBOSE
};

struct HighsIo {
  FILE* logging_file;
  bool* output_flag;
  bool* log_to_console;
  int* output_dev;
};

/**
 * @brief For _single-line_ user logging with message type notification
 */
// Printing format: must contain exactly one "\n" at end of format
void highsOutputUser(const HighsIo& io, const HighsMessageType type,
                     const char* format, ...);

/**
 * @brief For development logging
 */
void highsOutputDev(const HighsIo& io, const HighsMessageType type,
                    const char* format, ...);

/*
 * @brief sets the callbacks used to print output and and log
 *
 * Set to NULL to reset to default, which is to print to logfile and output file
 */
void highsSetMessageCallback(
    void (*printmsgcb_)(int level, const char* msg, void* msgcb_data),
    void (*logmsgcb_)(HighsMessageType type, const char* msg, void* msgcb_data),
    void* msgcb_data_);

/*
 * @brief sets callbacks from options
 */
void highsSetMessageCallback(HighsOptions& options  //!< the options
);

void highsReportIoOptions(const HighsIo& io);

#endif
