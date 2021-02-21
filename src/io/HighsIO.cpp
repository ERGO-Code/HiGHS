/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.cpp
 * @brief IO methods for HiGHS - currently just print/log messages
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HighsIO.h"

#include <cstdarg>
#include <cstdio>
//#include <ctime>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

void (*printmsgcb)(int, const char*, void*) = NULL;
void (*logmsgcb)(HighsLogType, const char*, void*) = NULL;
void* msgcb_data = NULL;

char msgbuffer[65536];

void highsLogUser(const HighsLogOptions& log_options, const HighsLogType type,
                     const char* format, ...) {
  if (!*log_options.output_flag || (log_options.log_file_stream == NULL && !*log_options.log_to_console))
    return;
  // highsLogUser should not be passed HighsLogType::DETAILED or HighsLogType::VERBOSE
  assert(type != HighsLogType::DETAILED);
  assert(type != HighsLogType::VERBOSE);
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    fprintf(log_options.log_file_stream, "%-9s", HighsLogTypeTag[(int)type]);
    vfprintf(log_options.log_file_stream, format, argptr);
    if (*log_options.log_to_console) {
      fprintf(log_options.log_file_stream, "%-9s", HighsLogTypeTag[(int)type]);
      va_start(argptr, format);
      vfprintf(stdout, format, argptr);
    }
  } else {
    int len;
    len = snprintf(msgbuffer, sizeof(msgbuffer), "%-9s",
                   HighsLogTypeTag[(int)type]);
    if (len < (int)sizeof(msgbuffer))
      len +=
          vsnprintf(msgbuffer + len, sizeof(msgbuffer) - len, format, argptr);
    if (len >= (int)sizeof(msgbuffer)) {
      // Output was truncated: for now just ensure string is null-terminated
      msgbuffer[sizeof(msgbuffer) - 1] = '\0';
    }
    logmsgcb(type, msgbuffer, msgcb_data);
  }
  va_end(argptr);
}

void highsLogDev(const HighsLogOptions& log_options, const HighsLogType type,
                    const char* format, ...) {
  if (log_options.log_file_stream == NULL || !*log_options.log_dev_level) return;
  // Always report HighsLogType INFO, WARNING or ERROR 
  // Report HighsLogType DETAILED if *log_options.log_dev_level >= LOG_DEV_LEVEL_DETAILED
  // Report HighsLogType VERBOSE if *log_options.log_dev_level >= LOG_DEV_LEVEL_VERBOSE
  if (type == HighsLogType::DETAILED && *log_options.log_dev_level < LOG_DEV_LEVEL_DETAILED) return;
  if (type == HighsLogType::VERBOSE && *log_options.log_dev_level < LOG_DEV_LEVEL_VERBOSE) return;
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    vfprintf(log_options.log_file_stream, format, argptr);
    if (*log_options.log_to_console) {
      va_start(argptr, format);
      vfprintf(stdout, format, argptr);
    }
  } else {
    int len;
    len = vsnprintf(msgbuffer, sizeof(msgbuffer), format, argptr);
    if (len >= (int)sizeof(msgbuffer)) {
      // Output was truncated: for now just ensure string is null-terminated
      msgbuffer[sizeof(msgbuffer) - 1] = '\0';
    }
    logmsgcb(type, msgbuffer, msgcb_data);
  }
  va_end(argptr);
}

void highsSetLogCallback(
    void (*printmsgcb_)(int level, const char* msg, void* msgcb_data),
    void (*logmsgcb_)(HighsLogType type, const char* msg, void* msgcb_data),
    void* msgcb_data_) {
  printmsgcb = printmsgcb_;
  logmsgcb = logmsgcb_;
  msgcb_data = msgcb_data_;
}

void highsSetLogCallback(HighsOptions& options) {
  printmsgcb = options.printmsgcb;
  logmsgcb = options.logmsgcb;
  msgcb_data = options.msgcb_data;
}

void highsReportLogOptions(const HighsLogOptions& log_options) {
  printf("\nHighs IO settings\n");
  if (log_options.log_file_stream == NULL) {
    printf("   log_file_stream = NULL\n");
  } else {
    printf("   log_file_stream = Not NULL\n");
  }
  printf("   output_flag = %d\n", *log_options.output_flag);
  printf("   log_to_console = %d\n", *log_options.log_to_console);
  printf("   log_dev_level = %d\n\n", *log_options.log_dev_level);
}
