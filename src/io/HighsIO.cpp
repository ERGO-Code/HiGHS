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

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

void (*printmsgcb)(int, const char*, void*) = NULL;
void (*logmsgcb)(HighsLogType, const char*, void*) = NULL;
void* msgcb_data = NULL;

char msgbuffer[65536];

void highsLogUser(const HighsLogOptions& log_options, const HighsLogType type,
                  const char* format, ...) {
  if (!*log_options.output_flag ||
      (log_options.log_file_stream == NULL && !*log_options.log_to_console))
    return;
  // highsLogUser should not be passed HighsLogType::DETAILED or
  // HighsLogType::VERBOSE
  assert(type != HighsLogType::DETAILED);
  assert(type != HighsLogType::VERBOSE);
  const bool prefix =
      type == HighsLogType::WARNING || type == HighsLogType::ERROR;
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    if (log_options.log_file_stream != NULL) {
      // Write to log file stream
      if (prefix)
        fprintf(log_options.log_file_stream, "%-9s",
                HighsLogTypeTag[(int)type]);
      vfprintf(log_options.log_file_stream, format, argptr);
      va_start(argptr, format);
    }
    if (*log_options.log_to_console && log_options.log_file_stream != stdout) {
      // Write to stdout unless log file stream is stdout
      if (prefix) fprintf(stdout, "%-9s", HighsLogTypeTag[(int)type]);
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
  if (log_options.log_file_stream == NULL || !*log_options.log_dev_level)
    return;
  // Always report HighsLogType INFO, WARNING or ERROR
  //
  // Report HighsLogType DETAILED if *log_options.log_dev_level >=
  // LOG_DEV_LEVEL_DETAILED
  //
  // Report HighsLogType VERBOSE if *log_options.log_dev_level >=
  // LOG_DEV_LEVEL_VERBOSE
  if (type == HighsLogType::DETAILED &&
      *log_options.log_dev_level < LOG_DEV_LEVEL_DETAILED)
    return;
  if (type == HighsLogType::VERBOSE &&
      *log_options.log_dev_level < LOG_DEV_LEVEL_VERBOSE)
    return;
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    if (log_options.log_file_stream != NULL) {
      // Write to log file stream
      vfprintf(log_options.log_file_stream, format, argptr);
      va_start(argptr, format);
    }
    if (*log_options.log_to_console && log_options.log_file_stream != stdout) {
      // Write to stdout unless log file stream is stdout
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

void highsSetLogCallback(void (*printmsgcb_)(int level, const char* msg,
                                             void* msgcb_data),
                         void (*logmsgcb_)(HighsLogType type, const char* msg,
                                           void* msgcb_data),
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

void highsSetLogOptions(HighsLogOptions& log_options, const bool* output_flag_,
                        FILE* log_file_stream_, const bool* log_to_console_,
                        const int* log_dev_level_) {
  bool output_flag = true;
  if (output_flag_ != NULL) output_flag = *output_flag_;
  log_options.output_flag = &output_flag;
  log_options.log_file_stream = log_file_stream_;
  bool log_to_console = true;
  if (log_to_console_ != NULL) log_to_console = *log_to_console_;
  log_options.log_to_console = &log_to_console;

  int log_dev_level = LOG_DEV_LEVEL_NONE;
  if (log_dev_level_ != NULL) log_dev_level = *log_dev_level_;
  log_options.log_dev_level = &log_dev_level;
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

std::string highsFormatToString(const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  int len = vsnprintf(msgbuffer, sizeof(msgbuffer), format, argptr);
  if (len >= (int)sizeof(msgbuffer)) {
    // Output was truncated: for now just ensure string is null-terminated
    msgbuffer[sizeof(msgbuffer) - 1] = '\0';
  }
  va_end(argptr);
  std::string local_string(msgbuffer);
  return local_string;
}
