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
void (*logmsgcb)(HighsMessageType, const char*, void*) = NULL;
void* msgcb_data = NULL;

char msgbuffer[65536];

void highsOutputUser(const HighsIo& io, const HighsMessageType type,
                     const char* format, ...) {
  if (!*io.output_flag || (io.logging_file == NULL && !*io.log_to_console))
    return;
  // highsOutputUser should not be passed HighsMessageType::DETAILED or HighsMessageType::VERBOSE
  assert(type != HighsMessageType::DETAILED);
  assert(type != HighsMessageType::VERBOSE);
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    fprintf(io.logging_file, "%-9s", HighsMessageTypeTag[(int)type]);
    vfprintf(io.logging_file, format, argptr);
    if (*io.log_to_console) {
      fprintf(io.logging_file, "%-9s", HighsMessageTypeTag[(int)type]);
      va_start(argptr, format);
      vfprintf(stdout, format, argptr);
    }
  } else {
    int len;
    len = snprintf(msgbuffer, sizeof(msgbuffer), "%-9s",
                   HighsMessageTypeTag[(int)type]);
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

void highsOutputDev(const HighsIo& io, const HighsMessageType type,
                    const char* format, ...) {
  if (io.logging_file == NULL || !*io.output_dev) return;
  // Always report HighsMessageType INFO, WARNING or ERROR 
  // Report HighsMessageType DETAILED if *io.output_dev >= OUTPUT_DEV_DETAILED
  // Report HighsMessageType VERBOSE if *io.output_dev >= OUTPUT_DEV_VERBOSE
  if (type == HighsMessageType::DETAILED && *io.output_dev < OUTPUT_DEV_DETAILED) return;
  if (type == HighsMessageType::VERBOSE && *io.output_dev < OUTPUT_DEV_VERBOSE) return;
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    vfprintf(io.logging_file, format, argptr);
    if (*io.log_to_console) {
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

void HighsSetMessageCallback(
    void (*printmsgcb_)(int level, const char* msg, void* msgcb_data),
    void (*logmsgcb_)(HighsMessageType type, const char* msg, void* msgcb_data),
    void* msgcb_data_) {
  printmsgcb = printmsgcb_;
  logmsgcb = logmsgcb_;
  msgcb_data = msgcb_data_;
}

void HighsSetIO(HighsOptions& options) {
  printmsgcb = options.printmsgcb;
  logmsgcb = options.logmsgcb;
  msgcb_data = options.msgcb_data;
}

void highsReportIo(const HighsIo& io) {
  printf("\nHighs IO settings\n");
  if (io.logging_file == NULL) {
    printf("   logging_file = NULL\n");
  } else {
    printf("   logging_file = Not NULL\n");
  }
  
  bool output_flag = *io.output_flag;
  bool log_to_console = *io.log_to_console;
  int output_dev = *io.output_dev;
  printf("   output_flag = %d\n", output_flag);
  printf("   log_to_console = %d\n", log_to_console);
  printf("   output_dev = %d\n\n", output_dev);
}
