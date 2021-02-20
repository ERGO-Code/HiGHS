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
  if (!io.output_flag || (io.logging_file == NULL && !io.log_to_console))
    return;
  // highsOutputUser should not be passed HighsMessageType::VERBOSE
  assert(type != HighsMessageType::VERBOSE);
  va_list argptr;
  va_start(argptr, format);
  if (logmsgcb == NULL) {
    fprintf(io.logging_file, "%-9s", HighsMessageTypeTag[(int)type]);
    vfprintf(io.logging_file, format, argptr);
    if (io.log_to_console) {
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
  if (io.logging_file == NULL || !io.output_dev) return;
  va_list argptr;
  va_start(argptr, format);
  vfprintf(io.logging_file, format, argptr);
  if (io.log_to_console) {
    va_start(argptr, format);
    vfprintf(stdout, format, argptr);
  }
}

void HighsPrintMessage(FILE* pass_output, const int pass_message_level,
                       const int level, const char* format, ...) {
  if (pass_output == NULL) return;
  if (pass_message_level & level) {
    va_list argptr;
    va_start(argptr, format);
    if (printmsgcb == NULL)
      vfprintf(pass_output, format, argptr);
    else {
      int len;
      len = vsnprintf(msgbuffer, sizeof(msgbuffer), format, argptr);
      if (len >= (int)sizeof(msgbuffer)) {
        // Output was truncated: for now just ensure string is null-terminated
        msgbuffer[sizeof(msgbuffer) - 1] = '\0';
      }
      printmsgcb(level, msgbuffer, msgcb_data);
    }
    va_end(argptr);
  }
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
