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

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

FILE* logfile = stdout;
FILE* output = stdout;
int message_level = ML_MINIMAL;
void (*printmsgcb)(int, const char*, void*) = NULL;
void (*logmsgcb)(HighsMessageType, const char*, void*) = NULL;
void* msgcb_data = NULL;

char msgbuffer[65536];

void HighsPrintMessage(int level, const char* format, ...) {
  if (message_level & level) {
    va_list argptr;
    va_start(argptr, format);
    if (printmsgcb == NULL)
      vfprintf(output, format, argptr);
    else {
      int len;
      len = vsnprintf(msgbuffer, sizeof(msgbuffer), format, argptr);
      if (len >= (int)sizeof(msgbuffer)) {
        /* output was truncated: for now just ensure string is null-terminated
         */
        msgbuffer[sizeof(msgbuffer) - 1] = '\0';
      }
      printmsgcb(level, msgbuffer, msgcb_data);
    }
    va_end(argptr);
  }
}

void HighsLogMessage(HighsMessageType type, const char* format, ...) {
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  va_list argptr;
  va_start(argptr, format);

  if (logmsgcb == NULL) {
    fprintf(logfile, "%02d:%02d:%02d [%-7s] ", timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec, HighsMessageTypeTag[(int)type]);
    vfprintf(logfile, format, argptr);
    fprintf(logfile, "\n");
  } else {
    int len;
    len = snprintf(msgbuffer, sizeof(msgbuffer), "%02d:%02d:%02d [%-7s] ",
                   timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec,
                   HighsMessageTypeTag[(int)type]);
    if (len < (int)sizeof(msgbuffer))
      len +=
          vsnprintf(msgbuffer + len, sizeof(msgbuffer) - len, format, argptr);
    if (len < (int)sizeof(msgbuffer) - 1) {
      msgbuffer[len] = '\n';
      ++len;
      msgbuffer[len] = '\0';
    } else
      msgbuffer[sizeof(msgbuffer) - 1] = '\0';
    logmsgcb(type, msgbuffer, msgcb_data);
  }

  va_end(argptr);
}

void HighsSetLogfile(FILE* lf) { logfile = lf; }

void HighsSetOutput(FILE* op) { output = op; }

void HighsSetMessagelevel(int level) { message_level = level; }

void HighsSetMessageCallback(
    void (*printmsgcb_)(int level, const char* msg, void* msgcb_data),
    void (*logmsgcb_)(HighsMessageType type, const char* msg, void* msgcb_data),
    void* msgcb_data_) {
  printmsgcb = printmsgcb_;
  logmsgcb = logmsgcb_;
  msgcb_data = msgcb_data_;
}

void HighsSetIO(HighsOptions& options) {
  logfile = options.logfile;
  output = options.output;
  message_level = options.message_level;
  printmsgcb = options.printmsgcb;
  logmsgcb = options.logmsgcb;
  msgcb_data = options.msgcb_data;
}
