#include "HighsIO.h"

#include <stdarg.h>
#include <stdio.h>
#include <time.h>

void HighsPrintMessage(HighsMessageType type, const char* format, ...) {
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  va_list argptr;
  va_start(argptr, format);

  // TODO: read from options what files the output should be written to
  // TODO: read from options whether timestamp should be printed
  if (type == HighsMessageType::INFO) {
    fprintf(stdout, "%02d:%02d:%02d [INFO] ", timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);
    vfprintf(stdout, format, argptr);
  } else if (type == HighsMessageType::DEBUG) {
#ifdef HiGHSDEV
    fprintf(stdout, "%02d:%02d:%02d [DEBUG] ", timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);
    vfprintf(stdout, format, argptr);
#endif
  } else if (type == HighsMessageType::WARNING) {
    fprintf(stderr, "%02d:%02d:%02d [WARNING] ", timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);
    vfprintf(stderr, format, argptr);
  } else if (type == HighsMessageType::ERROR) {
    fprintf(stderr, "%02d:%02d:%02d [ERROR] ", timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);
    vfprintf(stderr, format, argptr);
  }

  va_end(argptr);
}