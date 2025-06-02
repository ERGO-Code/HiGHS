#ifndef HIGHSPM_LOG_H
#define HIGHSPM_LOG_H

#include <sstream>

#include "ipm/hpm/auxiliary/IntConfig.h"
#include "io/HighsIO.h"

// Interface to Highs logging.
// Call Log::setOptions to set the HighsLogOptions.
// If log_options_ is null, nothing is printed.

namespace highspm {

class Log {
  static const HighsLogOptions* log_options_;

  // Private ctor and dtor
  Log();
  ~Log() = default;

 public:
  static void setOptions(const HighsLogOptions& log_options);
  static bool debug(Int level);

  // Logging normal, warning and error, with stream
  static void print(std::stringstream& ss);
  static void printw(std::stringstream& ss);
  static void printe(std::stringstream& ss);

  // Logging normal, warning and error, formatted
  template <typename... Args>
  static void printf(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kInfo, format, args...);
  }
  template <typename... Args>
  static void printw(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kWarning, format, args...);
  }
  template <typename... Args>
  static void printe(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kError, format, args...);
  }

  // Dev logging info, detailed and verbose, with stream
  static void printDevInfo(std::stringstream& ss);
  static void printDevDetailed(std::stringstream& ss);
  static void printDevVerbose(std::stringstream& ss);

  // Dev logging info, detailed and verbose, formatted
  template <typename... Args>
  static void printDevInfo(const char* format, Args... args) {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kInfo, format, args...);
  }
  template <typename... Args>
  static void printDevDetailed(const char* format, Args... args) {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kDetailed, format, args...);
  }
  template <typename... Args>
  static void printDevVerbose(const char* format, Args... args) {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kVerbose, format, args...);
  }
};

// Functions to print using streams, taken from IPX.
std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield);
std::string integer(Int i, Int width = 0);
inline std::string sci(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::scientific);
}
inline std::string fix(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::fixed);
}
template <typename T>
std::string textline(const T& text) {
  std::ostringstream s;
  s.setf(std::ios_base::left);
  s.width(32);
  s << text;
  return s.str();
}

}  // namespace highspm

#endif