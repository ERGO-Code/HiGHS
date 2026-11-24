#ifndef HIPO_LOGGER_H
#define HIPO_LOGGER_H

#include <sstream>

#include "io/HighsIO.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

// Class for logging.
// Messages are printed through the HiGHS logging functions, using either the
// log options passed with setOptions, or the log options owned by the Logger,
// which print only user-level messages to console.
// If no log options are in use, nothing is printed.

namespace hipo {

class Logger {
  // Log options owned by the Logger, used to print to console outside of
  // HiGHS. The pointer fields of own_log_options_ point to the own_* members
  // declared before it.
  bool own_output_flag_ = true;
  bool own_log_to_console_ = true;
  HighsInt own_log_dev_level_ = kHighsLogDevLevelNone;
  HighsLogOptions own_log_options_;

  // Log options in use: either &own_log_options_ or options passed with
  // setOptions. nullptr if logging is off.
  const HighsLogOptions* log_options_ = nullptr;

  template <typename... Args>
  void log(bool dev, HighsLogType type, const char* format,
           Args... args) const {
    if (!log_options_) return;
    if (dev)
      highsLogDev(*log_options_, type, format, args...);
    else
      highsLogUser(*log_options_, type, format, args...);
  }

 public:
  // If to_console is true, print to console using the log options owned by
  // the Logger. Otherwise, nothing is printed until setOptions is called.
  Logger(bool to_console = false);

  // Not copyable, since the pointer fields of own_log_options_ point to the
  // Logger's own members.
  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  // Use log_options to print. Pass nullptr to switch logging off.
  void setOptions(const HighsLogOptions* log_options);

  bool debug(Int level) const;

  // Shared loggers that never print, or that print to console.
  static const Logger& silent();
  static const Logger& console();

  template <typename... Args>
  void print(const char* format, Args... args) const {
    log(false, HighsLogType::kInfo, format, args...);
  }

  template <typename... Args>
  void printw(const char* format, Args... args) const {
    log(false, HighsLogType::kWarning, format, args...);
  }

  template <typename... Args>
  void printe(const char* format, Args... args) const {
    log(false, HighsLogType::kError, format, args...);
  }

  template <typename... Args>
  void printInfo(const char* format, Args... args) const {
    log(true, HighsLogType::kInfo, format, args...);
  }

  template <typename... Args>
  void printDetailed(const char* format, Args... args) const {
    log(true, HighsLogType::kDetailed, format, args...);
  }

  template <typename... Args>
  void printVerbose(const char* format, Args... args) const {
    log(true, HighsLogType::kVerbose, format, args...);
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

}  // namespace hipo
#endif