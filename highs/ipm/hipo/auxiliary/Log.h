#ifndef HIPO_LOG_H
#define HIPO_LOG_H

#include <sstream>

#include "io/HighsIO.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

// Base class for logging.
// Use print for normal logging, printw for warnings, printe for errors.
// dev_level = 0 does not print any debug information.
// dev_level = 1 prints debug level Info.
// dev_level = 2 prints debug levels Info, Detailed.
// dev_level = 3 prints debug levels Info, Detailed, Verbose.

namespace hipo {

class Log {
  Int dev_level_ = 0;

 public:
  Log(Int level = 0);

  virtual void print(std::stringstream& ss) const;
  virtual void printw(std::stringstream& ss) const;
  virtual void printe(std::stringstream& ss) const;
  virtual void print(const char* c) const;
  virtual void printw(const char* c) const;
  virtual void printe(const char* c) const;

  virtual void printDevInfo(std::stringstream& ss) const;
  virtual void printDevDetailed(std::stringstream& ss) const;
  virtual void printDevVerbose(std::stringstream& ss) const;
  virtual void printDevInfo(const char* c) const;
  virtual void printDevDetailed(const char* c) const;
  virtual void printDevVerbose(const char* c) const;
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