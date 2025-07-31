#include "Log.h"

namespace hipo {

Log::Log(Int level) : dev_level_{level} {}

void Log::print(std::stringstream& ss) const {
  printf("%s", ss.str().c_str());
}
void Log::printw(std::stringstream& ss) const {
  printf("WARNING: %s", ss.str().c_str());
}
void Log::printe(std::stringstream& ss) const {
  printf("ERROR: %s", ss.str().c_str());
}
void Log::print(const char* c) const { printf("%s", c); }
void Log::printw(const char* c) const { printf("WARNING: %s", c); }
void Log::printe(const char* c) const { printf("ERROR: %s", c); }
void Log::printDevInfo(std::stringstream& ss) const {
  if (dev_level_ >= 1) printf("%s", ss.str().c_str());
}
void Log::printDevDetailed(std::stringstream& ss) const {
  if (dev_level_ >= 2) printf("%s", ss.str().c_str());
}
void Log::printDevVerbose(std::stringstream& ss) const {
  if (dev_level_ >= 3) printf("%s", ss.str().c_str());
}
void Log::printDevInfo(const char* c) const {
  if (dev_level_ >= 1) printf("%s", c);
}
void Log::printDevDetailed(const char* c) const {
  if (dev_level_ >= 2) printf("%s", c);
}
void Log::printDevVerbose(const char* c) const {
  if (dev_level_ >= 3) printf("%s", c);
}

std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield) {
  std::ostringstream s;
  s.precision(prec);
  s.width(width);
  s.setf(floatfield, std::ios_base::floatfield);
  s << d;
  return s.str();
}

std::string integer(Int i, Int width) {
  std::ostringstream s;
  s.width(width);
  s << i;
  return s.str();
}

}  // namespace hipo
