#include "HpmLog.h"

namespace hipo {

const HighsLogOptions* Log::log_options_ = nullptr;

void Log::setOptions(const HighsLogOptions& log_options) {
  log_options_ = &log_options;
}

bool Log::debug(Int level) {
  // Return true if level agrees with log_dev_level

  if (log_options_ && log_options_->log_dev_level) {
    if (*log_options_->log_dev_level == kHighsLogDevLevelInfo)
      return level == 1;

    if (*log_options_->log_dev_level == kHighsLogDevLevelDetailed)
      return level == 1 || level == 2;

    if (*log_options_->log_dev_level == kHighsLogDevLevelVerbose)
      return level == 1 || level == 2 || level == 3;
  }
  return false;
}

void Log::print(std::stringstream& ss) {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kInfo, "%s", ss.str().c_str());
}
void Log::printw(std::stringstream& ss) {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kWarning, "%s", ss.str().c_str());
}
void Log::printe(std::stringstream& ss) {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kError, "%s", ss.str().c_str());
}

void Log::printDevInfo(std::stringstream& ss) {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kInfo, "%s", ss.str().c_str());
}
void Log::printDevDetailed(std::stringstream& ss) {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kDetailed, "%s", ss.str().c_str());
}
void Log::printDevVerbose(std::stringstream& ss) {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kVerbose, "%s", ss.str().c_str());
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
