#include "LogHighs.h"

namespace hipo {

void LogHighs::setOptions(const HighsLogOptions& log_options) {
  log_options_ = &log_options;
}

bool LogHighs::debug(Int level) const {
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

void LogHighs::print(std::stringstream& ss) const {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kInfo, "%s", ss.str().c_str());
}
void LogHighs::printw(std::stringstream& ss) const {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kWarning, "%s", ss.str().c_str());
}
void LogHighs::printe(std::stringstream& ss) const {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kError, "%s", ss.str().c_str());
}
void LogHighs::print(const char* c) const {
  if (log_options_) highsLogUser(*log_options_, HighsLogType::kInfo, "%s", c);
}
void LogHighs::printw(const char* c) const {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kWarning, "%s", c);
}
void LogHighs::printe(const char* c) const {
  if (log_options_) highsLogUser(*log_options_, HighsLogType::kError, "%s", c);
}

void LogHighs::printDevInfo(std::stringstream& ss) const {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kInfo, "%s", ss.str().c_str());
}
void LogHighs::printDevDetailed(std::stringstream& ss) const {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kDetailed, "%s", ss.str().c_str());
}
void LogHighs::printDevVerbose(std::stringstream& ss) const {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kVerbose, "%s", ss.str().c_str());
}
void LogHighs::printDevInfo(const char* c) const {
  if (log_options_) highsLogDev(*log_options_, HighsLogType::kInfo, "%s", c);
}
void LogHighs::printDevDetailed(const char* c) const {
  if (log_options_)
    highsLogDev(*log_options_, HighsLogType::kDetailed, "%s", c);
}
void LogHighs::printDevVerbose(const char* c) const {
  if (log_options_) highsLogDev(*log_options_, HighsLogType::kVerbose, "%s", c);
}

}  // namespace hipo
