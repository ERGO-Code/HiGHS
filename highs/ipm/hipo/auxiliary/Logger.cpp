#include "Logger.h"

namespace hipo {

Logger::Logger(bool to_console) {
  own_log_options_.log_stream = nullptr;
  own_log_options_.output_flag = &own_output_flag_;
  own_log_options_.log_to_console = &own_log_to_console_;
  own_log_options_.log_dev_level = &own_log_dev_level_;
  if (to_console) log_options_ = &own_log_options_;
}

void Logger::setOptions(const HighsLogOptions* log_options) {
  log_options_ = log_options;
}

bool Logger::debug(Int level) const {
  // Return true if level is between 1 and log_dev_level
  if (!log_options_ || !log_options_->log_dev_level) return false;
  return level >= 1 && level <= *log_options_->log_dev_level;
}

const Logger& Logger::silent() {
  static const Logger logger(false);
  return logger;
}

const Logger& Logger::console() {
  static const Logger logger(true);
  return logger;
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
