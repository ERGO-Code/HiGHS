#ifndef HIPO_LOG_HIGHS_H
#define HIPO_LOG_HIGHS_H

#include "ipm/hipo/auxiliary/Log.h"

// Class for Highs logging
// Use setOptions to pass the log_options.

namespace hipo {

class LogHighs : public Log {
  const HighsLogOptions* log_options_;

 public:
  void setOptions(const HighsLogOptions& log_options);
  bool debug(Int level) const;

  void print(std::stringstream& ss) const override;
  void printw(std::stringstream& ss) const override;
  void printe(std::stringstream& ss) const override;
  void print(const char* c) const override;
  void printw(const char* c) const override;
  void printe(const char* c) const override;

  void printDevInfo(std::stringstream& ss) const override;
  void printDevDetailed(std::stringstream& ss) const override;
  void printDevVerbose(std::stringstream& ss) const override;
  void printDevInfo(const char* c) const override;
  void printDevDetailed(const char* c) const override;
  void printDevVerbose(const char* c) const override;
};
}  // namespace hipo
#endif