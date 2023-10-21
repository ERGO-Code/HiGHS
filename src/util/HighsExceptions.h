#pragma once
#include <stdexcept>

class DataStackOverflow : public std::runtime_error {
 public:
  explicit DataStackOverflow(const std::string& msg)
      : std::runtime_error(msg) {}
};
