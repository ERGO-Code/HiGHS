#pragma once
#include <stdexcept>

class DataStackOverflow : public std::runtime_error {
public:
    explicit DataStackOverflow(const std::string& msg)
        : std::runtime_error(msg) {}
};

class ProblemTooLargeException : public std::runtime_error {
public:
    explicit ProblemTooLargeException(const std::string& msg)
        : std::runtime_error(msg) {}
};


class PresolveTooLarge : public std::runtime_error {
public:
    explicit PresolveTooLarge(const std::string& msg)
        : std::runtime_error(msg) {}
};
