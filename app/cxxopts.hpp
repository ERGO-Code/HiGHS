/*

Copyright (c) 2014 - 2021 Jarryd Beck
Copyright (c) 2021 - 2023 Pavel Artemkin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

#ifndef CXXOPTS_HPP_INCLUDED
#define CXXOPTS_HPP_INCLUDED

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef __has_include
# if __has_include(<optional>)
#  include <optional>
#  ifdef __cpp_lib_optional
#   define CXXOPTS_HAS_OPTIONAL
#  endif
# endif
#endif

#ifdef CXXOPTS_USE_UNICODE
# include <unicode/unistr.h>
#endif

#if __cplusplus >= 200809L
# define CXXOPTS_NORETURN [[noreturn]]
#else
# define CXXOPTS_NORETURN
#endif

#if __cplusplus >= 201603L
# define CXXOPTS_NODISCARD [[nodiscard]]
#else
# define CXXOPTS_NODISCARD
#endif

#if __cplusplus >= 202002L
# define CXXOPTS_CONSTEXPR constexpr
#else
# define CXXOPTS_CONSTEXPR
#endif

// Disable exceptions if the specific compiler flags are set.
#if !(defined(__cpp_exceptions) || defined(__EXCEPTIONS) || defined(_CPPUNWIND))
# define CXXOPTS_NO_EXCEPTIONS
#endif

#ifdef CXXOPTS_NO_EXCEPTIONS
# include <iostream>
#endif

#ifndef CXXOPTS_VECTOR_DELIMITER
# define CXXOPTS_VECTOR_DELIMITER ','
#endif

#define CXXOPTS__VERSION_MAJOR 5
#define CXXOPTS__VERSION_MINOR 3
#define CXXOPTS__VERSION_PATCH 1

namespace cxxopts {

static constexpr struct {
  uint8_t major, minor, patch;
} version = {CXXOPTS__VERSION_MAJOR, CXXOPTS__VERSION_MINOR,
             CXXOPTS__VERSION_PATCH};

#ifdef _WIN32
static const std::string LQUOTE("\'");
static const std::string RQUOTE("\'");
#else
static const std::string LQUOTE("‘");
static const std::string RQUOTE("’");
#endif

static constexpr std::size_t OPTION_LONGEST = 30;
static constexpr std::size_t OPTION_DESC_GAP = 2;
static constexpr std::size_t OPTION_TAB_SIZE = 8;

#ifdef CXXOPTS_USE_UNICODE
using cxx_string = icu::UnicodeString;
#else
using cxx_string = std::string;
#endif

} // namespace cxxopts

// when we ask cxxopts to use Unicode, help strings are processed using ICU,
// which results in the correct lengths being computed for strings when they
// are formatted for the help output
// it is necessary to make sure that <unicode/unistr.h> can be found by the
// compiler, and that icu-uc is linked in to the binary.

#ifdef CXXOPTS_USE_UNICODE

namespace cxxopts {

static inline cxx_string to_local_string(std::string s) {
  return icu::UnicodeString::fromUTF8(std::move(s));
}

class unicode_string_iterator {
public:
  using value_type = int32_t;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;
  using pointer = value_type*;
  using reference = value_type&;

  unicode_string_iterator(const icu::UnicodeString* string, int32_t pos)
    : s(string)
    , i(pos) {
  }

  value_type operator*() const {
    return s->char32At(i);
  }

  bool operator==(const unicode_string_iterator& rhs) const {
    return s == rhs.s && i == rhs.i;
  }

  bool operator!=(const unicode_string_iterator& rhs) const {
    return !(*this == rhs);
  }

  unicode_string_iterator& operator++() {
    ++i;
    return *this;
  }

  unicode_string_iterator operator+(int32_t v) {
    return unicode_string_iterator(s, i + v);
  }

private:
  const icu::UnicodeString* s;
  int32_t i;
};

static inline cxx_string& string_append(cxx_string& s, cxx_string a) {
  return s.append(std::move(a));
}

static inline cxx_string& string_append(cxx_string& s,
                                        std::size_t n,
                                        UChar32 c) {
  for (std::size_t i = 0; i != n; ++i) {
    s.append(c);
  }

  return s;
}

template <typename Iterator>
static inline cxx_string& string_append(cxx_string& s,
                                        Iterator begin,
                                        Iterator end) {
  while (begin != end) {
    s.append(*begin);
    ++begin;
  }

  return s;
}

static inline std::size_t string_length(const cxx_string& s) {
  return s.length();
}

static inline std::string to_utf8_string(const cxx_string& s) {
  std::string result;
  s.toUTF8String(result);

  return result;
}

static inline bool empty(const cxx_string& s) {
  return s.isEmpty();
}

} // namespace cxxopts

namespace std {

cxxopts::unicode_string_iterator begin(const icu::UnicodeString& s) {
  return cxxopts::unicode_string_iterator(&s, 0);
}

cxxopts::unicode_string_iterator end(const icu::UnicodeString& s) {
  return cxxopts::unicode_string_iterator(&s, s.length());
}

} // namespace std

#else // ifdef CXXOPTS_USE_UNICODE

namespace cxxopts {

template <typename T>
static inline T to_local_string(T&& t) {
  return std::forward<T>(t);
}

CXXOPTS_CONSTEXPR
static inline size_t string_length(const cxx_string& s) noexcept {
  return s.length();
}

CXXOPTS_CONSTEXPR
static inline cxx_string& string_append(cxx_string& s, const cxx_string& a) {
  return s.append(a);
}

CXXOPTS_CONSTEXPR
static inline cxx_string& string_append(cxx_string& s, std::size_t n, char c) {
  return s.append(n, c);
}

template <typename Iterator>
CXXOPTS_CONSTEXPR static inline cxx_string& string_append(cxx_string& s,
                                                          Iterator begin,
                                                          Iterator end) {
  return s.append(begin, end);
}

template <typename T>
static inline std::string to_utf8_string(T&& t) {
  return std::forward<T>(t);
}

CXXOPTS_CONSTEXPR
static inline bool empty(const std::string& s) noexcept {
  return s.empty();
}

} // namespace cxxopts

#endif // ifdef CXXOPTS_USE_UNICODE

/**
 * \defgroup Exceptions
 * @{
 */

namespace cxxopts {

class option_error : public std::runtime_error {
public:
  explicit option_error(const std::string& what_arg)
    : std::runtime_error(what_arg) {
  }
};

class parse_error : public option_error {
public:
  explicit parse_error(const std::string& what_arg)
    : option_error(what_arg) {
  }
};

class spec_error : public option_error {
public:
  explicit spec_error(const std::string& what_arg)
    : option_error(what_arg) {
  }
};

class option_exists_error : public spec_error {
public:
  explicit option_exists_error(const std::string& option)
    : spec_error("Option " + LQUOTE + option + RQUOTE + " already exists") {
  }
};

class invalid_option_format_error : public spec_error {
public:
  explicit invalid_option_format_error(const std::string& format)
    : spec_error("Invalid option format " + LQUOTE + format + RQUOTE) {
  }
};

class option_syntax_error : public parse_error {
public:
  explicit option_syntax_error(const std::string& text)
    : parse_error("Argument " + LQUOTE + text + RQUOTE +
                  " starts with '-' but has incorrect syntax") {
  }
};

class option_not_exists_error : public parse_error {
public:
  explicit option_not_exists_error(const std::string& option)
    : parse_error("Option " + LQUOTE + option + RQUOTE + " does not exist") {
  }
};

class missing_argument_error : public parse_error {
public:
  explicit missing_argument_error(const std::string& option)
    : parse_error("Option " + LQUOTE + option + RQUOTE +
                  " is missing an argument") {
  }
};

class option_requires_argument_error : public parse_error {
public:
  explicit option_requires_argument_error(const std::string& option)
    : parse_error("Option " + LQUOTE + option + RQUOTE +
                  " requires an argument") {
  }
};

class option_not_present_error : public parse_error {
public:
  explicit option_not_present_error(const std::string& option)
    : parse_error("Option " + LQUOTE + option + RQUOTE + " not present") {
  }
};

class argument_incorrect_type : public parse_error {
public:
  explicit argument_incorrect_type(const std::string& arg,
                                   const std::string& type = {})
    : parse_error(
        "Argument " + LQUOTE + arg + RQUOTE + " failed to parse" +
        (type.empty() ? std::string() : (": " + type + " expected"))) {
  }
};

class option_has_no_value_error : public option_error {
public:
  explicit option_has_no_value_error(const std::string& name)
    : option_error(name.empty()
                     ? "Option has no value"
                     : "Option " + LQUOTE + name + RQUOTE + " has no value") {
  }
};

namespace detail {

template <typename T, typename... Args>
CXXOPTS_NORETURN void throw_or_mimic(Args&&... args) {
  static_assert(std::is_base_of<std::exception, T>::value,
                "throw_or_mimic only works on std::exception and "
                "deriving classes");

#ifndef CXXOPTS_NO_EXCEPTIONS
  // If CXXOPTS_NO_EXCEPTIONS is not defined, just throw
  throw T{std::forward<Args>(args)...};
#else
  // Otherwise manually instantiate the exception, print what() to stderr,
  // and terminate.
  T exception{std::forward<Args>(args)...};
  std::cerr << exception.what() << std::endl;
  std::terminate();
#endif
}

} // namespace detail
} // namespace cxxopts

/**@}*/

/**
 * \defgroup Value parsing
 * @{
 */

namespace cxxopts {
namespace detail {

template <typename T, bool B>
struct signed_check;

template <typename T>
struct signed_check<T, true> {
  template <typename U>
  CXXOPTS_CONSTEXPR bool operator()(const U u,
                                    const bool negative) const noexcept {
    return ((static_cast<U>(std::numeric_limits<T>::min()) >= u) && negative) ||
           (static_cast<U>(std::numeric_limits<T>::max()) >= u);
  }
};

template <typename T>
struct signed_check<T, false> {
  template <typename U>
  CXXOPTS_CONSTEXPR bool operator()(const U, const bool) const noexcept {
    return true;
  }
};

template <typename T, typename U>
CXXOPTS_CONSTEXPR inline bool check_signed_range(const U value,
                                                 const bool negative) noexcept {
  return signed_check<T, std::numeric_limits<T>::is_signed>()(value, negative);
}

template <typename R, typename T>
CXXOPTS_CONSTEXPR inline bool checked_negate(R& r,
                                             T&& t,
                                             std::true_type) noexcept {
  // if we got to here, then `t` is a positive number that fits into
  // `R`. So to avoid MSVC C4146, we first cast it to `R`.
  // See https://github.com/jarro2783/cxxopts/issues/62 for more details.
  r = static_cast<R>(-static_cast<R>(t - 1) - 1);
  return true;
}

template <typename R, typename T>
CXXOPTS_CONSTEXPR inline bool checked_negate(R&,
                                             T&&,
                                             std::false_type) noexcept {
  return false;
}

inline bool parse_uint64(const std::string& text,
                         uint64_t& value,
                         bool& negative) noexcept {
  const char* p = text.c_str();
  // String should not be empty.
  if (*p == 0) {
    return false;
  }
  // Parse sign.
  if (*p == '+') {
    ++p;
  } else if (*p == '-') {
    negative = true;
    ++p;
  }
  // Not an integer value.
  if (*p == 0) {
    return false;
  } else {
    value = 0;
  }
  // Hex number.
  if (*p == '0' && *(p + 1) == 'x') {
    p += 2;
    if (*p == 0) {
      return false;
    }
    for (; *p; ++p) {
      uint64_t digit = 0;

      if (*p >= '0' && *p <= '9') {
        digit = static_cast<uint64_t>(*p - '0');
      } else if (*p >= 'a' && *p <= 'f') {
        digit = static_cast<uint64_t>(*p - 'a' + 10);
      } else if (*p >= 'A' && *p <= 'F') {
        digit = static_cast<uint64_t>(*p - 'A' + 10);
      } else {
        return false;
      }

      const uint64_t next = value * 16u + digit;
      if (value > next) {
        return false;
      } else {
        value = next;
      }
    }
    // Decimal number.
  } else {
    for (; *p; ++p) {
      uint64_t digit = 0;

      if (*p >= '0' && *p <= '9') {
        digit = static_cast<uint64_t>(*p - '0');
      } else {
        return false;
      }
      if ((value > std::numeric_limits<uint64_t>::max() / 10u) ||
          (value == std::numeric_limits<uint64_t>::max() / 10u && digit > 5))
      {
        return false;
      } else {
        value = value * 10u + digit;
      }
    }
  }
  return true;
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
inline void parse_value(const std::string& text, T& value) {
  using US = typename std::make_unsigned<T>::type;

  uint64_t u64_result{0};
  US result{0};
  bool negative{false};

  // Parse text to the uint64_t value.
  if (!parse_uint64(text, u64_result, negative)) {
    throw_or_mimic<argument_incorrect_type>(text, "integer");
  }
  // Check unsigned overflow.
  if (u64_result > std::numeric_limits<US>::max()) {
    throw_or_mimic<argument_incorrect_type>(text, "integer");
  } else {
    result = static_cast<US>(u64_result);
  }
  // Check signed overflow.
  if (!check_signed_range<T>(result, negative)) {
    throw_or_mimic<argument_incorrect_type>(text, "integer");
  }
  // Negate value.
  if (negative) {
    if (!checked_negate<T>(
          value, result,
          std::integral_constant<bool, std::numeric_limits<T>::is_signed>()))
    {
      throw_or_mimic<argument_incorrect_type>(text, "integer");
    }
  } else {
    value = static_cast<T>(result);
  }
}

inline void parse_value(const std::string& text, float& value) {
  value = std::stof(text);
}

inline void parse_value(const std::string& text, double& value) {
  value = std::stod(text);
}

inline void parse_value(const std::string& text, long double& value) {
  value = std::stold(text);
}

inline void parse_value(const std::string& text, bool& value) {
  switch (text.size()) {
    case 1: {
      const char ch = text[0];
      if (ch == '1' || ch == 't' || ch == 'T') {
        value = true;
        return;
      }
      if (ch == '0' || ch == 'f' || ch == 'F') {
        value = false;
        return;
      }
      break;
    }
    case 4:
      if ((text[0] == 't' || text[0] == 'T') &&
          (text[1] == 'r' && text[2] == 'u' && text[3] == 'e'))
      {
        value = true;
        return;
      }
      break;
    case 5:
      if ((text[0] == 'f' || text[0] == 'F') &&
          (text[1] == 'a' && text[2] == 'l' && text[3] == 's' &&
           text[4] == 'e'))
      {
        value = false;
        return;
      }
      break;
  }
  throw_or_mimic<argument_incorrect_type>(text, "bool");
}

inline void parse_value(const std::string& text, char& c) {
  if (text.length() != 1) {
    throw_or_mimic<argument_incorrect_type>(text, "char");
  }

  c = text[0];
}

inline void parse_value(const std::string& text, std::string& value) {
  value = text;
}

// The fallback parser. It uses the stringstream parser to parse all types
// that have not been overloaded explicitly.  It has to be placed in the
// source code before all other more specialized templates.
template <typename T,
          typename std::enable_if<!std::is_integral<T>::value>::type* = nullptr>
void parse_value(const std::string& text, T& value) {
  std::istringstream in{text};
  in >> value;
  if (!in) {
    throw_or_mimic<argument_incorrect_type>(text);
  }
}

#ifdef CXXOPTS_HAS_OPTIONAL
template <typename T>
void parse_value(const std::string& text, std::optional<T>& value) {
  T result;
  parse_value(text, result);
  value = std::move(result);
}
#endif

} // namespace detail

/**
 * Settings for customizing parser behaviour.
 */
struct parse_context {
  char delimiter{CXXOPTS_VECTOR_DELIMITER};
};

/**
 * A parser for values of type T.
 */
template <typename T>
struct value_parser {
  using value_type = T;
  /// By default, value cannot act as a container.
  static constexpr bool is_container = false;

  void parse(const parse_context&, const std::string& text, T& value) {
    detail::parse_value(text, value);
  }
};

template <typename T>
struct value_parser<std::vector<T>> {
  using value_type = T;
  /// Value of type std::vector<T> can act as container.
  static constexpr bool is_container = true;

  void parse(const parse_context& ctx,
             const std::string& text,
             std::vector<T>& value) {
    using parser_type = value_parser<T>;

    static_assert(
      !parser_type::is_container ||
        !value_parser<typename parser_type::value_type>::is_container,
      "dimensions of a container type should not exceed 2");

    auto parse_item = [&ctx](const std::string& txt) {
      T v;
      parser_type().parse(ctx, txt, v);
      return v;
    };

    if (text.empty() || parser_type::is_container) {
      value.push_back(parse_item(text));
    } else {
      std::istringstream in{text};
      std::string token;
      while (!in.eof() && std::getline(in, token, ctx.delimiter)) {
        value.push_back(parse_item(token));
      }
    }
  }
};

} // namespace cxxopts

/**@}*/

/**
 * \defgroup Value setup
 * @{
 */

namespace cxxopts {
namespace detail {

#if defined(__GNUC__)
// GNU GCC with -Weffc++ will issue a warning regarding the upcoming class, we
// want to silence it: warning: base class 'class
// std::enable_shared_from_this<cxxopts::value_base>' has accessible non-virtual
// destructor
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
// This will be ignored under other compilers like LLVM clang.
#endif
class value_base : public std::enable_shared_from_this<value_base> {
public:
  value_base() = default;

  virtual ~value_base() = default;

  /** Returns whether the default value was set. */
  CXXOPTS_NODISCARD
  bool has_default() const noexcept {
    return default_;
  }

  /** Returns whether the env variable was set. */
  CXXOPTS_NODISCARD
  bool has_env() const noexcept {
    return env_;
  }

  /** Returns whether the implicit value was set. */
  CXXOPTS_NODISCARD
  bool has_implicit() const noexcept {
    return implicit_;
  }

  /** Returns default value. */
  CXXOPTS_NODISCARD
  const std::string& get_default_value() const noexcept {
    return default_value_;
  }

  /** Returns env variable. */
  CXXOPTS_NODISCARD
  const std::string& get_env_var() const noexcept {
    return env_var_;
  }

  /** Returns implicit value. */
  CXXOPTS_NODISCARD
  const std::string& get_implicit_value() const noexcept {
    return implicit_value_;
  }

  CXXOPTS_NODISCARD
  bool get_no_value() const noexcept {
    return no_value_;
  }

  /**
   * Sets default value.
   *
   * Templated version is used to avoid implicit conversion 0u or nullptr
   * to a string, that leads to runtime error.
   */
  template <typename T>
  typename std::enable_if<
    !std::is_same<std::nullptr_t, typename std::remove_cv<T>::type>::value,
    std::shared_ptr<value_base>>::type
  default_value(T&& value) {
    default_ = true;
    default_value_.assign(std::forward<T>(value));
    return shared_from_this();
  }

  /** Sets delimiter for list values. */
  std::shared_ptr<value_base> delimiter(const char del) {
    parse_ctx_.delimiter = del;
    return shared_from_this();
  }

  /** Sets env variable. */
  template <typename T>
  typename std::enable_if<
    !std::is_same<std::nullptr_t, typename std::remove_cv<T>::type>::value,
    std::shared_ptr<value_base>>::type
  env(T&& var) {
    env_ = true;
    env_var_.assign(std::forward<T>(var));
    return shared_from_this();
  }

  /**
   * Sets implicit value.
   *
   * Templated version is used to avoid implicit conversion 0u or nullptr
   * to a string, that leads to runtime error.
   */
  template <typename T>
  typename std::enable_if<
    !std::is_same<std::nullptr_t, typename std::remove_cv<T>::type>::value,
    std::shared_ptr<value_base>>::type
  implicit_value(T&& value) {
    implicit_ = true;
    implicit_value_.assign(std::forward<T>(value));
    return shared_from_this();
  }

  /** Clears implicit value. */
  std::shared_ptr<value_base> no_implicit_value() {
    no_value_ = false;
    implicit_ = false;
    implicit_value_.clear();
    return shared_from_this();
  }

  /** Sets no-value field. */
  std::shared_ptr<value_base> no_value(const bool on = true) {
    no_value_ = on;
    return shared_from_this();
  }

  /** Returns whether the type of the value is boolean. */
  bool is_boolean() const noexcept {
    return do_is_boolean();
  }

  /** Returns whether the type of the value is container. */
  bool is_container() const noexcept {
    return do_is_container();
  }

  /** Parses the given text into the value. */
  void parse(const std::string& text) {
    return do_parse(parse_ctx_, text);
  }

  /** Parses the default value. */
  void parse() {
    return do_parse(parse_ctx_, default_value_);
  }

protected:
  virtual bool do_is_boolean() const noexcept = 0;

  virtual bool do_is_container() const noexcept = 0;

  virtual void do_parse(const parse_context& ctx, const std::string& text) = 0;

  void set_default_and_implicit(const bool set_default) {
    if (is_boolean()) {
      if (set_default) {
        default_ = true;
        default_value_ = "false";
      }
      implicit_ = true;
      implicit_value_ = "true";
      no_value_ = true;
    }
  }

private:
  /// A default value for the option.
  std::string default_value_{};
  /// A name of an environment variable from which a default
  /// value will be read.
  std::string env_var_{};
  /// An implicit value for the option.
  std::string implicit_value_{};
  /// Configuration of the value parser.
  parse_context parse_ctx_{};

  /// The default value has been set.
  bool default_{false};
  /// The environment variable has been set.
  bool env_{false};
  /// The implicit value has been set.
  bool implicit_{false};
  /// There should be no value for the option.
  bool no_value_{false};
};
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

template <typename T>
class basic_value : public value_base {
  using parser_type = value_parser<T>;

public:
  basic_value()
    : result_(new T{})
    , store_(result_.get()) {
    set_default_and_implicit(true);
  }

  explicit basic_value(T* const t)
    : store_(t) {
    set_default_and_implicit(false);
  }

  const T& get() const noexcept {
    return *store_;
  }

protected:
  bool do_is_boolean() const noexcept final override {
    return std::is_same<T, bool>::value;
  }

  bool do_is_container() const noexcept final override {
    return parser_type::is_container;
  }

  void do_parse(const parse_context& ctx, const std::string& text) override {
    parser_type().parse(ctx, text, *store_);
  }

private:
  basic_value(const basic_value& rhs) = delete;
  basic_value& operator=(const basic_value& rhs) = delete;

private:
  std::unique_ptr<T> result_{};
  T* store_{};
};

} // namespace detail

/**
 * Creates value holder for the specific type.
 */
template <typename T>
std::shared_ptr<detail::basic_value<T>> inline value() {
  return std::make_shared<detail::basic_value<T>>();
}

/**
 * Creates value holder for the specific type.
 */
template <typename T>
std::shared_ptr<detail::basic_value<T>> inline value(T& t) {
  return std::make_shared<detail::basic_value<T>>(&t);
}

} // namespace cxxopts

/**@}*/

namespace cxxopts {

class option_details {
public:
  option_details(std::string short_name,
                 std::string long_name,
                 std::string arg_help,
                 cxx_string desc,
                 std::shared_ptr<detail::value_base> val)
    : short_(std::move(short_name))
    , long_(std::move(long_name))
    , arg_help_(std::move(arg_help))
    , desc_(std::move(desc))
    , hash_(std::hash<std::string>{}(long_ + short_))
    , value_(std::move(val)) {
  }

  CXXOPTS_NODISCARD
  const std::string& arg_help() const noexcept {
    return arg_help_;
  }

  CXXOPTS_NODISCARD
  const cxx_string& description() const noexcept {
    return desc_;
  }

  CXXOPTS_NODISCARD
  std::shared_ptr<detail::value_base> value() const {
    return value_;
  }

  CXXOPTS_NODISCARD
  const std::string& canonical_name() const noexcept {
    return long_.empty() ? short_ : long_;
  }

  CXXOPTS_NODISCARD
  const std::string& short_name() const noexcept {
    return short_;
  }

  CXXOPTS_NODISCARD
  const std::string& long_name() const noexcept {
    return long_;
  }

  CXXOPTS_NODISCARD
  std::size_t hash() const noexcept {
    return hash_;
  }

  CXXOPTS_NODISCARD
  const std::string& default_value() const noexcept {
    return value_->get_default_value();
  }

  CXXOPTS_NODISCARD
  const std::string& implicit_value() const noexcept {
    return value_->get_implicit_value();
  }

  CXXOPTS_NODISCARD
  bool has_default() const noexcept {
    return value_->has_default();
  }

  CXXOPTS_NODISCARD
  bool has_implicit() const noexcept {
    return value_->has_implicit();
  }

  CXXOPTS_NODISCARD
  bool is_container() const noexcept {
    return value_->is_container();
  }

  CXXOPTS_NODISCARD
  bool is_boolean() const noexcept {
    return value_->is_boolean();
  }

private:
  /// Short name of the option.
  std::string short_;
  /// Long name of the option.
  std::string long_;
  std::string arg_help_;
  /// Description of the option.
  cxx_string desc_;
  std::size_t hash_;
  std::shared_ptr<detail::value_base> value_;
};

/**
 * Parsed value of an option.
 */
class option_value {
public:
  /**
   * A number of occurrences of the option value in
   * the command line arguments.
   */
  CXXOPTS_NODISCARD
  std::size_t count() const noexcept {
    return count_;
  }

  // TODO: maybe default options should count towards the number of arguments
  CXXOPTS_NODISCARD
  bool has_default() const noexcept {
    return default_;
  }

  /**
   * Returns true if there was a value for the option in the
   * command line arguments.
   */
  CXXOPTS_NODISCARD
  bool has_value() const noexcept {
    return value_ != nullptr;
  }

  /**
   * Casts option value to the specific type.
   */
  template <typename T>
  const T& as() const {
    if (!has_value()) {
      detail::throw_or_mimic<option_has_no_value_error>(long_name_);
    }
#ifdef CXXOPTS_NO_RTTI
    return static_cast<const detail::basic_value<T>&>(*value_).get();
#else
    return dynamic_cast<const detail::basic_value<T>&>(*value_).get();
#endif
  }

public:
  /**
   * Parses option value from the given text.
   */
  void parse(const option_details& details, const std::string& text) {
    ensure_value(details);
    ++count_;
    value_->parse(text);
    long_name_ = details.long_name();
  }

  /**
   * Parses option value from the default value.
   */
  void parse_default(const option_details& details) {
    ensure_value(details);
    default_ = true;
    long_name_ = details.long_name();
    value_->parse();
  }

  void parse_no_value(const option_details& details) {
    long_name_ = details.long_name();
  }

private:
  void ensure_value(const option_details& details) {
    if (value_ == nullptr) {
      value_ = details.value();
    }
  }

  std::string long_name_{};
  // Holding this pointer is safe, since option_value's only exist
  // in key-value pairs, where the key has the string we point to.
  std::shared_ptr<detail::value_base> value_{};
  std::size_t count_{0};
  bool default_{false};
};

/**
 * Provides the result of parsing of the command line arguments.
 */
class parse_result {
public:
  /// Maps option name to hash of the name.
  using name_hash_map = std::unordered_map<std::string, std::size_t>;
  /// Maps hash of an option name to the option value.
  using parsed_hash_map = std::unordered_map<std::size_t, option_value>;

  class key_value {
  public:
    key_value(std::string key, std::string value) noexcept
      : key_(std::move(key))
      , value_(std::move(value)) {
    }

    CXXOPTS_NODISCARD
    const std::string& key() const noexcept {
      return key_;
    }

    CXXOPTS_NODISCARD
    const std::string& value() const noexcept {
      return value_;
    }

    /**
     * Parses the value to a variable of the specific type.
     */
    template <typename T>
    T as() const {
      T result;
      value_parser<T>().parse(parse_context(), value_, result);
      return result;
    }

  private:
    const std::string key_;
    const std::string value_;
  };

public:
  parse_result() = default;
  parse_result(const parse_result&) = default;
  parse_result(parse_result&&) = default;
  parse_result(name_hash_map&& keys,
               parsed_hash_map&& values,
               std::vector<key_value>&& sequential,
               std::vector<std::string>&& unmatched_args,
               std::size_t consumed)
    : keys_(std::move(keys))
    , values_(std::move(values))
    , sequential_(std::move(sequential))
    , unmatched_(std::move(unmatched_args))
    , consumed_arguments_(consumed) {
  }

  parse_result& operator=(const parse_result&) = default;
  parse_result& operator=(parse_result&&) = default;

  /**
   * Returns a number of occurrences of the option in
   * the command line arguments.
   */
  CXXOPTS_NODISCARD
  std::size_t count(const std::string& name) const {
    const auto ki = keys_.find(name);
    if (ki == keys_.end()) {
      return 0;
    }

    const auto vi = values_.find(ki->second);
    if (vi == values_.end()) {
      return 0;
    }

    return vi->second.count();
  }

  CXXOPTS_NODISCARD
  bool has(const std::string& name) const {
    return count(name) != 0;
  }

  const option_value& operator[](const std::string& name) const {
    const auto ki = keys_.find(name);
    if (ki == keys_.end()) {
      detail::throw_or_mimic<option_not_present_error>(name);
    }

    const auto vi = values_.find(ki->second);
    if (vi == values_.end()) {
      detail::throw_or_mimic<option_not_present_error>(name);
    }

    return vi->second;
  }

  /**
   * Returns list of recognized options with non empty value.
   */
  CXXOPTS_NODISCARD
  const std::vector<key_value>& arguments() const noexcept {
    return sequential_;
  }

  /**
   * Returns number of consumed command line arguments.
   */
  CXXOPTS_NODISCARD
  std::size_t consumed() const noexcept {
    return consumed_arguments_;
  }

  /**
   * Returns list of unmatched arguments.
   */
  CXXOPTS_NODISCARD
  const std::vector<std::string>& unmatched() const noexcept {
    return unmatched_;
  }

private:
  name_hash_map keys_{};
  parsed_hash_map values_{};
  std::vector<key_value> sequential_{};
  /// List of arguments that did not match to any defined option.
  std::vector<std::string> unmatched_{};
  /// Number of consument command line arguments.
  std::size_t consumed_arguments_{0};
};

namespace detail {

class option_parser {
  using option_map =
    std::unordered_map<std::string, std::shared_ptr<option_details>>;
  using positional_list = std::vector<std::string>;
  using positional_list_iterator = positional_list::const_iterator;

  struct option_data {
    std::string name{};
    std::string value{};
    bool is_long{false};
    bool has_value{false};
  };

public:
  option_parser(const option_map& options,
                const positional_list& positional,
                bool allow_unrecognised,
                bool stop_on_positional)
    : options_(options)
    , positional_(positional)
    , allow_unrecognised_(allow_unrecognised)
    , stop_on_positional_(stop_on_positional) {
  }

  parse_result parse(const int argc, const char* const* argv) {
    int current = 1;
    auto next_positional = positional_.begin();
    std::vector<std::string> unmatched;

    while (current < argc) {
      if (is_dash_dash(argv[current])) {
        // Skip dash-dash argument.
        ++current;
        if (stop_on_positional_) {
          break;
        }
        // Try to consume all remaining arguments as positional.
        for (; current < argc; ++current) {
          if (!consume_positional(argv[current], next_positional)) {
            break;
          }
        }
        // Adjust argv for any that couldn't be swallowed.
        for (; current != argc; ++current) {
          unmatched.emplace_back(argv[current]);
        }
        break;
      }

      option_data result;

      if (!parse_argument(argv[current], result)) {
        // Not a flag.
        // But if it starts with a `-`, then it's an error.
        if (argv[current][0] == '-' && argv[current][1] != '\0') {
          if (!allow_unrecognised_) {
            detail::throw_or_mimic<option_syntax_error>(argv[current]);
          }
        }
        if (stop_on_positional_) {
          break;
        }
        // If true is returned here then it was consumed, otherwise it
        // is ignored.
        if (!consume_positional(argv[current], next_positional)) {
          unmatched.emplace_back(argv[current]);
        }
        // If we return from here then it was parsed successfully, so
        // continue.
      } else if (result.is_long) {
        // Long option.
        const std::string& name = result.name;
        const auto oi = options_.find(name);

        if (oi == options_.end()) {
          if (allow_unrecognised_) {
            // Keep unrecognised options in argument list,
            // skip to next argument.
            unmatched.emplace_back(argv[current]);
            ++current;
            continue;
          }
          // Error.
          detail::throw_or_mimic<option_not_exists_error>(name);
        }

        const auto& opt = oi->second;
        // Equal sign provided for the long option?
        if (result.has_value) {
          // Parse the option given.
          parse_option(opt, result.value);
        } else {
          // Parse the next argument.
          checked_parse_arg(argc, argv, current, opt, name);
        }
      } else {
        // Single short option or a group of short options.
        const std::string& seq = result.name;
        // Iterate over the sequence of short options.
        for (std::size_t i = 0; i != seq.size(); ++i) {
          const std::string name(1, seq[i]);
          const auto oi = options_.find(name);

          if (oi == options_.end()) {
            if (allow_unrecognised_) {
              unmatched.push_back(std::string("-") + seq[i]);
              continue;
            }
            // Error.
            detail::throw_or_mimic<option_not_exists_error>(name);
          }

          const auto& opt = oi->second;
          if (i + 1 == seq.size()) {
            // It must be the last argument.
            checked_parse_arg(argc, argv, current, opt, name);
          } else if (opt->has_implicit()) {
            parse_option(opt, opt->implicit_value());
          } else {
            parse_option(opt, seq.substr(i + 1));
            break;
          }
        }
      }

      ++current;
    }

    // Setup default or env values.
    for (auto& opt : options_) {
      auto& detail = opt.second;
      auto& store = parsed_[detail->hash()];
      const auto& value = detail->value();

      // Skip options with parsed values.
      if (store.count() || store.has_default()) {
        continue;
      }
      // Try to setup env value.
      if (value->has_env()) {
        if (const char* env = std::getenv(value->get_env_var().c_str())) {
          store.parse(*detail, std::string(env));
          continue;
        }
      }
      // Try to setup default value.
      if (value->has_default()) {
        store.parse_default(*detail);
      } else {
        store.parse_no_value(*detail);
      }
    }

    parse_result::name_hash_map keys;
    // Finalize aliases.
    for (const auto& option : options_) {
      const auto& detail = option.second;
      const auto hash = detail->hash();

      if (detail->short_name().size()) {
        keys[detail->short_name()] = hash;
      }
      if (detail->long_name().size()) {
        keys[detail->long_name()] = hash;
      }
    }

    assert(stop_on_positional_ || argc == current || argc == 0);

    return parse_result(std::move(keys), std::move(parsed_),
                        std::move(sequential_), std::move(unmatched), current);
  }

private:
  bool consume_positional(const std::string& arg,
                          positional_list_iterator& next) {
    for (; next != positional_.end(); ++next) {
      const auto oi = options_.find(*next);
      if (oi == options_.end()) {
        detail::throw_or_mimic<option_not_exists_error>(*next);
      }
      if (oi->second->value()->is_container()) {
        parse_option(oi->second, arg);
        return true;
      }
      if (parsed_[oi->second->hash()].count() == 0) {
        parse_option(oi->second, arg);
        ++next;
        return true;
      }
    }
    return false;
  }

  bool is_dash_dash(const char* str) const noexcept {
    return (str[0] != 0 && str[0] == '-') && (str[1] != 0 && str[1] == '-') &&
           (str[2] == 0);
  }

  bool is_dash_dash_or_option_name(const char* const arg) const {
    // The dash-dash symbol has a special meaning and cannot
    // be interpreted as an option value.
    if (is_dash_dash(arg)) {
      return true;
    }

    option_data result;
    // The argument does not match an option format
    // so that it can be safely consumed as a value.
    if (!parse_argument(arg, result)) {
      return false;
    }

    auto check_name = [this](const std::string& opt) {
      return options_.find(opt) != options_.end();
    };
    // Check that the argument does not match any
    // existing option.
    if (result.is_long) {
      return check_name(result.name);
    } else {
      return check_name(result.name.substr(0, 1));
    }

    return false;
  }

  void checked_parse_arg(const int argc,
                         const char* const* argv,
                         int& current,
                         const std::shared_ptr<option_details>& value,
                         const std::string& name) {
    auto parse_implicit = [&]() {
      if (value->has_implicit()) {
        parse_option(value, value->implicit_value());
      } else {
        detail::throw_or_mimic<missing_argument_error>(name);
      }
    };

    if (current + 1 == argc || value->value()->get_no_value()) {
      // Last argument or the option without value.
      parse_implicit();
    } else {
      const char* const arg = argv[current + 1];
      // Check that we do not silently consume any option as a value
      // of another option.
      if (arg[0] == '-' && is_dash_dash_or_option_name(arg)) {
        parse_implicit();
      } else {
        // Parse argument as a value for the option.
        parse_option(value, arg);
        ++current;
      }
    }
  }

  bool parse_argument(const std::string& text, option_data& data) const {
    const char* p = text.c_str();
    const char* end = text.c_str() + text.size();
    // The string should be at least two character long and starts with '-'.
    if (*p == 0 || *(p + 1) == 0 || *p != '-') {
      return false;
    }
    // Skip the '-'.
    ++p;
    // Long option starts with '--'.
    if (*p == '-') {
      ++p;
      if (std::isalnum(*p) && *(p + 1) != 0) {
        data.is_long = true;
        data.name += *p;
        ++p;
      } else {
        return false;
      }
      for (; *p; ++p) {
        if (*p == '=') {
          ++p;
          data.has_value = true;
          data.value.assign(p, end);
          break;
        }
        if (*p == '-' || *p == '_' || *p == '.' || std::isalnum(*p)) {
          data.name += *p;
        } else {
          return false;
        }
      }
      return data.name.size() > 1;
    } else {
      // Single char short option should start with an alnum or
      // be a question mark.
      if (!(std::isalnum(*p) || (*p == '?' && *(p + 1) == 0))) {
        return false;
      }
      // Copy the whole string and interpret it later as
      // a group of short options.
      data.name.assign(p, end);
      return true;
    }
    return false;
  }

  void parse_option(const std::shared_ptr<option_details>& details,
                    const std::string& arg) {
    parsed_[details->hash()].parse(*details, arg);
    sequential_.emplace_back(details->canonical_name(), arg);
  }

private:
  const option_map& options_;
  const positional_list& positional_;
  const bool allow_unrecognised_;
  const bool stop_on_positional_;

  std::vector<parse_result::key_value> sequential_{};
  parse_result::parsed_hash_map parsed_{};
};

} // namespace detail

class options;

class option {
public:
  option(std::string opts,
         std::string desc,
         std::shared_ptr<detail::value_base> value = ::cxxopts::value<bool>(),
         std::string arg_help = {}) noexcept
    : opts_(std::move(opts))
    , desc_(std::move(desc))
    , value_(std::move(value))
    , arg_help_(std::move(arg_help)) {
  }

private:
  friend class ::cxxopts::options;

  /// Short and long names of the option.
  const std::string opts_;
  /// Description of the option.
  const std::string desc_;
  const std::shared_ptr<detail::value_base> value_;
  const std::string arg_help_;
};

/**
 * Specification of command line options.
 */
class options {
public:
  struct help_group_details {
    std::string name{};
    std::vector<std::shared_ptr<option_details>> options{};
  };

  class option_adder {
  public:
    option_adder(std::string group, options& options) noexcept
      : group_(std::move(group))
      , options_(options) {
    }

    option_adder& operator()(const std::string& opts,
                             const std::string& desc,
                             const std::shared_ptr<detail::value_base>& value =
                               ::cxxopts::value<bool>(),
                             const std::string arg_help = {}) {
      std::string s;
      std::string l;
      if (parse_option_specifier(opts, s, l)) {
        assert(s.empty() || s.size() == 1);
        assert(l.empty() || l.size() > 1);

        options_.add_option(group_, std::move(s), std::move(l), desc, value,
                            std::move(arg_help));
      } else {
        detail::throw_or_mimic<invalid_option_format_error>(opts);
      }
      return *this;
    }

  private:
    bool parse_option_specifier(const std::string& text,
                                std::string& s,
                                std::string& l) const {
      const char* p = text.c_str();
      if (*p == 0) {
        return false;
      } else {
        s.clear();
        l.clear();
      }
      // Short option.
      if (*(p + 1) == 0 || *(p + 1) == ',') {
        if (*p == '?' || std::isalnum(*p)) {
          s = *p;
          ++p;
        } else {
          return false;
        }
      }
      // Skip comma.
      if (*p == ',') {
        if (s.empty()) {
          return false;
        }
        ++p;
      }
      // Skip spaces.
      while (*p && *p == ' ') {
        ++p;
      }
      // Valid specifier without long option.
      if (*p == 0) {
        return true;
      } else {
        l.reserve((text.c_str() + text.size()) - p);
      }
      // First char of an option name should be alnum.
      if (std::isalnum(*p)) {
        l += *p;
        ++p;
      }
      for (; *p; ++p) {
        if (*p == '-' || *p == '_' || *p == '.' || std::isalnum(*p)) {
          l += *p;
        } else {
          return false;
        }
      }
      return l.size() > 1;
    }

  private:
    const std::string group_;
    options& options_;
  };

public:
  explicit options(std::string program, std::string help_string = {})
    : program_(std::move(program))
    , help_string_(to_local_string(std::move(help_string)))
    , custom_help_("[OPTION...]")
    , positional_help_("positional parameters") {
  }

  /**
   * Adds list of options to the specific group.
   */
  void add_options(const std::string& group,
                   std::initializer_list<option> opts) {
    option_adder adder(group, *this);
    for (const auto& opt : opts) {
      adder(opt.opts_, opt.desc_, opt.value_, opt.arg_help_);
    }
  }

  /**
   * Adds an option to the specific group.
   */
  void add_option(const std::string& group, const option& opt) {
    add_options(group, {opt});
  }

  option_adder add_options(std::string group = {}) {
    return option_adder(std::move(group), *this);
  }

  options& allow_unrecognised_options(const bool value = true) noexcept {
    allow_unrecognised_ = value;
    return *this;
  }

  options& custom_help(std::string help_text) noexcept {
    custom_help_ = std::move(help_text);
    return *this;
  }

  options& footer(std::string text) noexcept {
    footer_ = std::move(text);
    return *this;
  }

  template <typename... Args>
  void parse_positional(Args&&... args) {
    parse_positional(std::vector<std::string>{std::forward<Args>(args)...});
  }

  template <typename I,
            typename std::enable_if<
              !std::is_same<typename std::iterator_traits<I>::value_type,
                            void>::value>::type>
  void parse_positional(const I begin, const I end) {
    parse_positional(std::vector<std::string>(begin, end));
  }

  void parse_positional(std::vector<std::string> opts) {
    positional_ = std::move(opts);
    positional_set_ =
      std::unordered_set<std::string>(positional_.begin(), positional_.end());
  }

  options& positional_help(std::string help_text) noexcept {
    positional_help_ = std::move(help_text);
    return *this;
  }

  options& set_tab_expansion(bool expansion = true) noexcept {
    tab_expansion_ = expansion;
    return *this;
  }

  options& set_width(std::size_t width) noexcept {
    width_ = width;
    return *this;
  }

  options& show_positional_help(const bool value = true) noexcept {
    show_positional_ = value;
    return *this;
  }

  /**
   * Stop parsing at first positional argument.
   */
  options& stop_on_positional(const bool value = true) noexcept {
    stop_on_positional_ = value;
    return *this;
  }

public:
  /**
   * Parses the command line arguments according to the current specification.
   */
  parse_result parse(int argc, const char* const* argv) const {
    return detail::option_parser(options_, positional_, allow_unrecognised_,
                                 stop_on_positional_)
      .parse(argc, argv);
  }

  /**
   * Generates help for the options.
   */
  std::string help(const std::vector<std::string>& help_groups = {},
                   const bool print_usage = true) const {
    cxx_string result;

    if (!empty(help_string_)) {
      result += wrap_string(help_string_, 0, width_);
      result += '\n';
    }

    if (print_usage) {
      result += "usage: ";
      result += to_local_string(program_);
      result += " ";
      result += to_local_string(custom_help_);
    }

    if (!positional_.empty() && !positional_help_.empty()) {
      result += " ";
      result += to_local_string(positional_help_);
    }

    result += "\n\n";

    if (help_groups.empty()) {
      generate_all_groups_help(result);
    } else {
      generate_group_help(result, help_groups);
    }

    if (!empty(footer_)) {
      result += "\n";
      result += wrap_string(to_local_string(footer_), 0, width_);
    }

    return to_utf8_string(result);
  }

  /**
   * Returns list of the defined groups.
   */
  std::vector<std::string> groups() const {
    std::vector<std::string> names;

    names.reserve(help_.size());

    for (auto hi = help_.cbegin(); hi != help_.cend(); ++hi) {
      names.push_back(hi->first);
    }

    return names;
  }

  const help_group_details& group_help(const std::string& group) const {
    return help_.at(group);
  }

  const std::string& program() const noexcept {
    return program_;
  }

private:
  void add_option(const std::string& group,
                  const std::string& s,
                  const std::string& l,
                  std::string desc,
                  const std::shared_ptr<detail::value_base>& value,
                  std::string arg_help) {
    auto details = std::make_shared<option_details>(
      s, l, std::move(arg_help), to_local_string(std::move(desc)), value);

    if (!s.empty()) {
      add_one_option(s, details);
    }
    if (!l.empty()) {
      add_one_option(l, details);
    }

    // Add the help details.
    help_[group].options.push_back(std::move(details));
  }

  void add_one_option(const std::string& name,
                      const std::shared_ptr<option_details>& details) {
    const auto in = options_.emplace(name, details);

    if (!in.second) {
      detail::throw_or_mimic<option_exists_error>(name);
    }
  }

  cxx_string format_option(const option_details& o) const {
    const auto& s = o.short_name();
    const auto& l = o.long_name();

    cxx_string result = "  ";

    if (!s.empty()) {
      result += "-";
      result += to_local_string(s);
      if (!l.empty()) {
        result += ",";
      }
    } else {
      result += "   ";
    }

    if (!l.empty()) {
      result += " --";
      result += to_local_string(l);
    }

    if (!o.is_boolean()) {
      const auto arg =
        !o.arg_help().empty() ? to_local_string(o.arg_help()) : "arg";

      if (o.has_implicit()) {
        result += " [=";
        result += arg;
        result += "(=";
        result += to_local_string(o.implicit_value());
        result += ")]";
      } else {
        result += " ";
        result += arg;
      }
    }

    return result;
  }

  cxx_string format_description(const option_details& o,
                                std::size_t start,
                                std::size_t allowed,
                                bool tab_expansion) const {
    cxx_string desc = o.description();

    if (o.has_default() && (!o.is_boolean() || o.default_value() != "false")) {
      if (!o.default_value().empty()) {
        desc += to_local_string(" (default: " + o.default_value() + ")");
      } else {
        desc += to_local_string(" (default: \"\")");
      }
    }

    if (tab_expansion) {
      desc = expand_tab_character(desc);
    }

    return wrap_string(std::move(desc), start, allowed);
  }

  cxx_string help_one_group(const std::string& group_name) const {
    using option_help =
      std::vector<std::pair<cxx_string, const option_details*>>;

    const auto gi = help_.find(group_name);
    if (gi == help_.end()) {
      return cxx_string();
    }

    option_help format;
    std::size_t longest = 0;
    cxx_string result;

    if (!group_name.empty()) {
      result += to_local_string(group_name);
      result += '\n';
    }
    // Preallocate buffer for list of options.
    format.reserve(gi->second.options.size() -
                   (show_positional_ ? 0 : positional_set_.size()));

    for (const auto& o : gi->second.options) {
      if (!show_positional_ &&
          positional_set_.find(o->long_name()) != positional_set_.end())
      {
        continue;
      }

      cxx_string s = format_option(*o);
      longest = std::max(longest, string_length(s));
      format.emplace_back(std::move(s), o.get());
    }
    longest = std::min(longest, OPTION_LONGEST);

    // widest allowed description -- min 10 chars for helptext/line
    std::size_t allowed = 10;
    if (width_ > allowed + longest + OPTION_DESC_GAP) {
      allowed = width_ - longest - OPTION_DESC_GAP;
    }

    for (auto fi = std::begin(format); fi != std::end(format); ++fi) {
      const auto& d = format_description(*fi->second, longest + OPTION_DESC_GAP,
                                         allowed, tab_expansion_);

      result += fi->first;
      if (string_length(fi->first) > longest) {
        result += '\n';
        result += to_local_string(std::string(longest + OPTION_DESC_GAP, ' '));
      } else {
        result += to_local_string(std::string(
          longest + OPTION_DESC_GAP - string_length(fi->first), ' '));
      }
      result += d;
      result += '\n';
    }

    return result;
  }

  void generate_group_help(cxx_string& result,
                           const std::vector<std::string>& print_groups) const {
    for (std::size_t i = 0; i != print_groups.size(); ++i) {
      const cxx_string& group_help_text = help_one_group(print_groups[i]);
      if (empty(group_help_text)) {
        continue;
      }
      result += group_help_text;
      if (i < print_groups.size() - 1) {
        result += '\n';
      }
    }
  }

  void generate_all_groups_help(cxx_string& result) const {
    generate_group_help(result, groups());
  }

  cxx_string expand_tab_character(const cxx_string& text) const {
    cxx_string result;
    cxx_string::const_iterator pi = std::begin(text);
    std::size_t size = 0;
    // Preallocate result buffer.
    result.reserve(text.size());
    // Process source string.
    for (auto ci = std::begin(text); ci != std::end(text); ++ci) {
      if (*ci == '\n') {
        size = 0;
      } else if (*ci == '\t') {
        const std::size_t skip = OPTION_TAB_SIZE - size % OPTION_TAB_SIZE;
        result.append(pi, ci);
        string_append(result, skip, ' ');
        size += skip;
        pi = ci + 1;
      } else {
        ++size;
      }
    }
    // Append rest of the source string.
    result.append(pi, std::end(text));

    return result;
  }

  cxx_string wrap_string(const cxx_string& desc,
                         const std::size_t start,
                         const std::size_t allowed) const {
    cxx_string result;

    // Nothing to wrap for empty string.
    if (std::begin(desc) == std::end(desc)) {
      return result;
    }

    auto current = std::begin(desc);
    auto previous = current;
    auto start_line = current;
    auto last_space = current;
    auto size = std::size_t{};
    bool only_whitespace = true;

    result.reserve(desc.size());

    for (; current != std::end(desc); ++current) {
      bool append_new_line = false;

      if (std::isblank(*previous)) {
        last_space = current;
      }

      if (!std::isblank(*current)) {
        only_whitespace = false;
      }
      // Skip all line feed characters.
      if (*current == '\n') {
        append_new_line = true;
        do {
          previous = current;
          ++current;
        } while (current != std::end(desc) && *current == '\n');
      }

      if (!append_new_line && size >= allowed) {
        if (last_space != start_line) {
          current = last_space;
          previous = current;
        }
        append_new_line = true;
      }

      if (append_new_line) {
        string_append(result, start_line, current);
        start_line = current;
        last_space = current;

        if (*previous != '\n') {
          string_append(result, "\n");
        }

        string_append(result, start, ' ');

        if (*previous != '\n') {
          string_append(result, last_space, current);
        }

        only_whitespace = true;
        size = 0;
      }

      if (current == std::end(desc)) {
        break;
      }

      previous = current;
      ++size;
    }

    // Append whatever is left but ignore whitespace.
    if (!only_whitespace) {
      string_append(result, start_line, current);
    }

    return result;
  }

private:
  using option_map =
    std::unordered_map<std::string, std::shared_ptr<option_details>>;
  using positional_list = std::vector<std::string>;

  std::string program_;
  cxx_string help_string_;
  std::string custom_help_;
  std::string positional_help_;
  std::string footer_{};
  std::size_t width_{76};
  /// Allow consume unrecognized options
  /// instead of throwing an error.
  bool allow_unrecognised_{false};
  /// Show help for options bonded to positional arguments.
  bool show_positional_{false};
  /// Stop parsing at first positional argument.
  bool stop_on_positional_{false};
  /// Replace tab with spaces.
  bool tab_expansion_{false};

  /// Named options.
  /// Short and long names exist as separate entries but
  /// point to the same object.
  option_map options_{};
  /// List of named positional arguments.
  positional_list positional_{};
  std::unordered_set<std::string> positional_set_{};
  /// Mapping from groups to help options.
  std::map<std::string, help_group_details> help_{};
};

} // namespace cxxopts

#endif // CXXOPTS_HPP_INCLUDED