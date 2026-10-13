#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;

const std::string kLpFile = "lp-file-format.lp";

// Write `content` to an LP file and read it with `highs`. Everything that
// HiGHS logs while reading the file is stored in `log`.
static HighsStatus readLp(Highs& highs, const std::string& content,
                          std::string& log) {
  FILE* file = fopen(kLpFile.c_str(), "wb");
  fwrite(content.data(), 1, content.size(), file);
  fclose(file);
  log.clear();
  highs.setOptionValue("output_flag", true);
  highs.setOptionValue("log_to_console", dev_run);
  highs.setCallback(
      [](int, const std::string& message, const HighsCallbackOutput*,
         HighsCallbackInput*, void* user_data) {
        *static_cast<std::string*>(user_data) += message;
      },
      &log);
  highs.startCallback(kCallbackLogging);
  const HighsStatus status = highs.readModel(kLpFile);
  highs.stopCallback(kCallbackLogging);
  std::remove(kLpFile.c_str());
  if (dev_run) printf("%s", log.c_str());
  return status;
}

static HighsStatus readLp(Highs& highs, const std::string& content) {
  std::string log;
  return readLp(highs, content, log);
}

struct DiagnosticCase {
  std::string name;
  std::string content;
  std::string diagnostic;
};

TEST_CASE("lp-file-format-quad-no-space", "[LpFileFormat]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/qcqp.lp";
  // HiGHS cannot handle quadratic constraints as there are in qcqp.lp
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kError);
}

TEST_CASE("lp-file-format-errors", "[LpFileFormat]") {
  // Each file is invalid, and reading it logs the diagnostic
  const std::vector<DiagnosticCase> cases = {
    {"control-character",
     "min\n"
     " x\n"
     "\x01\n"
     "end\n",
     "ERROR:   unexpected control character 0x01\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 |  \n"
     "  | ^ this character is not valid in an LP file\n"},
    {"missing-end",
     "min\n"
     " x\n"
     "st\n"
     " c: x >= 1\n"
     "\n"
     "  \n",
     "ERROR:   expected `end`, found end of file\n"
     " --> lp-file-format.lp:4:11\n"
     "  |\n"
     "4 |  c: x >= 1\n"
     "  |           ^ expected `end` after this\n"
     "  |\n"
     "  = help: an LP file must finish with the keyword `end`\n"},
    {"no-section-keyword",
     "x + y\n"
     "min\n"
     " x\n"
     "end\n",
     "ERROR:   expected a section keyword, found `x`\n"
     " --> lp-file-format.lp:1:1\n"
     "  |\n"
     "1 | x + y\n"
     "  | ^ expected a keyword like `minimize` or `maximize`\n"
     "  |\n"
     "  = help: an LP file must start with a section keyword\n"},
    {"subject-to-split",
     "subject\n"
     "to\n"
     " c: x >= 1\n"
     "end\n",
     "ERROR:   expected a section keyword, found `subject`\n"
     " --> lp-file-format.lp:1:1\n"
     "  |\n"
     "1 | subject\n"
     "  | ^^^^^^^ expected a keyword like `minimize` or `maximize`\n"
     "  |\n"
     "  = help: an LP file must start with a section keyword\n"},
    {"duplicate-objective",
     "min\n"
     " x\n"
     "max\n"
     " x\n"
     "end\n",
     "ERROR:   the objective is defined more than once\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | max\n"
     "  | ^^^ second objective section\n"
     "  |\n"
     "  = help: the first objective section is on line 1\n"},
    {"sos-section",
     "min\n"
     " x\n"
     "sos\n"
     " s1: S1:: x:1 y:2\n"
     "end\n",
     "ERROR:   SOS constraints are not supported by HiGHS\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | sos\n"
     "  | ^^^ SOS section\n"},
    {"content-after-end",
     "min\n"
     " x\n"
     "end\n"
     "x\n",
     "ERROR:   expected end of file, found `x`\n"
     " --> lp-file-format.lp:4:1\n"
     "  |\n"
     "4 | x\n"
     "  | ^ unexpected content after `end`\n"},
    {"missing-sign",
     "min\n"
     " obj: x y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x y\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"missing-sign-keyword",
     "Minimize a subject to a >= 1\n",
     "ERROR:   expected `+` or `-`, found `subject`\n"
     " --> lp-file-format.lp:1:12\n"
     "  |\n"
     "1 | Minimize a subject to a >= 1\n"
     "  |            ^^^^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: section keywords like `subject` must be at the start of a line\n"},
    {"missing-sign-keyword-semi",
     "min\n"
     " obj: x semi\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `semi`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x semi\n"
     "  |         ^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: section keywords like `semi` must be at the start of a line\n"},
    {"missing-sign-keyword-such",
     "min\n"
     " obj: x such that y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `such`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x such that y\n"
     "  |         ^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: section keywords like `such` must be at the start of a line\n"},
    {"missing-sign-digit-name",
     "min\n"
     " obj: x + 0.2 5C0ST\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `5`\n"
     " --> lp-file-format.lp:2:15\n"
     "  |\n"
     "2 |  obj: x + 0.2 5C0ST\n"
     "  |               ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: `5C0ST` looks like a name, but names cannot start with a digit\n"},
    {"missing-sign-constant",
     "min\n"
     " obj: 2 3 x\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `3`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: 2 3 x\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"coefficient-bracket",
     "min\n"
     " obj: x + 2 [ x * y ]\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `[`\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: x + 2 [ x * y ]\n"
     "  |             ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: a coefficient cannot multiply `[ ]`; move it inside the brackets\n"},
    {"expected-term",
     "min\n"
     " obj: >= 1\n"
     "end\n",
     "ERROR:   expected a term, found `>=`\n"
     " --> lp-file-format.lp:2:7\n"
     "  |\n"
     "2 |  obj: >= 1\n"
     "  |       ^^ expected a number, a variable, or `[`\n"},
    {"expected-term-after-sign",
     "min\n"
     " obj: x + :\n"
     "end\n",
     "ERROR:   expected a term, found `:`\n"
     " --> lp-file-format.lp:2:11\n"
     "  |\n"
     "2 |  obj: x + :\n"
     "  |           ^ expected a number, a variable, or `[` after the sign\n"},
    {"expected-term-keyword",
     "min\n"
     " obj: x +\n"
     "st\n"
     "end\n",
     "ERROR:   expected a term, found `st`\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | st\n"
     "  | ^^ expected a number, a variable, or `[` after the sign\n"
     "  |\n"
     "  = help: `st` is a section keyword because it is at the start of a line\n"},
    {"expected-variable-after-times",
     "min\n"
     " obj: 2 * 3\n"
     "end\n",
     "ERROR:   expected a variable, found `3`\n"
     " --> lp-file-format.lp:2:11\n"
     "  |\n"
     "2 |  obj: 2 * 3\n"
     "  |           ^ expected a variable after `*`\n"},
    {"times-outside-brackets",
     "min\n"
     " obj: x * y\n"
     "end\n",
     "ERROR:   quadratic terms must be inside `[` and `]`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x * y\n"
     "  |         ^ unexpected `*` outside `[ ]`\n"},
    {"power-outside-brackets",
     "min\n"
     " obj: 2 x ^ 2\n"
     "end\n",
     "ERROR:   quadratic terms must be inside `[` and `]`\n"
     " --> lp-file-format.lp:2:11\n"
     "  |\n"
     "2 |  obj: 2 x ^ 2\n"
     "  |           ^ unexpected `^` outside `[ ]`\n"},
    {"unclosed-bracket-eof",
     "min\n"
     " obj: [ x * y\n",
     "ERROR:   expected `]`, found end of file\n"
     " --> lp-file-format.lp:2:14\n"
     "  |\n"
     "2 |  obj: [ x * y\n"
     "  |              ^ expected `]` before this\n"
     "  |\n"
     "  = help: the `[` on line 2 is never closed\n"},
    {"unclosed-bracket-keyword",
     "min\n"
     " obj: [ x * y\n"
     "end\n",
     "ERROR:   expected `]`, found `end`\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | end\n"
     "  | ^^^ expected `]` before this\n"
     "  |\n"
     "  = help: the `[` on line 2 is never closed\n"},
    {"quadratic-missing-sign",
     "min\n"
     " obj: [ x ^ 2 y ^ 2 ] / 2\n"
     "end\n",
     "ERROR:   expected `+`, `-`, or `]`, found `y`\n"
     " --> lp-file-format.lp:2:15\n"
     "  |\n"
     "2 |  obj: [ x ^ 2 y ^ 2 ] / 2\n"
     "  |               ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"quadratic-bad-divisor",
     "min\n"
     " obj: [ x ^ 2 ] / 3\n"
     "end\n",
     "ERROR:   expected `2`, found `3`\n"
     " --> lp-file-format.lp:2:19\n"
     "  |\n"
     "2 |  obj: [ x ^ 2 ] / 3\n"
     "  |                   ^ expected `2` after `/`\n"
     "  |\n"
     "  = help: the quadratic part of the objective can only be divided by 2\n"},
    {"quadratic-missing-divisor",
     "min\n"
     " obj: [ x ^ 2 ] /\n"
     "end\n",
     "ERROR:   expected `2`, found `end`\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | end\n"
     "  | ^^^ expected `2` after `/`\n"
     "  |\n"
     "  = help: the quadratic part of the objective can only be divided by 2\n"},
    {"quadratic-bad-power",
     "min\n"
     " obj: [ x ^ 3 ]\n"
     "end\n",
     "ERROR:   expected `2`, found `3`\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: [ x ^ 3 ]\n"
     "  |             ^ expected `2` after `^`\n"
     "  |\n"
     "  = help: only squared terms like `x ^ 2` are supported\n"},
    {"quadratic-missing-power",
     "min\n"
     " obj: [ x ^ y ]\n"
     "end\n",
     "ERROR:   expected `2`, found `y`\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: [ x ^ y ]\n"
     "  |             ^ expected `2` after `^`\n"
     "  |\n"
     "  = help: only squared terms like `x ^ 2` are supported\n"},
    {"quadratic-linear-term",
     "min\n"
     " obj: [ x ]\n"
     "end\n",
     "ERROR:   expected `^` or `*`, found `]`\n"
     " --> lp-file-format.lp:2:11\n"
     "  |\n"
     "2 |  obj: [ x ]\n"
     "  |           ^ expected `^ 2` or `* <variable>`\n"
     "  |\n"
     "  = help: linear terms must be outside `[` and `]`\n"},
    {"quadratic-expected-variable",
     "min\n"
     " obj: [ 2 * 3 ]\n"
     "end\n",
     "ERROR:   expected a variable, found `3`\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: [ 2 * 3 ]\n"
     "  |             ^ expected a variable in this term\n"},
    {"quadratic-expected-variable-after-times",
     "min\n"
     " obj: [ x * 2 ]\n"
     "end\n",
     "ERROR:   expected a variable, found `2`\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: [ x * 2 ]\n"
     "  |             ^ expected a variable after `*`\n"},
    {"quadratic-constraint",
     "min\n"
     " x\n"
     "st\n"
     " c: x + [ x * y ] >= 1\n"
     "end\n",
     "ERROR:   quadratic constraints are not supported by HiGHS\n"
     " --> lp-file-format.lp:4:9\n"
     "  |\n"
     "4 |  c: x + [ x * y ] >= 1\n"
     "  |         ^ quadratic term in a constraint\n"},
    {"sos-constraint",
     "min\n"
     " x\n"
     "st\n"
     " S1:: x:1 y:2\n"
     "end\n",
     "ERROR:   SOS constraints are not supported by HiGHS\n"
     " --> lp-file-format.lp:4:2\n"
     "  |\n"
     "4 |  S1:: x:1 y:2\n"
     "  |  ^^ SOS constraint\n"},
    {"sos-constraint-named",
     "min\n"
     " x\n"
     "st\n"
     " c: S1:: x:1 y:2\n"
     "end\n",
     "ERROR:   SOS constraints are not supported by HiGHS\n"
     " --> lp-file-format.lp:4:5\n"
     "  |\n"
     "4 |  c: S1:: x:1 y:2\n"
     "  |     ^^ SOS constraint\n"},
    {"indicator-constraint",
     "min\n"
     " x\n"
     "st\n"
     " c: z = 1 -> x >= 1\n"
     "end\n",
     "ERROR:   indicator constraints are not supported by HiGHS\n"
     " --> lp-file-format.lp:4:11\n"
     "  |\n"
     "4 |  c: z = 1 -> x >= 1\n"
     "  |           ^^ indicator constraint\n"},
    {"missing-comparison-keyword",
     "min\n"
     " x\n"
     "st\n"
     " c: x\n"
     "bounds\n"
     "end\n",
     "ERROR:   expected a comparison, found `bounds`\n"
     " --> lp-file-format.lp:5:1\n"
     "  |\n"
     "5 | bounds\n"
     "  | ^^^^^^ expected `<=`, `>=`, or `=` before this\n"
     "  |\n"
     "  = help: `bounds` is a section keyword because it is at the start of a line\n"},
    {"missing-comparison-eof",
     "min\n"
     " x\n"
     "st\n"
     " c: x\n",
     "ERROR:   expected a comparison, found end of file\n"
     " --> lp-file-format.lp:4:6\n"
     "  |\n"
     "4 |  c: x\n"
     "  |      ^ expected `<=`, `>=`, or `=` before this\n"},
    {"missing-right-hand-side",
     "min\n"
     " x\n"
     "st\n"
     " c: x >= y\n"
     "end\n",
     "ERROR:   expected a number, found `y`\n"
     " --> lp-file-format.lp:4:10\n"
     "  |\n"
     "4 |  c: x >= y\n"
     "  |          ^ expected the right-hand side\n"},
    {"right-hand-side-keyword",
     "min\n"
     " x\n"
     "st\n"
     " c: x >= end\n",
     "ERROR:   expected a number, found `end`\n"
     " --> lp-file-format.lp:4:10\n"
     "  |\n"
     "4 |  c: x >= end\n"
     "  |          ^^^ expected the right-hand side\n"
     "  |\n"
     "  = help: section keywords like `end` must be at the start of a line\n"},
    {"right-hand-side-variable",
     "min\n"
     " x\n"
     "st\n"
     " c: 3 x >= 2 y\n"
     "end\n",
     "ERROR:   expected a new line, found `y`\n"
     " --> lp-file-format.lp:4:14\n"
     "  |\n"
     "4 |  c: 3 x >= 2 y\n"
     "  |              ^ unexpected variable on the right-hand side\n"
     "  |\n"
     "  = help: move the variables to the left-hand side\n"},
    {"constraint-trailing-content",
     "min\n"
     " x\n"
     "st\n"
     " c: x + y >= 1 + 1\n"
     "end\n",
     "ERROR:   expected a new line, found `+`\n"
     " --> lp-file-format.lp:4:16\n"
     "  |\n"
     "4 |  c: x + y >= 1 + 1\n"
     "  |                ^ expected the constraint to end before this\n"
     "  |\n"
     "  = help: a constraint ends with a comparison and a number; the next constraint must start on a new line\n"},
    {"bound-missing-comparison",
     "min\n"
     " x\n"
     "bounds\n"
     " x 1\n"
     "end\n",
     "ERROR:   expected `free` or a comparison, found `1`\n"
     " --> lp-file-format.lp:4:4\n"
     "  |\n"
     "4 |  x 1\n"
     "  |    ^ expected `free`, `<=`, `>=`, or `=`\n"},
    {"bound-bad-start",
     "min\n"
     " x\n"
     "bounds\n"
     " <= x\n"
     "end\n",
     "ERROR:   expected a bound, found `<=`\n"
     " --> lp-file-format.lp:4:2\n"
     "  |\n"
     "4 |  <= x\n"
     "  |  ^^ expected a variable or a number\n"},
    {"bound-number-missing-comparison",
     "min\n"
     " x\n"
     "bounds\n"
     " 1 x\n"
     "end\n",
     "ERROR:   expected a comparison, found `x`\n"
     " --> lp-file-format.lp:4:4\n"
     "  |\n"
     "4 |  1 x\n"
     "  |    ^ expected `<=`, `>=`, or `=`\n"},
    {"bound-expected-variable",
     "min\n"
     " x\n"
     "bounds\n"
     " 1 <= 2\n"
     "end\n",
     "ERROR:   expected a variable, found `2`\n"
     " --> lp-file-format.lp:4:7\n"
     "  |\n"
     "4 |  1 <= 2\n"
     "  |       ^ expected a variable\n"},
    {"bound-missing-value",
     "min\n"
     " x\n"
     "bounds\n"
     " x <= y\n"
     "end\n",
     "ERROR:   expected a number, found `y`\n"
     " --> lp-file-format.lp:4:7\n"
     "  |\n"
     "4 |  x <= y\n"
     "  |       ^ expected the value of the bound\n"},
    {"bound-double-equal",
     "min\n"
     " x\n"
     "bounds\n"
     " 1 = x <= 2\n"
     "end\n",
     "ERROR:   a bound on both sides of a variable cannot use `=`\n"
     " --> lp-file-format.lp:4:4\n"
     "  |\n"
     "4 |  1 = x <= 2\n"
     "  |    ^ expected `<=` or `>=`\n"},
    {"bound-double-mismatch-less",
     "min\n"
     " x\n"
     "bounds\n"
     " 1 <= x >= 2\n"
     "end\n",
     "ERROR:   the comparisons in a bound on both sides of a variable must match\n"
     " --> lp-file-format.lp:4:9\n"
     "  |\n"
     "4 |  1 <= x >= 2\n"
     "  |         ^^ expected `<=`\n"},
    {"bound-double-mismatch-greater",
     "min\n"
     " x\n"
     "bounds\n"
     " 2 >= x <= 1\n"
     "end\n",
     "ERROR:   the comparisons in a bound on both sides of a variable must match\n"
     " --> lp-file-format.lp:4:9\n"
     "  |\n"
     "4 |  2 >= x <= 1\n"
     "  |         ^^ expected `>=`\n"},
    {"bound-trailing-content",
     "min\n"
     " x\n"
     "bounds\n"
     " x <= 1 y <= 2\n"
     "end\n",
     "ERROR:   expected a new line, found `y`\n"
     " --> lp-file-format.lp:4:9\n"
     "  |\n"
     "4 |  x <= 1 y <= 2\n"
     "  |         ^ expected the bound to end before this\n"
     "  |\n"
     "  = help: each bound must be on its own line\n"},
    {"general-not-variable",
     "min\n"
     " x\n"
     "general\n"
     " x 1\n"
     "end\n",
     "ERROR:   expected a variable, found `1`\n"
     " --> lp-file-format.lp:4:4\n"
     "  |\n"
     "4 |  x 1\n"
     "  |    ^ expected a variable\n"},
    {"render-gutter",
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "\\ comment\n"
     "min\n"
     " obj: x y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     "  --> lp-file-format.lp:11:9\n"
     "   |\n"
     "11 |  obj: x y\n"
     "   |         ^ expected `+` or `-` before this\n"
     "   |\n"
     "   = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-crlf",
     "min\r\n"
     " obj: x y\r\n"
     "end\r\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x y\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-tab",
     "min\n"
     "\tobj:\tx\ty\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x y\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-long-line",
     "min\n"
     " obj: 1 x1 + 2 x2 + 3 x3 + 4 x4 + 5 x5 + 6 x6 + 7 x7 + 8 x8 + 9 x9 + 10 x10 + 11 x11 + 12 x12 + 13 x13 + 14 x14 y 1 z1 +  2 z2 +  3 z3 +  4 z4 +  5 z5 +  6 z6 +  7 z7 +  8 z8 +  9 z9 +  10 z10 +  11 z11 +  12 z12 +  13 z13 +  14 z14\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:113\n"
     "  |\n"
     "2 | ...x10 + 11 x11 + 12 x12 + 13 x13 + 14 x14 y 1 z1 +  2 z2 +  3 z3 +  4 z4 +  5 z5 + ...\n"
     "  |                                            ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-long-token",
     "min\n"
     " obj: x aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa...`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n"
     "  |         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-utf8",
     "min\n"
     " obj: caf\xC3\xA9 \xC3\xA9t\xC3\xA9\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `\xC3\xA9t\xC3\xA9`\n"
     " --> lp-file-format.lp:2:12\n"
     "  |\n"
     "2 |  obj: caf\xC3\xA9 \xC3\xA9t\xC3\xA9\n"
     "  |            ^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-utf8-truncate",
     "min\n"
     " obj: x a\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `a\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9...`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x a\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\n"
     "  |         ^^^^^^^^^^^^^^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-utf8-window",
     "min\n"
     " obj: x\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 + z\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9...`\n"
     " --> lp-file-format.lp:2:39\n"
     "  |\n"
     "2 | ...\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 + z\n"
     "  |                        ^^^^^^^^^^^^^^^^^^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-utf8-window-left",
     "min\n"
     " obj: a\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 + z\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9x y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:74\n"
     "  |\n"
     "2 | ...\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9x y\n"
     "  |                         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"missing-end-no-newline",
     "min\n"
     " x",
     "ERROR:   expected `end`, found end of file\n"
     " --> lp-file-format.lp:2:3\n"
     "  |\n"
     "2 |  x\n"
     "  |   ^ expected `end` after this\n"
     "  |\n"
     "  = help: an LP file must finish with the keyword `end`\n"},
    // The tokens of "semi-continuous" must be on one line
    {"semi-dash-next-line",
     "min\n"
     " x\n"
     "semi\n"
     "- continuous\n"
     "end\n",
     "ERROR:   expected a variable, found `-`\n"
     " --> lp-file-format.lp:4:1\n"
     "  |\n"
     "4 | - continuous\n"
     "  | ^ expected a variable\n"},
    {"semi-continuous-next-line",
     "min\n"
     " x\n"
     "semi -\n"
     "continuous\n"
     "end\n",
     "ERROR:   expected a variable, found `-`\n"
     " --> lp-file-format.lp:3:6\n"
     "  |\n"
     "3 | semi -\n"
     "  |      ^ expected a variable\n"},
    {"render-utf8-window-right",
     "min\n"
     " obj: x y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9 + z\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9...`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x y\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9...\n"
     "  |         ^^^^^^^^^^^^^^^^^^^^^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"render-gutter-crlf-eof",
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "\\ comment\r\n"
     "min\r\n"
     " obj: x\r\n",
     "ERROR:   expected `end`, found end of file\n"
     "  --> lp-file-format.lp:11:8\n"
     "   |\n"
     "11 |  obj: x\n"
     "   |        ^ expected `end` after this\n"
     "   |\n"
     "   = help: an LP file must finish with the keyword `end`\n"},
    {"delete-character",
     "min\n"
     " x\x7F\n"
     "end\n",
     "ERROR:   unexpected control character 0x7F\n"
     " --> lp-file-format.lp:2:3\n"
     "  |\n"
     "2 |  x\x7F\n"
     "  |   ^ this character is not valid in an LP file\n"},
    {"missing-end-after-number",
     "min\n"
     " obj: x + 1",
     "ERROR:   expected `end`, found end of file\n"
     " --> lp-file-format.lp:2:12\n"
     "  |\n"
     "2 |  obj: x + 1\n"
     "  |            ^ expected `end` after this\n"
     "  |\n"
     "  = help: an LP file must finish with the keyword `end`\n"},
    {"missing-end-after-exponent",
     "min\n"
     " obj: x + 2e",
     "ERROR:   expected `end`, found end of file\n"
     " --> lp-file-format.lp:2:13\n"
     "  |\n"
     "2 |  obj: x + 2e\n"
     "  |             ^ expected `end` after this\n"
     "  |\n"
     "  = help: an LP file must finish with the keyword `end`\n"},
    {"missing-sign-new-line",
     "min\n"
     " obj: x\n"
     " y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:3:2\n"
     "  |\n"
     "3 |  y\n"
     "  |  ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"missing-sign-constant-then-sign",
     "min\n"
     " obj: x 3 + y\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `3`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x 3 + y\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
    {"expected-variable-keyword",
     "min\n"
     " obj: 2 *\n"
     "end\n",
     "ERROR:   expected a variable, found `end`\n"
     " --> lp-file-format.lp:3:1\n"
     "  |\n"
     "3 | end\n"
     "  | ^^^ expected a variable after `*`\n"},
    {"right-hand-side-symbol",
     "min\n"
     " x\n"
     "st\n"
     " c: x >= [\n"
     "end\n",
     "ERROR:   expected a number, found `[`\n"
     " --> lp-file-format.lp:4:10\n"
     "  |\n"
     "4 |  c: x >= [\n"
     "  |          ^ expected the right-hand side\n"},
    {"constraint-equal-variable",
     "min\n"
     " x\n"
     "st\n"
     " c: x = y\n"
     "end\n",
     "ERROR:   expected a number, found `y`\n"
     " --> lp-file-format.lp:4:9\n"
     "  |\n"
     "4 |  c: x = y\n"
     "  |         ^ expected the right-hand side\n"},
    {"bound-not-free",
     "min\n"
     " x\n"
     "bounds\n"
     " x y\n"
     "end\n",
     "ERROR:   expected `free` or a comparison, found `y`\n"
     " --> lp-file-format.lp:4:4\n"
     "  |\n"
     "4 |  x y\n"
     "  |    ^ expected `free`, `<=`, `>=`, or `=`\n"},
    {"semi-dash-number",
     "min\n"
     " x\n"
     "semi-2\n"
     "end\n",
     "ERROR:   expected a variable, found `-`\n"
     " --> lp-file-format.lp:3:5\n"
     "  |\n"
     "3 | semi-2\n"
     "  |     ^ expected a variable\n"},
    {"semi-dash-other",
     "min\n"
     " x\n"
     "semi-cont\n"
     "end\n",
     "ERROR:   expected a variable, found `-`\n"
     " --> lp-file-format.lp:3:5\n"
     "  |\n"
     "3 | semi-cont\n"
     "  |     ^ expected a variable\n"},
    {"render-utf8-window-right-cut",
     "min\n"
     " obj: x y \xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\n"
     "end\n",
     "ERROR:   expected `+` or `-`, found `y`\n"
     " --> lp-file-format.lp:2:9\n"
     "  |\n"
     "2 |  obj: x y \xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9\xC3\xA9...\n"
     "  |         ^ expected `+` or `-` before this\n"
     "  |\n"
     "  = help: terms in an expression must be separated by `+` or `-`\n"},
  };
  for (const DiagnosticCase& c : cases) {
    INFO(c.name);
    Highs highs;
    std::string log;
    REQUIRE(readLp(highs, c.content, log) == HighsStatus::kError);
    CAPTURE(log);
    REQUIRE(log.find(c.diagnostic) != std::string::npos);
    REQUIRE(log.find("Parser error reading " + kLpFile) != std::string::npos);
  }
}

TEST_CASE("lp-file-format-warnings", "[LpFileFormat]") {
  // Each file is valid, but reading it logs warnings
  const std::vector<DiagnosticCase> cases = {
    {"constant",
     "min\n"
     " x\n"
     "st\n"
     " c: x - 1 >= 1\n"
     "end\n",
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:4:7\n"
     "  |\n"
     "4 |  c: x - 1 >= 1\n"
     "  |       ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"},
    {"duplicate",
     "min\n"
     " x\n"
     "st\n"
     " c: x + y + 2 x - x >= 1\n"
     "end\n",
     "WARNING: variable `x` appears more than once in this constraint\n"
     " --> lp-file-format.lp:4:15\n"
     "  |\n"
     "4 |  c: x + y + 2 x - x >= 1\n"
     "  |               ^ repeated here\n"
     "  |\n"
     "  = help: the coefficients are summed\n"},
    {"eleven",
     "min\n"
     " x\n"
     "st\n"
     " c1: x + 1 >= 1\n"
     " c2: x + 2 >= 2\n"
     " c3: x + 3 >= 3\n"
     " c4: x + 4 >= 4\n"
     " c5: x + 5 >= 5\n"
     " c6: x + 6 >= 6\n"
     " c7: x + 7 >= 7\n"
     " c8: x + 8 >= 8\n"
     " c9: x + 9 >= 9\n"
     " c10: x + 10 >= 10\n"
     " c11: x + 11 >= 11\n"
     "end\n",
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:4:8\n"
     "  |\n"
     "4 |  c1: x + 1 >= 1\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:5:8\n"
     "  |\n"
     "5 |  c2: x + 2 >= 2\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:6:8\n"
     "  |\n"
     "6 |  c3: x + 3 >= 3\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:7:8\n"
     "  |\n"
     "7 |  c4: x + 4 >= 4\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:8:8\n"
     "  |\n"
     "8 |  c5: x + 5 >= 5\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:9:8\n"
     "  |\n"
     "9 |  c6: x + 6 >= 6\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:10:8\n"
     "   |\n"
     "10 |  c7: x + 7 >= 7\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:11:8\n"
     "   |\n"
     "11 |  c8: x + 8 >= 8\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:12:8\n"
     "   |\n"
     "12 |  c9: x + 9 >= 9\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:13:9\n"
     "   |\n"
     "13 |  c10: x + 10 >= 10\n"
     "   |         ^^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: 1 more warning was not shown\n"},
    {"twelve",
     "min\n"
     " x\n"
     "st\n"
     " c1: x + 1 >= 1\n"
     " c2: x + 2 >= 2\n"
     " c3: x + 3 >= 3\n"
     " c4: x + 4 >= 4\n"
     " c5: x + 5 >= 5\n"
     " c6: x + 6 >= 6\n"
     " c7: x + 7 >= 7\n"
     " c8: x + 8 >= 8\n"
     " c9: x + 9 >= 9\n"
     " c10: x + 10 >= 10\n"
     " c11: x + 11 >= 11\n"
     " c12: x + 12 >= 12\n"
     "end\n",
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:4:8\n"
     "  |\n"
     "4 |  c1: x + 1 >= 1\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:5:8\n"
     "  |\n"
     "5 |  c2: x + 2 >= 2\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:6:8\n"
     "  |\n"
     "6 |  c3: x + 3 >= 3\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:7:8\n"
     "  |\n"
     "7 |  c4: x + 4 >= 4\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:8:8\n"
     "  |\n"
     "8 |  c5: x + 5 >= 5\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:9:8\n"
     "  |\n"
     "9 |  c6: x + 6 >= 6\n"
     "  |        ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:10:8\n"
     "   |\n"
     "10 |  c7: x + 7 >= 7\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:11:8\n"
     "   |\n"
     "11 |  c8: x + 8 >= 8\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:12:8\n"
     "   |\n"
     "12 |  c9: x + 9 >= 9\n"
     "   |        ^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     "  --> lp-file-format.lp:13:9\n"
     "   |\n"
     "13 |  c10: x + 10 >= 10\n"
     "   |         ^^^^ this constant is ignored\n"
     "   |\n"
     "   = help: move the constant to the right-hand side\n"
     "WARNING: 2 more warnings were not shown\n"},
    {"two-constants",
     "min\n"
     " x\n"
     "st\n"
     " c: x - 1 + 2 >= 1\n"
     "end\n",
     "WARNING: a constant on the left-hand side of a constraint is ignored\n"
     " --> lp-file-format.lp:4:7\n"
     "  |\n"
     "4 |  c: x - 1 + 2 >= 1\n"
     "  |       ^^^ this constant is ignored\n"
     "  |\n"
     "  = help: move the constant to the right-hand side\n"},
  };
  for (const DiagnosticCase& c : cases) {
    INFO(c.name);
    Highs highs;
    std::string log;
    REQUIRE(readLp(highs, c.content, log) == HighsStatus::kWarning);
    CAPTURE(log);
    REQUIRE(log.find(c.diagnostic) != std::string::npos);
  }
}

TEST_CASE("lp-file-format-warning-semantics", "[LpFileFormat]") {
  Highs highs;
  // The constant is ignored, and the repeated variables are summed
  REQUIRE(readLp(highs,
                 "min\n x\nst\n c: 2 x - 1 + y + x - 2 y >= 1\nend\n") ==
          HighsStatus::kWarning);
  const HighsLp& lp = highs.getLp();
  REQUIRE(lp.row_lower_ == std::vector<double>{1});
  REQUIRE(lp.a_matrix_.index_ == std::vector<HighsInt>{0, 0});
  REQUIRE(lp.a_matrix_.value_ == std::vector<double>{3, -1});
  // A zero constant is not worth a warning
  REQUIRE(readLp(highs, "min\n x\nst\n c: x + 0 >= 1\nend\n") ==
          HighsStatus::kOk);
}

TEST_CASE("lp-file-format-model", "[LpFileFormat]") {
  const std::string content =
      "\\ A model that uses most of the LP file format\n"
      "MAXIMIZE\n"
      " profit: 2 x + 3.5 y - z + 1.5 x + 10 - 2.5 + .5e1 w\n"
      "SUBJECT TO\n"
      " c1: x + y < 4\n"
      " c2: x - y > -1\n"
      " 55: x + w = 2\n"
      " 9.9: - - y + + z =< 5\n"
      " c5: 2 * x + 0 z => 1\n"
      " x + y + z == 3\n"
      " c7: >= -2\n"
      " c8: x\n"
      "   + y\n"
      "   <=\n"
      "   1E1\n"
      "bounds\n"
      " x <= 10\n"
      " y >= -5\n"
      " -1 <= z <= 1\n"
      " 3 >= w >= 1\n"
      " v Free\n"
      " u = 4\n"
      " t >= -INF\n"
      " -infinity <= s\n"
      " 7 >= r\n"
      " 2 = q\n"
      " p <= +inf\n"
      " 2 <= sc <= 10\n"
      " b1 <= 5\n"
      "general\n"
      " r q\n"
      "binary\n"
      " b1 b2\n"
      "semi-continuous\n"
      " sc\n"
      "integers\n"
      " sc\n"
      "end\n";
  Highs highs;
  REQUIRE(readLp(highs, content) == HighsStatus::kOk);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = 14;
  REQUIRE(lp.num_col_ == num_col);
  REQUIRE(lp.num_row_ == 8);
  REQUIRE(lp.col_names_ ==
          std::vector<std::string>{"x", "y", "z", "w", "v", "u", "t", "s",
                                   "r", "q", "p", "sc", "b1", "b2"});
  REQUIRE(lp.row_names_ == std::vector<std::string>{"c1", "c2", "55", "9.9",
                                                    "c5", "HiGHS_R5", "c7",
                                                    "c8"});
  REQUIRE(lp.objective_name_ == "profit");
  REQUIRE(lp.sense_ == ObjSense::kMaximize);
  REQUIRE(lp.offset_ == 7.5);
  std::vector<double> cost(num_col, 0);
  cost[0] = 3.5;
  cost[1] = 3.5;
  cost[2] = -1;
  cost[3] = 5;
  REQUIRE(lp.col_cost_ == cost);
  REQUIRE(lp.col_lower_ == std::vector<double>{0, -5, -1, 1, -inf, 4, -inf,
                                               -inf, 0, 2, 0, 2, 0, 0});
  REQUIRE(lp.col_upper_ == std::vector<double>{10, inf, 1, 3, inf, 4, inf,
                                               inf, 7, 2, inf, 10, 5, 1});
  std::vector<HighsVarType> integrality(num_col, HighsVarType::kContinuous);
  integrality[8] = HighsVarType::kInteger;
  integrality[9] = HighsVarType::kInteger;
  integrality[11] = HighsVarType::kSemiInteger;
  integrality[12] = HighsVarType::kInteger;
  integrality[13] = HighsVarType::kInteger;
  REQUIRE(lp.integrality_ == integrality);
  REQUIRE(lp.row_lower_ ==
          std::vector<double>{-inf, -1, 2, -inf, 1, 3, -2, -inf});
  REQUIRE(lp.row_upper_ == std::vector<double>{4, inf, 2, 5, inf, 3, inf, 10});
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  REQUIRE(matrix.isColwise());
  std::vector<HighsInt> start(num_col + 1, 14);
  start[0] = 0;
  start[1] = 6;
  start[2] = 11;
  start[3] = 13;
  REQUIRE(matrix.start_ == start);
  REQUIRE(matrix.index_ ==
          std::vector<HighsInt>{0, 1, 2, 4, 5, 7, 0, 1, 3, 5, 7, 3, 5, 2});
  REQUIRE(matrix.value_ == std::vector<double>{1, 1, 1, 2, 1, 1, 1, -1, 1, 1,
                                               1, 1, 1, 1});
  REQUIRE(!highs.getModel().isQp());
}

TEST_CASE("lp-file-format-objective-and-bounds", "[LpFileFormat]") {
  Highs highs;
  // An objective without a name, a constraint on one variable, and bounds
  // that start with infinity
  REQUIRE(readLp(highs,
                 "min\n - x + 1\nst\n c: x = 1\nbounds\n inf >= x\n"
                 " infinity >= y\n x >= -1\nend\n") == HighsStatus::kOk);
  const HighsLp& lp = highs.getLp();
  REQUIRE(lp.objective_name_.empty());
  REQUIRE(lp.col_cost_ == std::vector<double>{-1, 0});
  REQUIRE(lp.offset_ == 1);
  REQUIRE(lp.row_lower_ == std::vector<double>{1});
  REQUIRE(lp.row_upper_ == std::vector<double>{1});
  REQUIRE(lp.col_lower_ == std::vector<double>{-1, 0});
  REQUIRE(lp.col_upper_ == std::vector<double>{inf, inf});
}

TEST_CASE("lp-file-format-hessian", "[LpFileFormat]") {
  Highs highs;
  // With "/ 2" the bracketed terms are 0.5 x'Qx, and without it they are
  // x'Qx. The terms in w cancel, so they are not in the Hessian
  REQUIRE(readLp(highs,
                 "min\n"
                 " obj: x + [ 2 x ^ 2 + 4 x * y - y * x + 3 y ^ 2 ] / 2\n"
                 " + [ z * z ] - [ w ^ 2 ] + [ w ^ 2 + 2 * w * x - 2 x * w ]\n"
                 "end\n") == HighsStatus::kOk);
  const HighsHessian& hessian = highs.getModel().hessian_;
  REQUIRE(hessian.dim_ == 4);
  // HiGHS stores the lower triangle, with an explicit zero for w
  REQUIRE(hessian.start_ == std::vector<HighsInt>{0, 2, 3, 4, 5});
  REQUIRE(hessian.index_ == std::vector<HighsInt>{0, 1, 1, 2, 3});
  REQUIRE(hessian.value_ == std::vector<double>{2, 1.5, 3, 2, 0});
  REQUIRE(highs.getLp().col_cost_ == std::vector<double>{1, 0, 0, 0});

  // When all of the quadratic terms cancel, there is no Hessian
  REQUIRE(readLp(highs, "min\n obj: x + [ x ^ 2 - x * x ] / 2\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(!highs.getModel().isQp());

  // Brackets can be empty
  REQUIRE(readLp(highs, "min\n obj: x + [ ] / 2\nend\n") == HighsStatus::kOk);
  REQUIRE(!highs.getModel().isQp());
}

TEST_CASE("lp-file-format-keywords", "[LpFileFormat]") {
  Highs highs;
  for (const std::string& keyword :
       {"minimize", "Minimise", "MINIMUM", "min"}) {
    REQUIRE(readLp(highs, keyword + "\n obj: x\nend\n") == HighsStatus::kOk);
    REQUIRE(highs.getLp().sense_ == ObjSense::kMinimize);
  }
  for (const std::string& keyword :
       {"maximize", "Maximise", "MAXIMUM", "max"}) {
    REQUIRE(readLp(highs, keyword + "\n obj: x\nend\n") == HighsStatus::kOk);
    REQUIRE(highs.getLp().sense_ == ObjSense::kMaximize);
  }
  for (const std::string& keyword :
       {"subject to", "Such That", "st", "S.T.", "st."}) {
    REQUIRE(readLp(highs, "min\n obj: x\n" + keyword + "\n c: x >= 1\nend\n") ==
            HighsStatus::kOk);
    REQUIRE(highs.getLp().num_row_ == 1);
  }
  for (const std::string& keyword : {"bounds", "BOUND"}) {
    REQUIRE(readLp(highs, "min\n obj: x\n" + keyword + "\n x <= 1\nend\n") ==
            HighsStatus::kOk);
    REQUIRE(highs.getLp().col_upper_[0] == 1);
  }
  for (const std::string& keyword :
       {"general", "Generals", "GEN", "integer", "integers"}) {
    REQUIRE(readLp(highs, "min\n obj: x\n" + keyword + "\n x\nend\n") ==
            HighsStatus::kOk);
    REQUIRE(highs.getLp().integrality_[0] == HighsVarType::kInteger);
  }
  for (const std::string& keyword : {"binary", "Binaries", "BIN"}) {
    REQUIRE(readLp(highs, "min\n obj: x\n" + keyword + "\n x\nend\n") ==
            HighsStatus::kOk);
    REQUIRE(highs.getLp().integrality_[0] == HighsVarType::kInteger);
    REQUIRE(highs.getLp().col_upper_[0] == 1);
  }
  for (const std::string& keyword :
       {"semi-continuous", "Semi-Continuous", "semi", "SEMIS",
        "semi - continuous", "semi -continuous"}) {
    REQUIRE(readLp(highs, "min\n obj: x\nbounds\n 1 <= x <= 2\n" + keyword +
                              "\n x\nend\n") == HighsStatus::kOk);
    REQUIRE(highs.getLp().integrality_[0] == HighsVarType::kSemiContinuous);
  }
  // Keywords followed by ":" are names, so this is a constraint named "end"
  REQUIRE(readLp(highs, "min\n obj: x\nst\n end: x >= 1\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().row_names_ == std::vector<std::string>{"end"});
  // "subject" and "such" are only keywords when followed by "to" and "that"
  REQUIRE(readLp(highs,
                 "min\n obj: x\nst\n subject: x >= 1\n such: x >= 2\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().row_names_ ==
          std::vector<std::string>{"subject", "such"});
  REQUIRE(readLp(highs,
                 "min\n obj: x\nst\nsubject + x >= 1\nsuch + x >= 2\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().col_names_ ==
          std::vector<std::string>{"x", "subject", "such"});
  // A keyword that isn't at the start of a line is a name
  REQUIRE(readLp(highs, "min\n obj: x + end + bounds\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().col_names_ ==
          std::vector<std::string>{"x", "end", "bounds"});
  // Sections can appear in any order, and be repeated, apart from the
  // objective
  REQUIRE(readLp(highs,
                 "bounds\n x <= 1\nst\n c1: x + y >= 1\nmin\n obj: x\n"
                 "st\n c2: y >= 0\nbounds\n y <= 2\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().num_row_ == 2);
  REQUIRE(highs.getLp().col_upper_ == std::vector<double>{1, 2});
}

TEST_CASE("lp-file-format-lexer", "[LpFileFormat]") {
  Highs highs;
  // Comments, blank lines, tabs, CRLF line endings, and unusual names
  REQUIRE(readLp(highs,
                 "\\ comment\r\n"
                 "\r\n"
                 "min \\ comment after a keyword\r\n"
                 "\tobj:\t...100 + x.y_z(1){2}!\"#$%&;,?@`'|~ + caf\xC3\xA9\r\n"
                 "st\r\n"
                 " c: 2e3 x1 + 2E-1 x2 + 1. x3 + 2e + x4 >= 1e+2\r\n"
                 "end") == HighsStatus::kOk);
  const HighsLp& lp = highs.getLp();
  REQUIRE(lp.col_names_ ==
          std::vector<std::string>{"...100", "x.y_z(1){2}!\"#$%&;,?@`'|~",
                                   "caf\xC3\xA9", "x1", "x2", "x3", "e", "x4"});
  REQUIRE(lp.a_matrix_.value_ == std::vector<double>{2000, 0.2, 1, 2, 1});
  REQUIRE(lp.row_lower_ == std::vector<double>{100});
  // An empty file and a file with only comments are empty models
  for (const std::string& content : {"", " \n\t\n", "\\ comment"}) {
    REQUIRE(readLp(highs, content) == HighsStatus::kOk);
    REQUIRE(highs.getLp().num_col_ == 0);
    REQUIRE(highs.getLp().num_row_ == 0);
  }
}

TEST_CASE("lp-file-format-row-names", "[LpFileFormat]") {
  Highs highs;
  // A name beginning HiGHS_R is fine if there are no unnamed rows
  REQUIRE(readLp(highs, "min\n obj: x\nst\n HiGHS_R7: x >= 1\nend\n") ==
          HighsStatus::kOk);
  REQUIRE(highs.getLp().row_names_ == std::vector<std::string>{"HiGHS_R7"});
  // Otherwise the generated names might clash, so all row names are removed
  std::string log;
  REQUIRE(readLp(highs,
                 "min\n obj: x\nst\n HiGHS_R7: x >= 1\n x >= 2\nend\n",
                 log) == HighsStatus::kOk);
  REQUIRE(log.find("WARNING: Cannot create row name beginning \"HiGHS_R\"") !=
          std::string::npos);
  REQUIRE(highs.getLp().num_row_ == 2);
  REQUIRE(highs.getLp().row_names_.empty());
}

TEST_CASE("lp-file-format-file-not-found", "[LpFileFormat]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel("does-not-exist.lp") == HighsStatus::kError);
}

TEST_CASE("lp-file-format-integrality", "[LpFileFormat]") {
  // Semi-continuous variables that are also general or binary are
  // semi-integer, whatever order the sections are in
  Highs highs;
  REQUIRE(readLp(highs,
                 "min\n obj: a + b + c\nbounds\n 1 <= a <= 2\n 1 <= b <= 2\n"
                 " 1 <= c <= 2\ngeneral\n a\nsemi\n a b c\ngeneral\n b c\n"
                 "binary\n b\nsemi\n c\nend\n") == HighsStatus::kOk);
  REQUIRE(highs.getLp().integrality_ ==
          std::vector<HighsVarType>(3, HighsVarType::kSemiInteger));
  REQUIRE(highs.getLp().col_upper_ == std::vector<double>{2, 2, 2});
}

// A file that is larger than the chunks that the reader reads, with tokens
// that cross the boundaries between chunks
static std::string largeLpFile(const HighsInt num_row,
                               const std::string& extra) {
  std::string content = "min\n obj: y\nst\n";
  for (HighsInt i = 0; i < num_row; i++)
    content += " c" + std::to_string(i) + ": " + std::to_string(i + 1) + " x" +
               std::to_string(i % 1000) + " + y >= " + std::to_string(i) +
               "\n";
  return content + extra + "end\n";
}

TEST_CASE("lp-file-format-large-file", "[LpFileFormat]") {
  const HighsInt num_row = 100000;
  Highs highs;
  REQUIRE(readLp(highs, largeLpFile(num_row, "")) == HighsStatus::kOk);
  const HighsLp& lp = highs.getLp();
  REQUIRE(lp.num_row_ == num_row);
  REQUIRE(lp.num_col_ == 1001);
  // Check every value, by summing them
  double sum_value = 0;
  for (const double value : lp.a_matrix_.value_) sum_value += value;
  REQUIRE(sum_value == double(num_row) * (num_row + 1) / 2 + num_row);
  double sum_lower = 0;
  for (const double lower : lp.row_lower_) sum_lower += lower;
  REQUIRE(sum_lower == double(num_row) * (num_row - 1) / 2);
  REQUIRE(lp.row_names_.back() == "c" + std::to_string(num_row - 1));

  // Diagnostics near the end of the file show the right line
  std::string log;
  REQUIRE(readLp(highs, largeLpFile(num_row, " c: x y >= 1\n"), log) ==
          HighsStatus::kError);
  REQUIRE(log.find("ERROR:   expected `+` or `-`, found `y`\n"
                   "      --> lp-file-format.lp:100004:7\n"
                   "       |\n"
                   "100004 |  c: x y >= 1\n"
                   "       |       ^ expected `+` or `-` before this\n") !=
          std::string::npos);
  REQUIRE(readLp(highs, largeLpFile(num_row, " d: x0 + 1 >= 1\n"), log) ==
          HighsStatus::kWarning);
  REQUIRE(log.find("WARNING: a constant on the left-hand side of a "
                   "constraint is ignored\n"
                   "      --> lp-file-format.lp:100004:8\n"
                   "       |\n"
                   "100004 |  d: x0 + 1 >= 1\n"
                   "       |        ^^^ this constant is ignored\n") !=
          std::string::npos);

  // A name that is longer than a chunk
  const std::string name(1500000, 'x');
  REQUIRE(readLp(highs, "min\n obj: " + name + "\nend\n") == HighsStatus::kOk);
  REQUIRE(highs.getLp().col_names_ == std::vector<std::string>{name});
}
