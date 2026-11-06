// Runs as its own process because HighsExternalApi::tryLoad only acts once.
// On Windows, argv[1] is a folder holding planted_extras.dll, a name that
// nothing else in the build uses: the folder of this executable is searched.
#include <algorithm>
#include <iostream>
#include <string>

#include "HighsExternalApi.h"

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

static int failures = 0;

static void check(const bool ok, const std::string& what) {
  std::cout << (ok ? "ok: " : "FAILED: ") << what << std::endl;
  if (!ok) ++failures;
}

#if defined(_WIN32) || defined(_WIN64)
static std::wstring toWide(const std::string& utf8) {
  const int len = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, nullptr, 0);
  std::wstring wide(len, L'\0');
  MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, &wide[0], len);
  return wide;
}

// shortcut: argv[1] is in the ANSI code page but used as UTF-8, so this fails
// in a non-ASCII build folder; use GetCommandLineW if that ever matters
static void checkWindowsSearch(const std::string& planted_dir) {
  const std::string name = "planted_extras.dll";
  check(SetCurrentDirectoryA(planted_dir.c_str()) != 0, "setup: working dir");
  check(SetEnvironmentVariableA("PATH", planted_dir.c_str()) != 0,
        "setup: PATH");

  HighsDynamicLibrary library;
  check(!library.load(name), "bare name ignores working directory and PATH");
  check(library.load(name, planted_dir), "absolute path loads");
  library.unload();

  std::string forward_slash_dir = planted_dir;
  std::replace(forward_slash_dir.begin(), forward_slash_dir.end(), '\\', '/');
  check(library.load(name, forward_slash_dir), "forward slashes load");
  library.unload();

  // e-acute, euro sign: the path decodes differently in UTF-8 and ANSI
  const std::string utf8_dir = planted_dir + "\\extras-\xC3\xA9\xE2\x82\xAC";
  CreateDirectoryW(toWide(utf8_dir).c_str(), nullptr);
  CopyFileW(toWide(planted_dir + "\\" + name).c_str(),
            toWide(utf8_dir + "\\" + name).c_str(), FALSE);
  check(library.load(name, utf8_dir), "non-ASCII UTF-8 path loads");
}
#endif

int main(int argc, char** argv) {
#if defined(_WIN32) || defined(_WIN64)
  if (argc != 2) return 2;
  checkWindowsSearch(argv[1]);
#else
  (void)argc;
  (void)argv;
#endif

  // Like highspy when highspy_extras is absent, then a later feature check.
  // Any attempt to load a library would give a different status.
  check(!HighsExternalApi::tryLoad("", false), "no-search load is unavailable");
  check(!HighsExternalApi::tryLoad(), "later default load is unavailable");
  const std::string status = HighsExternalApi::getLoadStatus();
  check(status == "Extras: no extras package found",
        "no library load was attempted, status: " + status);

  return failures;
}
