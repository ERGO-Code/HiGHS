/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsDynamicLibrary.cpp
 * @brief Lightweight cross-platform runtime shared library loader
 */

#include "util/HighsDynamicLibrary.h"

#include <algorithm>

#include "util/stringutil.h"

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

HighsDynamicLibrary::~HighsDynamicLibrary() { unload(); }

bool HighsDynamicLibrary::load(const std::string& filename,
                               const std::string& path) {
  unload();

#if defined(_WIN32) || defined(_WIN64)
  std::string full_path = path + (path.empty() ? "" : "\\") + filename;
  // Python can give a mix of separators, and LoadLibraryExW needs backslashes
  std::replace(full_path.begin(), full_path.end(), '/', '\\');
  const int wide_length =
      MultiByteToWideChar(CP_UTF8, 0, full_path.c_str(), -1, nullptr, 0);
  std::wstring wide_path(wide_length, L'\0');
  MultiByteToWideChar(CP_UTF8, 0, full_path.c_str(), -1, &wide_path[0],
                      wide_length);
  // Never search the working directory or PATH, which an attacker may control.
  // DLL_LOAD_DIR finds dependencies beside a DLL given by absolute path; it is
  // rejected for a bare name or a relative path.
  const DWORD search_flags =
      LOAD_LIBRARY_SEARCH_DEFAULT_DIRS |
      (path.empty() ? 0 : LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR);
  handle_ = static_cast<void*>(
      LoadLibraryExW(wide_path.c_str(), nullptr, search_flags));
  if (!handle_) {
    DWORD error = GetLastError();
    char* message = nullptr;
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                       FORMAT_MESSAGE_IGNORE_INSERTS,
                   nullptr, error, 0, reinterpret_cast<LPSTR>(&message), 0,
                   nullptr);
    status_ = "Failed to load " + filename + ": " +
              (message ? message : "Unknown error");
    rtrim(status_);  // FormatMessage ends its text with a line break
    if (message) LocalFree(message);
    return false;
  }
#else
#if defined(__APPLE__)
  std::string full_path =
      path + (path.empty() ? "@loader_path/" : "/") + filename;
#else
  std::string full_path = path + (path.empty() ? "" : "/") + filename;
#endif

  handle_ = dlopen(full_path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!handle_) {
    const char* error = dlerror();
    status_ =
        "Failed to load " + filename + ": " + (error ? error : "Unknown error");
    return false;
  }
#endif

  status_ = "Successfully loaded " + filename;
  return true;
}

void HighsDynamicLibrary::unload() {
  if (!handle_) return;

#if defined(_WIN32) || defined(_WIN64)
  FreeLibrary(static_cast<HMODULE>(handle_));
#else
  dlclose(handle_);
#endif

  handle_ = nullptr;
}

void* HighsDynamicLibrary::resolveRaw(const char* name) const {
#if defined(_WIN32) || defined(_WIN64)
  return reinterpret_cast<void*>(
      GetProcAddress(static_cast<HMODULE>(handle_), name));
#else
  return dlsym(handle_, name);
#endif
}