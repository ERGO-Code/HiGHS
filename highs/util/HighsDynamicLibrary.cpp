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

#include "util/stringutil.h"

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOGDI
#define NOGDI
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
// any address inside this module, to find the module's directory
static void moduleAnchor() {}
#endif

HighsDynamicLibrary::~HighsDynamicLibrary() { unload(); }

bool HighsDynamicLibrary::load(const std::string& filename,
                               const std::string& path) {
  unload();

#if defined(_WIN32) || defined(_WIN64)
  std::string absolute_filename = path + (path.empty() ? "" : "\\") + filename;

  // allow application directory if no path is specified, otherwise use DLL load
  // directory and absolute path
  DWORD flags = (path.empty() ? LOAD_LIBRARY_SEARCH_APPLICATION_DIR
                              : LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR);
  handle_ =
      static_cast<void*>(LoadLibraryExA(absolute_filename.c_str(), nullptr,
                                        LOAD_LIBRARY_SEARCH_SYSTEM32 | flags));

  // then the directory of the module holding HiGHS, like @loader_path on
  // macOS, e.g. runtimes/<rid>/native for a .NET app
  HMODULE self = nullptr;
  char self_path[MAX_PATH];
  if (!handle_ && path.empty() &&
      GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                             GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                         reinterpret_cast<LPCSTR>(&moduleAnchor), &self)) {
    const DWORD len = GetModuleFileNameA(self, self_path, MAX_PATH);
    if (len > 0 && len < MAX_PATH) {
      std::string self_dir(self_path, len);
      self_dir.resize(self_dir.find_last_of('\\') + 1);
      handle_ = static_cast<void*>(
          LoadLibraryExA((self_dir + filename).c_str(), nullptr,
                         LOAD_LIBRARY_SEARCH_SYSTEM32 |
                             LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR));
    }
  }

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