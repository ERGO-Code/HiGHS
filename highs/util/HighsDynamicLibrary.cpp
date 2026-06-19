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

#include "HConfig.h"

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
  handle_ = static_cast<void*>(LoadLibraryA(full_path.c_str()));
  if (!handle_) {
    DWORD error = GetLastError();
    char* message = nullptr;
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                       FORMAT_MESSAGE_IGNORE_INSERTS,
                   nullptr, error, 0, reinterpret_cast<LPSTR>(&message), 0,
                   nullptr);
    status_ = "Failed to load " + filename + ": " +
              (message ? message : "Unknown error");
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

#ifdef HIPO_USES_OPENBLAS
  // Openblas creates a thread pool that may still exist after dlclose-ing the
  // library, leading to seg fault. blas_shutdown should prevent this, but its
  // symbol is not always exposed. If the symbol exists, use it and then call
  // dlclose. If the symbol does not exist, avoid calling dlclose.

  using shutdown_t = void (*)();
  shutdown_t shutdown_fn = nullptr;
  if (auto p = dlsym(handle_, "blas_shutdown"))
    shutdown_fn = reinterpret_cast<shutdown_t>(p);
  else if (auto p = dlsym(handle_, "openblas_shutdown"))
    shutdown_fn = reinterpret_cast<shutdown_t>(p);

  if (shutdown_fn) {
    shutdown_fn();
#if defined(_WIN32) || defined(_WIN64)
    FreeLibrary(static_cast<HMODULE>(handle_));
#else
    dlclose(handle_);
#endif
  }

#else

#if defined(_WIN32) || defined(_WIN64)
  FreeLibrary(static_cast<HMODULE>(handle_));
#else
  dlclose(handle_);
#endif
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