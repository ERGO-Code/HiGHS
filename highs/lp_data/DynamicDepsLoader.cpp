/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/DynamicDepsLoader.cpp
 * @brief Dynamic loader for the optional HiPO library
 */

#include "DynamicDepsLoader.h"

#include <cstring>
#include <vector>

// Platform-specific includes for dynamic loading
#if defined(_WIN32) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define PATH_SEPARATOR "\\"
#else
#include <dlfcn.h>
#define PATH_SEPARATOR "/"
#endif

DynamicDepsLoader& DynamicDepsLoader::instance() {
  static DynamicDepsLoader loader;
  return loader;
}

DynamicDepsLoader::DynamicDepsLoader() = default;

DynamicDepsLoader::~DynamicDepsLoader() { unloadLibrary(); }

bool DynamicDepsLoader::isAvailable() {
  return available_;
}

std::string DynamicDepsLoader::getVersion() const { return version_; }

std::string DynamicDepsLoader::getLibraryFilename() const {
#if defined(_WIN32) || defined(_WIN64)
  return "highs_hipo.dll";
#elif defined(__APPLE__)
  return "libhighs_hipo.dylib";
#else
  return "libhighs_hipo.so";
#endif
}

bool DynamicDepsLoader::loadLibrary(const std::string& path) {
#if defined(_WIN32) || defined(_WIN64)
  lib_handle_ = static_cast<void*>(LoadLibraryA(path.c_str()));
  if (!lib_handle_) {
    DWORD error = GetLastError();
    char* msg = nullptr;
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
                   nullptr, error, 0, reinterpret_cast<LPSTR>(&msg), 0,
                   nullptr);
    last_error_ = "Failed to load " + path + ": " + (msg ? msg : "Unknown error");
    if (msg) LocalFree(msg);
    return false;
  }
#else
  lib_handle_ = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!lib_handle_) {
    const char* err = dlerror();
    last_error_ = "Failed to load " + path + ": " + (err ? err : "Unknown error");
    return false;
  }
#endif
  return true;
}

void DynamicDepsLoader::unloadLibrary() {
  if (lib_handle_) {
#if defined(_WIN32) || defined(_WIN64)
    FreeLibrary(static_cast<HMODULE>(lib_handle_));
#else
    dlclose(lib_handle_);
#endif
    lib_handle_ = nullptr;
  }
  fn_get_abi_version_ = nullptr;
  fn_get_version_ = nullptr;
  fn_solve_lp_ = nullptr;
}

void* DynamicDepsLoader::resolveSymbol(const char* name) {
#if defined(_WIN32) || defined(_WIN64)
  return reinterpret_cast<void*>(
      GetProcAddress(static_cast<HMODULE>(lib_handle_), name));
#else
  return dlsym(lib_handle_, name);
#endif
}

bool DynamicDepsLoader::resolveFunctions() {
  if (!lib_handle_) return false;

  fn_get_abi_version_ =
      reinterpret_cast<hipo_get_abi_version_t>(resolveSymbol("hipo_get_abi_version"));
  fn_get_version_ =
      reinterpret_cast<hipo_get_version_t>(resolveSymbol("hipo_get_version"));
  fn_solve_lp_ =
      reinterpret_cast<hipo_solve_lp_t>(resolveSymbol("hipo_solve_lp"));

  if (!fn_get_abi_version_ || !fn_get_version_ || !fn_solve_lp_) {
    last_error_ = "Failed to resolve required HiPO functions";
    return false;
  }

  return true;
}

bool DynamicDepsLoader::tryLoad(const std::string path) {
  if (path.empty()) {
    last_error_ = "HiPO not available. Install with: pip install highspy[hipo]";
    return false;
  }

  if (!initialized_) {
    initialized_ = true;
    available_ = false;

    if (loadLibrary(path + PATH_SEPARATOR + getLibraryFilename())) {
      if (resolveFunctions()) {
        // Check ABI compatibility
        int loaded_abi_version = fn_get_abi_version_();
        if (loaded_abi_version != kHipoAbiVersion) {
          last_error_ =
              "HiPO ABI version mismatch: expected " +
              std::to_string(kHipoAbiVersion) + ", got " +
              std::to_string(loaded_abi_version) +
              ". Please reinstall: pip install --force-reinstall highspy[hipo]";
          unloadLibrary();
          return false;
        }

        // Get version string
        const char* ver = fn_get_version_();
        version_ = ver ? ver : "";
        available_ = true;
        return true;
      }
      unloadLibrary();
    }

    if (last_error_.empty()) {
      last_error_ =
          "HiPO not available. Install with: pip install highspy[hipo]";
    }
  } else
    return available_;

  return false;
}

HighsStatus DynamicDepsLoader::solveLp(HighsLpSolverObject& solver_object) {
  std::cout << "Using HiPO version: " << getVersion() << std::endl;

  if (!isAvailable()) {
    return HighsStatus::kError;
  }

  // Call the dynamically loaded solve function with actual types
  return fn_solve_lp_(
      solver_object.options_,
      solver_object.timer_,
      solver_object.lp_,
      solver_object.basis_,
      solver_object.solution_,
      solver_object.model_status_,
      solver_object.highs_info_,
      solver_object.callback_);
}
