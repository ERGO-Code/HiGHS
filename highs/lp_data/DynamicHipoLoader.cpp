/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/DynamicHipoLoader.cpp
 * @brief Dynamic loader for the optional HiPO library
 */

#include "lp_data/DynamicHipoLoader.h"

#include <cstring>
#include <vector>

#include "lp_data/HighsLpSolverObject.h"

// Platform-specific includes for dynamic loading
#if defined(_WIN32) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define PATH_SEPARATOR "\\"
#else
#include <dlfcn.h>
#define PATH_SEPARATOR "/"
#endif

DynamicHipoLoader& DynamicHipoLoader::instance() {
  static DynamicHipoLoader loader;
  return loader;
}

DynamicHipoLoader::DynamicHipoLoader() = default;

DynamicHipoLoader::~DynamicHipoLoader() { unloadLibrary(); }

bool DynamicHipoLoader::isAvailable() {
  if (!initialized_) {
    initialized_ = true;
    available_ = tryLoad();
  }
  return available_;
}

std::string DynamicHipoLoader::getVersion() const { return version_; }

std::string DynamicHipoLoader::getLibraryFilename() const {
#if defined(_WIN32) || defined(_WIN64)
  return "highs_hipo.dll";
#elif defined(__APPLE__)
  return "libhighs_hipo.dylib";
#else
  return "libhighs_hipo.so";
#endif
}

std::vector<std::string> DynamicHipoLoader::getSearchPaths() const {
  std::vector<std::string> paths;
  const std::string lib_name = getLibraryFilename();

  // 1. Explicit path via environment variable (for testing/advanced users)
  const char* explicit_path = std::getenv("HIGHS_HIPO_LIBRARY");
  if (explicit_path && std::strlen(explicit_path) > 0) {
    paths.push_back(std::string(explicit_path));
  }

  // 2. Python highspy_hipo package location
  //    Set by highspy_hipo.__init__ when the package is imported
  const char* pkg_path = std::getenv("HIGHSPY_HIPO_LIBRARY_PATH");
  if (pkg_path && std::strlen(pkg_path) > 0) {
    paths.push_back(std::string(pkg_path) + PATH_SEPARATOR + lib_name);
  }

  return paths;
}

bool DynamicHipoLoader::loadLibrary(const std::string& path) {
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

void DynamicHipoLoader::unloadLibrary() {
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

void* DynamicHipoLoader::resolveSymbol(const char* name) {
#if defined(_WIN32) || defined(_WIN64)
  return reinterpret_cast<void*>(
      GetProcAddress(static_cast<HMODULE>(lib_handle_), name));
#else
  return dlsym(lib_handle_, name);
#endif
}

bool DynamicHipoLoader::resolveFunctions() {
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

bool DynamicHipoLoader::tryLoad() {
  const auto paths = getSearchPaths();

  if (paths.empty()) {
    last_error_ = "HiPO not available. Install with: pip install highspy[hipo]";
    return false;
  }

  for (const auto& path : paths) {
    if (loadLibrary(path)) {
      if (resolveFunctions()) {
        // Check ABI compatibility
        int loaded_abi_version = fn_get_abi_version_();
        if (loaded_abi_version != kHipoAbiVersion) {
          last_error_ = "HiPO ABI version mismatch: expected " +
                        std::to_string(kHipoAbiVersion) + ", got " +
                        std::to_string(loaded_abi_version) +
                        ". Please reinstall: pip install --force-reinstall highspy[hipo]";
          unloadLibrary();
          continue;
        }

        // Get version string
        const char* ver = fn_get_version_();
        version_ = ver ? ver : "";

        return true;
      }
      unloadLibrary();
    }
  }

  if (last_error_.empty()) {
    last_error_ = "HiPO not available. Install with: pip install highspy[hipo]";
  }
  return false;
}

HighsStatus DynamicHipoLoader::solveLp(HighsLpSolverObject& solver_object) {
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
