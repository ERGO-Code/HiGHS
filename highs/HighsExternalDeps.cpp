/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalDeps.cpp
 * @brief Manages access (dynamic or static) to optional external dependencies
 */

#include "HighsExternalDeps.h"
#include "HConfig.h"

// Platform-specific includes for dynamic loading
#if defined(_WIN32) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define PATH_SEPARATOR "\\"
#else
#include <dlfcn.h>
#define PATH_SEPARATOR "/"
#endif

std::string getLibraryFilename() {
#if defined(_WIN32) || defined(_WIN64)
  return "highs_extras.dll";
#elif defined(__APPLE__)
  return "libhighs_extras.dylib";
#else
  return "libhighs_extras.so";
#endif
}

template <typename FuncType>
bool resolveSymbol(void* handle, FuncType& target, const char* name) {
#if defined(_WIN32) || defined(_WIN64)
  target = reinterpret_cast<FuncType>(
      GetProcAddress(static_cast<HMODULE>(handle), name));
#else
  target = reinterpret_cast<FuncType>(dlsym(handle, name));
#endif

  return target != nullptr;
}

// c++11 does not support inline static definition
HighsExternalDeps::amd HighsExternalDeps::amd_;
HighsExternalDeps::blas HighsExternalDeps::blas_;
HighsExternalDeps::metis HighsExternalDeps::metis_;
HighsExternalDeps::rcm HighsExternalDeps::rcm_;

HighsExternalDeps& HighsExternalDeps::instance() {
  static HighsExternalDeps _instance;
  return _instance;
}

#ifdef HIPO
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
void HighsExternalDeps::clear() {
  amd_ = amd{};
  blas_ = blas{};
  metis_ = metis{};
  rcm_ = rcm{};
}

void HighsExternalDeps::unload() {
  HighsExternalDeps& inst = instance();

  if (inst.lib_handle_) {
#if defined(_WIN32) || defined(_WIN64)
    FreeLibrary(static_cast<HMODULE>(inst.lib_handle_));
#else
    dlclose(inst.lib_handle_);
#endif
    inst.lib_handle_ = nullptr;
  }
  inst.clear();
  inst.available_ = false;
}

#define STRINGFY(s) STRINGFY0(s)
#define STRINGFY0(s) #s

bool HighsExternalDeps::tryLoad(const std::string& path) {
  HighsExternalDeps& inst = instance();

  // Allow multiple attempts.
  if (inst.available_)
      return true;

  // printf("Attempting to load HiGHS Extras from: %s\n", path.c_str());
  inst.initialized_ = true;
  inst.available_ = false;

  // Load library
  const std::string full_path = path + PATH_SEPARATOR + getLibraryFilename();

#if defined(_WIN32) || defined(_WIN64)
  inst.lib_handle_ = static_cast<void*>(LoadLibraryA(full_path.c_str()));
  if (!inst.lib_handle_) {
    DWORD error = GetLastError();
    char* msg = nullptr;
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
                   nullptr, error, 0, reinterpret_cast<LPSTR>(&msg), 0,
                   nullptr);
    inst.last_error_ =
        "Failed to load " + full_path + ": " + (msg ? msg : "Unknown error");
    if (msg) LocalFree(msg);
    return false;
  }
#else
  inst.lib_handle_ = dlopen(full_path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!inst.lib_handle_) {
    const char* err = dlerror();
    inst.last_error_ =
        "Failed to load " + full_path + ": " + (err ? err : "Unknown error");
    return false;
  }
#endif

  // Resolve all function pointers
  void* h = inst.lib_handle_;
  bool ok = true;

  // AMD
  ok &= resolveSymbol(h, amd_.defaults_,"highs_extras_amd_defaults");
  ok &= resolveSymbol(h, amd_.order_, "highs_extras_amd_order");

  // BLAS
  ok &= resolveSymbol(h, blas_.daxpy_, "highs_extras_daxpy");
  ok &= resolveSymbol(h, blas_.dcopy_, "highs_extras_dcopy");
  ok &= resolveSymbol(h, blas_.dscal_, "highs_extras_dscal");
  ok &= resolveSymbol(h, blas_.dswap_, "highs_extras_dswap");
  ok &= resolveSymbol(h, blas_.dgemv_, "highs_extras_dgemv");
  ok &= resolveSymbol(h, blas_.dtpsv_, "highs_extras_dtpsv");
  ok &= resolveSymbol(h, blas_.dtrsv_, "highs_extras_dtrsv");
  ok &= resolveSymbol(h, blas_.dger_, "highs_extras_dger");
  ok &= resolveSymbol(h, blas_.dgemm_, "highs_extras_dgemm");
  ok &= resolveSymbol(h, blas_.dsyrk_, "highs_extras_dsyrk");
  ok &= resolveSymbol(h, blas_.dtrsm_, "highs_extras_dtrsm");
  ok &= resolveSymbol(h, blas_.set_num_threads_, "highs_extras_openblas_set_num_threads");
  ok &= resolveSymbol(h, blas_.library_, "highs_extras_blas_library");

  // METIS
  ok &= resolveSymbol(h, metis_.set_default_options_,
                     "highs_extras_metis_set_default_options");
  ok &= resolveSymbol(h, metis_.nodend_, "highs_extras_metis_nodend");

  //// RCM
  //ok &= resolveSymbol(h, rcm.genrcm_, "highs_extras_genrcm");

  // Check ABI compatibility
  highs_extras_api::core::get_version_t get_version = nullptr;
  ok &= resolveSymbol(h, get_version, "highs_extras_get_version");

  ok &= resolveSymbol(h, instance().get_copyright_, "highs_extras_get_copyright");

  if (!ok) {
    inst.last_error_ = "Failed to resolve required external functions";
    inst.unload();
    return false;
  }

  std::string highs_version = STRINGFY(HIGHS_VERSION_MAJOR) "." STRINGFY(
      HIGHS_VERSION_MINOR) "." STRINGFY(HIGHS_VERSION_PATCH);
  std::string extras_version = get_version();
  if (extras_version != highs_version) {
    inst.last_error_ =
        "HiGHS Extras ABI version mismatch: expected " + highs_version +
        ", got " + extras_version +
        ". Please reinstall: pip install --force-reinstall highspy[extras]";
    inst.unload();
    return false;
  }

  inst.available_ = true;
  return true;
}
#endif
#endif
