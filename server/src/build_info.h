// build_info.h
// Reports GPU capability at runtime. Determined by the HIGHS_SERVER_CUPDLP_GPU
// macro injected at compile time by CMake.
#pragma once

namespace highs_server {

// 0 = CPU build (CUPDLP_GPU=OFF), 1 = GPU build (CUPDLP_GPU=ON)
constexpr bool kGpuBuilt =
#if defined(HIGHS_SERVER_CUPDLP_GPU) && HIGHS_SERVER_CUPDLP_GPU
    true;
#else
    false;
#endif

inline const char* GpuBuildString() {
  return kGpuBuilt ? "GPU-enabled (CUPDLP_GPU=ON)"
                   : "CPU-only (CUPDLP_GPU=OFF)";
}

}  // namespace highs_server
