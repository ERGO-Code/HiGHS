// build_info.h
// 运行时上报编译期 GPU 能力。由 CMake 注入的 HIGHS_SERVER_CUPDLP_GPU 宏决定。
#pragma once

namespace highs_server {

// 0 = CPU 构建 (CUPDLP_GPU=OFF), 1 = GPU 构建 (CUPDLP_GPU=ON)
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
