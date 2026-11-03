# x64-windows-static, but with every port built by clang-cl, so that OpenBLAS
# gets its optimised kernels (MSVC can only build the generic ones)
set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE static)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_CHAINLOAD_TOOLCHAIN_FILE "${CMAKE_CURRENT_LIST_DIR}/clang-cl-toolchain.cmake")
# Still use the MSVC environment (linker, Windows SDK); this is off by default
# with a chainloaded toolchain
set(VCPKG_LOAD_VCVARS_ENV ON)

if(PORT STREQUAL "openblas")
    # Select the kernels for the CPU at runtime
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS -DDYNAMIC_ARCH=ON)
    # Recent clang reports incompatible pointer types in some kernels as errors
    set(VCPKG_C_FLAGS "-Wno-error=incompatible-pointer-types -w")
    set(VCPKG_CXX_FLAGS "")
endif()
