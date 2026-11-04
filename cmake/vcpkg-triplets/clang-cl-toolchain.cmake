# Makes the vcpkg build of each port use clang-cl as its compiler.
include("${CMAKE_CURRENT_LIST_DIR}/find-clang-cl.cmake")
highs_find_clang_cl(CLANG_CL_EXECUTABLE)
if(NOT CLANG_CL_EXECUTABLE)
    message(FATAL_ERROR "clang-cl not found")
endif()
set(CMAKE_C_COMPILER "${CLANG_CL_EXECUTABLE}")
set(CMAKE_CXX_COMPILER "${CLANG_CL_EXECUTABLE}")
# Ninja would pass rc.exe flags to clang-cl when scanning .rc files
set(CMAKE_NINJA_CMCLDEPS_RC OFF)

# Keep vcpkg's usual Windows compiler settings
list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES _VCPKG_ROOT_DIR)
include("${_VCPKG_ROOT_DIR}/scripts/toolchains/windows.cmake")
