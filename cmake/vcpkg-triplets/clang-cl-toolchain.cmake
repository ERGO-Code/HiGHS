# Makes the vcpkg build of each port use clang-cl as its compiler. clang-cl is
# searched for on PATH and in the default LLVM install location.
find_program(CLANG_CL_EXECUTABLE clang-cl
    PATHS "$ENV{ProgramFiles}/LLVM/bin" REQUIRED)
set(CMAKE_C_COMPILER "${CLANG_CL_EXECUTABLE}")
set(CMAKE_CXX_COMPILER "${CLANG_CL_EXECUTABLE}")
# Ninja would pass rc.exe flags to clang-cl when scanning .rc files
set(CMAKE_NINJA_CMCLDEPS_RC OFF)

# Keep vcpkg's usual Windows compiler settings
list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES _VCPKG_ROOT_DIR)
include("${_VCPKG_ROOT_DIR}/scripts/toolchains/windows.cmake")
