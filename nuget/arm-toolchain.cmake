# Toolchain file for cross-compiling for ARM

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

# Specify the cross-compiler paths
find_program(CMAKE_C_COMPILER NAMES aarch64-linux-gnu-gcc)
find_program(CMAKE_CXX_COMPILER NAMES aarch64-linux-gnu-g++)

# Compiler flags
set(CMAKE_C_FLAGS_INIT "-O3" CACHE STRING "")
set(CMAKE_CXX_FLAGS_INIT "-O3" CACHE STRING "")

# Set this to true so CMake knows it's cross-compiling
set(CMAKE_CROSSCOMPILING TRUE)