# CMakeLists.txt for highspy-hipo wheel
# Builds the HiPO shared library using the main HiGHS CMake infrastructure

cmake_minimum_required(VERSION 3.15...3.27)

# Use parent directory as source - forward to main CMakeLists with HIPO_PYTHON_BUILD enabled
set(HIGHS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Use VCPKG_ROOT environment variable for dependencies
if(DEFINED ENV{VCPKG_ROOT})
    set(VCPKG_INSTALLED "$ENV{VCPKG_ROOT}/installed/x64-windows-static")
    set(ZLIB_LIBRARY "${VCPKG_INSTALLED}/lib/zlib.lib" CACHE STRING "Root directory of zlib")
    set(ZLIB_INCLUDE_DIR "${VCPKG_INSTALLED}/include" CACHE STRING "Include directory for zlib")
    set(BLAS_ROOT "${VCPKG_INSTALLED}/share/openblas" CACHE STRING "Root directory of BLAS/LAPACK")
    set(BLAS_LIBRARIES "${VCPKG_INSTALLED}/lib/openblas.lib" CACHE STRING "BLAS/LAPACK libraries to link against")
endif()

set(HIPO_PYTHON_BUILD ON CACHE BOOL "" FORCE)
add_subdirectory(${HIGHS_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/highs)
