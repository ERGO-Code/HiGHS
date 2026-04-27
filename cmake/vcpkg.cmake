# cmake/vcpkg.cmake
#
# Shared vcpkg integration for HiGHS-related projects.
#
# Usage:
#   cmake_minimum_required(VERSION 3.15)
#   list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
#   include(vcpkg)        # MUST be before the first project() call
#   project(my_project LANGUAGES C CXX)
#
# Behavior:
#   * If $ENV{VCPKG_ROOT} is set and CMAKE_TOOLCHAIN_FILE is not, point
#     CMAKE_TOOLCHAIN_FILE at vcpkg's toolchain so manifest mode kicks in.
#   * Auto-pick a sensible VCPKG_TARGET_TRIPLET unless the caller (or env
#     var VCPKG_DEFAULT_TRIPLET) already set one. Windows defaults to the
#     *-windows-static triplet so wheels stay self-contained.
#
# Notes on cross-compilation:
#   * Visual Studio generators expose the target arch via
#     CMAKE_GENERATOR_PLATFORM (-A Win32/x64/ARM64) — handled below.
#   * For Ninja/Makefile generators we cannot tell host vs. target before
#     project(); we fall back to the host arch (queried via
#     cmake_host_system_information OS_PLATFORM, which works pre-project()
#     unlike CMAKE_HOST_SYSTEM_PROCESSOR which is empty until project()
#     runs on CMake < 3.25). When cross-compiling with such a generator,
#     set VCPKG_TARGET_TRIPLET (or the VCPKG_DEFAULT_TRIPLET env var)
#     explicitly, e.g. -DVCPKG_TARGET_TRIPLET=arm64-windows-static.
#
# Skip when consumed as a vcpkg port:
#   When HiGHS is built as a port (e.g. installed via `vcpkg install highs`),
#   the outer vcpkg invocation already owns the toolchain, triplet, and
#   manifest settings. Detecting that situation and bailing out early
#   prevents this module from forcing a *-windows-static triplet, hijacking
#   VCPKG_MANIFEST_DIR to extern/vcpkg.json, or otherwise interfering with
#   the port build. Indicators that we're inside a vcpkg-driven build:
#     * CMAKE_TOOLCHAIN_FILE is already set (any toolchain in effect)
#     * Z_VCPKG_ROOT_DIR exported by vcpkg's portfile

if(DEFINED CMAKE_TOOLCHAIN_FILE
   OR DEFINED ENV{Z_VCPKG_ROOT_DIR})
  return()
endif()

if (DEFINED ENV{VCPKG_INSTALLATION_ROOT} AND NOT DEFINED ENV{VCPKG_ROOT})
  set(ENV{VCPKG_ROOT} $ENV{VCPKG_INSTALLATION_ROOT})
endif()

# Treat an empty VCPKG_ROOT the same as undefined (allows overrides to
# disable vcpkg by setting VCPKG_ROOT="").
if(DEFINED ENV{VCPKG_ROOT} AND NOT "$ENV{VCPKG_ROOT}" STREQUAL "")
  message(STATUS "Using VCPKG_ROOT: $ENV{VCPKG_ROOT}")
  file(TO_CMAKE_PATH "$ENV{VCPKG_ROOT}" VCPKG_ROOT_CMAKE)

  set(CMAKE_TOOLCHAIN_FILE
    "${VCPKG_ROOT_CMAKE}/scripts/buildsystems/vcpkg.cmake"
    CACHE FILEPATH "vcpkg toolchain file")

  # The manifest lives in extern/ alongside the code that consumes the
  # dependencies. vcpkg only auto-discovers vcpkg.json in CMAKE_SOURCE_DIR,
  # so point VCPKG_MANIFEST_DIR at extern/ explicitly. This makes manifest
  # mode work whether the top-level project is highspy-extras/, extern/
  # itself, or any future consumer.
  if(NOT DEFINED VCPKG_MANIFEST_DIR
     AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/../extern/vcpkg.json")
    get_filename_component(VCPKG_MANIFEST_DIR
      "${CMAKE_CURRENT_LIST_DIR}/../extern" ABSOLUTE)
    set(VCPKG_MANIFEST_DIR "${VCPKG_MANIFEST_DIR}"
      CACHE PATH "vcpkg manifest directory")
  endif()

  if(NOT DEFINED VCPKG_TARGET_TRIPLET AND NOT DEFINED ENV{VCPKG_DEFAULT_TRIPLET})
    # CMAKE_HOST_SYSTEM_PROCESSOR isn't reliably populated before project()
    # on CMake < 3.25, so query the host arch explicitly.
    if(NOT _vcpkg_host_arch)
      cmake_host_system_information(RESULT _vcpkg_host_arch
        QUERY OS_PLATFORM)
    endif()

    if(CMAKE_GENERATOR_PLATFORM STREQUAL "Win32")
      set(_vcpkg_arch x86)
    elseif(CMAKE_GENERATOR_PLATFORM STREQUAL "ARM64")
      set(_vcpkg_arch arm64)
    elseif(CMAKE_GENERATOR_PLATFORM STREQUAL "x64")
      set(_vcpkg_arch x64)
    elseif(_vcpkg_host_arch MATCHES "^(x86_64|AMD64|x64)$")
      set(_vcpkg_arch x64)
    elseif(_vcpkg_host_arch MATCHES "^(aarch64|arm64|ARM64)$")
      set(_vcpkg_arch arm64)
    elseif(_vcpkg_host_arch MATCHES "^(i.86|x86|X86)$")
      set(_vcpkg_arch x86)
    endif()

    message(STATUS "Host arch (OS_PLATFORM): ${_vcpkg_host_arch}")
    message(STATUS "Auto-detected vcpkg architecture: ${_vcpkg_arch}")

    if(WIN32 AND DEFINED _vcpkg_arch)
      set(VCPKG_TARGET_TRIPLET "${_vcpkg_arch}-windows-static"
        CACHE STRING "vcpkg triplet")
      set(VCPKG_HOST_TRIPLET "${_vcpkg_arch}-windows-static"
        CACHE STRING "vcpkg host triplet")

      # Match the static CRT used by the *-windows-static triplet to
      # avoid LNK4098 "defaultlib 'LIBCMT' conflicts with ..." warnings.
      set(CMAKE_MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>"
        CACHE STRING "MSVC runtime library")
        
    elseif(APPLE AND DEFINED _vcpkg_arch)
      set(VCPKG_TARGET_TRIPLET "${_vcpkg_arch}-osx"
        CACHE STRING "vcpkg triplet")
    elseif(UNIX AND DEFINED _vcpkg_arch)
      set(VCPKG_TARGET_TRIPLET "${_vcpkg_arch}-linux"
        CACHE STRING "vcpkg triplet")
    endif()

    unset(_vcpkg_arch)
    unset(_vcpkg_host_arch)
  endif()

  # Fall back to the env-var default if nothing else picked a triplet.
  if(NOT DEFINED VCPKG_TARGET_TRIPLET AND DEFINED ENV{VCPKG_DEFAULT_TRIPLET})
    set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_DEFAULT_TRIPLET}"
      CACHE STRING "vcpkg triplet")
  endif()

  if(DEFINED VCPKG_TARGET_TRIPLET)
    message(STATUS "Using vcpkg triplet: ${VCPKG_TARGET_TRIPLET}")
  else()
    message(WARNING
      "VCPKG_TARGET_TRIPLET not set; vcpkg will use its built-in default. "
      "Set -DVCPKG_TARGET_TRIPLET=... when cross-compiling.")
  endif()
endif()
