include_guard(GLOBAL)

# Fetch OpenBLAS
message(STATUS "Running FindHipoDeps, BUILD_OPENBLAS=${BUILD_OPENBLAS}")

if(APPLE)
    set(BUILD_OPENBLAS OFF CACHE BOOL "OpenBLAS not required on macOS" FORCE)
    set(BUILD_OPENBLAS OFF)
endif()

set(BLAS_ROOT "" CACHE STRING "Root directory of BLAS or OpenBLAS")

# set appropriate BLAS linking information
function(highs_configure_blas_target)
    set(HIGHS_BLAS_TARGET "" PARENT_SCOPE)
    set(HIGHS_BLAS_INCLUDE_DIRS "" PARENT_SCOPE)

    if(OpenBLAS_FOUND AND TARGET OpenBLAS::OpenBLAS)
        set(HIGHS_BLAS_TARGET OpenBLAS::OpenBLAS PARENT_SCOPE)
    elseif(OPENBLAS_LIB)
        set(HIGHS_BLAS_TARGET "${OPENBLAS_LIB}" PARENT_SCOPE)
        set(HIGHS_BLAS_INCLUDE_DIRS "${OPENBLAS_INCLUDE_DIR}" PARENT_SCOPE)
    elseif(TARGET BLAS::BLAS)
        set(HIGHS_BLAS_TARGET BLAS::BLAS PARENT_SCOPE)
    elseif(BLAS_LIBRARIES)
        set(HIGHS_BLAS_TARGET "${BLAS_LIBRARIES}" PARENT_SCOPE)
    endif()
endfunction()

# set appropriate BLAS dependency licensing metadata
function(highs_configure_blas_metadata)
    set(HIGHS_BLAS_COMPILE_DEFINITION "" PARENT_SCOPE)

    if(OpenBLAS_FOUND OR OPENBLAS_LIB OR BLA_VENDOR MATCHES "OpenBLAS")
        set(HIGHS_BLAS_VENDOR OpenBLAS PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${OpenBLAS_VERSION}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE BSD-3-Clause PARENT_SCOPE)
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_OPENBLAS PARENT_SCOPE)
    elseif(BLA_VENDOR MATCHES "Intel|MKL")
        set(HIGHS_BLAS_VENDOR "Intel oneMKL" PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${MKL_VERSION}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE "Intel Simplified Software License" PARENT_SCOPE)
    elseif(APPLE AND BLA_VENDOR STREQUAL "Apple")
        set(HIGHS_BLAS_VENDOR "Apple Accelerate" PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${CMAKE_OSX_DEPLOYMENT_TARGET}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE "Apple SDK license" PARENT_SCOPE)
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_APPLE_BLAS PARENT_SCOPE)
    else()
        set(HIGHS_BLAS_VENDOR "${BLA_VENDOR}" PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION unknown PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE unknown PARENT_SCOPE)
    endif()
endfunction()

# set cmake variables for BLAS for parent scope
macro(highs_export_blas_state)
    foreach(_highs_blas_var
            HIGHS_BLAS_CONFIGURED
            HIGHS_BLAS_TARGET
            HIGHS_BLAS_INCLUDE_DIRS
            HIGHS_BLAS_COMPILE_DEFINITION
            HIGHS_BLAS_VENDOR
            HIGHS_BLAS_VERSION
            HIGHS_BLAS_LICENSE)
        if(DEFINED ${_highs_blas_var})
            set(${_highs_blas_var} "${${_highs_blas_var}}" PARENT_SCOPE)
        else()
            unset(${_highs_blas_var} PARENT_SCOPE)
        endif()
    endforeach()
endmacro()

# discover BLAS (build if necessary) and set cmake variables for linking and metadata
function(highs_configure_blas)
    if(HIGHS_BLAS_CONFIGURED)
        return()
    endif()

    if(NOT BLAS_ROOT STREQUAL "")
        message(STATUS "BLAS_ROOT is ${BLAS_ROOT}")
    endif()

    set(USE_CMAKE_FIND_BLAS ON)

    if(BUILD_OPENBLAS)
        set(_highs_openblas_version "0.3.30")

        include(FetchContent)
        set(FETCHCONTENT_QUIET OFF)
        set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
        set(CMAKE_POSITION_INDEPENDENT_CODE ON)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(CMAKE_Fortran_COMPILER OFF)

        # NOTE: FetchContent_MakeAvailable() populates OpenBLAS via add_subdirectory(),
        # not a separate `cmake` invocation, so FetchContent_Declare's CMAKE_ARGS has no
        # effect here (that option only applies to ExternalProject-style population) -
        # it was silently a no-op. Every OpenBLAS build option below is instead a real
        # CACHE variable set in this (shared) scope, so OpenBLAS's own option()/
        # if(DEFINED ...) checks in its CMakeLists.txt/cmake/system.cmake actually see it.

        # Exclude components not used by HiGHS
        set(ONLY_CBLAS ON CACHE BOOL "" FORCE)
        set(NO_LAPACK ON CACHE BOOL "" FORCE)
        set(NO_LAPACKE ON CACHE BOOL "" FORCE)
        set(NO_COMPLEX ON CACHE BOOL "" FORCE)
        set(NO_COMPLEX16 ON CACHE BOOL "" FORCE)
        set(NO_DOUBLE_COMPLEX ON CACHE BOOL "" FORCE)
        set(NO_SINGLE ON CACHE BOOL "" FORCE)

        if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|armv8|arm")
            if(CMAKE_SIZEOF_VOID_P EQUAL 4)
                message(FATAL_ERROR "The HiGHS build with OpenBLAS does not yet support 32-bit ARM architectures. \
                        You could try to compile OpenBLAS separately on your machine, see https://github.com/OpenMathLib/OpenBLAS. \
                        Then link with HiGHS by passing the path to the OpenBLAS installation via BLAS_ROOT. \
                        Please don't hesitate to get in touch with us with details about your related issues.")
            endif()
        endif()

        if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
            # OpenBLAS' DYNAMIC_ARCH kernel objects (and some of its generic,
            # non-per-arch fallback kernels, e.g. kernel/x86_64/dscal.c and
            # common_arm64.h) rely on GNU inline assembly (__asm__ __volatile__,
            # __attribute__) that plain MSVC (cl.exe) cannot parse at all - unlike
            # clang-cl, which accepts it under its GCC-compatibility mode. Rather
            # than silently building a single-target library (which is what
            # happened before this file started actually forwarding these options
            # to OpenBLAS - see the FetchContent_Declare CMAKE_ARGS note above),
            # explicitly pin the portable GENERIC target with DYNAMIC_ARCH off, so
            # the result stays correct (if unoptimized) on any x86_64/ARM64 host
            # regardless of what CPU it happened to be built on.
            message(STATUS "MSVC detected: building OpenBLAS as a single portable GENERIC target (no DYNAMIC_ARCH) since MSVC cannot compile OpenBLAS's GNU-inline-asm kernel objects.")
            set(DYNAMIC_ARCH OFF CACHE BOOL "" FORCE)
            set(TARGET GENERIC CACHE STRING "" FORCE)
        else()
            message(STATUS "Enabling DYNAMIC_ARCH for runtime CPU detection.")
            set(DYNAMIC_ARCH ON CACHE BOOL "" FORCE)
        endif()

        # CMAKE_SIZEOF_VOID_P is 4 for 32-bit and 8 for 64-bit
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            message(STATUS "32-bit target detected. Applying 32-bit configuration flags for OpenBLAS.")

            if(UNIX AND NOT APPLE)
                # OpenBLAS' build-time cpuid probe (getarch) reports the host's real
                # (64-bit-capable) microarchitecture (e.g. ZEN), which has no
                # 32-bit-compatible kernel and leaves PREFETCH-related macros undefined,
                # breaking assembly of the 32-bit kernels (e.g. gemm_kernel_2x4_sse3.S /
                # gemm_kernel_4x4_sse3.S). Pin a generic 32-bit-safe default target, and
                # set BINARY=32 so OpenBLAS's own cmake/system.cmake also passes
                # -DNO_AVX{,2,512} to the getarch probe itself (its CMAKE_SIZEOF_VOID_P
                # fallback check for this is broken - see system.cmake's GETARCH_FLAGS
                # setup). DYNAMIC_ARCH still dispatches to better kernels at runtime.

                # message(STATUS "Pinning OpenBLAS TARGET=PRESCOTT for 32-bit Linux build to avoid cpuid misdetection.")
                # set(TARGET PRESCOTT CACHE STRING "" FORCE)

                set(BINARY 32 CACHE STRING "" FORCE)
            endif()

            set(INTERFACE64 0 CACHE STRING "" FORCE)
        endif()

        # OpenBLAS' DYNAMIC_ARCH CMake build repeatedly mis-scopes AVX512-gated macros
        # (e.g. HAVE_AVX512F leaking into generic kernel objects, like kernel/arm/sum.c's
        # ssum_k/dsum_k, that are compiled without a matching -mavx512f flag), which both
        # clang and gcc reject as a hard error. This reproduces across unrelated hosts/
        # toolchains (Skylake register spills were only one manifestation), so always
        # disable AVX512 for Linux builds of OpenBLAS rather than special-casing Skylake.
        if(UNIX AND NOT APPLE)
            message(STATUS "Disabling AVX512 in OpenBLAS (unreliable with DYNAMIC_ARCH's CMake build).")
            set(NO_AVX512 ON CACHE BOOL "Build OpenBLAS without AVX512" FORCE)
        endif()

        set(OPENBLAS_BUILD_TYPE "Release" CACHE STRING "Build type for OpenBLAS" FORCE)

        if(NOT DEBUG_MEMORY STREQUAL "Off")
            # OpenBLAS's own CMake build (getarch's compiler probe in particular,
            # see cmake/prebuild.cmake) compiles a throwaway test program with
            # whatever CMAKE_<LANG>_FLAGS_<CONFIG> is in scope, but does not
            # forward the matching sanitizer runtime to the link step it uses to
            # link that program, so the probe fails with undefined references to
            # e.g. __tsan_func_entry. OpenBLAS is vendored/third-party code we
            # don't sanitize anyway, so strip the sanitizer flags HiGHS added to
            # these variables just for fetching/building it here; this is a
            # function-local shadow (no CACHE), so it reverts automatically once
            # highs_configure_blas() returns and has no effect on HiGHS's own
            # targets.
            message(STATUS "DEBUG_MEMORY=${DEBUG_MEMORY}: building OpenBLAS without sanitizer instrumentation")
            foreach(_highs_blas_flags_var
                    CMAKE_C_FLAGS_DEBUG
                    CMAKE_C_FLAGS_RELWITHDEBINFO
                    CMAKE_CXX_FLAGS_DEBUG
                    CMAKE_CXX_FLAGS_RELWITHDEBINFO
                    CMAKE_EXE_LINKER_FLAGS_DEBUG
                    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
                    CMAKE_SHARED_LINKER_FLAGS_DEBUG
                    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO)
                string(REGEX REPLACE "-f(no-)?sanitize[-=][A-Za-z0-9,-]+" "" ${_highs_blas_flags_var} "${${_highs_blas_flags_var}}")
                set(${_highs_blas_flags_var} "${${_highs_blas_flags_var}}")
            endforeach()
        endif()

        if(DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
            set(_highs_blas_ipo_backup "${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
        endif()
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF CACHE BOOL "" FORCE)

        message(CHECK_START "Fetching OpenBLAS")
        list(APPEND CMAKE_MESSAGE_INDENT "  ")
        FetchContent_Declare(
                openblas
                GIT_REPOSITORY "https://github.com/OpenMathLib/OpenBLAS.git"
                GIT_TAG        "v${_highs_openblas_version}"
                GIT_SHALLOW TRUE
                UPDATE_COMMAND git reset --hard
        )
        FetchContent_MakeAvailable(openblas)
        get_property(all_targets DIRECTORY ${openblas_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
        message(STATUS "OpenBLAS targets: ${all_targets}")

        foreach(_lapacke_target LAPACKE genlapacke)
            if(TARGET ${_lapacke_target})
                set_target_properties(${_lapacke_target} PROPERTIES EXCLUDE_FROM_ALL TRUE)
            endif()
        endforeach()

        set_property(DIRECTORY ${openblas_SOURCE_DIR}
                PROPERTY CTEST_EXCLUDE_FROM_MAIN TRUE)

        if(ALL_TESTS)
            set(BUILD_TESTING ON CACHE BOOL "" FORCE)
        endif()

        if(DEFINED _highs_blas_ipo_backup)
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION "${_highs_blas_ipo_backup}" CACHE BOOL "" FORCE)
        endif()

        list(POP_BACK CMAKE_MESSAGE_INDENT)
        message(CHECK_PASS "fetched")

        if(TARGET openblas)
            get_target_property(_openblas_aliased openblas ALIASED_TARGET)
            if(_openblas_aliased)
                set(_openblas_target ${_openblas_aliased})
                message(STATUS "OpenBLAS is an alias for: ${_openblas_target}")
            else()
                set(_openblas_target openblas)
            endif()
        elseif(TARGET openblas_static)
            set(_openblas_target openblas_static)
        elseif(TARGET openblas_shared)
            set(_openblas_target openblas_shared)
        else()
            message(FATAL_ERROR "OpenBLAS target not found")
        endif()

        if(NOT DEFINED OpenBLAS_VERSION OR OpenBLAS_VERSION STREQUAL "")
            message(STATUS "OpenBLAS_VERSION not reported; using ${_highs_openblas_version}")
            set(OpenBLAS_VERSION "${_highs_openblas_version}")
        endif()

        set(OpenBLAS_FOUND TRUE)
        set(HIGHS_BLAS_TARGET ${_openblas_target})
        set(HIGHS_BLAS_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/openblas-src/include")
        highs_configure_blas_metadata()

        set(HIGHS_BLAS_CONFIGURED TRUE)
        highs_export_blas_state()
        return()
    endif()

    if(NOT USE_CMAKE_FIND_BLAS)
        if(WIN32)
            if(NOT (BLAS_ROOT STREQUAL ""))
                message(STATUS "Looking for blas in ${BLAS_ROOT}")
                set(OpenBLAS_ROOT ${BLAS_ROOT})
                message(STATUS "OpenBLAS_ROOT is ${OpenBLAS_ROOT}")
                find_package(OpenBLAS CONFIG NO_DEFAULT_PATH)

                if(OpenBLAS_FOUND)
                    message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
                else()
                    message(STATUS "OpenBLAS not found in ${BLAS_ROOT}")
                endif()
            endif()

            if((BLAS_ROOT STREQUAL "") OR (NOT OpenBLAS_FOUND))
                message(STATUS "Looking for blas")
                find_package(OpenBLAS REQUIRED)

                if(OpenBLAS_FOUND)
                    if(TARGET OpenBLAS::OpenBLAS)
                        message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
                    elseif(OPENBLAS_LIB)
                        message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
                    endif()
                else()
                    message(FATAL_ERROR "No BLAS library found")
                endif()
            endif()
        elseif(NOT APPLE)
            # If a BLAS install was specified, try to find it first
            if(NOT (BLAS_ROOT STREQUAL ""))
                message(STATUS "Looking for blas in ${BLAS_ROOT}")

                find_library(OPENBLAS_LIB
                        NAMES openblas
                        HINTS "${BLAS_ROOT}/lib"
                        NO_DEFAULT_PATH)

                if(OPENBLAS_LIB)
                    message(STATUS "Found OpenBLAS library at ${OPENBLAS_LIB}")
                else()
                    find_library(BLAS_LIB
                            NAMES blas
                            HINTS "${BLAS_ROOT}/lib"
                            NO_DEFAULT_PATH)

                    if(BLAS_LIB)
                        message(STATUS "Found BLAS library at ${BLAS_LIB}")
                    else()
                        message(STATUS "Did not find blas library at ${BLAS_ROOT}")
                        message(STATUS "Attempting default locations search")
                    endif()
                endif()
            endif()

            if((BLAS_ROOT STREQUAL "") OR (NOT OPENBLAS_LIB AND NOT BLAS_LIB))
                find_library(OPENBLAS_LIB
                        NAMES openblas
                        HINTS "${BLAS_ROOT}/lib")

                if(OPENBLAS_LIB)
                    message(STATUS "Found OpenBLAS library at ${OPENBLAS_LIB}")
                else()
                    find_library(BLAS_LIB
                            NAMES blas
                            HINTS "${BLAS_ROOT}/lib")

                    if(BLAS_LIB)
                        message(STATUS "Found BLAS library at ${BLAS_LIB}")
                    else()
                        message(FATAL_ERROR "No BLAS library found")
                    endif()
                endif()
            endif()
        endif()
    else()
        if(NOT BLAS_LIBRARIES AND NOT BLA_VENDOR)
            find_package(OpenBLAS CONFIG)
            if(OpenBLAS_FOUND)
                message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")

                # Ubuntu's OpenBLASConfig.cmake sets OpenBLAS_LIBRARIES, not OPENBLAS_LIB
                if(NOT OPENBLAS_LIB AND OpenBLAS_LIBRARIES)
                    set(OPENBLAS_LIB ${OpenBLAS_LIBRARIES})
                    set(OPENBLAS_INCLUDE_DIR ${OpenBLAS_INCLUDE_DIRS})
                endif()
            endif()
        endif()

        if(NOT OpenBLAS_FOUND)
            if(NOT BLA_VENDOR AND NOT BLAS_LIBRARIES)
                if(APPLE)
                    set(BLA_VENDOR Apple)
                elseif(LINUX)
                    set(BLA_VENDOR OpenBLAS)
                elseif(WIN32)
                    set(BLA_VENDOR OpenBLAS)
                endif()

                find_package(BLAS QUIET)
                if(BLAS_FOUND)
                    message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                    if(BLAS_INCLUDE_DIRS)
                        message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                    endif()
                else()
                    unset(BLA_VENDOR)
                endif()
            else()
                if(BLA_VENDOR)
                    message(STATUS "Specified BLA_VENDOR: ${BLA_VENDOR}")
                endif()
                if(BLAS_LIBRARIES)
                    message(STATUS "Specified BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
                endif()

                if(BLA_VENDOR MATCHES "Intel|MKL")
                    find_package(MKL CONFIG REQUIRED)
                endif()
            endif()

            # try libblas on linux
            if(LINUX AND NOT BLAS_FOUND)
                find_package(BLAS QUIET)
                if(BLAS_FOUND)
                    message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                    if(BLAS_INCLUDE_DIRS)
                        message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                    endif()
                endif()
            endif()

            if(NOT BLAS_FOUND)
                find_package(BLAS REQUIRED)

                if(BLAS_FOUND)
                    message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                    if(BLAS_INCLUDE_DIRS)
                        message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                    endif()
                else()
                    message(FATAL_ERROR "No BLAS library found!")
                endif()
            endif()
        endif()
    endif()

    set(HIGHS_BLAS_CONFIGURED TRUE)
    highs_configure_blas_target()
    highs_configure_blas_metadata()
    highs_export_blas_state()
endfunction()

# configure the given target for linking to BLAS
function(highs_link_blas target_name)
    if(NOT HIGHS_BLAS_CONFIGURED)
        message(FATAL_ERROR "Ensure highs_configure_blas() called before highs_link_blas(${target_name}).")
    endif()

    target_compile_definitions(${target_name} PRIVATE ${HIGHS_BLAS_COMPILE_DEFINITION})
    target_compile_definitions(${target_name} PRIVATE HIGHS_BLAS_VENDOR="${HIGHS_BLAS_VENDOR}")
    target_compile_definitions(${target_name} PRIVATE HIGHS_BLAS_VERSION="${HIGHS_BLAS_VERSION}")
    target_compile_definitions(${target_name} PRIVATE HIGHS_BLAS_LICENSE="${HIGHS_BLAS_LICENSE}")

    if(NOT HIGHS_BLAS_TARGET)
        message(FATAL_ERROR "BLAS not found for ${target_name}.")
    endif()

    if(BUILD_OPENBLAS AND NOT TARGET ${HIGHS_BLAS_TARGET})
        message(FATAL_ERROR "OpenBLAS target not found for ${target_name}.")
    endif()

    target_link_libraries(${target_name} PUBLIC ${HIGHS_BLAS_TARGET})

    if(HIGHS_BLAS_INCLUDE_DIRS)
        target_include_directories(${target_name} PUBLIC
            $<BUILD_INTERFACE:${HIGHS_BLAS_INCLUDE_DIRS}>)
    endif()
endfunction()