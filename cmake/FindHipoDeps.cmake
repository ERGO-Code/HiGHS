# Fetch OpenBLAS
if (BUILD_OPENBLAS)
    include(FetchContent)
    set(FETCHCONTENT_QUIET OFF)
    set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
    # set(BUILD_SHARED_LIBS ON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(BUILD_TESTING OFF)
    set(CMAKE_Fortran_COMPILER OFF)

    # Define the size-minimizing flags as a list
    set(OPENBLAS_MINIMAL_FLAGS
        # Exclude components not used by HiGHS
        -DONLY_CBLAS:BOOL=ON
        -DNO_LAPACK:BOOL=ON
        -DNO_LAPACKE:BOOL=ON
        -DNO_COMPLEX:BOOL=ON
        -DNO_COMPLEX16:BOOL=ON
        -DNO_DOUBLE_COMPLEX:BOOL=ON
        -DNO_SINGLE:BOOL=ON
    )

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|armv8|arm")
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            message(FATAL_ERROR "The HiGHS build with OpenBLAS does not yet support 32-bit ARM architectures. \
            You could try to compile OpenBLAS separately on your machine, see https://github.com/OpenMathLib/OpenBLAS. \
            Then link with HiGHS by passing the path to the OpenBLAS installation via BLAS_ROOT. \
            Please don't hesitate to get in touch with us with details about your related issues.")

            # Unreachable, revisit later. Could not get it to work on the CI, -DNOASM=1 is not being respected and openblas
            # keeps trying to use 64bit registers which are not available. Works fine on 32bit amd and 64bit arm.
            message(STATUS "ARM architecture detected. 32bit.")

            # list(APPEND OPENBLAS_MINIMAL_FLAGS -DARMV7:BOOL=ON)
             # Set environment variable to disable assembly
            # set(ENV{NOASM} "1")
            # set(NOASM 1)

            list(APPEND OPENBLAS_MINIMAL_FLAGS
                -DTARGET=GENERIC
                -DBINARY=32
                -DNOASM=1
                -DDYNAMIC_ARCH:BOOL=OFF
                -DUSE_THREAD:BOOL=OFF
                # Aggressively disable complex operations
                -DNO_CGEMM:BOOL=ON
                -DNO_ZGEMM:BOOL=ON
                -DNO_CTRMM:BOOL=ON
                -DNO_ZTRMM:BOOL=ON
                -DNO_CTRSM:BOOL=ON
                -DNO_ZTRSM:BOOL=ON
                # Disable all Level 3 BLAS (includes TRMM, TRSM, etc.)
                -DNO_LEVEL3:BOOL=ON
                -DCMAKE_C_FLAGS="-march=armv7-a -mfpu=vfpv3-d16 -mfloat-abi=softfp"
                -DCMAKE_ASM_FLAGS="-march=armv7-a -mfpu=vfpv3-d16 -mfloat-abi=softfp"
                -DCMAKE_CXX_FLAGS="-march=armv7-a -mfpu=vfpv3-d16 -mfloat-abi=softfp"
                # -DARM_SOFTFP_ABI=1
                # -DCMAKE_ASM_FLAGS="-mfpu=vfpv3-d16"
                # -DCMAKE_C_FLAGS="-march=armv7-a -mfpu=vfpv3-d16"
                # -DCMAKE_ASM_FLAGS="-march=armv7-a -mfpu=vfpv3-d16"
                # -DCMAKE_CXX_FLAGS="-march=armv7-a -mfpu=vfpv3-d16"
            )
            # list(APPEND OPENBLAS_MINIMAL_FLAGS -DTARGET=GENERIC)
            # list(APPEND OPENBLAS_MINIMAL_FLAGS
            #     -DDYNAMIC_ARCH:BOOL=OFF
            #     -DUSE_THREAD:BOOL=OFF        # Simplify build
            #     -DNO_WARMUP:BOOL=ON          # Skip warmup routine
            #     # -DNO_GETARCH:BOOL=ON
            #     # -DUSE_VFPV3:BOOL=ON
            #     # -DUSE_VFPV3_D32:BOOL=OFF   # crucial: only use d0–d15
            #     # -DNO_TRMM:BOOL=ON
            #     # -DNO_TRSM:BOOL=ON
            #     -DNO_L3:BOOL=ON               # skip complex Level-3 kernels
            #     # -DCMAKE_ASM_FLAGS="-mfpu=vfpv3-d16"
            #     # -DUSE_GENERIC:BOOL=ON
            # )
            # # Explicitly disable assembly

            # set(CMAKE_ASM_COMPILER "")
            # set(NOASM 1)

            # set(SKIP_PARSE_GETARCH TRUE)
        else()
            message(STATUS "ARM architecture detected. Applying -DTARGET=ARMV8.")
            list(APPEND OPENBLAS_MINIMAL_FLAGS -DTARGET=ARMV8)
            # list(APPEND OPENBLAS_MINIMAL_FLAGS -DONLY_BLAS=ON -DNO_LAPACK=ON -DNO_LAPACKE=ON)
        endif()
    # else()
        # list(APPEND OPENBLAS_MINIMAL_FLAGS -DONLY_BLAS=ON -DNO_LAPACK=ON -DNO_LAPACKE=ON)
    endif()

    # CMAKE_SIZEOF_VOID_P is 4 for 32-bit builds, 8 for 64-bit builds.
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        message(STATUS "32-bit target detected. Applying 32-bit configuration flags for OpenBLAS.")

        if (WIN32)
            list(APPEND OPENBLAS_MINIMAL_FLAGS -DCMAKE_GENERATOR_PLATFORM=Win32)
        endif()

        # Crucial for static linking: Force OpenBLAS to use the static runtime
        # if (NOT BUILD_SHARED_LIBS)
        #     set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
        # endif()
        # now global

        # list(APPEND OPENBLAS_MINIMAL_FLAGS -DUSE_THREAD=OFF)
        list(APPEND OPENBLAS_MINIMAL_FLAGS -DINTERFACE64=0)

        # If the MSVC runtime library issue persists, you can try this flag as well,
        # though CMAKE_GENERATOR_PLATFORM should usually be sufficient.
        # list(APPEND OPENBLAS_MINIMAL_FLAGS -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL)
    endif()

    if(UNIX AND NOT APPLE)
        execute_process(
            COMMAND bash -c "grep -m1 'model name' /proc/cpuinfo | grep -i skylake"
            RESULT_VARIABLE SKYLAKE_CHECK
            OUTPUT_QUIET
            ERROR_QUIET
        )

        if(SKYLAKE_CHECK EQUAL 0)
            message(STATUS "Skylake detected - adjusting OpenBLAS target to avoid register spills")
            set(OPENBLAS_TARGET "HASWELL" CACHE STRING "OpenBLAS target architecture" FORCE)
            set(NO_AVX512 ON CACHE BOOL "Disable AVX512" FORCE)
            # set(CMAKE_C_FLAGS_OPENBLAS "-DTARGET=HASWELL -DNO_AVX512=1")
        else()
            message(STATUS "NOT Skylake")
        endif()

        if(NO_AVX512)
            message(STATUS "NO_AVX512 - adjusting OpenBLAS possibly for valgrind")
            set(NO_AVX512 ON CACHE BOOL "Disable AVX512" FORCE)
            # set(CMAKE_C_FLAGS_OPENBLAS "-DTARGET=HASWELL -DNO_AVX512=1")
        endif()
    endif()

    set(OPENBLAS_BUILD_TYPE "Release" CACHE STRING "Build type for OpenBLAS" FORCE)

    # Override CMAKE_BUILD_TYPE for OpenBLAS subdirectory
    set(CMAKE_BUILD_TYPE_BACKUP ${CMAKE_BUILD_TYPE})
    set(CMAKE_BUILD_TYPE Release)

    message(CHECK_START "Fetching OpenBLAS")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    FetchContent_Declare(
        openblas
        GIT_REPOSITORY "https://github.com/OpenMathLib/OpenBLAS.git"
        GIT_TAG        "v0.3.30"
        GIT_SHALLOW TRUE
        UPDATE_COMMAND git reset --hard
        CMAKE_ARGS
            ${OPENBLAS_MINIMAL_FLAGS}
            # Force optimization even in Debug builds to avoid register spills
            # Force high optimization and strip debug symbols for the kernels
            # -DCMAKE_BUILD_TYPE=Release
            # -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS_OPENBLAS}
            # -DCMAKE_C_FLAGS="-O3 -fomit-frame-pointer"
            # -DCMAKE_C_FLAGS_RELEASE="-O3 -fomit-frame-pointer"
            # -DCMAKE_C_FLAGS_DEBUG="-O2 -fomit-frame-pointer"
            # -DCORE_OPTIMIZATION="-O3"
    )
    FetchContent_MakeAvailable(openblas)

    if (ALL_TESTS)
        set(BUILD_TESTING ON)
    endif()

    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_BACKUP})

    list(POP_BACK CMAKE_MESSAGE_INDENT)
    message(CHECK_PASS "fetched")

    if (TARGET openblas)
        get_target_property(_openblas_aliased openblas ALIASED_TARGET)
        if(_openblas_aliased)
            set(_openblas_target ${_openblas_aliased})
            message(STATUS "OpenBLAS is an alias for: ${_openblas_target}")
        else()
            set(_openblas_target openblas)
        endif()
    elseif (TARGET openblas_static)
        set(_openblas_target openblas_static)
    elseif (TARGET openblas_shared)
        set(_openblas_target openblas_shared)
    else()
        message(FATAL_ERROR "OpenBLAS target not found")
    endif()
    message(STATUS "OpenBLAS target: ${_openblas_target}")

    return()
endif()

# Find BLAS
set(BLAS_ROOT "" CACHE STRING "Root directory of BLAS or OpenBLAS")
if (NOT BLAS_ROOT STREQUAL "")
    message(STATUS "BLAS_ROOT is " ${BLAS_ROOT})
endif()

set(USE_CMAKE_FIND_BLAS ON)

# Optionally set the vendor:
# set(BLA_VENDOR libblastrampoline)

if (NOT USE_CMAKE_FIND_BLAS)
    if (WIN32)
        if (NOT (BLAS_ROOT STREQUAL ""))
            message(STATUS "Looking for blas in " ${BLAS_ROOT})
            set(OpenBLAS_ROOT ${BLAS_ROOT})
            message(STATUS "OpenBLAS_ROOT is ${OpenBLAS_ROOT} ")
            find_package(OpenBLAS CONFIG NO_DEFAULT_PATH)

            if(OpenBLAS_FOUND)
                message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
            else()
                message(STATUS "OpenBLAS not found in ${BLAS_ROOT}")
            endif()
        endif()

        if ((BLAS_ROOT STREQUAL "") OR (NOT OpenBLAS_FOUND))
            # (NOT OpenBLAS_FOUND AND NOT BLAS_FOUND))
            message(STATUS "Looking for blas")

            find_package(OpenBLAS REQUIRED)

            if(OpenBLAS_FOUND)
                if(TARGET OpenBLAS::OpenBLAS)
                    message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
                elseif(OPENBLAS_LIB)
                    message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
                else()
                    # try blas
                    # find_package(BLAS)
                    # if (BLAS_FOUND)
                    #     message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                    #     message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                    # else()
                    #     message(STATUS "OpenBLAS found but no target?")
                    # endif()
                endif()
            else()
                message(FATAL_ERROR "No BLAS library found")
            endif()
        endif()

    elseif(NOT APPLE)
        # LINUX

        # If a BLAS install was specified try to use it first.
        if (NOT (BLAS_ROOT STREQUAL ""))
            message(STATUS "Looking for blas in " ${BLAS_ROOT})

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
        if ((BLAS_ROOT STREQUAL "") OR (NOT OPENBLAS_LIB AND NOT BLAS_LIB))

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

    if (WIN32 AND NOT BLAS_LIBRARIES AND NOT BLA_VENDOR)
        find_package(OpenBLAS CONFIG)
        if(OpenBLAS_FOUND)
            message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
        endif()
    endif()

    if (NOT OpenBLAS_FOUND)
        if (NOT BLA_VENDOR)
            if (APPLE)
                set (BLA_VENDOR Apple)
            elseif(LINUX)
                set (BLA_VENDOR OpenBLAS)
            elseif(WIN32)
                set (BLA_VENDOR OpenBLAS)
            endif()

            find_package(BLAS QUIET)
            if (BLAS_FOUND)
                message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                if (BLAS_INCLUDE_DIRS)
                    message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                endif()
            else()
                unset(BLA_VENDOR)
            endif()
        else()
            message(STATUS "Specified BLA_VENDOR: ${BLA_VENDOR}")
        endif()

        if (NOT BLAS_FOUND)
            find_package(BLAS REQUIRED)
            if (BLAS_FOUND)
                message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
                if (BLAS_INCLUDE_DIRS)
                    message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
                endif()
            else()
                message(FATAL_ERROR "No BLAS library found!")
            endif()
        endif()
    endif()
endif()
