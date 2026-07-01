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
    set(HIGHS_BLAS_COMPILE_DEFINITION "" PARENT_SCOPE)

    if(APPLE)
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_APPLE_BLAS PARENT_SCOPE)
    elseif(OpenBLAS_FOUND AND TARGET OpenBLAS::OpenBLAS)
        set(HIGHS_BLAS_TARGET OpenBLAS::OpenBLAS PARENT_SCOPE)
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_OPENBLAS PARENT_SCOPE)
    elseif(OPENBLAS_LIB)
        set(HIGHS_BLAS_TARGET "${OPENBLAS_LIB}" PARENT_SCOPE)
        set(HIGHS_BLAS_INCLUDE_DIRS "${OPENBLAS_INCLUDE_DIR}" PARENT_SCOPE)
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_OPENBLAS PARENT_SCOPE)
    elseif(TARGET BLAS::BLAS)
        set(HIGHS_BLAS_TARGET BLAS::BLAS PARENT_SCOPE)
    endif()
endfunction()

# set appropriate BLAS dependency licensing metadata
function(highs_configure_blas_metadata)
    if(APPLE)
        set(HIGHS_BLAS_VENDOR "Apple Accelerate" PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${CMAKE_OSX_DEPLOYMENT_TARGET}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE "Apple SDK license" PARENT_SCOPE)
    elseif(OpenBLAS_FOUND OR OPENBLAS_LIB)
        set(HIGHS_BLAS_VENDOR OpenBLAS PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${OpenBLAS_VERSION}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE BSD-3-Clause PARENT_SCOPE)
    elseif(BLA_VENDOR MATCHES "Intel|MKL")
        set(HIGHS_BLAS_VENDOR "Intel oneMKL" PARENT_SCOPE)
        set(HIGHS_BLAS_VERSION "${MKL_VERSION}" PARENT_SCOPE)
        set(HIGHS_BLAS_LICENSE "Intel Simplified Software License" PARENT_SCOPE)
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

        # Exclude components not used by HiGHS
        set(OPENBLAS_MINIMAL_FLAGS
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
            endif()
        endif()

        message(STATUS "Enabling DYNAMIC_ARCH for runtime CPU detection.")
        list(APPEND OPENBLAS_MINIMAL_FLAGS -DDYNAMIC_ARCH=ON)

        # CMAKE_SIZEOF_VOID_P is 4 for 32-bit and 8 for 64-bit
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            message(STATUS "32-bit target detected. Applying 32-bit configuration flags for OpenBLAS.")

            if(WIN32)
                list(APPEND OPENBLAS_MINIMAL_FLAGS -DCMAKE_GENERATOR_PLATFORM=Win32)
            endif()

            list(APPEND OPENBLAS_MINIMAL_FLAGS -DINTERFACE64=0)
        endif()

        # TODO: potentially improve (not great for cross-compilation)
        # can use cmake to read /proc/cpuinfo instead of using bash
        if(UNIX AND NOT APPLE)
            execute_process(
                    COMMAND bash -c "grep -m1 'model name' /proc/cpuinfo | grep -i skylake"
                    RESULT_VARIABLE SKYLAKE_CHECK
                    OUTPUT_QUIET
                    ERROR_QUIET
            )

            if(SKYLAKE_CHECK EQUAL 0)
                message(STATUS "Skylake detected - disabling AVX512 to avoid register spills")
                list(APPEND OPENBLAS_MINIMAL_FLAGS -DNO_AVX512=ON)
            else()
                message(STATUS "NOT Skylake")
            endif()

            if(NO_AVX512)
                message(STATUS "NO_AVX512 set - disabling AVX512 in OpenBLAS")
                list(APPEND OPENBLAS_MINIMAL_FLAGS -DNO_AVX512=ON)
            endif()
        endif()

        set(OPENBLAS_BUILD_TYPE "Release" CACHE STRING "Build type for OpenBLAS" FORCE)

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
                CMAKE_ARGS
                        ${OPENBLAS_MINIMAL_FLAGS}
        )
        set(NO_LAPACKE ON CACHE BOOL "" FORCE)
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
        set(HIGHS_BLAS_COMPILE_DEFINITION HIPO_USES_OPENBLAS)
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
            if(NOT BLA_VENDOR)
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
                message(STATUS "Specified BLA_VENDOR: ${BLA_VENDOR}")

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

    if(APPLE)
        target_link_libraries(${target_name} PRIVATE "-framework Accelerate")
    else()
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
    endif()
endfunction()