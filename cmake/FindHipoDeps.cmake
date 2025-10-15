# BLAS
set(BLAS_ROOT "" CACHE STRING "Root directory of BLAS or OpenBLAS")
message(STATUS "BLAS_ROOT is " ${BLAS_ROOT})

if (WIN32)
    find_package(OpenBLAS CONFIG REQUIRED)
    message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
elseif(NOT APPLE)
    # LINUX
    find_library(OPENBLAS_LIB
        NAMES openblas
        HINTS "${BLAS_ROOT}/lib")

    if(OPENBLAS_LIB)
        message("Found OpenBLAS library at ${OPENBLAS_LIB}")

    else(OPENBLAS_LIB)
        find_library(BLAS_LIB
            NAMES blas HINTS
            "${BLAS_ROOT}/lib")

        if(NOT BLAS_LIB)
            message(FATAL_ERROR "No BLAS library found")
        endif(NOT BLAS_LIB)
        message("Found BLAS library at ${BLAS_LIB}")
    endif(OPENBLAS_LIB)
endif()

# METIS
set(METIS_ROOT "" CACHE STRING "Root directory of METIS")
message(STATUS "METIS_ROOT is " ${METIS_ROOT})

find_package(metis CONFIG)

if(metis_FOUND)
    message(STATUS "metis CMake config path: ${metis_DIR}")
else()
    find_path(METIS_PATH
        NAMES "metis.h"
        REQUIRED
        PATHS "${METIS_ROOT}/include"
        NO_DEFAULT_PATH)

    message("Found Metis header at ${METIS_PATH}")

    find_library(METIS_LIB
        NAMES metis libmetis
        REQUIRED
        PATHS "${METIS_ROOT}/lib" "${METIS_ROOT}/bin"
        NO_DEFAULT_PATH)

    message("Found Metis library at ${METIS_LIB}")
endif()

# GKlib optional for newer versions on ubuntu and macos
set(GKLIB_ROOT "" CACHE STRING "Root directory of GKlib")
if (NOT (GKLIB_ROOT STREQUAL ""))
    message(STATUS "GKLIB_ROOT is " ${GKLIB_ROOT})

    find_package(GKlib CONFIG)

    if(GKlib_FOUND)
        message(STATUS "gklib CMake config path: ${GKlib_DIR}")

        # get_cmake_property(_vars VARIABLES)
        # foreach(_v IN LISTS _vars)
        #     if(_v MATCHES "GKlib")
        #     message(STATUS "${_v} = ${${_v}}")
        #     endif()
        # endforeach()

        # get_property(_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY IMPORTED_TARGETS)
        # foreach(_t IN LISTS _targets)
        # if(_t MATCHES "GKlib")
        #     message(STATUS "GKlib exported target: ${_t}")
        # endif()
        # endforeach()

    else()
        find_path(GKLIB_PATH
            NAMES "gklib.h"
            REQUIRED
            PATHS "${GKLIB_ROOT}/include"
            NO_DEFAULT_PATH)

        message("Found GKlib header at ${GKLIB_PATH}")

        find_library(GKLIB_LIB
            NAMES GKlib libGKlib
            REQUIRED
            PATHS "${GKLIB_ROOT}/lib"
            NO_DEFAULT_PATH)

        message("Found GKlib library at ${GKLIB_LIB}")
    endif()
endif()