# BLAS
set(BLAS_ROOT "" CACHE STRING "Root directory of BLAS or OpenBLAS")
message(STATUS "BLAS_ROOT is " ${BLAS_ROOT})

if (WIN32)
    find_package(OpenBLAS CONFIG REQUIRED)
    message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")
elseif(NOT APPLE)
    # LINUX

    # If a BLAS install was specified try to use it first.
    if (NOT (BLAS_ROOT STREQUAL ""))
        message(STATUS "Looking for blas CMake targets file in " ${BLAS_ROOT})

        find_library(OPENBLAS_LIB
            NAMES openblas
            HINTS "${BLAS_ROOT}/lib"
            NO_DEFAULT_PATH)

        if(OPENBLAS_LIB)
            message("Found OpenBLAS library at ${OPENBLAS_LIB}")
        else()
            find_library(BLAS_LIB
                NAMES blas
                HINTS "${BLAS_ROOT}/lib"
                NO_DEFAULT_PATH)

            if(BLAS_LIB)
                message("Found BLAS library at ${BLAS_LIB}")
            else()
                message("Did not find blas library at ${BLAS_ROOT}")
                message("Attempting default locations search")
            endif()
    endif()
    if ((BLAS_ROOT STREQUAL "") OR 
        (NOT OPENBLAS_LIB and NOT BLAS_LIB))

        find_library(OPENBLAS_LIB
            NAMES openblas
            HINTS "${BLAS_ROOT}/lib")

        if(OPENBLAS_LIB)
            message("Found OpenBLAS library at ${OPENBLAS_LIB}")
        else()
            find_library(BLAS_LIB
                NAMES blas
                HINTS "${BLAS_ROOT}/lib")

            if(BLAS_LIB)
                message("Found BLAS library at ${BLAS_LIB}")
            else()
                message(FATAL_ERROR "No BLAS library found")
            endif()
        endif()
    endif()
endif()

# METIS
set(METIS_ROOT "" CACHE STRING "Root directory of METIS")
message(STATUS "METIS_ROOT is " ${METIS_ROOT})

# If a METIS install was specified try to use it first.
if (NOT (METIS_ROOT STREQUAL ""))
    message(STATUS "Looking for METIS CMake targets file in " ${METIS_ROOT})
    find_package(metis CONFIG NO_DEFAULT_PATH)
else()
    find_package(metis CONFIG)
endif()

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

    if(METIS_LIB)
        message("Found Metis library at ${METIS_LIB}")
    else()
        # METIS_ROOT was not successful
        message("Metis not found in METIS_PATH, fallback to default search.")
        if (NOT (METIS_ROOT STREQUAL ""))
            find_package(metis CONFIG)

            if (metis_FOUND)
                message(STATUS "metis CMake config path: ${metis_DIR}")
            else()
                message(FATAL_ERROR "No Metis library found")
            endif()
        endif()
    endif()
endif()

# GKlib optional for newer versions on ubuntu and macos
set(GKLIB_ROOT "" CACHE STRING "Root directory of GKlib")
if (NOT (GKLIB_ROOT STREQUAL ""))
    message(STATUS "GKLIB_ROOT is " ${GKLIB_ROOT})

    find_package(GKlib CONFIG NO_DEFAULT_PATH)

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
            NAMES "gklib.h" "GKlib.h"
            REQUIRED
            PATHS "${GKLIB_ROOT}/include"
            NO_DEFAULT_PATH)

        message("Found GKlib header at ${GKLIB_PATH}")

        find_library(GKLIB_LIB
            NAMES GKlib libGKlib
            REQUIRED
            PATHS "${GKLIB_ROOT}/lib"
            NO_DEFAULT_PATH)

        if(GKLIB_LIB)
            message("Found GKlib library at ${GKLIB_LIB}")
        else()
            # GKLIB_ROOT was not successful
            message("GKlib not found in GKLIB_PATH, fallback to default search.")
            find_package(GKlib CONFIG)

            if (GKlib_FOUND)
                message(STATUS "GKlib CMake config path: ${GKlib_DIR}")
            else()
                message(FATAL_ERROR "No GKLib library found")
            endif()
        endif()
    endif()
endif()