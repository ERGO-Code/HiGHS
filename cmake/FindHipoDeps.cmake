# BLAS
option(BLAS_ROOT "Root directory of BLAS or OpenBLAS" "")
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
option(METIS_ROOT "Root directory of METIS" "")
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