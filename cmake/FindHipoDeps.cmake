# BLAS
option(BLAS_ROOT "Root directory of BLAS or OpenBLAS" "")
message(STATUS "BLAS_ROOT is " ${BLAS_ROOT})

if (WIN32)

    find_package(OpenBLAS CONFIG REQUIRED) 
    message(STATUS "OpenBLAS CMake config path: ${OpenBLAS_DIR}")

    find_package(metis CONFIG REQUIRED)
    message(STATUS "metis CMake config path: ${metis_DIR}")

    find_package(gklib CONFIG REQUIRED)
    message(STATUS "gklib CMake config path: ${gklib_DIR}")
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

if (NOT WIN32)
    # Find_package works if deps are installed with vcpkg.

    # METIS
    option(METIS_ROOT "Root directory of METIS" "")
    message(STATUS "METIS_ROOT is " ${METIS_ROOT})

    find_path(METIS_PATH 
        NAMES "metis.h" 
        REQUIRED
        HINTS "${METIS_ROOT}/include")

    message("Found Metis header at ${METIS_PATH}")

    find_library(METIS_LIB 
        NAMES metis 
        REQUIRED
        HINTS "${METIS_ROOT}/lib")
    message("Found Metis library at ${METIS_LIB}")

    # GKlib
    option(GKLIB_ROOT "Root directory of GKlib" "")
    message(STATUS "GKLIB_ROOT is " ${GKLIB_ROOT})

    find_path(GKLIB_PATH 
        NAMES "GKlib.h" REQUIRED
        HINTS "${GKLIB_ROOT}/include")

    message("Found GKlib header at ${GKLIB_PATH}")

    find_library(GKLIB_LIB 
        NAMES GKlib
        REQUIRED
        HINTS "${GKLIB_ROOT}/lib")

    message("Found GKlib library at ${GKLIB_LIB}")
endif()