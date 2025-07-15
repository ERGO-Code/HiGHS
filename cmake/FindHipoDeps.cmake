# BLAS
option(BLAS_ROOT "Root directory of BLAS or OpenBLAS" "")
message(STATUS "BLAS_ROOT is " ${BLAS_ROOT})

if(NOT APPLE)
    find_library(BLAS_LIB 
        NAMES blas HINTS 
        "${BLAS_ROOT}/lib")

    if(BLAS_LIB)
        message("Found BLAS library at ${BLAS_LIB}")
    else(BLAS_LIB)
        find_library(OPENBLAS_LIB 
            NAMES openblas 
            HINTS "${BLAS_ROOT}/lib")
        if(NOT OPENBLAS_LIB)
            message(FATAL_ERROR "No BLAS library found")
        endif(NOT OPENBLAS_LIB)
        message("Found OpenBLAS library at ${OPENBLAS_LIB}")
    endif(BLAS_LIB)
    
endif(NOT APPLE)

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
message(STATUS "GKlib_ROOT is " ${GKLIB_ROOT})

find_path(GKLIB_PATH 
    NAMES "GKlib.h" REQUIRED
    HINTS "${GKLIB_ROOT}/include")

message("Found GKlib header at ${GKLIB_PATH}")

find_library(GKLIB_LIB 
    NAMES GKlib
    REQUIRED
    HINTS "${GKLIB_ROOT}/lib")

message("Found GKlib library at ${GKLIB_LIB}")
