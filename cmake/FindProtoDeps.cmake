
# BLAS
if(NOT APPLE)
    find_library(BLAS_LIB NAMES blas PATH /opt/intel/oneapi/mkl/latest/lib/)
    if(BLAS_LIB)
        message("Found BLAS library at ${BLAS_LIB}")
    else(BLAS_LIB)
        find_library(OPENBLAS_LIB NAMES openblas)
        if(NOT OPENBLAS_LIB)
            message(FATAL_ERROR "No BLAS library found")
        endif(NOT OPENBLAS_LIB)
        message("Found OpenBLAS library at ${OPENBLAS_LIB}")
    endif(BLAS_LIB)
endif(NOT APPLE)

# Metis
find_path(METIS_PATH NAMES "metis.h" REQUIRED HINTS "${METIS_DIR}/include")
message("Found Metis header at ${METIS_PATH}")

find_library(METIS_LIB NAMES metis REQUIRED HINTS "${METIS_PATH}/../lib")
message("Found Metis library at ${METIS_LIB}")

# GKlib
find_path(GKLIB_PATH NAMES "GKlib.h" REQUIRED HINTS "${GKLIB_DIR}/include" HINTS "${METIS_DIR}/include")
message("Found GKlib header at ${GKLIB_PATH}")

find_library(GKLIB_LIB NAMES GKlib REQUIRED HINTS "${GKLIB_PATH}/../lib")
message("Found GKlib library at ${GKLIB_LIB}")
