# BLAS

find_package(BLAS REQUIRED)

if(BLAS_FOUND)
    message(STATUS "Using BLAS library: ${BLAS_LIBRARIES}")
    message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIRS}")
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
