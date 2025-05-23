set(CUDA_LIBRARY-NOTFOUND, OFF)
message(NOTICE "Finding CUDA environment")
message(NOTICE "    - CUDA Home detected at $ENV{CUDA_HOME}")
set(CMAKE_CUDA_ARCHITECTURES "all")

# On Windows users should set -DCMAKE_CUDA_PATH="..." when configuring CMake. 
# For all test setups CUPDLP_FIND_CUDA was not required on Windows.
if (NOT WIN32)
        set(CMAKE_CUDA_PATH "$ENV{CUDA_HOME}")
endif()
# For local testing.
# set(CMAKE_CUDA_PATH "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.9")

# message(STATUS "CMAKE_CUDA_PATH  ${CMAKE_CUDA_PATH}")

set(CMAKE_CUDA_COMPILER "${CMAKE_CUDA_PATH}/bin/nvcc")

enable_language(CUDA)

find_library(CUDA_LIBRARY_ART
        NAMES cudart
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        HINTS "${CMAKE_CUDA_PATH}/lib/x64"
        REQUIRED
)
find_library(CUDA_LIBRARY_SPS
        NAMES cusparse
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        HINTS "${CMAKE_CUDA_PATH}/lib/x64/"
        REQUIRED
)
find_library(CUDA_LIBRARY_BLS
        NAMES cublas
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        HINTS "${CMAKE_CUDA_PATH}/lib/x64/"
        REQUIRED
)
if (${CUDA_LIBRARY-NOTFOUND})
    message(WARNING "    - CUDA Libraries not detected at $ENV{CUDA_HOME}")
else ()
    message(NOTICE "    - CUDA Libraries detected at $ENV{CUDA_HOME}")
    set(CUDA_LIBRARY ${CUDA_LIBRARY_ART} ${CUDA_LIBRARY_SPS} ${CUDA_LIBRARY_BLS})
    message(NOTICE "    -   :${CUDA_LIBRARY}")
endif ()