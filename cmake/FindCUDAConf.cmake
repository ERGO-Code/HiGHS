set(CUDA_LIBRARY-NOTFOUND, OFF)
message(NOTICE "Finding CUDA environment")
message(NOTICE "    - CUDA Home detected at $ENV{CUDA_HOME}")
set(CMAKE_CUDA_ARCHITECTURES "all")
set(CMAKE_CUDA_PATH "$ENV{CUDA_HOME}")
set(CMAKE_CUDA_COMPILER "${CMAKE_CUDA_PATH}/bin/nvcc")

enable_language(CUDA)

find_library(CUDA_LIBRARY_ART
        NAMES cudart
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        REQUIRED
)
find_library(CUDA_LIBRARY_SPS
        NAMES cusparse
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        REQUIRED
)
find_library(CUDA_LIBRARY_BLS
        NAMES cublas
        HINTS "${CMAKE_CUDA_PATH}/lib64/"
        REQUIRED
)
if (${CUDA_LIBRARY-NOTFOUND})
    message(WARNING "    - CUDA Libraries not detected at $ENV{CUDA_HOME}")
else ()
    message(NOTICE "    - CUDA Libraries detected at $ENV{CUDA_HOME}")
    set(CUDA_LIBRARY ${CUDA_LIBRARY_ART} ${CUDA_LIBRARY_SPS} ${CUDA_LIBRARY_BLS})
    message(NOTICE "    -   :${CUDA_LIBRARY}")
endif ()