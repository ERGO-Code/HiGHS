enable_language(CXX CUDA)

add_library(cudalin SHARED
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cuda_kernels.cu
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cuda_kernels.cuh
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cudalinalg.cuh
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cudalinalg.cu
)
set_target_properties(cudalin PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# set_target_properties(cudalin PROPERTIES CUDA_ARCHITECTURES native)

target_include_directories(cudalin PUBLIC "${CMAKE_CUDA_PATH}/include")


if (WIN32)
  target_link_libraries(cudalin ${CUDA_LIBRARY})
else()
  target_compile_definitions(cudalin
      PUBLIC
      # If the debug configuration pass the DEBUG define to the compiler
      $<$<CONFIG:Debug>:-DCUPDLP_DEBUG=1>
    )
  target_link_libraries(cudalin ${CUDA_LIBRARY} m)
endif()

# add a test
add_executable(testcudalin test_cuda_linalg.c)
add_executable(testcublas test_cublas.c)

set_target_properties(testcudalin PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#target_include_directories(cudalinalg PRIVATE ${CUPDLP_INCLUDE_DIR}/cuda)
target_link_libraries(testcudalin PRIVATE cudalin ${CUDA_LIBRARY})

set_target_properties(testcublas PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(testcublas PRIVATE cudalin ${CUDA_LIBRARY})

