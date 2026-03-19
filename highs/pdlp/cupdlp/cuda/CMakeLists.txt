if (HIGHS_GPU_LIB)
  return()
endif()

# separate libs for now

enable_language(CXX CUDA)

add_library(cudalin SHARED
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cuda_kernels.cu
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cuda_kernels.cuh
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cudalinalg.cuh
       ${CUPDLP_INCLUDE_DIR}/cuda/cupdlp_cudalinalg.cu
       ../../hipdlp/pdhg.cu
)

# latest has
# set_target_properties(cudalin PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_target_properties(cudalin PROPERTIES
CUDA_SEPARABLE_COMPILATION ON
#CUDA_ARCHITECTURES "60;70;75;80;86" )
#Deleted 60;70 due to errors
#nvcc fatal   : Unsupported gpu architecture 'compute_60'
#nvcc fatal   : Unsupported gpu architecture 'compute_70'
CUDA_ARCHITECTURES "75;80;86" )


# set_target_properties(cudalin PROPERTIES CUDA_ARCHITECTURES native)

target_include_directories(cudalin PUBLIC "$<BUILD_INTERFACE:${CMAKE_CUDA_PATH}/include>")

if (WIN32)
  target_link_libraries(cudalin ${CUDA_LIBRARY})  # linking happens here
  # else()
  #   # Stub build - no toolkit needed
  #   # Link stubs instead of real libs - they satisfy the linker but aren't shipped
  #   target_include_directories(cudalin PRIVATE
  #       ${CMAKE_SOURCE_DIR}/extern/cuda_stubs/include
  #   )
  #   target_link_directories(cudalin PRIVATE
  #       ${CMAKE_SOURCE_DIR}/extern/cuda_stubs/lib/Windows
  #   )
  #   target_link_libraries(cudalin PRIVATE cuda cublas cusparse cudart)

  # endif()
else() # LINUX
  # target_compile_definitions(cudalin
  #     PUBLIC
  #     # If the debug configuration pass the DEBUG define to the compiler
  #     $<$<CONFIG:Debug>:-DCUPDLP_DEBUG=1>
  #   )
  # if (NOT GPU_STUBS)
    target_link_libraries(cudalin ${CUDA_LIBRARY} m)
  # else()

    # target_link_options(cudalin PRIVATE
    #         -L${CMAKE_SOURCE_DIR}/extern/cuda_stubs/lib/Linux
    #         -lcuda -lcublas -lcusparse -lcudart
    #         -Wl,--allow-shlib-undefined
    # )

  # endif()
endif()

if (NOT GPU_STUBS)
  # add a test
  add_executable(testcudalin test_cuda_linalg.c)
  add_executable(testcublas test_cublas.c)

  set_target_properties(testcudalin PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  #target_include_directories(cudalinalg PRIVATE ${CUPDLP_INCLUDE_DIR}/cuda)
  if (NOT GPU_STUBS)
    target_link_libraries(testcudalin PRIVATE cudalin ${CUDA_LIBRARY})
  endif()

  set_target_properties(testcublas PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  if (NOT GPU_STUBS)
    target_link_libraries(testcublas PRIVATE cudalin ${CUDA_LIBRARY})

    target_include_directories(cudalin PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CUDA_PATH}/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../hipdlp>"  # Add this for pdhg_kernels.hpp
    )
  endif()

  if (GPU_STUBS)
    target_include_directories(cudalin PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../hipdlp>"  # Add this for pdhg_kernels.hpp
      "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/cuda_stubs>"  # Add this for stubs
    )
  endif()
endif()