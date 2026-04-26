# Hot-kernel microbenchmark suite.
# Built only when -DHIGHS_PERF=ON is passed to CMake.
add_executable(highs_perf
  perf_main.cpp
  perf_fixtures.cpp
  perf_factor.cpp
  perf_spmv.cpp
  perf_saxpy.cpp)

if (FAST_BUILD)
  target_link_libraries(highs_perf PRIVATE highs)
else()
  target_link_libraries(highs_perf PRIVATE libhighs)
endif()

target_include_directories(highs_perf PRIVATE
  ${PROJECT_SOURCE_DIR}/highs
  ${PROJECT_SOURCE_DIR}/extern
  ${HIGHS_BINARY_DIR})

target_compile_definitions(highs_perf PRIVATE
  HIGHS_DIR="${PROJECT_SOURCE_DIR}")

set_target_properties(highs_perf PROPERTIES UNITY_BUILD OFF)
