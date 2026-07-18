# highs-server: standalone gRPC service for HiGHS
#
# This is an external interface to HiGHS (like highspy / C API / Fortran),
# NOT a subdirectory integrated into HiGHS main CMake. It links against an
# already-installed HiGHS via find_package(HiGHS).
#
# Build:
#   1. Install HiGHS first (with CUPDLP_GPU=ON if GPU desired):
#        cmake -S. -Bbuild-highs -DCMAKE_INSTALL_PREFIX=/tmp/highs-install -DCUPDLP_GPU=ON
#        cmake --build build-highs --target install
#   2. Build highs-server pointing to the install prefix:
#        cmake -Sserver -Bbuild-server -DCMAKE_PREFIX_PATH=/tmp/highs-install
#        cmake --build build-server --target highs_grpc_server
#
# Or use the provided configure.sh which automates both steps.
cmake_minimum_required(VERSION 3.24)
project(highs_server LANGUAGES CXX)

# --- 1. Find dependencies ---
# CUDAToolkit: required when HiGHS was built with CUPDLP_GPU (highs::cudalin
# depends on CUDA::cudart). Find it before HiGHS so the target is available.
find_package(CUDAToolkit QUIET)

# HiGHS (must be installed first; provides highs::highs target)
find_package(HiGHS REQUIRED)
message(STATUS "HiGHS found: ${HIGHS_VERSION}")

# gRPC + Protobuf (prefer system/conda config; fall back to FetchContent)
find_package(gRPC CONFIG QUIET)
if(NOT gRPC_FOUND)
  message(STATUS "gRPC not found via config; falling back to FetchContent")
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
  include(FetchContent)
  set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
  FetchContent_Declare(
    grpc
    GIT_REPOSITORY https://github.com/grpc/grpc.git
    GIT_TAG v1.65.0
  )
  FetchContent_MakeAvailable(grpc)
else()
  message(STATUS "gRPC found via config: ${gRPC_VERSION}")
  if(NOT TARGET protobuf::protoc)
    find_package(Protobuf CONFIG REQUIRED)
  endif()
endif()

# SQLite3 (job persistence)
find_package(SQLite3 QUIET)
if(NOT SQLite3_FOUND)
  message(STATUS "SQLite3 cmake config not found, falling back to find_library")
  find_library(SQLITE3_LIB NAMES sqlite3 REQUIRED)
endif()

# --- 2. Determine GPU build flag ---
# HiGHS installs HConfig.h which defines (or #undef) CUPDLP_GPU. We inspect the
# installed header to set HIGHS_SERVER_CUPDLP_GPU automatically. If that fails,
# fall back to an explicit cache variable the user can set.
if(NOT DEFINED HIGHS_SERVER_CUPDLP_GPU)
  # Try to locate HConfig.h from the imported target's INTERFACE_INCLUDE_DIRECTORIES
  get_target_property(_highs_inc highs::highs INTERFACE_INCLUDE_DIRECTORIES)
  if(_highs_inc)
    foreach(_dir ${_highs_inc})
      if(EXISTS "${_dir}/HConfig.h")
        file(STRINGS "${_dir}/HConfig.h" _cfg REGEX "define CUPDLP_GPU")
        if(_cfg)
          set(HIGHS_SERVER_CUPDLP_GPU ON)
        else()
          set(HIGHS_SERVER_CUPDLP_GPU OFF)
        endif()
        break()
      endif()
    endforeach()
  endif()
endif()
if(HIGHS_SERVER_CUPDLP_GPU)
  message(STATUS "highs-server: GPU-enabled (HiGHS was built with CUPDLP_GPU=ON)")
else()
  message(STATUS "highs-server: CPU-only (HiGHS was built with CUPDLP_GPU=OFF)")
endif()

# --- 3. Generate proto sources ---
set(PROTO_FILE "${CMAKE_CURRENT_SOURCE_DIR}/protos/solver.proto")
set(GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${GENERATED_DIR}")

set(proto_srcs "${GENERATED_DIR}/solver.pb.cc")
set(proto_hdrs "${GENERATED_DIR}/solver.pb.h")
set(grpc_srcs  "${GENERATED_DIR}/solver.grpc.pb.cc")
set(grpc_hdrs  "${GENERATED_DIR}/solver.grpc.pb.h")

add_custom_command(
  OUTPUT  "${proto_srcs}" "${proto_hdrs}" "${grpc_srcs}" "${grpc_hdrs}"
  COMMAND $<TARGET_FILE:protobuf::protoc>
          --cpp_out="${GENERATED_DIR}"
          --grpc_out="${GENERATED_DIR}"
          -I "${CMAKE_CURRENT_SOURCE_DIR}/protos"
          --plugin=protoc-gen-grpc=$<TARGET_FILE:gRPC::grpc_cpp_plugin>
          "${PROTO_FILE}"
  DEPENDS "${PROTO_FILE}" protobuf::protoc gRPC::grpc_cpp_plugin
  COMMENT "Generating gRPC + protobuf sources for solver.proto"
)

# --- 4. Build server executable ---
add_executable(highs_grpc_server
  src/main.cpp
  src/service_impl.cpp
  src/job_store.cpp
  src/model_converter.cpp
  src/validator.cpp
  ${proto_srcs}
  ${grpc_srcs}
)
target_include_directories(highs_grpc_server PRIVATE
  "${GENERATED_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(highs_grpc_server PRIVATE
  highs::highs
  gRPC::grpc++
  protobuf::libprotobuf
)
if(SQLite3_FOUND)
  target_link_libraries(highs_grpc_server PRIVATE SQLite::SQLite3)
else()
  target_link_libraries(highs_grpc_server PRIVATE ${SQLITE3_LIB})
endif()
# Reflection: target name varies by gRPC version; link by existence for compatibility
if(TARGET gRPC::grpcpp_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpcpp_reflection)
elseif(TARGET gRPC::grpc++_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpc++_reflection)
endif()
target_compile_features(highs_grpc_server PRIVATE cxx_std_17)
target_compile_definitions(highs_grpc_server PRIVATE
  HIGHS_SERVER_CUPDLP_GPU=$<BOOL:${HIGHS_SERVER_CUPDLP_GPU}>
)
