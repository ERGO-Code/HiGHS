# highs-server gRPC subproject CMake
# Fixes v1 issues: (1) no self-referencing add_subdirectory; (2) no hardcoded CUDA flag
# (CUPDLP_GPU passed adaptively by configure.sh); (3) protoc DEPENDS to prevent parallel race.
# Key: server source does not depend on any CUDA macro; CPU/GPU binaries from same source.

# 1. Find gRPC (its config brings Protobuf/absl/c-ares automatically, no separate find)
#    Prefer system/conda preinstalled config (CMAKE_PREFIX_PATH injected by configure.sh);
#    fall back to FetchContent (with CMAKE_POLICY_VERSION_MINIMUM for cmake>=4 compat)
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
  # Protobuf brought in by gRPC; ensure protoc/plugin targets available
  if(NOT TARGET protobuf::protoc)
    find_package(Protobuf CONFIG REQUIRED)
  endif()
endif()

# 2. Generate proto (add DEPENDS to prevent parallel build race)
set(PROTO_FILE "${CMAKE_CURRENT_SOURCE_DIR}/protos/solver.proto")
set(GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${GENERATED_DIR}")

set(proto_srcs  "${GENERATED_DIR}/solver.pb.cc")
set(proto_hdrs  "${GENERATED_DIR}/solver.pb.h")
set(grpc_srcs   "${GENERATED_DIR}/solver.grpc.pb.cc")
set(grpc_hdrs   "${GENERATED_DIR}/solver.grpc.pb.h")

set(PROTOC       $<TARGET_FILE:protobuf::protoc>)
set(GRPC_PLUGIN  $<TARGET_FILE:gRPC::grpc_cpp_plugin>)

add_custom_command(
  OUTPUT  "${proto_srcs}" "${proto_hdrs}" "${grpc_srcs}" "${grpc_hdrs}"
  COMMAND ${PROTOC}
          --cpp_out="${GENERATED_DIR}"
          --grpc_out="${GENERATED_DIR}"
          -I "${CMAKE_CURRENT_SOURCE_DIR}/protos"
          --plugin=protoc-gen-grpc="${GRPC_PLUGIN}"
          "${PROTO_FILE}"
  DEPENDS "${PROTO_FILE}" protobuf::protoc gRPC::grpc_cpp_plugin
  COMMENT "Generating gRPC + protobuf sources for solver.proto"
)

# 3. server executable
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
  "${CMAKE_SOURCE_DIR}"        # HiGHS root (for Highs.h)
  "${CMAKE_BINARY_DIR}"        # HConfig.h (contains CUPDLP_GPU macro)
)
# SQLite3 (job persistence): prefer system, then conda env
find_package(SQLite3 QUIET)
if(SQLite3_FOUND)
  target_link_libraries(highs_grpc_server PRIVATE SQLite::SQLite3)
else()
  message(STATUS "SQLite3 cmake config not found, falling back to find_library")
  find_library(SQLITE3_LIB NAMES sqlite3)
  if(SQLITE3_LIB)
    target_link_libraries(highs_grpc_server PRIVATE ${SQLITE3_LIB})
  else()
    message(FATAL_ERROR "SQLite3 not found; install libsqlite3-dev or sqlite in conda env")
  endif()
endif()

target_link_libraries(highs_grpc_server PRIVATE
  highs
  gRPC::grpc++
  protobuf::libprotobuf
)
# Health check API in gRPC >=1.51 provided by grpc++ main lib + header, no separate target
# Reflection: target name varies by gRPC version; link by existence for compatibility
if(TARGET gRPC::grpcpp_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpcpp_reflection)
elseif(TARGET gRPC::grpc++_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpc++_reflection)
endif()
target_compile_features(highs_grpc_server PRIVATE cxx_std_17)

# 4. Runtime GPU build flag (inject macro from CUPDLP_GPU variable, no CUDA macro needed)
target_compile_definitions(highs_grpc_server PRIVATE
  HIGHS_SERVER_CUPDLP_GPU=$<BOOL:${CUPDLP_GPU}>
)
