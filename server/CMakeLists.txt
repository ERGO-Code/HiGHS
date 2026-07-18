# highs-server gRPC 子项目 CMake
# 修复 v1 三大致命问题：① 不再 add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}) 自引用；
# ② 不硬编码 CUDA 开关（由 configure.sh 自适应传入 CUPDLP_GPU）；
# ③ protoc 命令补 DEPENDS 防并行竞态。
# 关键：服务端源码本身不依赖任何 CUDA 宏，CPU/GPU 二进制由同一份代码产出。

# 1. 依赖 gRPC（gRPC 的 config 会自动带上 Protobuf/absl/c-ares 等，无需单独找）
#    优先用系统/conda 预装的 config（CMAKE_PREFIX_PATH 已由 configure.sh 注入），
#    失败才 FetchContent（兜底，并设 CMAKE_POLICY_VERSION_MINIMUM 兼容 cmake>=4）
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
  # Protobuf 由 gRPC 带入，确保 protoc/plugin target 可用
  if(NOT TARGET protobuf::protoc)
    find_package(Protobuf CONFIG REQUIRED)
  endif()
endif()

# 2. 生成 proto（补全 DEPENDS，避免并行构建竞态）
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

# 3. server 可执行
add_executable(highs_grpc_server
  src/main.cpp
  src/service_impl.cpp
  src/model_converter.cpp
  src/validator.cpp
  ${proto_srcs}
  ${grpc_srcs}
)
target_include_directories(highs_grpc_server PRIVATE
  "${GENERATED_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}"        # HiGHS 根（找 Highs.h）
  "${CMAKE_BINARY_DIR}"        # HConfig.h（含 CUPDLP_GPU 宏定义）
)
target_link_libraries(highs_grpc_server PRIVATE
  highs
  gRPC::grpc++
  protobuf::libprotobuf
)
# 健康检查 API 在 gRPC >=1.51 由 grpc++ 主库 + 头文件提供，无需独立 target
# 服务反射：target 名随 gRPC 版本变化，按存在性链接以兼容
if(TARGET gRPC::grpcpp_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpcpp_reflection)
elseif(TARGET gRPC::grpc++_reflection)
  target_link_libraries(highs_grpc_server PRIVATE gRPC::grpc++_reflection)
endif()
target_compile_features(highs_grpc_server PRIVATE cxx_std_17)

# 4. 运行时打印编译期 GPU 开关（读 CUPDLP_GPU 变量注入宏，无需 CUDA 宏）
target_compile_definitions(highs_grpc_server PRIVATE
  HIGHS_SERVER_CUPDLP_GPU=$<BOOL:${CUPDLP_GPU}>
)
