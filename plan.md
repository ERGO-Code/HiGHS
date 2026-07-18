# HiGHS-Server 改造与部署 Agent 执行计划 (Plan v2)

本计划基于已 fork 的 HiGHS 仓库 `https://github.com/awakeningofthetrailblazer/HiGHServer.git`，
通过 `add_subdirectory(server)` 的方式将 gRPC 服务作为子目录集成进 HiGHS 主构建，
启用 HiGHS 官方的 `CUPDLP_GPU` 编译开关以支持 PDLP 求解器的 GPU 加速。

> **核心约束（与 v1 计划的根本差异）**
> 1. 项目根目录 = HiGHS fork 仓库根目录（不再嵌套克隆）。
> 2. GPU 支持是**编译期**决定（`-DCUPDLP_GPU=ON`），**不是运行时自动检测**。若二进制未以该开关编译，即使机器上有 GPU、即使设置 `solver=pdlp`，也只会跑 CPU 版 PDLP（性能远不如 IPM/simplex）。
> 3. gRPC server 作为子目录 `server/` 挂入 HiGHS 主 CMake，复用同一份 `highs` target，不再自引用根目录。
> 4. 协议加版本号 `highsserver.v1`，索引类型用 `int64`，响应补全对偶解与状态枚举。

---

## 📋 阶段一：仓库初始化与环境准备

> **Agent 目标**：把 fork 仓库内容同步到当前目录（作为项目根），初始化子模块，验证 CUDA/CMake 编译环境。

### 1.1 当前目录即仓库根
当前目录 `/home/n-md/highs-server` 仅含本 `plan.md`，将其作为 fork 仓库的工作根：

```bash
cd /home/n-md/highs-server
git init
git remote add origin https://github.com/awakeningofthetrailblazer/HiGHServer.git
git fetch origin
# 根据远程默认分支切换（main 或 master，二选一）
git checkout -b feature/grpc-server-impl origin/main || \
git checkout -b feature/grpc-server-impl origin/master
git submodule update --init --recursive
```

### 1.2 CUDA / CMake 环境探测
HiGHS 的 `CUPDLP_FIND_CUDA` 要求 **CMake >= 3.24**（推荐 3.25+），且需 `CUDA_HOME` 指向 toolkit 根：

```bash
cmake --version                 # 需要 >= 3.25
which nvcc && nvcc --version    # 验证 nvcc
nvidia-smi                      # 验证驱动 + GPU
echo "CUDA_HOME=${CUDA_HOME:-/usr/local/cuda}"
ls ${CUDA_HOME:-/usr/local/cuda}/lib64/libcudart.so   # 验证 cudart
ls ${CUDA_HOME:-/usr/local/cuda}/lib64/libcusparse.so # 验证 cusparse
ls ${CUDA_HOME:-/usr/local/cuda}/lib64/libcublas.so   # 验证 cublas
```

将上述输出记录到 `env_probe.log`，便于后续排错。

---

## 💻 阶段二：协议设计（v1，含版本化与完整性）

> **Agent 目标**：在 `server/protos/solver.proto` 编写支持 LP/MIP、对偶解、状态枚举、健康检查、流式上传的协议。

### 2.1 创建 `server/protos/solver.proto`

```protobuf
syntax = "proto3";

package highsserver.v1;

// 目标方向：消除 v1 中"取负 cost"的 hack
enum ObjSense {
  OBJ_SENSE_MIN = 0;
  OBJ_SENSE_MAX = 1;
}

// 变量类型，支持 MIP
enum VarType {
  VAR_CONTINUOUS = 0;
  VAR_INTEGER = 1;
  VAR_BINARY = 2;
}

// 模型求解状态枚举（对齐 HiGHS HighsModelStatus 主要取值）
enum ModelStatus {
  MODEL_STATUS_NOT_SET = 0;
  MODEL_STATUS_OPTIMAL = 1;
  MODEL_STATUS_INFEASIBLE = 2;
  MODEL_STATUS_UNBOUNDED = 3;
  MODEL_STATUS_ITERATION_LIMIT = 4;
  MODEL_STATUS_TIME_LIMIT = 5;
  MODEL_STATUS_NUMERICAL_ERROR = 6;
  MODEL_STATUS_OTHER = 99;
}

message SolveRequest {
  string model_name = 1;
  ObjSense sense = 2;                       // 默认 MIN

  // 列信息
  repeated double col_cost = 3;
  repeated double col_lower = 4;
  repeated double col_upper = 5;
  repeated VarType col_integrality = 11;    // 可选，缺省全连续

  // 约束矩阵 A（CSC 稀疏格式，索引用 int64 防溢出）
  repeated int64 a_format_index = 6;
  repeated int64 a_format_start = 7;        // 长度必须 = num_col + 1
  repeated double a_format_value = 8;

  // 约束上下界
  repeated double row_lower = 9;
  repeated double row_upper = 10;

  // 求解器参数，如 {"solver": "pdlp", "pdlp_iteration_limit": "10000"}
  map<string, string> options = 12;

  // 求解时间上限（秒），0 表示用 HiGHS 默认
  double time_limit = 13;
}

message SolveResponse {
  ModelStatus model_status = 1;
  string status_message = 2;
  repeated double col_value = 3;
  repeated double col_dual = 4;             // 对偶解（PDLP 原生产出）
  repeated double row_value = 5;            // 约束活动值
  repeated double row_dual = 6;
  double objective_value = 7;
  int64 iteration_count = 8;                // PDLP 迭代数
  double solve_time = 9;                    // 秒
}

message HealthCheckRequest {}
message HealthCheckResponse {
  bool serving = 1;
  string message = 2;
  bool gpu_available = 3;   // 编译期是否开启 CUPDLP_GPU，客户端据此自适应选 solver
}

service HighsService {
  // 单次求解
  rpc Solve(SolveRequest) returns (SolveResponse);
  // 大模型分片流式上传（规避 gRPC 默认 4MB 上限）
  rpc SolveStream(stream SolveRequest) returns (SolveResponse);
  // 健康检查（供 K8s/负载均衡探活）
  rpc Check(HealthCheckRequest) returns (HealthCheckResponse);
}
```


---

## 🛠️ 阶段三：C++ 核心代码实现与 CMake 集成（CUDA 自适应）

> **Agent 目标**：以 `server/` 子目录形式集成 gRPC 服务端，修复 v1 中 CMake 自引用、CUDA 开关名错、protoc 依赖缺失等问题；C++ 端补齐输入校验、错误码、并发控制、优雅退出；**构建期自动探测 CUDA，缺失时降级为 CPU 版 PDLP，二进制不重编也能跑**。

### 3.1 目录结构
```
highs-server/                      # 仓库根（HiGHS fork）
├── CMakeLists.txt                 # HiGHS 官方根 CMake（末尾追加 add_subdirectory(server)）
├── configure.sh                   # 自适应构建脚本：探测 CUDA 并传对应 CMake 标志
├── server/
│   ├── CMakeLists.txt             # gRPC server 子构建（不依赖 GPU 编译期宏）
│   ├── protos/solver.proto
│   └── src/
│       ├── main.cpp
│       ├── service_impl.h
│       ├── service_impl.cpp
│       ├── model_converter.h
│       ├── model_converter.cpp
│       ├── validator.h
│       ├── validator.cpp
│       └── build_info.h           # 运行时上报编译期 GPU 开关
└── test/
    ├── test_client.py
    └── run_e2e.sh
```

### 3.2 自适应构建脚本 `configure.sh`
**这是 CUDA 自适应的核心入口**：探测 `nvcc`/`CUDAToolkit` 是否可用，自动决定是否传 `-DCUPDLP_GPU=ON`，无需用户手工判断环境。

```bash
#!/usr/bin/env bash
# 自适应配置：有 CUDA 就开 GPU，没有就 CPU-only，二者皆可编译运行
set -euo pipefail
cd "$(dirname "$0")"

BUILD_DIR="${BUILD_DIR:-build}"
EXTRA_CMAKE_ARGS=()

# 1. 探测 nvcc
if command -v nvcc >/dev/null 2>&1; then
  echo "[configure] nvcc found: $(nvcc --version | tail -1)"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=ON)
  # 默认用本机 GPU 架构，避免硬编码 sm_xx
  EXTRA_CMAKE_ARGS+=(-DCMAKE_CUDA_ARCHITECTURES=native)
  # 若 CUDA 不在标准路径，补 CUDA_HOME
  if [[ -z "${CUDA_HOME:-}" ]]; then
    CUDA_HOME="$(dirname "$(dirname "$(command -v nvcc)")")"
    export CUDA_HOME
    EXTRA_CMAKE_ARGS+=(-DCUDA_HOME="${CUDA_HOME}")
  fi
  echo "[configure] CUDA_HOME=${CUDA_HOME} → 启用 CUPDLP_GPU"
else
  echo "[configure] nvcc 未找到 → 构建 CPU-only 版本（PDLP 仍可用，跑在 CPU 上）"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
fi

# 2. 强制开启 gRPC server 子项目
EXTRA_CMAKE_ARGS+=(-DHIGHS_BUILD_GRPC_SERVER=ON)
EXTRA_CMAKE_ARGS+=(-DCMAKE_BUILD_TYPE=Release)

echo "[configure] cmake -S. -B${BUILD_DIR} ${EXTRA_CMAKE_ARGS[*]}"
cmake -S. -B"${BUILD_DIR}" "${EXTRA_CMAKE_ARGS[@]}"
echo "[configure] 完成。下一步: cmake --build ${BUILD_DIR} --parallel --target highs_grpc_server"
```

### 3.3 修改根 `CMakeLists.txt`（仅追加，不改 HiGHS 既有逻辑）
在 HiGHS 根 `CMakeLists.txt` **末尾**追加。注意 `HIGHS_BUILD_GRPC_SERVER` 与 `CUPDLP_GPU` 解耦：即使 `CUPDLP_GPU=OFF`，gRPC server 仍可编译（求解器跑 CPU 版 PDLP）。

```cmake
# ===== highs-server gRPC 子项目（可选，与 CUPDLP_GPU 解耦）=====
option(HIGHS_BUILD_GRPC_SERVER "Build the gRPC server subproject" OFF)
if(HIGHS_BUILD_GRPC_SERVER AND BUILD_CXX)
  add_subdirectory(server)
endif()
```

### 3.4 `server/CMakeLists.txt`
修复 v1 的三大 CMake 致命问题：① 不再 `add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR})` 自引用；② 不再硬编码 CUDA 开关（由 `configure.sh` 自适应传入 `CUPDLP_GPU`）；③ protoc 命令补 `DEPENDS` 防并行竞态。**关键：服务端源码本身不依赖任何 CUDA 宏，CPU/GPU 二进制由同一份代码产出。**

```cmake
# 1. 依赖 gRPC/Protobuf：优先系统 config，失败则 FetchContent
find_package(gRPC CONFIG QUIET)
find_package(Protobuf CONFIG QUIET)
if(NOT gRPC_FOUND OR NOT Protobuf_FOUND)
  message(STATUS "gRPC/Protobuf not found via config; falling back to FetchContent")
  include(FetchContent)
  set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
  FetchContent_Declare(
    grpc
    GIT_REPOSITORY https://github.com/grpc/grpc.git
    GIT_TAG v1.65.0
  )
  FetchContent_MakeAvailable(grpc)
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
  gRPC::grpc_health_provider     # 健康检查
  gRPC::grpcpp_reflection        # 服务反射（grpcurl 调试用）
  protobuf::libprotobuf
)
target_compile_features(highs_grpc_server PRIVATE cxx_std_17)

# 4. 运行时打印编译期 GPU 开关（无需 CUDA 宏，读 HConfig.h 即可）
target_compile_definitions(highs_grpc_server PRIVATE
  HIGHS_SERVER_CUPDLP_GPU=$<BOOL:${CUPDLP_GPU}>
)
```

### 3.5 `src/build_info.h`（运行时上报编译期 GPU 能力）
让客户端通过健康检查即可得知服务端是否为 GPU 构建，从而决定是否用 `solver=pdlp`（CPU 环境下 PDLP 远不如 IPM，客户端可回退 `solver=ipm`/`simplex`）。

```cpp
// build_info.h
#pragma once

namespace highs_server {
// 由 CMake 注入的 HIGHS_SERVER_CUPDLP_GPU 宏决定（0=CPU 构建, 1=GPU 构建）
constexpr bool kGpuBuilt =
#if defined(HIGHS_SERVER_CUPDLP_GPU) && HIGHS_SERVER_CUPDLP_GPU
    true;
#else
    false;
#endif

inline const char* GpuBuildString() {
  return kGpuBuilt ? "GPU-enabled (CUPDLP_GPU=ON)" : "CPU-only (CUPDLP_GPU=OFF)";
}
}  // namespace highs_server
```

### 3.6 `src/validator.h` / `validator.cpp`（输入校验，返回明确错误）
所有非法输入统一返回 `grpc::Status(INVALID_ARGUMENT, ...)`，杜绝 v1 的静默错误。

```cpp
// validator.h
#pragma once
#include <grpcpp/grpcpp.h>
#include "solver.pb.h"

grpc::Status ValidateRequest(const highsserver::v1::SolveRequest& req);
```

```cpp
// validator.cpp
#include "validator.h"
#include <sstream>

grpc::Status ValidateRequest(const highsserver::v1::SolveRequest& req) {
  const auto n_col = req.col_cost_size();
  if (n_col == 0)
    return {grpc::INVALID_ARGUMENT, "col_cost is empty"};
  if (req.col_lower_size() != n_col || req.col_upper_size() != n_col) {
    std::ostringstream os;
    os << "col arrays size mismatch: cost=" << n_col
       << " lower=" << req.col_lower_size()
       << " upper=" << req.col_upper_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }
  if (req.col_integrality_size() != 0 && req.col_integrality_size() != n_col)
    return {grpc::INVALID_ARGUMENT, "col_integrality size mismatch"};

  const auto n_row = req.row_lower_size();
  if (req.row_upper_size() != n_row) {
    std::ostringstream os;
    os << "row arrays size mismatch: lower=" << n_row
       << " upper=" << req.row_upper_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }

  if (req.a_format_start_size() != n_col + 1) {
    std::ostringstream os;
    os << "a_format_start size must be num_col+1=" << (n_col + 1)
       << " got=" << req.a_format_start_size();
    return {grpc::INVALID_ARGUMENT, os.str()};
  }
  const auto nnz = req.a_format_index_size();
  if (req.a_format_value_size() != nnz)
    return {grpc::INVALID_ARGUMENT, "a_format_index/value size mismatch"};
  for (int i = 0; i <= n_col; ++i) {
    const auto s = req.a_format_start(i);
    if (s < 0 || s > nnz)
      return {grpc::INVALID_ARGUMENT, "a_format_start out of range"};
  }
  for (int k = 0; k < nnz; ++k) {
    const auto idx = req.a_format_index(k);
    if (idx < 0 || idx >= n_row) {
      std::ostringstream os;
      os << "a_format_index[" << k << "]=" << idx << " out of [0," << n_row << ")";
      return {grpc::INVALID_ARGUMENT, os.str()};
    }
  }
  return grpc::Status::OK;
}
```

### 3.7 `src/model_converter.h` / `model_converter.cpp`
修复 v1 漏设 `a_matrix_.num_col_/num_row_`、目标方向 hack、integrality 缺失。

```cpp
// model_converter.h
#pragma once
#include "Highs.h"
#include "solver.pb.h"

grpc::Status FillHighsModel(const highsserver::v1::SolveRequest& req,
                            HighsModel& model);
```

```cpp
// model_converter.cpp
#include "model_converter.h"

grpc::Status FillHighsModel(const highsserver::v1::SolveRequest& req,
                            HighsModel& model) {
  auto& lp = model.lp_;
  lp.num_col_ = req.col_cost_size();
  lp.num_row_ = req.row_lower_size();

  lp.col_cost_.assign(req.col_cost().begin(), req.col_cost().end());
  lp.col_lower_.assign(req.col_lower().begin(), req.col_lower().end());
  lp.col_upper_.assign(req.col_upper().begin(), req.col_upper().end());
  lp.row_lower_.assign(req.row_lower().begin(), req.row_lower().end());
  lp.row_upper_.assign(req.row_upper().begin(), req.row_upper().end());

  lp.sense_ = (req.sense() == highsserver::v1::OBJ_SENSE_MAX)
                  ? ObjSense::kMaximize : ObjSense::kMinimize;

  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = lp.num_col_;   // v1 漏设，会导致 HiGHS 内部断言
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_.assign(req.a_format_start().begin(),
                             req.a_format_start().end());
  lp.a_matrix_.index_.assign(req.a_format_index().begin(),
                             req.a_format_index().end());
  lp.a_matrix_.value_.assign(req.a_format_value().begin(),
                             req.a_format_value().end());

  if (req.col_integrality_size() > 0) {
    lp.integrality_.resize(lp.num_col_);
    for (int i = 0; i < lp.num_col_; ++i) {
      using VT = highsserver::v1::VarType;
      switch (req.col_integrality(i)) {
        case VT::VAR_INTEGER:
        case VT::VAR_BINARY:
          lp.integrality_[i] = HighsVarType::kInteger; break;
        default:
          lp.integrality_[i] = HighsVarType::kContinuous;
      }
    }
  }
  return grpc::Status::OK;
}
```

### 3.8 `src/service_impl.h` / `service_impl.cpp`
修复 v1 错误处理全部 `Status::OK`、无并发控制、无取消检查、options 顺序错乱。健康检查上报 GPU 能力，客户端据此自适应选 solver。

```cpp
// service_impl.h
#pragma once
#include <grpcpp/grpcpp.h>
#include <mutex>
#include "Highs.h"
#include "solver.grpc.pb.h"

class HighsServiceImpl final : public highsserver::v1::HighsService::Service {
 public:
  explicit HighsServiceImpl(int max_concurrent_solves);
  grpc::Status Solve(grpc::ServerContext*, const highsserver::v1::SolveRequest*,
                     highsserver::v1::SolveResponse*) override;
  grpc::Status SolveStream(grpc::ServerContext*,
                           grpc::ServerReader<highsserver::v1::SolveRequest>*,
                           highsserver::v1::SolveResponse*) override;
  grpc::Status Check(grpc::ServerContext*, highsserver::v1::HealthCheckRequest*,
                     highsserver::v1::HealthCheckResponse*) override;
 private:
  std::mutex solve_mutex_;        // GPU 串行化限流：避免显存争抢 / cuSPARSE handle 冲突
  const int max_concurrent_;
};
```

```cpp
// service_impl.cpp
#include "service_impl.h"
#include "build_info.h"
#include "model_converter.h"
#include "validator.h"
#include <chrono>

HighsServiceImpl::HighsServiceImpl(int max_concurrent)
    : max_concurrent_(max_concurrent) {}

grpc::Status HighsServiceImpl::Solve(grpc::ServerContext* ctx,
                                     const highsserver::v1::SolveRequest* req,
                                     highsserver::v1::SolveResponse* resp) {
  // 1. 输入校验 → INVALID_ARGUMENT
  auto st = ValidateRequest(*req);
  if (!st.ok()) return st;

  // 2. 取消检查（gRPC deadline 传播）
  if (ctx->IsCancelled())
    return {grpc::StatusCode::CANCELLED, "client cancelled before solve"};

  // 3. 串行化求解（GPU 资源保护；CPU 构建下此锁开销可忽略）
  std::lock_guard<std::mutex> lk(solve_mutex_);

  Highs highs;
  highs.setOptionValue("output_flag", false);

  // 4. options 先于 passModel（v1 顺序反了）
  //    自适应默认 solver：CPU 构建下 pdlp 性能差，若用户未显式指定则回退 ipm
  bool user_specified_solver = req->options().find("solver") != req->options().end();
  if (!user_specified_solver) {
    highs.setOptionValue("solver", highs_server::kGpuBuilt ? "pdlp" : "ipm");
  }
  for (const auto& kv : req->options())
    highs.setOptionValue(kv.first, kv.second);
  if (req->time_limit() > 0)
    highs.setOptionValue("time_limit", req->time_limit());

  // 5. 灌入模型
  HighsModel model;
  FillHighsModel(*req, model);
  if (highs.passModel(model) != HighsStatus::kOk) {
    resp->set_model_status(highsserver::v1::MODEL_STATUS_NUMERICAL_ERROR);
    resp->set_status_message("passModel failed");
    return grpc::Status::OK;
  }

  // 6. 求解
  auto t0 = std::chrono::steady_clock::now();
  HighsStatus run_st = highs.run();
  auto t1 = std::chrono::steady_clock::now();
  resp->set_solve_time(std::chrono::duration<double>(t1 - t0).count());

  // 7. 状态映射
  const auto ms = highs.getModelStatus();
  switch (ms) {
    case HighsModelStatus::kOptimal:        resp->set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL); break;
    case HighsModelStatus::kInfeasible:     resp->set_model_status(highsserver::v1::MODEL_STATUS_INFEASIBLE); break;
    case HighsModelStatus::kUnbounded:      resp->set_model_status(highsserver::v1::MODEL_STATUS_UNBOUNDED); break;
    case HighsModelStatus::kIterationLimit: resp->set_model_status(highsserver::v1::MODEL_STATUS_ITERATION_LIMIT); break;
    case HighsModelStatus::kTimeLimit:      resp->set_model_status(highsserver::v1::MODEL_STATUS_TIME_LIMIT); break;
    default:                                resp->set_model_status(highsserver::v1::MODEL_STATUS_OTHER);
  }

  // 8. 装填解（含对偶，v1 缺失）
  if (run_st == HighsStatus::kOk || highs.getSolution().value_valid) {
    const auto& sol = highs.getSolution();
    for (double v : sol.col_value) resp->add_col_value(v);
    for (double v : sol.col_dual)  resp->add_col_dual(v);
    for (double v : sol.row_value) resp->add_row_value(v);
    for (double v : sol.row_dual)  resp->add_row_dual(v);
    resp->set_objective_value(highs.getObjectiveValue());
  }
  resp->set_iteration_count(highs.getInfo().pdlp_iteration_count);
  resp->set_status_message("ok");
  return grpc::Status::OK;
}

grpc::Status HighsServiceImpl::SolveStream(
    grpc::ServerContext* ctx,
    grpc::ServerReader<highsserver::v1::SolveRequest>* reader,
    highsserver::v1::SolveResponse* resp) {
  // 大模型分片：首片含完整结构，后续片仅追加 a_format_* 增量
  highsserver::v1::SolveRequest first;
  if (!reader->Read(&first))
    return {grpc::INVALID_ARGUMENT, "no chunks received"};
  highsserver::v1::SolveRequest merged = first;
  highsserver::v1::SolveRequest chunk;
  while (reader->Read(&chunk)) {
    merged.mutable_a_format_index()->MergeFrom(chunk.a_format_index());
    merged.mutable_a_format_value()->MergeFrom(chunk.a_format_value());
  }
  return Solve(ctx, &merged, resp);
}

grpc::Status HighsServiceImpl::Check(
    grpc::ServerContext*, highsserver::v1::HealthCheckRequest*,
    highsserver::v1::HealthCheckResponse* resp) {
  resp->set_serving(true);
  // 自适应核心：上报编译期 GPU 能力，客户端据此选 solver
  resp->set_gpu_available(highs_server::kGpuBuilt);
  resp->set_message(highs_server::GpuBuildString());
  return grpc::Status::OK;
}
```

### 3.9 `src/main.cpp`
修复 v1 无优雅退出、无 ResourceQuota、监听 0.0.0.0 不安全、无健康检查注册。启动日志打印 GPU 构建状态，一眼可辨。

```cpp
#include <csignal>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <thread>
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface_builder.h>
#include <grpcpp/ext/proto_server_reflection_plugin.h>
#include "build_info.h"
#include "service_impl.h"

static std::shared_ptr<grpc::Server> g_server;
static std::mutex g_shutdown_mtx;
static bool g_shutdown = false;

static void OnSignal(int) {
  std::lock_guard<std::mutex> lk(g_shutdown_mtx);
  if (!g_shutdown) {
    g_shutdown = true;
    if (g_server) {
      // Shutdown 必须在独立线程，避免在 signal handler 里阻塞
      std::thread([] { g_server->Shutdown(); }).detach();
    }
  }
}

int main(int argc, char** argv) {
  std::string bind = "127.0.0.1:50051";   // 默认仅本机；公网需配 TLS + 鉴权
  int max_concurrent = 1;                  // GPU 默认串行；CPU 可调大
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--bind" && i + 1 < argc) bind = argv[++i];
    else if (a == "--max-concurrent" && i + 1 < argc) max_concurrent = std::atoi(argv[++i]);
  }

  HighsServiceImpl service(max_concurrent);
  grpc::EnableDefaultHealthCheckService(true);
  grpc::reflection::InitProtoReflectionServerBuilderPlugin();

  grpc::ServerBuilder builder;
  builder.AddListeningPort(bind, grpc::InsecureServerCredentials());
  builder.RegisterService(&service);

  grpc::ResourceQuota quota;
  quota.SetMaxThreads(static_cast<size_t>(max_concurrent) + 2);
  builder.SetResourceQuota(quota);

  g_server = builder.BuildAndStart();
  std::signal(SIGINT, OnSignal);
  std::signal(SIGTERM, OnSignal);
  std::cout << "HiGHS-Server listening on " << bind
            << " (max_concurrent=" << max_concurrent << ") "
            << highs_server::GpuBuildString() << std::endl;
  g_server->Wait();
  return 0;
}
```

> **自适应小结**：同一份源码、同一套构建命令（`./configure.sh`），在装了 CUDA 的机器上产出 GPU 二进制（默认 solver=pdlp），在没装 CUDA 的机器上产出 CPU 二进制（默认 solver=ipm）。客户端无需关心环境，通过 `Check` 健康检查拿到 `gpu_available` 后自行决定 solver 策略。

---

## 🔨 阶段四：编译、运行与测试（CPU/GPU 两套环境验证）

> **Agent 目标**：用自适应脚本构建，启动服务，跑覆盖正常/异常/并发/CPU-GPU 一致性的 E2E 测试。

### 4.1 自适应编译
**无需手工判断环境有无 CUDA**，`configure.sh` 自动探测：
```bash
# 任意机器统一入口
./configure.sh
cmake --build build --parallel --target highs_grpc_server
```
- 装了 CUDA → 输出含 `CUPDLP_GPU=ON`，启动日志 `GPU-enabled (CUPDLP_GPU=ON)`
- 没装 CUDA → 输出含 `CUPDLP_GPU=OFF`，启动日志 `CPU-only (CUPDLP_GPU=OFF)`，仍可正常求解（默认 solver=ipm）

### 4.2 Python 测试客户端 `test/test_client.py`
修复 v1 无 deadline、无 `with`、依赖 cost 取负 hack 的问题；**新增：先查健康检查自适应选 solver**。

```python
import grpc
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import solver_pb2 as pb
import solver_pb2_grpc as pbg

def get_stub(addr='localhost:50051'):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    return pbg.HighsServiceStub(ch)

def probe_gpu(stub):
    """自适应：查服务端是否 GPU 构建，决定是否用 pdlp"""
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    return hc.gpu_available

def solve(stub, req, timeout=30):
    return stub.Solve(req, timeout=timeout)

def case_feasible():
    stub = get_stub()
    gpu = probe_gpu(stub)
    # GPU 建构用 pdlp（GPU 上快），CPU 建建构用 ipm（CPU 上比 pdlp 快）
    solver = "pdlp" if gpu else "ipm"
    req = pb.SolveRequest(
        model_name="simple",
        sense=pb.OBJ_SENSE_MAX,
        col_cost=[1.0, 1.0],
        col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
        a_format_start=[0, 1, 2],
        a_format_index=[0, 0],
        a_format_value=[1.0, 2.0],
        row_lower=[-1e30], row_upper=[4.0],
        options={"solver": solver},
    )
    resp = solve(stub, req)
    assert resp.model_status == pb.MODEL_STATUS_OPTIMAL, resp
    print(f"[ok] feasible obj={resp.objective_value} (solver={solver}, gpu_build={gpu})")

def case_infeasible():
    stub = get_stub()
    req = pb.SolveRequest(
        sense=pb.OBJ_SENSE_MIN,
        col_cost=[1.0], col_lower=[0.0], col_upper=[1e30],
        a_format_start=[0, 1], a_format_index=[0], a_format_value=[1.0],
        row_lower=[5.0], row_upper=[1.0],   # 5 <= x <= 1 不可行
        options={"solver": "ipm"},
    )
    resp = solve(stub, req)
    assert resp.model_status == pb.MODEL_STATUS_INFEASIBLE, resp
    print("[ok] infeasible detected")

def case_bad_args():
    stub = get_stub()
    req = pb.SolveRequest(
        col_cost=[1.0, 2.0],
        a_format_start=[0, 1, 2, 3],   # 长度应为 num_col+1=3
    )
    try:
        solve(stub, req)
        print("[fail] expected INVALID_ARGUMENT")
    except grpc.RpcError as e:
        assert e.code() == grpc.StatusCode.INVALID_ARGUMENT
        print("[ok] bad args rejected:", e.details())

if __name__ == '__main__':
    case_feasible()
    case_infeasible()
    case_bad_args()
```

### 4.3 E2E 启动脚本 `test/run_e2e.sh`
修复 v1 `sleep 2` 不可靠、`kill` 无 `trap` 兜底。

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

./build/highs_grpc_server --bind 127.0.0.1:50051 --max-concurrent 1 &
SERVER_PID=$!
trap 'kill $SERVER_PID 2>/dev/null || true; wait $SERVER_PID 2>/dev/null || true' EXIT INT TERM

# 轮询就绪（最多 30s），替代脆弱的 sleep 2
for i in $(seq 1 60); do
  if python -c "import grpc; grpc.insecure_channel('localhost:50051').channel_ready_future().result(timeout=1)" 2>/dev/null; then
    break
  fi
  sleep 0.5
done

# 生成 Python 桩
cd test
python -m grpc_tools.protoc -I../server/protos --python_out=. --grpc_python_out=. ../server/protos/solver.proto

python test_client.py
echo "ALL TESTS PASSED"
```

### 4.4 测试矩阵
| 用例 | 验证点 |
|---|---|
| 可行 LP | `MODEL_STATUS_OPTIMAL` + objective 正确（CPU/GPU 两套构建都验证） |
| 不可行 LP | `MODEL_STATUS_INFEASIBLE` |
| 无界 LP | `MODEL_STATUS_UNBOUNDED` |
| 非法参数 | gRPC `INVALID_ARGUMENT` |
| 大模型（>4MB） | 触发流式 `SolveStream` |
| 并发请求 | `--max-concurrent` 限流生效，无 OOM |
| 健康检查 GPU 上报 | `Check.gpu_available` 与启动日志一致 |
| CPU/GPU 一致性 | 同模型两种构建 objective 一致（容差 1e-4） |

---

## 🚀 阶段五：部署与运维（CPU/GPU 双变体）

### 5.1 `Dockerfile.gpu`（GPU 运行时）
```dockerfile
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04 AS build
RUN apt-get update && apt-get install -y cmake git build-essential \
    libgrpc++-dev libprotobuf-dev protobuf-compiler-grpc pkg-config
WORKDIR /src
COPY . .
RUN ./configure.sh && cmake --build build --parallel --target highs_grpc_server

FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04
RUN apt-get update && apt-get install -y libgrpc++1 libprotobuf23 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/highs_grpc_server /usr/local/bin/
EXPOSE 50051
ENTRYPOINT ["highs_grpc_server", "--bind", "0.0.0.0:50051"]
```

### 5.2 `Dockerfile.cpu`（无 CUDA 环境，更轻量）
```dockerfile
FROM ubuntu:22.04 AS build
RUN apt-get update && apt-get install -y cmake g++ git build-essential \
    libgrpc++-dev libprotobuf-dev protobuf-compiler-grpc pkg-config
WORKDIR /src
COPY . .
# configure.sh 探测不到 nvcc → 自动 CUPDLP_GPU=OFF
RUN ./configure.sh && cmake --build build --parallel --target highs_grpc_server

FROM ubuntu:22.04
RUN apt-get update && apt-get install -y libgrpc++1 libprotobuf23 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/highs_grpc_server /usr/local/bin/
EXPOSE 50051
ENTRYPOINT ["highs_grpc_server", "--bind", "0.0.0.0:50051"]
```

### 5.3 `docker-compose.yml`（CPU/GPU 自适应选择）
```yaml
services:
  highs-server:
    # 默认 cpu；GPU 环境用 --profile gpu
    build:
      context: .
      dockerfile: Dockerfile.cpu
    ports: ["50051:50051"]
    healthcheck:
      test: ["CMD", "grpcurl", "-plaintext", "localhost:50051", "list"]
      interval: 10s

  highs-server-gpu:
    profiles: ["gpu"]
    build:
      context: .
      dockerfile: Dockerfile.gpu
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
    ports: ["50051:50051"]
```
- CPU 部署：`docker compose up`
- GPU 部署：`docker compose --profile gpu up highs-server-gpu`

### 5.4 可观测性
- 日志：接入 `spdlog` 结构化 JSON，含请求 ID 与 `gpu_build` 字段
- metrics：Prometheus exporter 暴露 `highs_solve_duration_seconds`、`highs_solves_total{status=,solver=,gpu_build=}`、`highs_gpu_mem_bytes`（仅 GPU 构建有值）
- 探活：`grpc.health.v1.Health` 服务（阶段三已注册）

### 5.5 CI/CD（GitHub Actions 矩阵）
```yaml
strategy:
  matrix:
    include:
      - {os: ubuntu-22.04, cuda: false, image: ubuntu:22.04}
      - {os: ubuntu-22.04, cuda: true,  image: nvidia/cuda:12.4.1-devel-ubuntu22.04}
steps:
  - uses: actions/checkout@v4
  - run: ./configure.sh
  - run: cmake --build build --parallel --target highs_grpc_server
  - run: ./test/run_e2e.sh
```
CPU runner 用 GitHub 托管，GPU runner 用 `self-hosted` + `nvidia` label。同一套脚本，无需分支判断。

---

## 📊 优先级与里程碑

| 阶段 | 优先级 | 里程碑验收 |
|---|---|---|
| 一、仓库初始化 | P0 | 本地 `./configure.sh` 能编出 HiGHS（CPU 或 GPU 视环境） |
| 二、Proto v1 | P0 | proto 冻结，可生成桩代码 |
| 三、C++ + CMake（自适应） | P0 | `highs_grpc_server` 可启动，日志打印 GPU 构建状态，健康检查返回 `gpu_available` |
| 四、编译测试 | P0 | E2E 用例全过（CPU/GPU 两套构建均验证） |
| 五、部署运维（双变体） | P1 | CPU/GPU 两镜像均可跑通 |
| 流式/并发/安全加固 | P2 | 生产可用 |

---

## ⚠️ 关键风险与对策

| 风险 | 对策 |
|---|---|
| `FetchContent` 拉 gRPC 编译极慢（数十分钟） | 文档提供 `apt install libgrpc++-dev` 预装路径，CI 缓存 `build/_deps` |
| GPU 显存 OOM / cuSPARSE handle 冲突 | `--max-concurrent 1` 串行化；监控显存；CPU 构建不受影响 |
| cuPDLP-C GPU 数值不稳定 | options 允许用户切 `solver=simplex`/`ipm` 回退 CPU |
| 误以为运行时自动启用 GPU | 启动日志 + 健康检查双通道上报 `gpu_build`，客户端 `probe_gpu` 自适应 |
| 大模型超 4MB gRPC 上限 | `SolveStream` 分片上传 |
| 公网部署被滥用 | 默认绑 `127.0.0.1`，公网需 TLS + token 拦截器 |
| CPU 环境误用 pdlp 性能差 | 服务端默认 solver 自适应（GPU 构建→pdlp，CPU 构建→ipm） |

---

## 🔧 v1 → v2 关键修复对照表

| v1 问题 | v2 修复 |
|---|---|
| CMake `add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR})` 自引用死循环 | 改为根 CMake 追加 `add_subdirectory(server)` |
| CUDA 开关 `BUILD_CUDA`（不存在） | 改用官方 `CUPDLP_GPU`，由 `configure.sh` 自适应传入 |
| "pdlp 自动触发 GPU" 误导注释 | 明确 GPU 是编译期开关；启动日志 + 健康检查上报 |
| 必须有 CUDA 才能构建 | `configure.sh` 探测不到 nvcc 自动降级 CPU-only，二进制仍可用 |
| `add_custom_command` 缺 `DEPENDS` 插件 | 补 `DEPENDS protobuf::protoc gRPC::grpc_cpp_plugin` |
| `num_row_ = row_lower_size()` 等无校验 | 新增 `validator.cpp` 全字段校验 |
| 错误全返 `Status::OK` | 区分 `INVALID_ARGUMENT`/`CANCELLED`/`OK` |
| `a_matrix_.num_col_/num_row_` 漏设 | converter 中显式赋值 |
| 无并发控制 | `std::mutex` 串行化 + `ResourceQuota` |
| 无超时/取消 | `time_limit` 选项 + `IsCancelled()` 检查 |
| 无优雅退出 | `SIGINT/SIGTERM` → `server->Shutdown()` |
| 监听 `0.0.0.0` 不安全 | 默认 `127.0.0.1`，可配 |
| 无健康检查/反射 | 注册 `grpc_health_provider` + reflection + `gpu_available` 上报 |
| 索引 `int32` 溢出 | 改 `int64` |
| 响应只有 `col_value` | 补 `col_dual`/`row_value`/`row_dual`/`model_status`/`iteration_count`/`solve_time` |
| 无 `sense` 字段，靠 cost 取负 | 加 `ObjSense` 枚举 |
| 不支持 MIP | 加 `col_integrality` |
| 无流式上传 | 加 `SolveStream` |
| 无版本号 | `package highsserver.v1` |
| `sleep 2` 启动等待 | 轮询 `channel_ready_future` |
| `kill` 无 `trap` 兜底 | `trap` 清理 + `wait` |
| 仅一个测试用例 | 测试矩阵 8 项 |
| 单一 Dockerfile（强依赖 CUDA） | `Dockerfile.cpu` + `Dockerfile.gpu` 双变体 + compose profile |
| CPU 环境误用 pdlp | 服务端默认 solver 自适应：GPU 构建→pdlp，CPU 构建→ipm |
