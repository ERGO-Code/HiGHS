# HiGHS-Server

基于 [HiGHS](https://github.com/ERGO-Code/HiGHS) 的 gRPC 在线求解服务，支持 **CPU / GPU 自适应构建**，将 HiGHS 的 LP/MIP 求解能力通过 gRPC 暴露为远程调用服务。

> 本目录是 HiGHS fork 仓库的子项目，通过 `HIGHS_BUILD_GRPC_SERVER=ON` 开关集成进主构建，不改动 HiGHS 核心逻辑。

---

## ✨ 特性

- **自适应 GPU**：构建期自动探测 CUDA，有则启用 `CUPDLP_GPU`（PDLP on GPU），无则 CPU-only（默认 solver=ipm）。同一份源码、同一套命令，两种环境皆可编译运行。
- **能力上报**：服务端通过 gRPC 健康检查上报 `gpu_available`，客户端据此自适应选择 `pdlp`（GPU）或 `ipm`（CPU），无需人工判断环境。
- **完整协议**：支持 LP/MIP、对偶解、目标方向（Max/Min 无需取负 hack）、int64 索引防溢出、流式分片上传（绕过 4MB 限制）。
- **工程化**：输入校验（`INVALID_ARGUMENT`）、并发串行化（GPU 资源保护）、取消传播、优雅退出（SIGINT/SIGTERM）、ResourceQuota、健康检查、服务反射。
- **零侵入**：仅 `add_subdirectory(server)` 挂入 HiGHS 主 CMake，与上游解耦，便于跟踪上游更新。

---

## 📐 架构

```
┌──────────────┐   gRPC    ┌─────────────────────────────────┐
│  Client      │──────────▶│  highs_grpc_server              │
│ (Python/Go/  │           │  ┌───────────────────────────┐  │
│  Java/...)   │◀──────────│  │ HighsServiceImpl           │  │
└──────────────┘  Solve    │  │  ├─ ValidateRequest        │  │
                           │  │  ├─ FillHighsModel         │  │
                           │  │  ├─ highs.run()            │  │
                           │  │  │   ├─ CPU: ipm/simplex   │  │
                           │  │  │   └─ GPU: pdlp+CUDA     │  │
                           │  │  └─ 装填解(含对偶)          │  │
                           │  └───────────────────────────┘  │
                           │  健康检查: gpu_available 上报    │
                           └─────────────────────────────────┘
```

**自适应闭环**：
1. 构建期 `configure.sh` 探测 `nvcc` → 决定 `CUPDLP_GPU` 开关
2. 编译期 `HIGHS_SERVER_CUPDLP_GPU` 宏固化进二进制
3. 启动日志打印 `GPU-enabled` / `CPU-only`
4. 健康检查 RPC 上报 `gpu_available` 字段
5. 客户端 `probe_gpu()` 据此选 `solver=pdlp`（GPU）或 `solver=ipm`（CPU）

---

## 🚀 快速开始

### 环境要求

- C++17 编译器（gcc-11/12，**注意 nvcc 12.x 不支持 gcc-13+**）
- CMake ≥ 3.24
- gRPC + Protobuf（系统装或 conda 装，见下）
- **可选** CUDA Toolkit 12.x（启用 GPU 加速）

### 1. 准备依赖（推荐 conda，免 sudo）

```bash
# 创建专用环境
conda create -n highs-server python=3.11
conda activate highs-server

# C++ 构建依赖（清华镜像，按需替换）
conda install --channel https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ \
  cmake grpc protobuf gcc=12 gxx=12

# Python 测试依赖
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple grpcio grpcio-tools
```

系统装方式（需 sudo）：
```bash
sudo apt install cmake g++-11 libgrpc++-dev libprotobuf-dev protobuf-compiler-grpc \
  python3-venv python3-grpcio python3-grpc-tools
```

### 2. 配置 + 编译

```bash
# 自适应（自动探测 CUDA，推荐）
./configure.sh
cmake --build build --parallel --target highs_grpc_server

# 强制 CPU-only（即便有 CUDA）
./configure.sh --no-cuda
cmake --build build --parallel --target highs_grpc_server

# 强制 GPU（无 nvcc 则报错）
./configure.sh --cuda
cmake --build build --parallel --target highs_grpc_server
```

> 若系统 gcc-13 太新导致 nvcc 报错，追加：
> `./configure.sh --cuda -- -DCMAKE_CXX_COMPILER=/usr/bin/g++-11`

### 3. 启动服务

```bash
# 默认绑 127.0.0.1:50051（仅本机，安全）
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --max-concurrent 1
```

输出（GPU 构建）：
```
HiGHS-Server listening on 127.0.0.1:50051 (max_concurrent=1) GPU-enabled (CUPDLP_GPU=ON)
```

启动参数：
| 参数 | 默认 | 说明 |
|---|---|---|
| `--bind` | `127.0.0.1:50051` | 监听地址（公网需配 TLS） |
| `--max-concurrent` | `1` | 最大并发求解数（GPU 建议 1，CPU 可调大） |

### 4. 跑测试

```bash
./test/run_e2e.sh
```

预期输出：
```
[ok] health check: serving=True gpu_available=True msg='GPU-enabled (CUPDLP_GPU=ON)'
[ok] feasible obj=4.0 (solver=pdlp, gpu_build=True)
[ok] infeasible detected
[ok] bad args rejected: col arrays size mismatch
ALL TESTS PASSED
```

---

## 📡 协议

完整定义见 `protos/solver.proto`，要点：

### 服务方法
| RPC | 类型 | 用途 |
|---|---|---|
| `Solve` | unary | 单次求解（模型 < 4MB） |
| `SolveStream` | client stream | 分片上传大模型（> 4MB） |
| `Check` | unary | 健康检查 + GPU 能力上报 |

### 模型表示（CSC 稀疏格式）
```
SolveRequest {
  sense: MIN/MAX              // 目标方向，无需取负
  col_cost[], col_lower[], col_upper[]      // 列信息
  col_integrality[]           // 可选，变量类型（连续/整数/二值）→ MIP
  a_format_index[], a_format_start[], a_format_value[]  // CSC 约束矩阵
  row_lower[], row_upper[]    // 约束上下界
  options: map<string,string> // 求解器参数，如 {"solver":"pdlp","time_limit":60}
  time_limit: double          // 求解时间上限（秒）
}
```

### 响应
```
SolveResponse {
  model_status: enum           // OPTIMAL/INFEASIBLE/UNBOUNDED/...
  col_value[], col_dual[]      // 原始解 + 对偶解
  row_value[], row_dual[]      // 约束活动值 + 对偶
  objective_value              // 目标值
  iteration_count              // PDLP 迭代数
  solve_time                   // 求解耗时（秒）
}
```

---

## 🐍 Python 客户端示例

```python
import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg

# 生成桩代码：python -m grpc_tools.protoc -Iserver/protos --python_out=. --grpc_python_out=. server/protos/solver.proto

channel = grpc.insecure_channel('localhost:50051')
stub = pbg.HighsServiceStub(channel)

# 1. 探测 GPU 能力，自适应选 solver
hc = stub.Check(pb.HealthCheckRequest())
solver = "pdlp" if hc.gpu_available else "ipm"

# 2. 构造 LP: Maximize x1 + x2  s.t. x1 + 2*x2 <= 4, x >= 0
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[1.0, 1.0],
    col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
    a_format_start=[0, 1, 2],      # CSC: 长度 = num_col+1
    a_format_index=[0, 0],          # 行索引
    a_format_value=[1.0, 2.0],      # 非零值
    row_lower=[-1e30], row_upper=[4.0],
    options={"solver": solver},
)

resp = stub.Solve(req, timeout=60)
print(f"status={resp.model_status} obj={resp.objective_value}")
print(f"solution={list(resp.col_value)}")
```

---

## 🔧 构建选项

| CMake 选项 | 默认 | 说明 |
|---|---|---|
| `HIGHS_BUILD_GRPC_SERVER` | `OFF` | 启用 gRPC server 子项目 |
| `CUPDLP_GPU` | `OFF` | HiGHS 官方选项，启用 PDLP GPU 加速 |
| `CMAKE_PREFIX_PATH` | 自动 | conda env / venv 前缀（`configure.sh` 自动注入） |

---

## ⚠️ 注意事项

1. **GPU 支持是编译期决定**：`CUPDLP_GPU=ON` 编译的二进制才能用 GPU PDLP；`OFF` 编译的即便机器有 GPU 也只跑 CPU。启动日志会明确打印构建类型。
2. **nvcc 与 gcc 版本**：CUDA 12.x 的 nvcc 不支持 gcc-13+，需用 gcc-11 或 gcc-12（conda 装的 gcc-12 可用）。
3. **不可行模型构造**：HiGHS `passModel` 校验 `lower <= upper`，单行 `lower > upper` 会被拒。不可行性需用两行矛盾约束表达（如 `x>=5` 和 `x<=1`）。
4. **运行时库路径**：若用 conda 装依赖，运行前需 `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib`。
5. **公网部署安全**：默认绑 `127.0.0.1`，公网暴露需自行加 TLS + 鉴权拦截器。

---

## 📂 目录结构

```
server/
├── CMakeLists.txt          # 子项目构建
├── README.md               # 本文档
├── protos/
│   └── solver.proto        # gRPC 协议定义
└── src/
    ├── main.cpp            # 入口：优雅退出、ResourceQuota、健康检查注册
    ├── service_impl.{h,cpp}# gRPC 服务实现：Solve/SolveStream/Check
    ├── model_converter.{h,cpp} # proto → HighsModel 转换
    ├── validator.{h,cpp}   # 输入校验
    └── build_info.h        # 编译期 GPU 能力上报
test/
├── test_client.py          # E2E 测试客户端（4 用例）
└── run_e2e.sh              # E2E 启动脚本
configure.sh                # 自适应构建脚本
```

---

## 🧪 已验证环境

| 环境 | GPU | 状态 |
|---|---|---|
| Ubuntu 24.04 + CUDA 12.1 + RTX 4070 Ti + gcc-11 + conda gRPC 1.51 | ✅ | E2E 4/4 通过 |
| Ubuntu 24.04 + 无 CUDA + gcc-11 + conda gRPC 1.51 | ❌ | E2E 4/4 通过 |

---

## 🐳 Docker 部署

提供 CPU / GPU 双镜像，通过 docker compose profile 切换。

### 前置条件
- Docker + Docker Compose v2
- **GPU 版额外需要** [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)

### CPU 版（默认）
```bash
docker compose up -d --build
# 服务启动在 localhost:50051，CPU-only 模式，max-concurrent=4
```

### GPU 版
```bash
docker compose --profile gpu up -d --build highs-server-gpu
# 服务启动在 localhost:50051，GPU 模式，max-concurrent=1（串行避免显存争抢）
```

### 验证
```bash
# 查看启动日志确认 GPU 构建状态
docker compose logs highs-server-gpu | head
# 预期: HiGHS-Server listening on 0.0.0.0:50051 (max_concurrent=1) GPU-enabled (CUPDLP_GPU=ON)

# 健康检查（需 grpcurl）
grpcurl -plaintext localhost:50051 highsserver.v1.HighsService/Check
```

### 文件说明
| 文件 | 用途 |
|---|---|
| `Dockerfile.cpu` | CPU 镜像（ubuntu:22.04 基础，无 CUDA 依赖，镜像更小） |
| `Dockerfile.gpu` | GPU 镜像（nvidia/cuda:12.4.1 基础，含 PDLP GPU 加速） |
| `docker-compose.yml` | 编排文件，CPU 默认 / GPU 用 `--profile gpu` |
| `.dockerignore` | 减小构建上下文（排除 build 产物、.git、缓存） |

---

## 📖 相关文档

- [执行计划](../plan.md) — 完整设计与 v1→v2 修复对照
- [HiGHS 官方文档](https://github.com/ERGO-Code/HiGHS) — 求解器本体
- [HiGHS GPU 加速指南](https://github.com/ERGO-Code/HiGHS/blob/master/docs/src/guide/gpu.md)
- [cuPDLP-C](https://github.com/COPT-Public/cuPDLP-C) — GPU PDLP 实现

---

## 📄 许可证

继承 HiGHS 的 [MIT License](../LICENSE.txt)。
