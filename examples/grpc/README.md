# HiGHS-Server gRPC 示例

本目录提供调用 highs-server gRPC 服务的示例代码，涵盖从最简调用到带建模工具的完整流程。

## 📋 示例一览

| 示例 | 说明 | 依赖 |
|---|---|---|
| [`minimal_lp.py`](minimal_lp.py) | **同步**最小示例：手构造 `SolveRequest`，不依赖建模工具 | grpcio, grpcio-tools |
| [`pyomo_lp.py`](pyomo_lp.py) | **同步** Pyomo 建模 → 提取 CSC 矩阵 → 求解（含 MPS 文件路径） | pyomo, highspy, grpcio, grpcio-tools |
| [`async_job.py`](async_job.py) | **异步** job 模式：SubmitSolve → 轮询 GetResult → 取消 | grpcio, grpcio-tools |
| [`progress_demo.py`](progress_demo.py) | **进度上报**：高频轮询展示实时迭代数/目标值 | grpcio, grpcio-tools |

## 🚀 运行

### 前置
1. 启动 highs-server 服务：
   ```bash
   ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
   ```
2. 安装 Python 依赖：
   ```bash
   pip install grpcio grpcio-tools          # minimal_lp.py
   pip install grpcio grpcio-tools pyomo highspy  # pyomo_lp.py
   ```

> 示例脚本会自动从 `server/protos/solver.proto` 生成 gRPC Python 桩到 `_generated/`，无需手工生成。

### 跑示例
```bash
# 最小示例（快速验证）
python examples/grpc/minimal_lp.py

# Pyomo 示例（生产计划 LP，演示两条路径）
python examples/grpc/pyomo_lp.py
```

### 预期输出（GPU 构建）
```
server: GPU-enabled (CUPDLP_GPU=ON) → 使用 solver=pdlp
status: MODEL_STATUS_OPTIMAL
objective: 4.0  (期望 4.0)
✓ 验证通过
```

## 📐 关键概念

### CSC 稀疏矩阵格式
gRPC 请求用 CSC（Compressed Sparse Column）格式传输约束矩阵 $A$：

```
a_format_start: 长度 = num_col + 1，列指针
a_format_index: 非零元素的行索引
a_format_value: 非零元素的值
```

对 $A = \begin{pmatrix} 1 & 2 \end{pmatrix}$（1 行 2 列）：
```python
a_format_start = [0, 1, 2]   # 列0 从0开始, 列1 从1开始, 结束于2
a_format_index = [0, 0]      # 列0的非零在行0, 列1的非零在行0
a_format_value = [1.0, 2.0]
```

### 目标方向
直接用 `sense=OBJ_SENSE_MAX`，无需像某些接口那样取负 cost。服务端处理。

### 自适应 solver
示例先调 `Check` 健康检查拿到 `gpu_available`，据此选 `pdlp`（GPU）或 `ipm`（CPU），无需人工判断环境。

## 📂 目录结构
```
examples/grpc/
├── README.md          # 本文档
├── minimal_lp.py      # 最小示例
├── pyomo_lp.py        # Pyomo 建模示例
└── _generated/        # 自动生成的桩（gitignore，运行时生成）
    ├── solver_pb2.py
    └── solver_pb2_grpc.py
```

## 🔗 相关
- [协议定义](../../server/protos/solver.proto)
- [服务端文档](../../server/README.md)
- [E2E 测试](../../test/)
