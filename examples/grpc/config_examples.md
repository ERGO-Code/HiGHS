# 配置示例集

各种求解器/参数/场景的配置示例，通过 `SolveRequest.options` map 传入 HiGHS。

## 求解器选择

HiGHS 支持多种 LP/MIP 求解器，通过 `options["solver"]` 指定：

```python
# PDLP（一阶方法，GPU 加速，适合大规模稀疏 LP）
req.options["solver"] = "pdlp"

# IPM（内点法，CPU 上通常最快，适合中等规模 LP）
req.options["solver"] = "ipm"

# Simplex（单纯形法，适合需要基解的场景）
req.options["solver"] = "simplex"

# 不指定 → 服务端自适应：GPU 构建用 pdlp，CPU 构建用 ipm
```

## 场景 1：GPU 加速大规模 LP（PDLP）

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MIN,
    # ... 模型数据 ...
    options={
        "solver": "pdlp",
        "pdlp_iteration_limit": "10000",      # 最大迭代数
        "pdlp_optimality_tolerance": "1e-6",  # 最优性容差
        "time_limit": "300",                   # 5 分钟超时
    },
    time_limit=300,
)
```

## 场景 2：CPU 快速求解中等 LP（IPM）

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MIN,
    # ... 模型数据 ...
    options={
        "solver": "ipm",
        "run_crossover": "on",    # IPM 后跑 crossover 得基解
        "time_limit": "60",
    },
    time_limit=60,
)
```

## 场景 3：MIP 求解（整数规划）

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[3.0, 5.0, 2.0],
    col_lower=[0.0, 0.0, 0.0],
    col_upper=[10, 20, 15],
    col_integrality=[pb.VAR_INTEGER, pb.VAR_INTEGER, pb.VAR_CONTINUOUS],  # 前两个整数
    # ... 约束矩阵 ...
    options={
        "solver": "highs",         # MIP 用 highs 求解器
        "mip_rel_gap": "0.01",     # 相对间隙 1%
        "mip_max_nodes": "1000",   # 最大节点数
        "time_limit": "120",
    },
    time_limit=120,
)
```

## 场景 4：二值变量（0/1 整数）

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[1.0, 1.0, 1.0],
    col_lower=[0.0, 0.0, 0.0],
    col_upper=[1.0, 1.0, 1.0],     # 上下界 [0,1]
    col_integrality=[pb.VAR_BINARY, pb.VAR_BINARY, pb.VAR_BINARY],
    # ... 约束 ...
    options={"solver": "highs"},
)
```

## 场景 5：同步 vs 异步选择

```python
# 小任务（<10s）→ 同步 Solve，简单直接
resp = stub.Solve(req, timeout=60)

# 大任务（>10s）→ 异步 SubmitSolve，避免长连接
sub = stub.SubmitSolve(req, timeout=5)
# 轮询 GetResult...（见 async_job.py）
```

## 场景 6：带持久化的服务端启动

```bash
# 不持久化（默认，重启丢 job）
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2

# SQLite 持久化（重启后历史 job 可查）
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2 \
  --job-db /var/lib/highs-server/jobs.db

# CPU 多并发（生产环境）
./build/bin/highs_grpc_server --bind 0.0.0.0:50051 \
  --max-concurrent 8 --job-workers 8 \
  --job-db /var/lib/highs-server/jobs.db
```

## 场景 7：自适应 solver（推荐）

不指定 solver，让服务端根据构建类型自适应：

```python
# 客户端先查 GPU 能力
hc = stub.Check(pb.HealthCheckRequest())
# GPU 构建 → pdlp（GPU 加速），CPU 构建 → ipm（CPU 最优）
# 无需 options["solver"]，服务端自动选
req = pb.SolveRequest(..., time_limit=60)
resp = stub.Solve(req, timeout=60)
```

## 常用 HiGHS 选项速查

| 选项 | 值 | 说明 |
|---|---|---|
| `solver` | `pdlp`/`ipm`/`simplex`/`highs` | 求解器选择 |
| `time_limit` | 秒数 | 求解时间上限 |
| `pdlp_iteration_limit` | 整数 | PDLP 最大迭代 |
| `pdlp_optimality_tolerance` | 浮点 | PDLP 最优性容差 |
| `mip_rel_gap` | 浮点 | MIP 相对间隙 |
| `mip_max_nodes` | 整数 | MIP 最大节点 |
| `run_crossover` | `on`/`off` | IPM 后是否 crossover |
| `output_flag` | `true`/`false` | 求解器日志（服务端默认 false） |
| `threads` | 整数 | 并行线程数（simplex/IPM） |

> 完整选项见 [HiGHS options 文档](https://github.com/ERGO-Code/HiGHS#options)。
