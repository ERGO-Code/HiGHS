# HiGHS-Server

A gRPC solving service built on [HiGHS](https://github.com/ERGO-Code/HiGHS), with **CPU/GPU adaptive builds**, exposing HiGHS LP/MIP solving as a remote service.

> This directory is a subproject of the HiGHS fork, integrated via `HIGHS_BUILD_GRPC_SERVER=ON`. It does not modify HiGHS core logic.

---

## ✨ Features

- **Sync + Async dual mode**: small tasks use sync `Solve` (long connection, returns in seconds); large tasks use async `SubmitSolve`→`GetResult`→`CancelSolve` (short connections, avoids holding connection).
- **Progress reporting**: async jobs return real-time iteration count / objective / elapsed / MIP gap (HiGHS callback driven), so clients can show progress bars.
- **Job persistence**: optional SQLite backend (`--job-db`); historical jobs remain queryable after restart, unfinished jobs marked FAILED.
- **CPU/GPU adaptive**: build-time `nvcc` probe enables `CUPDLP_GPU` (PDLP on GPU) or falls back to CPU (default solver=ipm). Same source, same `configure.sh`, both environments build and run.
- **Capability reporting**: server reports `gpu_available` via health check; clients adaptively pick `pdlp` (GPU) or `ipm` (CPU) without manual env detection.
- **Engineering**: input validation (`INVALID_ARGUMENT`), concurrency control, cancellation propagation, graceful shutdown (SIGINT/SIGTERM), ResourceQuota, health check, server reflection.

---

## 📐 Architecture

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
                           │  │  └─ Fill solution+duals    │  │
                           │  └───────────────────────────┘  │
                           │  Health check: gpu_available    │
                           └─────────────────────────────────┘
```

**Adaptive loop**:
1. Build-time `configure.sh` probes `nvcc` → decides `CUPDLP_GPU`
2. Compile-time `HIGHS_SERVER_CUPDLP_GPU` macro baked into binary
3. Startup log prints `GPU-enabled` / `CPU-only`
4. Health check RPC reports `gpu_available`
5. Client `probe_gpu()` picks `solver=pdlp` (GPU) or `solver=ipm` (CPU)

---

## 🚀 Quick Start

### Requirements
- C++17 compiler (gcc-11/12; **nvcc 12.x does not support gcc-13+**)
- CMake ≥ 3.24
- gRPC + Protobuf (system or conda)
- **Optional** CUDA Toolkit 12.x (GPU acceleration)

### 1. Prepare dependencies (conda recommended, no sudo)

```bash
conda create -n highs-server python=3.11
conda activate highs-server

# C++ build deps (Tsinghua mirror, replace as needed)
conda install --channel https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ \
  cmake grpc protobuf gcc=12 gxx=12

# Python test deps
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple grpcio grpcio-tools
```

System install (needs sudo):
```bash
sudo apt install cmake g++-11 libgrpc++-dev libprotobuf-dev protobuf-compiler-grpc \
  python3-venv python3-grpcio python3-grpc-tools
```

### 2. Configure + build

```bash
# Adaptive (auto-detect CUDA, recommended)
./configure.sh
cmake --build build --parallel --target highs_grpc_server

# Force CPU-only (even if CUDA present)
./configure.sh --no-cuda

# Force GPU (error if no nvcc)
./configure.sh --cuda
```

> If system gcc-13 is too new for nvcc, append:
> `./configure.sh --cuda -- -DCMAKE_CXX_COMPILER=/usr/bin/g++-11`

### 3. Start service

```bash
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --max-concurrent 1 --job-workers 2
```

Output (GPU build):
```
HiGHS-Server listening on 127.0.0.1:50051 (max_concurrent=1, job_workers=2) GPU-enabled (CUPDLP_GPU=ON)
```

Startup flags:
| Flag | Default | Description |
|---|---|---|
| `--bind` | `127.0.0.1:50051` | Listen address (add TLS for public) |
| `--max-concurrent` | `1` | Sync `Solve` concurrency limit (GPU: 1, CPU: higher) |
| `--job-workers` | `1` | Async job worker count (GPU: 1, CPU: = core count) |
| `--job-db` | empty | Async job SQLite path; empty = no persistence (jobs lost on restart) |

### 4. Run tests

```bash
./server/test/run_e2e.sh
```

Expected:
```
[ok] health check: serving=True gpu_available=True msg='GPU-enabled (CUPDLP_GPU=ON)'
[ok] feasible obj=4.0 (solver=pdlp, gpu_build=True)
[ok] infeasible detected
[ok] bad args rejected
[ok] async job: obj=4.0 elapsed=0.002s
[ok] async cancel → CANCELLED
[ok] async not_found → NOT_FOUND
ALL TESTS PASSED
```

---

## 📡 Protocol

### Service methods

**Synchronous mode** (small tasks, long connection, returns in seconds)
| RPC | Type | Purpose |
|---|---|---|
| `Solve` | unary | Single solve (model < 4MB) |
| `SolveStream` | client stream | Sharded upload for large models (> 4MB), sync solve |

**Async job mode** (large tasks, short connections, avoids holding connection)
| RPC | Type | Purpose |
|---|---|---|
| `SubmitSolve` | unary | Submit task, returns `job_id` immediately (~ms) |
| `GetResult` | unary | Query / long-poll result (`wait`+`wait_timeout`), returns progress |
| `CancelSolve` | unary | Cancel job (only PENDING/RUNNING) |

**Common**
| RPC | Type | Purpose |
|---|---|---|
| `Check` | unary | Health check + GPU capability (`gpu_available`) |

### Model representation (CSC sparse format)
```
SolveRequest {
  sense: MIN/MAX              // objective sense, no need to negate cost
  col_cost[], col_lower[], col_upper[]      // column data
  col_integrality[]           // optional, var type (continuous/integer/binary) → MIP
  a_format_index[], a_format_start[], a_format_value[]  // CSC constraint matrix
  row_lower[], row_upper[]    // row bounds
  options: map<string,string> // solver params, e.g. {"solver":"pdlp","time_limit":60}
  time_limit: double          // solve time limit (seconds)
}
```

### Response

**Sync `SolveResponse`**
```
SolveResponse {
  model_status: enum           // OPTIMAL/INFEASIBLE/UNBOUNDED/...
  col_value[], col_dual[]      // primal + dual solution
  row_value[], row_dual[]      // constraint activity + dual
  objective_value              // objective
  iteration_count              // PDLP iteration count
  solve_time                   // solve elapsed (seconds)
}
```

**Async `GetResultResponse`** (with progress)
```
GetResultResponse {
  job_status: enum             // PENDING/RUNNING/SUCCEEDED/FAILED/CANCELLED
  status_message: string
  result: SolveResponse        // populated only when SUCCEEDED
  elapsed_time: double         // elapsed (seconds), available even when PENDING/RUNNING
  progress: SolveProgress {    // real-time progress while PENDING/RUNNING
    iteration_count            // current iterations (PDLP/simplex/ipm)
    objective_value            // current objective
    running_time               // elapsed seconds
    mip_gap                    // MIP gap (0 for non-MIP)
    mip_node_count             // MIP node count (0 for non-MIP)
  }
}
```

Progress is driven by HiGHS callbacks (`kCallbackLogging` + `kCallbackMipSolution`), covering PDLP/simplex/IPM/MIP.

---

## 🐍 Python Client Examples

### Sync Solve (small tasks)

```python
import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg

channel = grpc.insecure_channel('localhost:50051')
stub = pbg.HighsServiceStub(channel)

# Adaptive solver selection
hc = stub.Check(pb.HealthCheckRequest())
solver = "pdlp" if hc.gpu_available else "ipm"

req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[1.0, 1.0],
    col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
    a_format_start=[0, 1, 2], a_format_index=[0, 0], a_format_value=[1.0, 2.0],
    row_lower=[-1e30], row_upper=[4.0],
    options={"solver": solver},
)
resp = stub.Solve(req, timeout=60)
print(f"status={resp.model_status} obj={resp.objective_value}")
print(f"solution={list(resp.col_value)}")
```

### Async job mode (large tasks, avoid long connection)

```python
# 1. Submit, returns job_id immediately (~ms)
sub = stub.SubmitSolve(req, timeout=5)
job_id = sub.job_id

# 2. Long-poll for result (server blocks up to wait_timeout)
while True:
    gr = stub.GetResult(pb.GetResultRequest(job_id=job_id, wait=True, wait_timeout=2.0),
                        timeout=10)
    if gr.job_status in (pb.JOB_STATUS_SUCCEEDED, pb.JOB_STATUS_FAILED,
                         pb.JOB_STATUS_CANCELLED):
        break
    # Read progress while running
    print(f"running: iter={gr.progress.iteration_count} obj={gr.progress.objective_value}")

# 3. Fetch result
if gr.job_status == pb.JOB_STATUS_SUCCEEDED:
    print(f"obj={gr.result.objective_value}")

# 4. Cancel (optional)
stub.CancelSolve(pb.CancelRequest(job_id=job_id))
```

### Progress reporting

Async jobs return `progress` field while RUNNING:
```python
gr = stub.GetResult(pb.GetResultRequest(job_id=job_id, wait=False), timeout=5)
p = gr.progress
print(f"iter={p.iteration_count} obj={p.objective_value} "
      f"time={p.running_time}s mip_gap={p.mip_gap}")
```

More examples: [`examples/grpc/`](../examples/grpc/)

---

## 🔧 Build Options

| CMake option | Default | Description |
|---|---|---|
| `HIGHS_BUILD_GRPC_SERVER` | `OFF` | Enable gRPC server subproject |
| `CUPDLP_GPU` | `OFF` | HiGHS official option, enable PDLP GPU acceleration |
| `CMAKE_PREFIX_PATH` | auto | conda env / venv prefix (auto-injected by `configure.sh`) |

---

## ⚠️ Notes

1. **GPU support is compile-time**: `CUPDLP_GPU=ON` build uses GPU PDLP; `OFF` build runs CPU even if GPU present. Startup log clearly prints build type.
2. **nvcc vs gcc**: CUDA 12.x nvcc does not support gcc-13+; use gcc-11 or gcc-12 (conda gcc-12 works).
3. **Infeasible model construction**: HiGHS `passModel` validates `lower <= upper`; a single row with `lower > upper` is rejected. Express infeasibility via two contradictory constraints (e.g. `x>=5` and `x<=1`).
4. **Runtime library path**: if deps installed via conda, run `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib` before starting.
5. **Public deployment security**: default bind `127.0.0.1`; add TLS + auth interceptor for public exposure.
6. **Job persistence**: `--job-db` SQLite path makes jobs queryable across restarts; without it, jobs are in-memory only. Unfinished jobs marked FAILED on restart (cannot resume runtime state).
7. **Progress reporting limitation**: HiGHS `run()` callback frequency depends on solver (PDLP periodic logging); tiny problems may solve too fast to sample intermediate progress.

---

## 📂 Directory Structure

```
server/
├── CMakeLists.txt          # subproject build
├── README.md               # this doc
├── protos/
│   └── solver.proto        # gRPC protocol
└── src/
    ├── main.cpp            # entry: graceful shutdown, ResourceQuota, health check
    ├── service_impl.{h,cpp}# gRPC service: Solve/SolveStream/Submit/Get/Cancel/Check
    ├── job_store.{h,cpp}   # JobStore + WorkerPool + progress + SQLite
    ├── model_converter.{h,cpp} # proto -> HighsModel
    ├── validator.{h,cpp}   # input validation
    └── build_info.h        # compile-time GPU capability
test/
├── test_client.py          # E2E client (7 cases)
└── run_e2e.sh              # E2E runner
configure.sh                # adaptive build script
```

---

## 🧪 Verified Environments

| Environment | GPU | Status |
|---|---|---|
| Ubuntu 24.04 + CUDA 12.1 + RTX 4070 Ti + gcc-11 + conda gRPC 1.51 | ✅ | E2E 7/7 pass |
| Ubuntu 24.04 + no CUDA + gcc-11 + conda gRPC 1.51 | ❌ | E2E 7/7 pass |

---

## 🐳 Docker Deployment

See [root README](../README.md) and `Dockerfile.cpu` / `Dockerfile.gpu` / `docker-compose.yml`.

### CPU (default)
```bash
docker compose up -d --build
```

### GPU
```bash
docker compose --profile gpu up -d --build highs-server-gpu
```

---

## 📖 Related Docs

- [Execution plan](../plan.md)
- [Examples](../examples/grpc/)
- [HiGHS docs](https://github.com/ERGO-Code/HiGHS)
- [HiGHS GPU guide](https://github.com/ERGO-Code/HiGHS/blob/master/docs/src/guide/gpu.md)
- [cuPDLP-C](https://github.com/COPT-Public/cuPDLP-C)

---

## 📄 License

Inherits HiGHS [MIT License](../LICENSE.txt).
