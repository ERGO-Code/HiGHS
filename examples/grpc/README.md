# HiGHS-Server gRPC Examples

Client examples for calling the highs-server gRPC service, from minimal to full modeling-tool workflows.

## 📋 Examples

| Example | Description | Dependencies |
|---|---|---|
| [`minimal_lp.py`](minimal_lp.py) | **Sync** minimal: construct `SolveRequest` by hand, no modeling tool | grpcio, grpcio-tools |
| [`pyomo_lp.py`](pyomo_lp.py) | **Sync** Pyomo modeling -> extract CSC -> solve (incl. MPS file path) | pyomo, highspy, grpcio, grpcio-tools |
| [`pyomo_async.py`](pyomo_async.py) | **Async** Pyomo modeling -> SubmitSolve -> progress monitor -> batch submit | pyomo, grpcio, grpcio-tools |
| [`async_job.py`](async_job.py) | **Async** job mode: SubmitSolve -> poll GetResult -> cancel | grpcio, grpcio-tools |
| [`progress_demo.py`](progress_demo.py) | **Progress**: high-freq poll showing real-time iteration/objective | grpcio, grpcio-tools |
| [`config_examples.md`](config_examples.md) | Config reference: 7 scenarios + HiGHS options table | - |

## 🚀 Run

### Prerequisites
1. Start highs-server:
   ```bash
   ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
   ```
2. Install Python deps:
   ```bash
   pip install grpcio grpcio-tools                  # minimal_lp.py, async_job.py, progress_demo.py
   pip install grpcio grpcio-tools pyomo highspy    # pyomo_lp.py, pyomo_async.py
   ```

> Examples auto-generate gRPC Python stubs from `server/protos/solver.proto` into `_generated/`; no manual `grpc_tools.protoc` needed.

### Run an example
```bash
python examples/grpc/minimal_lp.py
python examples/grpc/pyomo_async.py
```

### Expected output (GPU build)
```
server: GPU-enabled (CUPDLP_GPU=ON) → solver=pdlp
status: MODEL_STATUS_OPTIMAL
objective: 4.0
✓ verified
```

## 📐 Key Concepts

### CSC sparse matrix format
The request uses CSC (Compressed Sparse Column) for constraint matrix $A$:

```
a_format_start: length = num_col + 1, column pointers
a_format_index: row indices of nonzeros
a_format_value: values of nonzeros
```

For $A = \begin{pmatrix} 1 & 2 \end{pmatrix}$ (1 row, 2 cols):
```python
a_format_start = [0, 1, 2]   # col0 starts at 0, col1 at 1, end at 2
a_format_index = [0, 0]      # col0 nonzero at row0; col1 nonzero at row0
a_format_value = [1.0, 2.0]
```

### Objective sense
Use `sense=OBJ_SENSE_MAX` directly; no need to negate cost like some APIs. Server handles it.

### Adaptive solver
Examples call `Check` health check to get `gpu_available`, then pick `pdlp` (GPU) or `ipm` (CPU) adaptively — no manual env detection needed.

## 📂 Directory Structure
```
examples/grpc/
├── README.md              # this doc
├── minimal_lp.py          # sync minimal
├── pyomo_lp.py            # sync Pyomo
├── pyomo_async.py         # async Pyomo + progress + batch
├── async_job.py           # async basic + cancel
├── progress_demo.py       # progress reporting
├── config_examples.md     # config reference
└── _generated/            # auto-generated stubs (gitignored, generated at runtime)
    ├── solver_pb2.py
    └── solver_pb2_grpc.py
```

## 🔗 Related
- [Protocol definition](../../server/protos/solver.proto)
- [Server docs](../../server/README.md)
- [E2E tests](../../server/test/)
