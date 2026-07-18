# Configuration Examples

Various solver/param/scenario configurations, passed to HiGHS via `SolveRequest.options` map.

## Solver Selection

HiGHS supports multiple LP/MIP solvers via `options["solver"]`:

```python
# PDLP (first-order, GPU-accelerated, for large sparse LP)
req.options["solver"] = "pdlp"

# IPM (interior point, usually fastest on CPU, for medium LP)
req.options["solver"] = "ipm"

# Simplex (for cases needing a basic solution)
req.options["solver"] = "simplex"

# Unspecified → server adaptive: GPU build uses pdlp, CPU build uses ipm
```

## Scenario 1: GPU-accelerated large LP (PDLP)

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MIN,
    # ... model data ...
    options={
        "solver": "pdlp",
        "pdlp_iteration_limit": "10000",
        "pdlp_optimality_tolerance": "1e-6",
        "time_limit": "300",
    },
    time_limit=300,
)
```

## Scenario 2: CPU fast medium LP (IPM)

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MIN,
    # ... model data ...
    options={
        "solver": "ipm",
        "run_crossover": "on",    # run crossover after IPM for basic solution
        "time_limit": "60",
    },
    time_limit=60,
)
```

## Scenario 3: MIP (integer programming)

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[3.0, 5.0, 2.0],
    col_lower=[0.0, 0.0, 0.0],
    col_upper=[10, 20, 15],
    col_integrality=[pb.VAR_INTEGER, pb.VAR_INTEGER, pb.VAR_CONTINUOUS],  # first two integer
    # ... constraint matrix ...
    options={
        "solver": "highs",         # MIP uses highs solver
        "mip_rel_gap": "0.01",     # 1% relative gap
        "mip_max_nodes": "1000",
        "time_limit": "120",
    },
    time_limit=120,
)
```

## Scenario 4: Binary variables (0/1 integer)

```python
req = pb.SolveRequest(
    sense=pb.OBJ_SENSE_MAX,
    col_cost=[1.0, 1.0, 1.0],
    col_lower=[0.0, 0.0, 0.0],
    col_upper=[1.0, 1.0, 1.0],     # bounds [0,1]
    col_integrality=[pb.VAR_BINARY, pb.VAR_BINARY, pb.VAR_BINARY],
    # ... constraints ...
    options={"solver": "highs"},
)
```

## Scenario 5: Sync vs Async

```python
# Small task (<10s) → sync Solve, simple
resp = stub.Solve(req, timeout=60)

# Large task (>10s) → async SubmitSolve, avoids long connection
sub = stub.SubmitSolve(req, timeout=5)
# poll GetResult... (see async_job.py)
```

## Scenario 6: Server startup with persistence

```bash
# No persistence (default, jobs lost on restart)
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2

# SQLite persistence (jobs queryable across restarts)
./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2 \
  --job-db /var/lib/highs-server/jobs.db

# CPU high concurrency (production)
./build/bin/highs_grpc_server --bind 0.0.0.0:50051 \
  --max-concurrent 8 --job-workers 8 \
  --job-db /var/lib/highs-server/jobs.db
```

## Scenario 7: Adaptive solver (recommended)

Don't specify solver; let server pick based on build type:

```python
# Client probes GPU capability first
hc = stub.Check(pb.HealthCheckRequest())
# GPU build → pdlp (GPU-accelerated), CPU build → ipm (CPU-optimal)
# No need for options["solver"]; server picks automatically
req = pb.SolveRequest(..., time_limit=60)
resp = stub.Solve(req, timeout=60)
```

## Common HiGHS Options Reference

| Option | Values | Description |
|---|---|---|
| `solver` | `pdlp`/`ipm`/`simplex`/`highs` | Solver selection |
| `time_limit` | seconds | Solve time limit |
| `pdlp_iteration_limit` | int | PDLP max iterations |
| `pdlp_optimality_tolerance` | float | PDLP optimality tolerance |
| `mip_rel_gap` | float | MIP relative gap |
| `mip_max_nodes` | int | MIP max nodes |
| `run_crossover` | `on`/`off` | Run crossover after IPM |
| `output_flag` | `true`/`false` | Solver logging (server default false) |
| `threads` | int | Parallel threads (simplex/IPM) |

> Full options: [HiGHS options docs](https://github.com/ERGO-Code/HiGHS#options).
