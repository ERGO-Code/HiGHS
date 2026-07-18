#!/usr/bin/env python3
"""
Progress reporting example: poll GetResult after submit, show real-time progress
(iteration count / objective / elapsed).

Suitable for: long solves where client needs a progress bar or monitoring.

Model: a 100-var random LP so PDLP runs enough iterations to observe progress.

Run:
  1. Start server: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2
  2. python examples/grpc/progress_demo.py
"""
import sys, os, subprocess, time, random

_EXAMPLE_DIR = os.path.dirname(os.path.abspath(__file__))
_PROTO_DIR = os.path.normpath(os.path.join(_EXAMPLE_DIR, "..", "..", "server", "protos"))
_GEN_DIR = os.path.join(_EXAMPLE_DIR, "_generated")
os.makedirs(_GEN_DIR, exist_ok=True)
if not os.path.exists(os.path.join(_GEN_DIR, "solver_pb2.py")):
    subprocess.check_call([
        sys.executable, "-m", "grpc_tools.protoc",
        f"-I{_PROTO_DIR}", f"--python_out={_GEN_DIR}", f"--grpc_python_out={_GEN_DIR}",
        os.path.join(_PROTO_DIR, "solver.proto"),
    ])
sys.path.insert(0, _GEN_DIR)

import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg


def build_random_lp(n=50):
    """Construct n-var random LP: min c·x s.t. Ax <= b, x>=0"""
    random.seed(42)
    m = n // 2  # number of constraints
    col_cost = [random.uniform(-1, 1) for _ in range(n)]
    col_lower = [0.0] * n
    col_upper = [1e30] * n
    row_lower = [-1e30] * m
    row_upper = [random.uniform(1, 10) for _ in range(m)]
    # CSC: ~m/2 nonzeros per column
    a_start = [0]
    a_index = []
    a_value = []
    for col in range(n):
        nnz = random.randint(1, max(1, m // 2))
        rows = sorted(random.sample(range(m), min(nnz, m)))
        for r in rows:
            a_index.append(r)
            a_value.append(random.uniform(-2, 2))
        a_start.append(len(a_index))
    return pb.SolveRequest(
        model_name=f"random_lp_{n}v{m}c",
        sense=pb.OBJ_SENSE_MIN,
        col_cost=col_cost, col_lower=col_lower, col_upper=col_upper,
        a_format_start=a_start, a_format_index=a_index, a_format_value=a_value,
        row_lower=row_lower, row_upper=row_upper,
    )


def main(addr="localhost:50051"):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)

    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    print(f"server: {hc.message}")

    # Build a larger LP to force PDLP to run many iterations
    req = build_random_lp(80)
    req.options["solver"] = "pdlp"
    req.options["pdlp_iteration_limit"] = "5000"  # ample iteration limit
    print(f"[model] {req.model_name}: {len(req.col_cost)} vars")

    # Submit
    sub = stub.SubmitSolve(req, timeout=5)
    print(f"[submit] job_id={sub.job_id[:12]}...")

    # High-freq poll for progress (wait=false returns current state immediately)
    print("\n[progress] real-time (sampled every 0.3s):")
    print(f"{'poll':>4} {'status':<20} {'iter':>6} {'objective':>14} {'time':>8}")
    poll = 0
    final = None
    while True:
        poll += 1
        gr = stub.GetResult(pb.GetResultRequest(job_id=sub.job_id, wait=False), timeout=5)
        p = gr.progress
        print(f"{poll:>4} {pb.JobStatus.Name(gr.job_status):<20} "
              f"{p.iteration_count:>6} {p.objective_value:>14.4f} {p.running_time:>8.3f}s")
        if gr.job_status in (pb.JOB_STATUS_SUCCEEDED, pb.JOB_STATUS_FAILED, pb.JOB_STATUS_CANCELLED):
            final = gr
            break
        time.sleep(0.3)

    print(f"\n[final] {pb.JobStatus.Name(final.job_status)} elapsed={final.elapsed_time:.3f}s")
    if final.job_status == pb.JOB_STATUS_SUCCEEDED:
        r = final.result
        print(f"[result] model_status={pb.ModelStatus.Name(r.model_status)} "
              f"obj={r.objective_value:.4f} iter={r.iteration_count} solve_time={r.solve_time:.3f}s")


if __name__ == "__main__":
    main()
