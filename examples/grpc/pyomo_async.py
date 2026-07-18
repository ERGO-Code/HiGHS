#!/usr/bin/env python3
"""
Pyomo modeling + async job mode example.

Most practical production scenario combo:
  Pyomo build large model -> extract CSC -> SubmitSolve async -> progress monitor -> fetch result

Reuses build_model and pyomo_model_to_request from pyomo_lp.py (DRY),
demonstrates migrating a sync Pyomo workflow to async job mode seamlessly.

Model: production planning LP (same as pyomo_lp.py), optimal obj=38

Run:
  1. Start server: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051 --job-workers 2
  2. python examples/grpc/pyomo_async.py
"""
import sys, os, subprocess, time

# --- Auto-generate gRPC Python stubs (self-contained) ---
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
sys.path.insert(0, _EXAMPLE_DIR)  # import pyomo_lp from same dir

import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg
# Reuse build_model + pyomo_model_to_request from pyomo_lp.py (DRY)
from pyomo_lp import build_model, pyomo_model_to_request


def submit_and_monitor(stub, req, poll_interval=0.3):
    """
    Async submit + progress monitor + fetch result
    
    Returns final GetResultResponse (terminal state)
    """
    # 1. Submit, returns job_id immediately
    t0 = time.time()
    sub = stub.SubmitSolve(req, timeout=5)
    job_id = sub.job_id
    print(f"[submit] job_id={job_id[:12]}... ({(time.time()-t0)*1000:.1f}ms)")

    # 2. Progress monitor loop
    #    wait=false returns current state + progress immediately, for progress bars
    print(f"\n{'poll':>4} {'status':<20} {'iter':>6} {'objective':>14} {'time':>8}")
    poll = 0
    while True:
        poll += 1
        gr = stub.GetResult(
            pb.GetResultRequest(job_id=job_id, wait=False),
            timeout=5,
        )
        p = gr.progress
        print(f"{poll:>4} {pb.JobStatus.Name(gr.job_status):<20} "
              f"{p.iteration_count:>6} {p.objective_value:>14.4f} {p.running_time:>8.3f}s")
        if gr.job_status in (pb.JOB_STATUS_SUCCEEDED, pb.JOB_STATUS_FAILED,
                             pb.JOB_STATUS_CANCELLED):
            return gr
        time.sleep(poll_interval)


def main(addr="localhost:50051"):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)

    # 0. Health check + adaptive solver selection
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    solver = "pdlp" if hc.gpu_available else "ipm"
    print(f"server: {hc.message} → solver={solver}")
    print(f"  Note: async job mode suits large models with long solves; this demo uses a small model\n")

    # 1. Pyomo modeling (reuse pyomo_lp.build_model)
    model = build_model()
    print(f"[model] {model.name}")

    # 2. Pyomo -> CSC matrix -> SolveRequest (reuse pyomo_lp.pyomo_model_to_request)
    req = pyomo_model_to_request(model)
    req.options["solver"] = solver

    # 3. Async submit + progress monitor
    print("\n=== Async solve ===")
    final = submit_and_monitor(stub, req)

    # 4. Output result
    print(f"\n=== Result ===")
    print(f"job_status: {pb.JobStatus.Name(final.job_status)}")
    print(f"elapsed: {final.elapsed_time:.3f}s")
    if final.job_status == pb.JOB_STATUS_SUCCEEDED:
        r = final.result
        print(f"model_status: {pb.ModelStatus.Name(r.model_status)}")
        print(f"objective: {r.objective_value}  (expected 38)")
        print(f"col_value: {list(r.col_value)}  (expected [6.0, 4.0])")
        print(f"solve_time: {r.solve_time:.4f}s, iter: {r.iteration_count}")
        assert abs(r.objective_value - 38) < 1e-4, "obj mismatch"
        print("\n✓ Pyomo + async job mode verified")
    else:
        print(f"msg: {final.status_message}")
        sys.exit(1)

    # 5. Demo: batch submit with different Pyomo models
    print("\n=== batch submit demo (max/min variants) ===")
    jobs = []
    for sense_name, sense in [("max", pb.OBJ_SENSE_MAX), ("min", pb.OBJ_SENSE_MIN)]:
        m = build_model()
        m.profit.sense = pyo_maximize if sense_name == "max" else pyo_minimize
        r = pyomo_model_to_request(m)
        r.options["solver"] = solver
        r.model_name = f"{m.name}_{sense_name}"
        sub = stub.SubmitSolve(r, timeout=5)
        jobs.append((sub.job_id, sense_name, r.model_name))
        print(f"  submitted {r.model_name}: job_id={sub.job_id[:12]}...")

    # Wait for all to complete
    print("\n  waiting for all jobs to complete...")
    for job_id, sense_name, name in jobs:
        gr = stub.GetResult(
            pb.GetResultRequest(job_id=job_id, wait=True, wait_timeout=10),
            timeout=15,
        )
        obj = gr.result.objective_value if gr.job_status == pb.JOB_STATUS_SUCCEEDED else None
        print(f"  {name}: {pb.JobStatus.Name(gr.job_status)} obj={obj}")


# Pyomo sense constants (build_model uses maximize; here for switching demo)
import pyomo.environ as pyo
pyo_maximize = pyo.maximize
pyo_minimize = pyo.minimize


if __name__ == "__main__":
    main()
