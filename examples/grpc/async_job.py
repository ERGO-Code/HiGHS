#!/usr/bin/env python3
"""
Async job mode example: SubmitSolve -> poll GetResult -> fetch result.

Suitable for: large models, long solves (>10s), avoids holding connection,
supports cancellation, resumable on disconnect.

Model: production planning LP (same as pyomo_lp.py), optimal obj=38

Flow:
  1. SubmitSolve(req) -> returns job_id immediately (connection ~1ms released)
  2. GetResult(job_id, wait=true, wait_timeout=2) -> short-poll until terminal
  3. Fetch result at terminal state

Run:
  1. Start server: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
  2. python examples/grpc/async_job.py
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

import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg


def build_production_lp():
    """Production planning LP: Max 3x1+5x2 s.t. 3 constraints, optimal x1=6,x2=4,obj=38"""
    return pb.SolveRequest(
        model_name="production_planning",
        sense=pb.OBJ_SENSE_MAX,
        col_cost=[3.0, 5.0],
        col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
        # 3 rows 2 cols CSC
        a_format_start=[0, 3, 5],       # col0 has 3 nonzeros, col1 has 2
        a_format_index=[0, 1, 2, 0, 2], # col0: rows 0,1,2; col1: rows 0,2
        a_format_value=[1.0, 1.0, 1.0, 2.0, 1.0],
        row_lower=[-1e30, -1e30, -1e30],
        row_upper=[14.0, 6.0, 10.0],
    )


def main(addr="localhost:50051"):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)

    # 0. Health check + adaptive solver selection
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    solver = "pdlp" if hc.gpu_available else "ipm"
    print(f"server: {hc.message} → solver={solver}")

    req = build_production_lp()
    req.options["solver"] = solver

    # 1. Submit task (returns job_id immediately)
    t0 = time.time()
    submit_resp = stub.SubmitSolve(req, timeout=5)
    job_id = submit_resp.job_id
    print(f"[submit] job_id={job_id} status={pb.JobStatus.Name(submit_resp.job_status)} "
          f"({(time.time()-t0)*1000:.1f}ms)")

    # 2. Poll / long-poll for result
    #    wait=true + wait_timeout=2: server blocks up to 2s, returns immediately if terminal
    #    Client loops until terminal (more efficient than pure client-side sleep polling)
    final = None
    poll_count = 0
    while True:
        poll_count += 1
        gr = stub.GetResult(pb.GetResultRequest(
            job_id=job_id, wait=True, wait_timeout=2.0
        ), timeout=10)
        print(f"[poll #{poll_count}] status={pb.JobStatus.Name(gr.job_status)} "
              f"elapsed={gr.elapsed_time:.3f}s")
        if gr.job_status in (pb.JOB_STATUS_SUCCEEDED, pb.JOB_STATUS_FAILED,
                             pb.JOB_STATUS_CANCELLED):
            final = gr
            break

    # 3. Output result
    if final.job_status == pb.JOB_STATUS_SUCCEEDED:
        r = final.result
        print(f"\n[result] model_status={pb.ModelStatus.Name(r.model_status)}")
        print(f"[result] objective={r.objective_value}  (expected 38)")
        print(f"[result] col_value={list(r.col_value)}  (expected [6.0, 4.0])")
        print(f"[result] solve_time={r.solve_time:.4f}s iter={r.iteration_count}")
        assert abs(r.objective_value - 38) < 1e-4, "obj mismatch"
        print("\n✓ async job mode verified")
    else:
        print(f"\n[!] job terminal but not SUCCEEDED: {pb.JobStatus.Name(final.job_status)} "
              f"msg={final.status_message}")
        sys.exit(1)

    # 4. Demo cancel (submit then cancel immediately)
    print("\n--- cancel demo ---")
    req2 = build_production_lp()
    req2.options["solver"] = solver
    req2.time_limit = 60  # long time limit for cancel demo
    sub2 = stub.SubmitSolve(req2, timeout=5)
    print(f"[submit] job_id={sub2.job_id}")
    cancel = stub.CancelSolve(pb.CancelRequest(job_id=sub2.job_id), timeout=5)
    print(f"[cancel] cancelled={cancel.cancelled} msg={cancel.message}")
    # Query final status
    gr2 = stub.GetResult(pb.GetResultRequest(job_id=sub2.job_id, wait=True, wait_timeout=2.0),
                         timeout=10)
    print(f"[final] status={pb.JobStatus.Name(gr2.job_status)} msg={gr2.status_message}")


if __name__ == "__main__":
    main()
