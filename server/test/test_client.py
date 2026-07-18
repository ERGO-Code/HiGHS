#!/usr/bin/env python3
"""E2E test client: covers feasible/infeasible/bad-args + async job/cancel/not-found."""
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
    """Adaptive: probe server GPU build to pick pdlp or ipm"""
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    return hc.gpu_available


def solve(stub, req, timeout=30):
    return stub.Solve(req, timeout=timeout)


def case_feasible():
    stub = get_stub()
    gpu = probe_gpu(stub)
    # GPU build uses pdlp (fast on GPU), CPU build uses ipm (faster than pdlp on CPU)
    solver = "pdlp" if gpu else "ipm"
    # Maximize x1 + x2 s.t. x1 + 2 x2 <= 4, x>=0  ->  optimal (4,0), obj=4
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
    assert resp.model_status == pb.MODEL_STATUS_OPTIMAL, \
        f"expected OPTIMAL, got {resp.model_status}: {resp.status_message}"
    assert abs(resp.objective_value - 4.0) < 1e-4, \
        f"expected obj=4.0, got {resp.objective_value}"
    print(f"[ok] feasible obj={resp.objective_value} "
          f"(solver={solver}, gpu_build={gpu})")


def case_infeasible():
    stub = get_stub()
    # Construct infeasible LP via two legal but contradictory constraints:
    #   row0: x >= 5    (row_lower=5,  row_upper=+inf)
    #   row1: x <= 1    (row_lower=-inf, row_upper=1)
    # i.e. 5 <= x <= 1, solver detects infeasible
    # Note: HiGHS passModel validates lower<=upper, cannot write lower>upper in one row
    req = pb.SolveRequest(
        sense=pb.OBJ_SENSE_MIN,
        col_cost=[1.0], col_lower=[0.0], col_upper=[1e30],
        # 2 rows 1 col, CSC: col0 has 2 nonzeros (row0 and row1)
        a_format_start=[0, 2], a_format_index=[0, 1], a_format_value=[1.0, 1.0],
        row_lower=[5.0, -1e30], row_upper=[1e30, 1.0],
        options={"solver": "ipm"},
    )
    resp = solve(stub, req)
    assert resp.model_status == pb.MODEL_STATUS_INFEASIBLE, \
        f"expected INFEASIBLE, got {resp.model_status}: {resp.status_message}"
    print("[ok] infeasible detected")


def case_bad_args():
    stub = get_stub()
    # a_format_start length should be num_col+1=3, here given 4
    req = pb.SolveRequest(
        col_cost=[1.0, 2.0],
        a_format_start=[0, 1, 2, 3],
    )
    try:
        solve(stub, req)
        print("[fail] expected INVALID_ARGUMENT")
        sys.exit(1)
    except grpc.RpcError as e:
        assert e.code() == grpc.StatusCode.INVALID_ARGUMENT, \
            f"expected INVALID_ARGUMENT, got {e.code()}"
        print(f"[ok] bad args rejected: {e.details()}")


def case_health_check():
    stub = get_stub()
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    assert hc.serving is True
    assert hc.gpu_available in (True, False)
    print(f"[ok] health check: serving={hc.serving} "
          f"gpu_available={hc.gpu_available} msg={hc.message!r}")


def case_async_job():
    """Async job mode: SubmitSolve -> GetResult poll -> verify result"""
    stub = get_stub()
    gpu = probe_gpu(stub)
    solver = "pdlp" if gpu else "ipm"
    req = pb.SolveRequest(
        sense=pb.OBJ_SENSE_MAX,
        col_cost=[1.0, 1.0],
        col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
        a_format_start=[0, 1, 2], a_format_index=[0, 0], a_format_value=[1.0, 2.0],
        row_lower=[-1e30], row_upper=[4.0],
        options={"solver": solver},
    )
    # Submit
    sub = stub.SubmitSolve(req, timeout=5)
    assert sub.job_id, "empty job_id"
    assert sub.job_status == pb.JOB_STATUS_PENDING
    # Poll until terminal
    import time
    t0 = time.time()
    while True:
        gr = stub.GetResult(pb.GetResultRequest(job_id=sub.job_id, wait=True, wait_timeout=2.0),
                            timeout=10)
        if gr.job_status in (pb.JOB_STATUS_SUCCEEDED, pb.JOB_STATUS_FAILED, pb.JOB_STATUS_CANCELLED):
            break
        assert time.time() - t0 < 30, "async job timeout"
    assert gr.job_status == pb.JOB_STATUS_SUCCEEDED,         f"expected SUCCEEDED, got {pb.JobStatus.Name(gr.job_status)}: {gr.status_message}"
    assert abs(gr.result.objective_value - 4.0) < 1e-4,         f"expected obj=4.0, got {gr.result.objective_value}"
    print(f"[ok] async job: {sub.job_id[:8]}... obj={gr.result.objective_value} "
          f"elapsed={gr.elapsed_time:.3f}s")


def case_async_cancel():
    """Async cancel: submit then cancel immediately, verify CANCELLED"""
    stub = get_stub()
    req = pb.SolveRequest(
        sense=pb.OBJ_SENSE_MIN,
        col_cost=[1.0], col_lower=[0.0], col_upper=[1e30],
        a_format_start=[0, 1], a_format_index=[0], a_format_value=[1.0],
        row_lower=[-1e30], row_upper=[1.0],
        options={"solver": "ipm"}, time_limit=60,
    )
    sub = stub.SubmitSolve(req, timeout=5)
    cancel = stub.CancelSolve(pb.CancelRequest(job_id=sub.job_id), timeout=5)
    # cancel.cancelled is True if job was PENDING/RUNNING; False if already
    # terminal. Both are valid for tiny LPs that finish in milliseconds.
    assert cancel.cancelled or "terminal" in cancel.message, \
        f"unexpected cancel response: cancelled={cancel.cancelled} msg={cancel.message}"
    gr = stub.GetResult(pb.GetResultRequest(job_id=sub.job_id, wait=True, wait_timeout=2.0),
                        timeout=10)
    assert gr.job_status in (pb.JOB_STATUS_CANCELLED, pb.JOB_STATUS_SUCCEEDED),         f"expected CANCELLED/SUCCEEDED, got {pb.JobStatus.Name(gr.job_status)}"
    print(f"[ok] async cancel: {sub.job_id[:8]}... → {pb.JobStatus.Name(gr.job_status)}")


def case_async_not_found():
    """Query non-existent job_id -> NOT_FOUND"""
    stub = get_stub()
    try:
        stub.GetResult(pb.GetResultRequest(job_id="nonexistent_job_id"), timeout=5)
        print("[fail] expected NOT_FOUND")
        sys.exit(1)
    except grpc.RpcError as e:
        assert e.code() == grpc.StatusCode.NOT_FOUND, f"expected NOT_FOUND, got {e.code()}"
        print(f"[ok] async not_found rejected: {e.details()}")


if __name__ == '__main__':
    case_health_check()
    case_feasible()
    case_infeasible()
    case_bad_args()
    case_async_job()
    case_async_cancel()
    case_async_not_found()
    print("ALL TESTS PASSED")
