#!/usr/bin/env python3
"""E2E 测试客户端：覆盖可行/不可行/非法参数，自适应探测 GPU 构建选 solver。"""
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
    """自适应：查服务端是否 GPU 构建，决定是否用 pdlp"""
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    return hc.gpu_available


def solve(stub, req, timeout=30):
    return stub.Solve(req, timeout=timeout)


def case_feasible():
    stub = get_stub()
    gpu = probe_gpu(stub)
    # GPU 构建用 pdlp（GPU 上快），CPU 构建用 ipm（CPU 上比 pdlp 快）
    solver = "pdlp" if gpu else "ipm"
    # Maximize x1 + x2 s.t. x1 + 2 x2 <= 4, x>=0  →  最优解 (4,0), obj=4
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
    # 用两个合法但矛盾的约束构造不可行 LP：
    #   row0: x >= 5    (row_lower=5,  row_upper=+inf)
    #   row1: x <= 1    (row_lower=-inf, row_upper=1)
    # 即 5 <= x <= 1，求解器判定 infeasible
    # 注意：HiGHS passModel 校验 lower<=upper，所以不能单行写 lower>upper
    req = pb.SolveRequest(
        sense=pb.OBJ_SENSE_MIN,
        col_cost=[1.0], col_lower=[0.0], col_upper=[1e30],
        # 2 行 1 列，CSC: 列0 有 2 个非零（行0和行1）
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
    # a_format_start 长度应为 num_col+1=3，这里给 4
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


if __name__ == '__main__':
    case_health_check()
    case_feasible()
    case_infeasible()
    case_bad_args()
    print("ALL TESTS PASSED")
