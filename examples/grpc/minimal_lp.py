#!/usr/bin/env python3
"""
Minimal example: construct SolveRequest by hand, no modeling tool needed.

Suitable for: quick connectivity check, users not using Pyomo, embedded calls.

Model: Maximize x1 + x2  s.t. x1 + 2*x2 <= 4, x >= 0  ->  optimal (4,0), obj=4

Run:
  1. Start server: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
  2. python examples/grpc/minimal_lp.py
"""
import sys, os, subprocess

# --- 自动生成 gRPC Python 桩（与 pyomo_lp.py 相同逻辑，自包含）---
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


def main(addr="localhost:50051"):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)

    # 1. 健康检查 + 自适应选 solver
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    solver = "pdlp" if hc.gpu_available else "ipm"
    print(f"server: {hc.message} → 使用 solver={solver}")

    # 2. 手构造 LP 模型（CSC 稀疏格式）
    #    Maximize x1 + x2
    #    s.t.  x1 + 2*x2 <= 4
    #          x1, x2 >= 0
    req = pb.SolveRequest(
        model_name="minimal_lp",
        sense=pb.OBJ_SENSE_MAX,
        col_cost=[1.0, 1.0],          # 目标系数
        col_lower=[0.0, 0.0],         # 变量下界
        col_upper=[1e30, 1e30],       # 变量上界（无界）
        # CSC 约束矩阵：2 列 → start 长度 3
        a_format_start=[0, 1, 2],     # 列0 非零起始0, 列1 起始1, 结束2
        a_format_index=[0, 0],        # 列0 的非零在行0; 列1 的非零在行0
        a_format_value=[1.0, 2.0],    # 列0 系数1; 列1 系数2
        row_lower=[-1e30],            # 约束下界（-inf）
        row_upper=[4.0],              # 约束上界 4
        options={"solver": solver},
    )

    # 3. 求解
    resp = stub.Solve(req, timeout=60)
    print(f"status: {pb.ModelStatus.Name(resp.model_status)}")
    print(f"objective: {resp.objective_value}  (期望 4.0)")
    print(f"col_value: {list(resp.col_value)}  (期望 [4.0, 0.0])")
    print(f"col_dual:  {list(resp.col_dual)}")
    print(f"solve_time: {resp.solve_time:.4f}s, iter: {resp.iteration_count}")

    assert abs(resp.objective_value - 4.0) < 1e-4, "obj 不符"
    print("\n✓ 验证通过")


if __name__ == "__main__":
    main()
