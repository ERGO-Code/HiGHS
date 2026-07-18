#!/usr/bin/env python3
"""
异步 job 模式示例：SubmitSolve → 轮询 GetResult → 取结果

适合：大模型、长求解（>10s），避免长连接占用、支持取消、断线可恢复。

模型：生产计划 LP（与 pyomo_lp.py 相同），最优解 obj=38

流程：
  1. SubmitSolve(req) → 立即返回 job_id（连接 ~1ms 释放）
  2. GetResult(job_id, wait=true, wait_timeout=2) → 短轮询等终态
  3. 终态后取 result

运行：
  1. 启动服务: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
  2. python examples/grpc/async_job.py
"""
import sys, os, subprocess, time

# --- 自动生成 gRPC Python 桩（自包含）---
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
    """生产计划 LP: Max 3x1+5x2 s.t. 3 约束，最优解 x1=6,x2=4,obj=38"""
    return pb.SolveRequest(
        model_name="production_planning",
        sense=pb.OBJ_SENSE_MAX,
        col_cost=[3.0, 5.0],
        col_lower=[0.0, 0.0], col_upper=[1e30, 1e30],
        # 3 行 2 列 CSC
        a_format_start=[0, 3, 5],       # 列0 有3个非零, 列1 有2个
        a_format_index=[0, 1, 2, 0, 2], # 列0: 行0,1,2; 列1: 行0,2
        a_format_value=[1.0, 1.0, 1.0, 2.0, 1.0],
        row_lower=[-1e30, -1e30, -1e30],
        row_upper=[14.0, 6.0, 10.0],
    )


def main(addr="localhost:50051"):
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)

    # 0. 健康检查 + 自适应选 solver
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    solver = "pdlp" if hc.gpu_available else "ipm"
    print(f"server: {hc.message} → solver={solver}")

    req = build_production_lp()
    req.options["solver"] = solver

    # 1. 提交任务（立即返回 job_id）
    t0 = time.time()
    submit_resp = stub.SubmitSolve(req, timeout=5)
    job_id = submit_resp.job_id
    print(f"[submit] job_id={job_id} status={pb.JobStatus.Name(submit_resp.job_status)} "
          f"({(time.time()-t0)*1000:.1f}ms)")

    # 2. 轮询/长轮询等结果
    #    wait=true + wait_timeout=2: 服务端阻塞最多 2s，终态则立即返回
    #    客户端循环直到终态（比纯客户端 sleep 轮询更高效）
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

    # 3. 输出结果
    if final.job_status == pb.JOB_STATUS_SUCCEEDED:
        r = final.result
        print(f"\n[result] model_status={pb.ModelStatus.Name(r.model_status)}")
        print(f"[result] objective={r.objective_value}  (期望 38)")
        print(f"[result] col_value={list(r.col_value)}  (期望 [6.0, 4.0])")
        print(f"[result] solve_time={r.solve_time:.4f}s iter={r.iteration_count}")
        assert abs(r.objective_value - 38) < 1e-4, "obj 不符"
        print("\n✓ 异步 job 模式验证通过")
    else:
        print(f"\n[!] job 终态非 SUCCEEDED: {pb.JobStatus.Name(final.job_status)} "
              f"msg={final.status_message}")
        sys.exit(1)

    # 4. 演示取消（提交后立即取消）
    print("\n--- 演示取消 ---")
    req2 = build_production_lp()
    req2.options["solver"] = solver
    req2.time_limit = 60  # 给个长时限，便于取消演示
    sub2 = stub.SubmitSolve(req2, timeout=5)
    print(f"[submit] job_id={sub2.job_id}")
    cancel = stub.CancelSolve(pb.CancelRequest(job_id=sub2.job_id), timeout=5)
    print(f"[cancel] cancelled={cancel.cancelled} msg={cancel.message}")
    # 查最终状态
    gr2 = stub.GetResult(pb.GetResultRequest(job_id=sub2.job_id, wait=True, wait_timeout=2.0),
                         timeout=10)
    print(f"[final] status={pb.JobStatus.Name(gr2.job_status)} msg={gr2.status_message}")


if __name__ == "__main__":
    main()
