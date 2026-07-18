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
sys.path.insert(0, _EXAMPLE_DIR)  # 导入同目录 pyomo_lp

import grpc
import solver_pb2 as pb
import solver_pb2_grpc as pbg
# 复用 pyomo_lp.py 的建模 + 矩阵提取逻辑（DRY）
from pyomo_lp import build_model, pyomo_model_to_request


def submit_and_monitor(stub, req, poll_interval=0.3):
    """
    异步提交 + 进度监控 + 取结果
    
    返回最终 GetResultResponse（终态）
    """
    # 1. 提交，立即返回 job_id
    t0 = time.time()
    sub = stub.SubmitSolve(req, timeout=5)
    job_id = sub.job_id
    print(f"[submit] job_id={job_id[:12]}... ({(time.time()-t0)*1000:.1f}ms)")

    # 2. 进度监控循环
    #    wait=false 立即返回当前状态 + 进度，便于展示进度条
    print(f"\n{'轮次':>4} {'状态':<20} {'迭代':>6} {'目标值':>14} {'耗时':>8}")
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

    # 0. 健康检查 + 自适应选 solver
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    solver = "pdlp" if hc.gpu_available else "ipm"
    print(f"server: {hc.message} → solver={solver}")
    print(f"  提示: 异步 job 模式适合大模型长求解，本例用小模型演示流程\n")

    # 1. Pyomo 建模（复用 pyomo_lp.build_model）
    model = build_model()
    print(f"[model] {model.name}")

    # 2. Pyomo → CSC 矩阵 → SolveRequest（复用 pyomo_lp.pyomo_model_to_request）
    req = pyomo_model_to_request(model)
    req.options["solver"] = solver

    # 3. 异步提交 + 进度监控
    print("\n=== 异步求解 ===")
    final = submit_and_monitor(stub, req)

    # 4. 输出结果
    print(f"\n=== 结果 ===")
    print(f"job_status: {pb.JobStatus.Name(final.job_status)}")
    print(f"elapsed: {final.elapsed_time:.3f}s")
    if final.job_status == pb.JOB_STATUS_SUCCEEDED:
        r = final.result
        print(f"model_status: {pb.ModelStatus.Name(r.model_status)}")
        print(f"objective: {r.objective_value}  (期望 38)")
        print(f"col_value: {list(r.col_value)}  (期望 [6.0, 4.0])")
        print(f"solve_time: {r.solve_time:.4f}s, iter: {r.iteration_count}")
        assert abs(r.objective_value - 38) < 1e-4, "obj 不符"
        print("\n✓ Pyomo + 异步 job 模式验证通过")
    else:
        print(f"msg: {final.status_message}")
        sys.exit(1)

    # 5. 演示：Pyomo 建不同模型批量提交
    print("\n=== 批量提交演示（3 个不同 sense 的变体）===")
    jobs = []
    for sense_name, sense in [("max", pb.OBJ_SENSE_MAX), ("min", pb.OBJ_SENSE_MIN)]:
        m = build_model()
        m.profit.sense = pyo_maximize if sense_name == "max" else pyo_minimize
        r = pyomo_model_to_request(m)
        r.options["solver"] = solver
        r.model_name = f"{m.name}_{sense_name}"
        sub = stub.SubmitSolve(r, timeout=5)
        jobs.append((sub.job_id, sense_name, r.model_name))
        print(f"  提交 {r.model_name}: job_id={sub.job_id[:12]}...")

    # 等所有完成
    print("\n  等待所有 job 完成...")
    for job_id, sense_name, name in jobs:
        gr = stub.GetResult(
            pb.GetResultRequest(job_id=job_id, wait=True, wait_timeout=10),
            timeout=15,
        )
        obj = gr.result.objective_value if gr.job_status == pb.JOB_STATUS_SUCCEEDED else None
        print(f"  {name}: {pb.JobStatus.Name(gr.job_status)} obj={obj}")


# pyomo sense 常量（build_model 用 maximize，这里便于切换演示）
import pyomo.environ as pyo
pyo_maximize = pyo.maximize
pyo_minimize = pyo.minimize


if __name__ == "__main__":
    main()
