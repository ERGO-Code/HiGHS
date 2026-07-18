#!/usr/bin/env python3
"""
端到端示例：用 Pyomo 建模 → 编码为 gRPC 请求 → 请求 highs-server 求解

演示两条路径：
  路径 A（推荐）：Pyomo 模型 → 直接提取 CSC 矩阵 → SolveRequest
  路径 B（文件）：Pyomo 模型 → 写 MPS 文件 → highspy 读取 → SolveRequest

模型：生产计划 LP
  Maximize  利润 = 3*x1 + 5*x2
  s.t.  x1 + 2*x2 <= 14   (原料 A)
        x1       <= 6      (原料 B)
        x1 +  x2 <= 10     (工时)
        x1, x2 >= 0
  最优解：x1=6, x2=4, obj=38

运行：
  1. 先启动服务: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
  2. python test/example_pyomo.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import grpc
import pyomo.environ as pyo
from pyomo.repn import generate_standard_repn
import solver_pb2 as pb
import solver_pb2_grpc as pbg


# ============================================================
# 1. 用 Pyomo 建模
# ============================================================
def build_model() -> pyo.ConcreteModel:
    """生产计划 LP 模型"""
    m = pyo.ConcreteModel(name="production_planning")

    # 决策变量（Pyomo 变量默认下界 0）
    m.x1 = pyo.Var(within=pyo.NonNegativeReals, name="product_1")
    m.x2 = pyo.Var(within=pyo.NonNegativeReals, name="product_2")

    # 目标：Maximize 利润
    m.profit = pyo.Objective(expr=3 * m.x1 + 5 * m.x2, sense=pyo.maximize)

    # 约束
    m.raw_a = pyo.Constraint(expr=m.x1 + 2 * m.x2 <= 14, name="raw_a")
    m.raw_b = pyo.Constraint(expr=m.x1 <= 6, name="raw_b")
    m.labor = pyo.Constraint(expr=m.x1 + m.x2 <= 10, name="labor")
    return m


# ============================================================
# 2A. Pyomo 模型 → CSC 矩阵 → SolveRequest（推荐路径）
# ============================================================
def pyomo_model_to_request(model: pyo.ConcreteModel) -> pb.SolveRequest:
    """
    把 Pyomo ConcreteModel 转换成 SolveRequest。

    关键点：
      - 用 generate_standard_repn 把每个约束的 body 拆成 (线性部分 + 常数)
      - 约束形式 lower <= body <= upper，移项后 row_lower = lower - constant,
        row_upper = upper - constant（constant 从 body 侧移到 bound 侧）
      - 矩阵按列压缩 (CSC)：start 长度 = num_col+1
    """
    # --- 收集变量（保持顺序） ---
    var_list = [v for v in model.component_data_objects(pyo.Var, active=True)]
    var_index = {id(v): i for i, v in enumerate(var_list)}
    num_col = len(var_list)

    # 列信息：cost / lower / upper
    col_cost = [0.0] * num_col
    col_lower = [0.0] * num_col
    col_upper = [1e30] * num_col
    for i, v in enumerate(var_list):
        col_lower[i] = float(v.lb) if v.lb is not None else -1e30
        col_upper[i] = float(v.ub) if v.ub is not None else 1e30

    # 目标函数系数
    obj = next(model.component_data_objects(pyo.Objective, active=True))
    obj_repn = generate_standard_repn(obj.expr)
    for var, coef in zip(obj_repn.linear_vars, obj_repn.linear_coefs):
        col_cost[var_index[id(var)]] = float(coef)
    sense = pb.OBJ_SENSE_MAX if obj.sense == pyo.maximize else pb.OBJ_SENSE_MIN

    # --- 收集约束 → CSR → CSC ---
    # 先按行收集 (row, col, value)，再按列排序转 CSC
    triples = []  # (row_idx, col_idx, value)
    row_lower = []
    row_upper = []
    row_names = []
    row_idx = 0
    for con in model.component_data_objects(pyo.Constraint, active=True):
        repn = generate_standard_repn(con.body)
        c = float(repn.constant or 0.0)
        # 约束的 lower/upper（Pyomo 里 lower/upper 可能为 None）
        cl = float(con.lower) if con.lower is not None else -1e30
        cu = float(con.upper) if con.upper is not None else 1e30
        # 移项：lower <= linear + c <= upper  =>  lower-c <= linear <= upper-c
        row_lower.append(cl - c)
        row_upper.append(cu - c)
        row_names.append(con.name)
        for var, coef in zip(repn.linear_vars, repn.linear_coefs):
            triples.append((row_idx, var_index[id(var)], float(coef)))
        row_idx += 1
    num_row = row_idx

    # CSR → CSC：按列排序，构造 start/index/value
    triples.sort(key=lambda t: (t[1], t[0]))  # 按列、再按行排序
    a_start = [0] * (num_col + 1)
    a_index = []
    a_value = []
    col_counts = [0] * num_col
    for (_, col, _) in triples:
        col_counts[col] += 1
    for col in range(num_col):
        a_start[col + 1] = a_start[col] + col_counts[col]
    for (row, col, val) in triples:
        a_index.append(row)
        a_value.append(val)

    req = pb.SolveRequest(
        model_name=model.name,
        sense=sense,
        col_cost=col_cost,
        col_lower=col_lower,
        col_upper=col_upper,
        a_format_start=a_start,
        a_format_index=a_index,
        a_format_value=a_value,
        row_lower=row_lower,
        row_upper=row_upper,
    )
    print(f"[model] {num_col} cols, {num_row} rows, {len(triples)} nnz")
    print(f"[model] vars: {[v.name for v in var_list]}")
    print(f"[model] cons: {row_names}")
    return req


# ============================================================
# 2B. Pyomo 模型 → MPS 文件 → highspy 读取 → SolveRequest（文件路径）
# ============================================================
def mps_file_to_request(mps_path: str) -> pb.SolveRequest:
    """从 MPS 文件读取模型（用 highspy），转成 SolveRequest。"""
    import highspy
    h = highspy.Highs()
    h.readModel(mps_path)
    lp = h.getLp()

    num_col = lp.num_col_
    num_row = lp.num_row_
    # 列信息
    col_cost = list(lp.col_cost_)
    col_lower = list(lp.col_lower_)
    col_upper = [u if u < 1e29 else 1e30 for u in lp.col_upper_]
    # 行界
    row_lower = [l if l > -1e29 else -1e30 for l in lp.row_lower_]
    row_upper = [u if u < 1e29 else 1e30 for u in lp.row_upper_]
    # 矩阵：HiGHS lp.a_matrix_ 是 CSC，直接取
    a_start = list(lp.a_matrix_.start_)
    a_index = list(lp.a_matrix_.index_)
    a_value = list(lp.a_matrix_.value_)
    # highspy 读 MPS 后保留原始 sense（kMaximize/kMinimize）和原始 cost（不取负）
    from highspy.highs import ObjSense
    sense = pb.OBJ_SENSE_MAX if lp.sense_ == ObjSense.kMaximize else pb.OBJ_SENSE_MIN

    req = pb.SolveRequest(
        model_name=os.path.basename(mps_path),
        sense=sense,
        col_cost=col_cost,
        col_lower=col_lower,
        col_upper=col_upper,
        a_format_start=a_start,
        a_format_index=a_index,
        a_format_value=a_value,
        row_lower=row_lower,
        row_upper=row_upper,
    )
    print(f"[mps] {num_col} cols, {num_row} rows, {len(a_index)} nnz (from {mps_path})")
    return req


# ============================================================
# 3. 调用 highs-server 求解
# ============================================================
def solve_via_grpc(req: pb.SolveRequest, addr="localhost:50051") -> pb.SolveResponse:
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)
    # 自适应选 solver
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    if "solver" not in req.options:
        req.options["solver"] = "pdlp" if hc.gpu_available else "ipm"
    print(f"[grpc] gpu_available={hc.gpu_available}, solver={req.options['solver']}")
    return stub.Solve(req, timeout=60)


# ============================================================
# 4. 主流程
# ============================================================
def main():
    model = build_model()

    # --- 路径 A：Pyomo 直接提取矩阵 ---
    print("=" * 60)
    print("路径 A: Pyomo 模型 → 直接提取 CSC → gRPC 求解")
    print("=" * 60)
    req_a = pyomo_model_to_request(model)
    resp_a = solve_via_grpc(req_a)
    print(f"[result] status={pb.ModelStatus.Name(resp_a.model_status)} "
          f"obj={resp_a.objective_value}")
    print(f"[result] col_value={list(resp_a.col_value)}")
    print(f"[result] solve_time={resp_a.solve_time:.4f}s, iter={resp_a.iteration_count}")

    # --- 路径 B：Pyomo → MPS 文件 → highspy 读 → gRPC ---
    print("\n" + "=" * 60)
    print("路径 B: Pyomo 模型 → MPS 文件 → highspy 读 → gRPC 求解")
    print("=" * 60)
    mps_file = "/tmp/production.mps"
    model.write(mps_file, format="mps")
    print(f"[file] wrote {mps_file}")
    req_b = mps_file_to_request(mps_file)
    resp_b = solve_via_grpc(req_b)
    print(f"[result] status={pb.ModelStatus.Name(resp_b.model_status)} "
          f"obj={resp_b.objective_value}")
    print(f"[result] col_value={list(resp_b.col_value)}")

    # --- 对比期望解 ---
    print("\n" + "=" * 60)
    print("期望解: x1=6, x2=4, obj=38")
    print("=" * 60)
    print(f"路径 A obj={resp_a.objective_value}  (期望 38)")
    print(f"路径 B obj={resp_b.objective_value}  (期望 38)")
    # 路径 B 因 MPS 把 max 转 min（cost 取负），obj 为 -38
    assert abs(resp_a.objective_value - 38) < 1e-4, "路径 A obj 不符"
    assert abs(resp_b.objective_value - 38) < 1e-4, "路径 B obj 不符"
    print("\n✓ 路径 A 和路径 B 均验证通过（obj≈38, x1≈6, x2≈4）")


if __name__ == "__main__":
    main()
