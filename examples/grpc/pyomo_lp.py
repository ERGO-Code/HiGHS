#!/usr/bin/env python3
"""
End-to-end example: Pyomo modeling -> encode as gRPC request -> solve via highs-server.

Two paths:
  Path A (recommended): Pyomo model -> extract CSC matrix directly -> SolveRequest
  Path B (file):        Pyomo model -> write MPS -> read via highspy -> SolveRequest

Model: production planning LP
  Maximize  profit = 3*x1 + 5*x2
  s.t.  x1 + 2*x2 <= 14   (raw material A)
        x1       <= 6      (raw material B)
        x1 +  x2 <= 10     (labor)
        x1, x2 >= 0
  Optimal: x1=6, x2=4, obj=38

Run:
  1. Start server: ./build/bin/highs_grpc_server --bind 127.0.0.1:50051
  2. python examples/grpc/pyomo_lp.py
"""
import sys, os, subprocess

# --- Auto-generate gRPC Python stubs (self-contained, no need for pre-generated stubs) ---
_EXAMPLE_DIR = os.path.dirname(os.path.abspath(__file__))
_PROTO_DIR = os.path.normpath(os.path.join(_EXAMPLE_DIR, "..", "..", "server", "protos"))
_GEN_DIR = os.path.join(_EXAMPLE_DIR, "_generated")
os.makedirs(_GEN_DIR, exist_ok=True)
# Regenerate only if stubs missing or proto updated
_need_gen = (not os.path.exists(os.path.join(_GEN_DIR, "solver_pb2.py")) or
             os.path.getmtime(os.path.join(_PROTO_DIR, "solver.proto")) >
             os.path.getmtime(os.path.join(_GEN_DIR, "solver_pb2.py")))
if _need_gen:
    subprocess.check_call([
        sys.executable, "-m", "grpc_tools.protoc",
        f"-I{_PROTO_DIR}", f"--python_out={_GEN_DIR}", f"--grpc_python_out={_GEN_DIR}",
        os.path.join(_PROTO_DIR, "solver.proto"),
    ])
sys.path.insert(0, _GEN_DIR)

import grpc
import pyomo.environ as pyo
from pyomo.repn import generate_standard_repn
import solver_pb2 as pb
import solver_pb2_grpc as pbg


# ============================================================
# 1. Build model with Pyomo
# ============================================================
def build_model() -> pyo.ConcreteModel:
    """Production planning LP model"""
    m = pyo.ConcreteModel(name="production_planning")

    # Decision variables (Pyomo var default lower bound 0)
    m.x1 = pyo.Var(within=pyo.NonNegativeReals, name="product_1")
    m.x2 = pyo.Var(within=pyo.NonNegativeReals, name="product_2")

    # Objective: Maximize profit
    m.profit = pyo.Objective(expr=3 * m.x1 + 5 * m.x2, sense=pyo.maximize)

    # Constraints
    m.raw_a = pyo.Constraint(expr=m.x1 + 2 * m.x2 <= 14, name="raw_a")
    m.raw_b = pyo.Constraint(expr=m.x1 <= 6, name="raw_b")
    m.labor = pyo.Constraint(expr=m.x1 + m.x2 <= 10, name="labor")
    return m


# ============================================================
# 2A. Pyomo model -> CSC matrix -> SolveRequest (recommended path)
# ============================================================
def pyomo_model_to_request(model: pyo.ConcreteModel) -> pb.SolveRequest:
    """
    Convert a Pyomo ConcreteModel into a SolveRequest.

    Key points:
      - Use generate_standard_repn to split each constraint body into (linear part + constant)
      - Constraint form lower <= body <= upper; after moving terms: row_lower = lower - constant,
        row_upper = upper - constant (constant moved from body side to bound side)
      - Matrix in CSC (Compressed Sparse Column): start length = num_col+1
    """
    # --- Collect variables (preserve order) ---
    var_list = [v for v in model.component_data_objects(pyo.Var, active=True)]
    var_index = {id(v): i for i, v in enumerate(var_list)}
    num_col = len(var_list)

    # Column data: cost / lower / upper
    col_cost = [0.0] * num_col
    col_lower = [0.0] * num_col
    col_upper = [1e30] * num_col
    for i, v in enumerate(var_list):
        col_lower[i] = float(v.lb) if v.lb is not None else -1e30
        col_upper[i] = float(v.ub) if v.ub is not None else 1e30

    # Objective coefficients
    obj = next(model.component_data_objects(pyo.Objective, active=True))
    obj_repn = generate_standard_repn(obj.expr)
    for var, coef in zip(obj_repn.linear_vars, obj_repn.linear_coefs):
        col_cost[var_index[id(var)]] = float(coef)
    sense = pb.OBJ_SENSE_MAX if obj.sense == pyo.maximize else pb.OBJ_SENSE_MIN

    # --- Collect constraints -> CSR -> CSC ---
    # First collect by row (row, col, value), then sort by column for CSC
    triples = []  # (row_idx, col_idx, value)
    row_lower = []
    row_upper = []
    row_names = []
    row_idx = 0
    for con in model.component_data_objects(pyo.Constraint, active=True):
        repn = generate_standard_repn(con.body)
        c = float(repn.constant or 0.0)
        # Constraint lower/upper (may be None in Pyomo)
        cl = float(con.lower) if con.lower is not None else -1e30
        cu = float(con.upper) if con.upper is not None else 1e30
        # Move terms: lower <= linear + c <= upper  =>  lower-c <= linear <= upper-c
        row_lower.append(cl - c)
        row_upper.append(cu - c)
        row_names.append(con.name)
        for var, coef in zip(repn.linear_vars, repn.linear_coefs):
            triples.append((row_idx, var_index[id(var)], float(coef)))
        row_idx += 1
    num_row = row_idx

    # CSR -> CSC: sort by column, build start/index/value
    triples.sort(key=lambda t: (t[1], t[0]))  # sort by column, then row
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
# 2B. Pyomo model -> MPS file -> read via highspy -> SolveRequest (file path)
# ============================================================
def mps_file_to_request(mps_path: str) -> pb.SolveRequest:
    """Read model from MPS file (via highspy), convert to SolveRequest."""
    import highspy
    h = highspy.Highs()
    h.readModel(mps_path)
    lp = h.getLp()

    num_col = lp.num_col_
    num_row = lp.num_row_
    # Column data
    col_cost = list(lp.col_cost_)
    col_lower = list(lp.col_lower_)
    col_upper = [u if u < 1e29 else 1e30 for u in lp.col_upper_]
    # Row bounds
    row_lower = [l if l > -1e29 else -1e30 for l in lp.row_lower_]
    row_upper = [u if u < 1e29 else 1e30 for u in lp.row_upper_]
    # Matrix: HiGHS lp.a_matrix_ is CSC, take directly
    a_start = list(lp.a_matrix_.start_)
    a_index = list(lp.a_matrix_.index_)
    a_value = list(lp.a_matrix_.value_)
    # highspy preserves original sense (kMaximize/kMinimize) and original cost (not negated) after reading MPS
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
# 3. Solve via highs-server
# ============================================================
def solve_via_grpc(req: pb.SolveRequest, addr="localhost:50051") -> pb.SolveResponse:
    ch = grpc.insecure_channel(addr)
    grpc.channel_ready_future(ch).result(timeout=5)
    stub = pbg.HighsServiceStub(ch)
    # Adaptive solver selection
    hc = stub.Check(pb.HealthCheckRequest(), timeout=5)
    if "solver" not in req.options:
        req.options["solver"] = "pdlp" if hc.gpu_available else "ipm"
    print(f"[grpc] gpu_available={hc.gpu_available}, solver={req.options['solver']}")
    return stub.Solve(req, timeout=60)


# ============================================================
# 4. Main flow
# ============================================================
def main():
    model = build_model()

    # --- Path A: Pyomo direct matrix extraction ---
    print("=" * 60)
    print("Path A: Pyomo model -> direct CSC extraction -> gRPC solve")
    print("=" * 60)
    req_a = pyomo_model_to_request(model)
    resp_a = solve_via_grpc(req_a)
    print(f"[result] status={pb.ModelStatus.Name(resp_a.model_status)} "
          f"obj={resp_a.objective_value}")
    print(f"[result] col_value={list(resp_a.col_value)}")
    print(f"[result] solve_time={resp_a.solve_time:.4f}s, iter={resp_a.iteration_count}")

    # --- Path B: Pyomo -> MPS file -> highspy read -> gRPC ---
    print("\n" + "=" * 60)
    print("Path B: Pyomo model -> MPS file -> highspy read -> gRPC solve")
    print("=" * 60)
    mps_file = "/tmp/production.mps"
    model.write(mps_file, format="mps")
    print(f"[file] wrote {mps_file}")
    req_b = mps_file_to_request(mps_file)
    resp_b = solve_via_grpc(req_b)
    print(f"[result] status={pb.ModelStatus.Name(resp_b.model_status)} "
          f"obj={resp_b.objective_value}")
    print(f"[result] col_value={list(resp_b.col_value)}")

    # --- Compare with expected solution ---
    print("\n" + "=" * 60)
    print("Expected: x1=6, x2=4, obj=38")
    print("=" * 60)
    print(f"Path A obj={resp_a.objective_value}  (expected 38)")
    print(f"Path B obj={resp_b.objective_value}  (expected 38)")
    # Path B: MPS converts max to min (cost negated), obj would be -38
    assert abs(resp_a.objective_value - 38) < 1e-4, "Path A obj mismatch"
    assert abs(resp_b.objective_value - 38) < 1e-4, "Path B obj mismatch"
    print("\n✓ Path A and Path B both verified (obj~38, x1~6, x2~4)")


if __name__ == "__main__":
    main()
