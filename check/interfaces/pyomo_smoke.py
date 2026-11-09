# Solve a small LP and MIP through Pyomo's highspy-based appsi_highs solver
import pyomo.environ as pyo

m = pyo.ConcreteModel()
m.x = pyo.Var(bounds=(0, 4.5))
m.y = pyo.Var(bounds=(1, 1))
m.z = pyo.Var(within=pyo.NonNegativeIntegers)
m.obj = pyo.Objective(expr=m.x + 4 * m.y + 9 * m.z, sense=pyo.minimize)
m.c1 = pyo.Constraint(expr=m.x + m.y <= 5.5)
m.c2 = pyo.Constraint(expr=m.x + m.z >= 10)
m.c3 = pyo.Constraint(expr=-m.y + m.z <= 7.5)

solver = pyo.SolverFactory("appsi_highs")
assert solver.available(), "appsi_highs (highspy) is not available"

# LP relaxation, then the MIP
for relax, expected in ((True, 58.0), (False, 62.0)):
    model = pyo.TransformationFactory("core.relax_integer_vars").create_using(m) if relax else m
    res = solver.solve(model)
    print(f"relax={relax}: {res.solver.termination_condition} obj={pyo.value(model.obj)}")
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert abs(pyo.value(model.obj) - expected) < 1e-6
