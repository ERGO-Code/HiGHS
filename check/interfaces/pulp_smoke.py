# Solve a small LP and MIP through PuLP's highspy-based HiGHS solver
import pulp

prob = pulp.LpProblem("smoke", pulp.LpMinimize)
x = prob.add_variable("x", 0, 4.5)
y = prob.add_variable("y", 1, 1)
z = prob.add_variable("z", 0, None, pulp.LpInteger)
prob += x + 4 * y + 9 * z
prob += x + y <= 5.5
prob += x + z >= 10
prob += -y + z <= 7.5

# LP relaxation, then the MIP
for mip, expected in ((False, 58.0), (True, 62.0)):
    stats = prob.solve(pulp.HiGHS(msg=False, mip=mip))
    print(f"mip={mip}: {stats.status_str} obj={pulp.value(prob.objective)}")
    assert stats.status_str == "Optimal"
    assert abs(pulp.value(prob.objective) - expected) < 1e-6
