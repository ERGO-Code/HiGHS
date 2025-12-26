import numpy as np
import highspy

hscb = highspy.cb 

h = highspy.Highs()
h.setOptionValue("output_flag", False);

inf = highspy.kHighsInf
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 3
lp.col_cost_ = np.array([0, 0], dtype=np.double)
lp.col_lower_ = np.array([0, 0], dtype=np.double)
lp.col_upper_ = np.array([inf, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, -inf, -inf], dtype=np.double)
lp.row_upper_ = np.array([18, 8, 14], dtype=np.double)
lp.a_matrix_.start_ = np.array([0, 3, 6])
lp.a_matrix_.index_ = np.array([0, 1, 2, 0, 1, 2])
lp.a_matrix_.value_ = np.array([3, 1, 1, 1, 1, 2], dtype=np.double)
h.passModel(lp)
# LP constraints are
#
# 3x +  y <= 18
#
#  x +  y <= 8
#
#  x + 2y <= 14
linear_objective0 = highspy.HighsLinearObjective()
linear_objective0.weight = -1
linear_objective0.offset = -1
linear_objective0.coefficients = np.array([1, 1], dtype=np.double)
linear_objective0.abs_tolerance = 0
linear_objective0.rel_tolerance = 0
linear_objective0.priority = 10

linear_objective1 = highspy.HighsLinearObjective()
linear_objective1.weight = 1e-4
linear_objective1.offset = 0
linear_objective1.coefficients = np.array([1, 0], dtype=np.double)
linear_objective1.abs_tolerance = -1
linear_objective1.rel_tolerance = -1
linear_objective1.priority = 0

h.addLinearObjective(linear_objective0)
h.addLinearObjective(linear_objective1)
# Objectives
#
# f0: x + y - 1 (weight -1)
#
# f1: x         (weight 1e-4)
#
# With objective min -f0 (since its weight is -1) the LP has nonunique
# optimal solutions on the line joining (2, 6) and (5, 3)
#
# Blending -f0 + 1e-4 f1 minimizes x along the line joining (2, 6) and
# (5, 3) to give optimal solution at (2, 6) with objective -7
h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}) for min -f0 + 1e-4 f1")

# Switch to lexicographic optimization 
h.setOptionValue("blend_multi_objectives", False);

linear_objective0.coefficients = np.array([1.0001, 1], dtype=np.double)
linear_objective0.abs_tolerance = 1e-5;
linear_objective0.rel_tolerance = 0.05;
linear_objective1.weight = 1e-3

h.clearLinearObjectives()
h.addLinearObjective(linear_objective0)
h.addLinearObjective(linear_objective1)
# Objectives
#
# f0 = 1.0001 x + y - 1 (with priority 10 and absolute tolerance 1e-5)
#
# f1 =        x         (with priority  0)
#
# Lexicographically: HiGHS 
#
# minimizes f0 (to give objective 7.0005) 
#
# adds a constraint that
#
# 1.0001 x + y - 1 >= 7.0005 - 1e-5 => 1.0001 x + y >= 8.00049
#
# to ensure that the initial objective is within 1e-5 of its optimal
# value
#
# minimizes f1 to give optimal solution at (4.90000, 3.10000) - where
# x + y = 8 and 1.0001 x + y = 8.00049

h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}) for min f1 and 1.0001 x + y >= 8.00049")

linear_objective0.abs_tolerance = -1;
h.clearLinearObjectives()
h.addLinearObjective(linear_objective0)
h.addLinearObjective(linear_objective1)
# Objectives as before, but absolute tolerence for f0 now -1 (negative
# => ignored) so relative tolerance of 0.05 is used
#
# Lexicographically: HiGHS 
#
# minimizes f0 (to give objective 7.0005) 
#
# adds a constraint that
#
# 1.0001 x + y - 1 >= 7.0005 * 0.95 => 1.0001 x + y >= 7.650475
#
# to ensure that the initial objective 95% of its optimal value
#
# minimizes f1 to give optimal solution at (1.30069, 6.34966) - where
# x + y = 8 and 1.0001 x + y = 7.650475
h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}) for min f1 and 1.0001 x + y = 7.650475")


