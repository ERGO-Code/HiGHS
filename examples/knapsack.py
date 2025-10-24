import numpy as np
import highspy

h = highspy.Highs()
h.setOptionValue("output_flag", False);
h.setOptionValue("presolve", "off");
inf = highspy.kHighsInf
lp = highspy.HighsLp()
lp.sense_ = highspy.ObjSense.kMaximize
lp.num_col_ = 5
lp.num_row_ = 1
lp.col_cost_ = np.array([8, 5, 3, 11, 7], dtype=np.double)
lp.col_lower_ = np.array([0, 0, 0, 0, 0], dtype=np.double)
lp.col_upper_ = np.array([1, 1, 1, 1, 1], dtype=np.double)
lp.row_lower_ = np.array([-inf], dtype=np.double)
lp.row_upper_ = np.array([11], dtype=np.double)
lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
lp.a_matrix_.start_ = np.array([0, 5])
lp.a_matrix_.index_ = np.array([0, 1, 2, 3, 4])
lp.a_matrix_.value_ = np.array([4, 3, 1, 5, 4], dtype=np.double)
lp.integrality_ = np.array([highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger])
h.passModel(lp)

h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}, {solution.col_value[2]}, {solution.col_value[3]}, {solution.col_value[4]})")

# Solution is [1, 0, 1, 1, 0]

# Illustrate setSolution

# First by passing back the optimal solution

h.clearSolver()
h.setSolution(solution)
h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}, {solution.col_value[2]}, {solution.col_value[3]}, {solution.col_value[4]})")

# Now passing back the optimal values of two variables as a sparse solution
h.clearSolver()
index = np.array([0, 3])
value = np.array([1, 1], dtype=np.double)
h.setSolution(2, index, value)
h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}, {solution.col_value[2]}, {solution.col_value[3]}, {solution.col_value[4]})")

# Test passing back the optimal value of one variable, and a non-optimal value of another, as a sparse solution in untyped array
h.clearSolver()
h.setSolution(2, np.array([0, 4]), np.array([1, 1]))
h.run()
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]}, {solution.col_value[2]}, {solution.col_value[3]}, {solution.col_value[4]})")

