# Simple knapsack example to illustrate setSolution
# max 8x1 + 5x2 + 3x3 + 11x4 + 7x5
#     4x1 + 3x2 + 1x3 +  5x4 + 4x5 <= 11

import numpy as np
import highspy

def print_solution(h, x):
    print(f"Solution: ({', '.join(map(str, abs(h.val(x))))})")

h = highspy.Highs()
h.silent()
h.setOptionValue("presolve", "off")

x = h.addBinaries(5)
h.addConstr((x * [4, 3, 1, 5, 4]).sum() <= 11)
h.maximize((x * [8, 5, 3, 11, 7]).sum())
    
# solution is [1, 0, 1, 1, 0]
print_solution(h, x)


# Illustrate setSolution
# First by passing back the optimal solution
h.clearSolver()
h.setSolution(h.getSolution())
h.run()
print_solution(h, x)

# Now passing back the optimal values of two variables as a sparse solution
h.clearSolver()
index = np.array([0, 3])
value = np.array([1, 1], dtype=np.float64)
h.setSolution(2, index, value)
h.run()
print_solution(h, x)

# Test passing back the optimal value of one variable, and a non-optimal value of another, as a sparse solution in untyped array
h.clearSolver()
h.setSolution(2, np.array([0, 4]), np.array([1, 1]))
h.run()
print_solution(h, x)

