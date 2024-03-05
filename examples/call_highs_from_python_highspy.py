
import highspy
import numpy as np

# Highs h
h = highspy.Highs()

# Set up problem
inf = highspy.kHighsInf
h.addVars(2, np.array([-inf, -inf]), np.array([inf, inf]))
h.changeColsCost(2, np.array([0, 1]), np.array([0, 1], dtype=np.double));
num_cons = 2
lower = np.array([2, 0], dtype=np.double)
upper = np.array([inf, inf], dtype=np.double)
num_new_nz = 4
starts = np.array([0, 2])
indices = np.array([0, 1, 0, 1])
values = np.array([-1, 1, 1, 1], dtype=np.double)
h.addRows(num_cons, lower, upper, num_new_nz, starts, indices, values)

# Access LP
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns',
      lp.num_row_, ' rows and ', num_nz, ' nonzeros')

print('Solving...')
# Disable output from HiGHS for very small LP
h.setOptionValue("log_to_console", False)

# Solve problem
h.run()

print('Problem solved.')
print()

# Print solution information
# solution = h.getSolution()
# basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()

print('Model status = ', h.modelStatusToString(model_status))
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ',
      h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ',
      h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))

h.clear()
