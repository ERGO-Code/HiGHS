
# The path to the MPS file instance assumes that this is run 
# from the root directory.

import highspy
import numpy as np

# Highs h
h = highspy.Highs()

# Solve from mps file
# Initialize an instance of Highs
# h = highspy.Highs()
# Here we are re-using the one from above.
h.readModel('check/instances/avgas.mps')

# Print
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns',
      lp.num_row_, ' rows and ', num_nz, ' nonzeros')

# Solve the problem
h.run()

# Print solution information
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print()
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ',
      h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ',
      h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))

# basis.col_status is already a list, but accessing values in
# solution.col_value directly is very inefficient, so convert it to a
# list
col_status = basis.col_status
row_status = basis.row_status
col_value = list(solution.col_value)
row_value = list(solution.row_value)

num_var = h.getNumCol()
num_row = h.getNumRow()
print("Variables")
for icol in range(num_var):
    print(icol, col_value[icol],
          h.basisStatusToString(col_status[icol]))
print("Constraints")
for irow in range(num_row):
    print(irow, row_value[irow],
          h.basisStatusToString(row_status[irow]))

