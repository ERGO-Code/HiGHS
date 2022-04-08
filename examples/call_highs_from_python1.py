import highspy
import numpy as np
inf = highspy.kHighsInf
h = highspy.Highs()
h.addVar(-inf, inf)
h.addVar(-inf, inf)
h.changeColsCost(2, np.array([0, 1]), np.array([0, 1], dtype=np.double))
num_cons = 2
lower = np.array([2, 0], dtype=np.double)
upper = np.array([inf, inf], dtype=np.double)
num_new_nz = 4
starts = np.array([0, 2])
indices = np.array([0, 1, 0, 1])
values = np.array([-1, 1, 1, 1], dtype=np.double)
h.addRows(num_cons, lower, upper, num_new_nz, starts, indices, values)
h.setOptionValue('log_to_console', True)
h.run()
num_var = h.getNumCol()
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
#model_status = h.getModelStatus()
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
for icol in range(num_var):
    print(icol, solution.col_value[icol])
h.readModel("check/instances/avgas.mps")
h.writeModel("ml.mps")
h.run()
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns', lp.num_row_, ' rows and ', num_nz, ' nonzeros')
lp.num_col_ = 99
num_col = h.getNumCol()
print('LP has ', num_col, ' columns')
