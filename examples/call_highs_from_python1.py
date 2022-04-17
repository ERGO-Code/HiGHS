import highspy
import numpy as np
inf = highspy.kHighsInf

h = highspy.Highs()
alt_inf = h.getInfinity()
print('highspy.kHighsInf = ', inf, '; h.getInfinity() = ', alt_inf)

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
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))
for icol in range(num_var):
    print(icol, solution.col_value[icol], h.basisStatusToString(basis.col_status[icol]))

# Read in and solve avgas
h.readModel("check/instances/avgas.mps")
#h.writeModel("ml.mps")
h.run()
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns', lp.num_row_, ' rows and ', num_nz, ' nonzeros')

# Clear so that incumbent model is empty
h.clear()

# Now define the blending model as a HighsLp instance
#
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 2
lp.sense_ = highspy.ObjSense.kMaximize
lp.col_cost_ = np.array([8, 10], dtype=np.double)
lp.col_lower_ = np.array([0, 0], dtype=np.double)
lp.col_upper_ = np.array([inf, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, -inf], dtype=np.double)
lp.row_upper_ = np.array([120, 210], dtype=np.double)
lp.a_matrix_.start_ = np.array([0, 2, 4])
lp.a_matrix_.index_ = np.array([0, 1, 0, 1])
lp.a_matrix_.value_ = np.array([0.3, 0.7, 0.5, 0.5], dtype=np.double)
h.passModel(lp)
h.run()
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))
num_var = h.getNumCol()
num_row = h.getNumRow()
print('Variables')
for icol in range(num_var):
    print(icol, solution.col_value[icol], h.basisStatusToString(basis.col_status[icol]))
print('Constraints')
for irow in range(num_row):
    print(irow, solution.row_value[irow], h.basisStatusToString(basis.row_status[irow]))

# Clear so that incumbent model is empty
h.clear()
# Now define the test-semi-definite0 model (from TestQpSolver.cpp) as a HighsModel instance
#
model = highspy.HighsModel()
model.lp_.model_name_ = "semi-definite"
model.lp_.num_col_ = 3
model.lp_.num_row_ = 1
model.lp_.col_cost_ = np.array([1.0, 1.0, 2.0], dtype=np.double)
model.lp_.col_lower_ = np.array([0, 0, 0], dtype=np.double)
model.lp_.col_upper_ = np.array([inf, inf, inf], dtype=np.double)
model.lp_.row_lower_ = np.array([2], dtype=np.double)
model.lp_.row_upper_ = np.array([inf], dtype=np.double)
model.lp_.a_matrix_.format_ = highspy.MatrixFormat.kColwise
model.lp_.a_matrix_.start_ = np.array([0, 1, 2, 3])
model.lp_.a_matrix_.index_ = np.array([0, 0, 0])
model.lp_.a_matrix_.value_ = np.array([1.0, 1.0, 1.0], dtype=np.double)
model.hessian_.dim_ = model.lp_.num_col_
model.hessian_.start_ = np.array([0, 2, 2, 3])
model.hessian_.index_ = np.array([0, 2, 2])
model.hessian_.value_ = np.array([2.0, -1.0, 1.0], dtype=np.double)

print('test-semi-definite0 as HighsModel')
h.passModel(model)
h.run()

# Clear so that incumbent model is empty
h.clear()
num_col = 3
num_row = 1
sense = highspy.ObjSense.kMinimize
offset = 0
col_cost = np.array([1.0, 1.0, 2.0], dtype=np.double)
col_lower = np.array([0, 0, 0], dtype=np.double)
col_upper = np.array([inf, inf, inf], dtype=np.double)
row_lower = np.array([2], dtype=np.double)
row_upper = np.array([inf], dtype=np.double)
a_matrix_format = highspy.MatrixFormat.kColwise
a_matrix_start = np.array([0, 1, 2, 3])
a_matrix_index = np.array([0, 0, 0])
a_matrix_value = np.array([1.0, 1.0, 1.0], dtype=np.double)
a_matrix_num_nz = a_matrix_start[num_col]
hessian_format = highspy.HessianFormat.kTriangular
hessian_start = np.array([0, 2, 2, 3])
hessian_index = np.array([0, 2, 2])
hessian_value = np.array([2.0, -1.0, 1.0], dtype=np.double)
hessian_num_nz = hessian_start[num_col]
integrality = np.array([0, 0, 0])

print('test-semi-definite0 as pointers')
h.passModel(num_col, num_row, a_matrix_num_nz, hessian_num_nz,
            a_matrix_format, hessian_format, sense, offset,
            col_cost, col_lower, col_upper, row_lower, row_upper,
            a_matrix_start, a_matrix_index, a_matrix_value,
            hessian_start, hessian_index, hessian_value,
            integrality)
h.run()
h.writeSolution("", 1)

# Clear so that incumbent model is empty
h.clear()
print('25fv47 as HighsModel')

h.readModel("../../../check/instances/25fv47.mps")
h.presolve()
presolved_lp = h.getPresolvedLp()
# Create a HiGHS instance to solve the presolved LP
print('\nCreate Highs instance to solve presolved LP')
h1 = highspy.Highs()
h1.passModel(presolved_lp)
options = h1.getOptions()
options.presolve = 'off'
options.solver = 'ipm'
h1.passOptions(options)
h1.writeOptions("", True)
h1.run()
solution = h1.getSolution()
basis = h1.getBasis()
print('\nCrossover then postsolve LP using solution and basis from other Highs instance')
h.postsolve(solution, basis)
info = h.getInfo()
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)

run_time = h.getRunTime()
print('Total HiGHS run time is ', run_time)
