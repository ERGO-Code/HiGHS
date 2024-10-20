# The paths to MPS file instances assumes that this is run in the
# root directory of HiGHS

import numpy as np
import highspy

hscb = highspy.cb 

h = highspy.Highs()

# h.setOptionValue("log_to_console", True)

inf = highspy.kHighsInf
alt_inf = h.getInfinity()
print('highspy.kHighsInf = ', inf,
      'h.getInfinity() = ', alt_inf)

# ~~~
# Read in and solve avgas
h.readModel("check/instances/avgas.mps")

# h.writeModel("ml.mps")

h.run()
lp = h.getLp()
num_nz = h.getNumNz()
print("LP has ", lp.num_col_,
      " columns", lp.num_row_,
      " rows and ", num_nz, " nonzeros.")

# ~~~
# Clear so that incumbent model is empty
h.clear()

# Now define the blending model as a HighsLp instance
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 2
lp.sense_ = highspy.HighsObjSense.kMaximize
lp.col_cost_ = np.array([8, 10], dtype=np.double)
lp.col_lower_ = np.array([0, 0], dtype=np.double)
lp.col_upper_ = np.array([inf, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, -inf], dtype=np.double)
lp.row_upper_ = np.array([120, 210], dtype=np.double)
lp.a_matrix_.start_ = np.array([0, 2, 4])
lp.a_matrix_.index_ = np.array([0, 1, 0, 1])
lp.a_matrix_.value_ = np.array([0.3, 0.7, 0.5, 0.5], dtype=np.double)
h.passModel(lp)

# Solve
h.run()

# Print solution
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)
print(
    "Primal solution status = ", h.solutionStatusToString(
        info.primal_solution_status)
)
print("Dual solution status = ",
      h.solutionStatusToString(info.dual_solution_status))
print("Basis validity = ", h.basisValidityToString(info.basis_validity))
num_var = h.getNumCol()
num_row = h.getNumRow()
print("Variables")
for icol in range(num_var):
    print(icol, solution.col_value[icol],
          h.basisStatusToString(basis.col_status[icol]))
print("Constraints")
for irow in range(num_row):
    print(irow, solution.row_value[irow],
          h.basisStatusToString(basis.row_status[irow]))

# ~~~
# Clear so that incumbent model is empty
h.clear()

# Now define the test-semi-definite0 model (from TestQpSolver.cpp) 
# as a HighsModel instance
model = highspy.HighsModel()
model.lp_.model_name_ = "semi-definite"
model.lp_.num_col_ = 3
model.lp_.num_row_ = 1
model.lp_.col_cost_ = np.array([1.0, 1.0, 2.0], dtype=np.double)
model.lp_.col_lower_ = np.array([0, 0, 0], dtype=np.double)
model.lp_.col_upper_ = np.array([inf, inf, inf], dtype=np.double)
model.lp_.row_lower_ = np.array([2], dtype=np.double)
model.lp_.row_upper_ = np.array([inf], dtype=np.double)
model.lp_.a_matrix_.format_ = highspy.HighsMatrixFormat.kColwise
model.lp_.a_matrix_.start_ = np.array([0, 1, 2, 3])
model.lp_.a_matrix_.index_ = np.array([0, 0, 0])
model.lp_.a_matrix_.value_ = np.array([1.0, 1.0, 1.0], dtype=np.double)
model.hessian_.dim_ = model.lp_.num_col_
model.hessian_.start_ = np.array([0, 2, 2, 3])
model.hessian_.index_ = np.array([0, 2, 2])
model.hessian_.value_ = np.array([2.0, -1.0, 1.0], dtype=np.double)

print("test-semi-definite0 as HighsModel")
h.passModel(model)
h.run()

# ~~~
# Clear so that incumbent model is empty
h.clear()
num_col = 3
num_row = 1
sense = highspy.HighsObjSense.kMinimize
offset = 0
col_cost = np.array([1.0, 1.0, 2.0], dtype=np.double)
col_lower = np.array([0, 0, 0], dtype=np.double)
col_upper = np.array([inf, inf, inf], dtype=np.double)
row_lower = np.array([2], dtype=np.double)
row_upper = np.array([inf], dtype=np.double)
a_matrix_format = highspy.HighsMatrixFormat.kColwise
a_matrix_start = np.array([0, 1, 2, 3])
a_matrix_index = np.array([0, 0, 0])
a_matrix_value = np.array([1.0, 1.0, 1.0], dtype=np.double)
a_matrix_num_nz = a_matrix_start[num_col]
hessian_format = highspy.HighsHessianFormat.kTriangular
hessian_start = np.array([0, 2, 2, 3])
hessian_index = np.array([0, 2, 2])
hessian_value = np.array([2.0, -1.0, 1.0], dtype=np.double)
hessian_num_nz = hessian_start[num_col]
integrality = np.array([0, 0, 0])

print("test-semi-definite0 as pointers")
h.passModel(
    num_col,
    num_row,
    a_matrix_num_nz,
    hessian_num_nz,
    a_matrix_format,
    hessian_format,
    sense,
    offset,
    col_cost,
    col_lower,
    col_upper,
    row_lower,
    row_upper,
    a_matrix_start,
    a_matrix_index,
    a_matrix_value,
    hessian_start,
    hessian_index,
    hessian_value,
    integrality,
)
h.run()
h.writeSolution("", 1)

# ~~~
# Clear so that incumbent model is empty
h.clear()
print("25fv47 as HighsModel")

h.readModel("check/instances/25fv47.mps")

h.presolve()
presolved_lp = h.getPresolvedLp()

# Create a HiGHS instance to solve the presolved LP
print('\nCreate Highs instance to solve presolved LP')
h1 = highspy.Highs()
h1.passModel(presolved_lp)

# Get and set options
options = h1.getOptions()
options.presolve = "off"
options.solver = "ipm"

h1.passOptions(options)

# can be used to check option values
# h1.writeOptions("")

h1.run()
solution = h1.getSolution()
basis = h1.getBasis()

print("Crossover, then postsolve using solution and basis from another instance")

h.postsolve(solution, basis)

# Get solution
info = h.getInfo()
model_status = h.getModelStatus()
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)

run_time = h.getRunTime()
print("Total HiGHS run time is ", run_time)

# Get an optimal basis

# Clear so that incumbent model is empty
h.clear()
print("25fv47 as HighsModel")

h.readModel("check/instances/25fv47.mps")

h.run()
simplex_iteration_count = h.getInfo().simplex_iteration_count
print("From initial basis, simplex iteration count =", simplex_iteration_count)
basis = h.getBasis()
h.clearSolver()

h.setBasis(basis)
h.run()
simplex_iteration_count = h.getInfo().simplex_iteration_count
print("From optimal basis, simplex iteration count =", simplex_iteration_count)
status = h.setBasis()
h.run()
simplex_iteration_count = h.getInfo().simplex_iteration_count
print("From logical basis, simplex iteration count =", simplex_iteration_count)


# Define a callback

def user_interrupt_callback(
    callback_type,
    message,
    data_out,
    data_in,
    user_callback_data
):
    # dev_run = True
    dev_run = False

    # Constants for iteration limits or objective targets, adjust as required
    SIMPLEX_ITERATION_LIMIT = 100
    IPM_ITERATION_LIMIT = 100
    EGOUT_OBJECTIVE_TARGET = 1.0

    # Callback for MIP Improving Solution
    if callback_type == hscb.HighsCallbackType.kCallbackMipImprovingSolution:
        # Assuming it is a list or array
        assert user_callback_data is not None, "User callback data is None!"
        local_callback_data = user_callback_data[0]

        if dev_run:
            print(f"userCallback(type {callback_type};")
            print(f"data {local_callback_data:.4g}): {message}")
            print(f"with objective {data_out.objective_function_value}")
            print(f"and solution[0] = {data_out.mip_solution[0]}")

        # Check and update the objective function value
        assert (
            local_callback_data >= data_out.objective_function_value
        ), "Objective function value is invalid!"
        user_callback_data[0] = data_out.objective_function_value

    else:
        # Various other callback types
        if callback_type == hscb.HighsCallbackType.kCallbackLogging:
            if dev_run:
                print(f"userInterruptCallback(type {callback_type}): {message}")

        elif callback_type == hscb.HighsCallbackType.kCallbackSimplexInterrupt:
            if dev_run:
                print(f"userInterruptCallback(type {callback_type}): {message}")
                print("with iteration count = ", 
                      data_out.simplex_iteration_count)

            data_in.user_interrupt = (
                data_out.simplex_iteration_count > SIMPLEX_ITERATION_LIMIT
            )

        elif callback_type == hscb.HighsCallbackType.kCallbackIpmInterrupt:
            if dev_run:
                print(f"userInterruptCallback(type {callback_type}): {message}")
                print(f"with iteration count = {data_out.ipm_iteration_count}")

            data_in.user_interrupt = (
                data_out.ipm_iteration_count > IPM_ITERATION_LIMIT
            )

        elif callback_type == hscb.HighsCallbackType.kCallbackMipInterrupt:
            if dev_run:
                print(f"userInterruptCallback(type {callback_type}): {message}")
                print(f"Dual bound = {data_out.mip_dual_bound:.4g}")
                print(f"Primal bound = {data_out.mip_primal_bound:.4g}")
                print(f"Gap = {data_out.mip_gap:.4g}")
                print(f"Objective = {data_out.objective_function_value:.4g}")

            data_in.user_interrupt = (
                data_out.objective_function_value < EGOUT_OBJECTIVE_TARGET
            )


# Define model
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


# Set callback and run
h.setCallback(user_interrupt_callback, None)
h.startCallback(hscb.HighsCallbackType.kCallbackLogging)

h.run()
h.stopCallback(hscb.HighsCallbackType.kCallbackLogging)

# Get solution
num_var = h.getNumCol()
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)
print(
    "Primal solution status = ", h.solutionStatusToString(
        info.primal_solution_status)
)
print("Dual solution status = ",
      h.solutionStatusToString(info.dual_solution_status))
print("Basis validity = ", h.basisValidityToString(info.basis_validity))
print("Variables:")
for icol in range(0, 5):
    print(icol, solution.col_value[icol],
          h.basisStatusToString(basis.col_status[icol]))
print("...")
for icol in range(num_var-2, num_var):
    print(icol, solution.col_value[icol],
          h.basisStatusToString(basis.col_status[icol]))

print("computing IIS for lp-incompatible-bounds")
"""
LP has row0 and col2 with inconsistent bounds.

When prioritising rows, row0 and its constituent columns (1, 2) should be found
When prioritising columns, col2 and its constituent rows (0, 1) should be found
"""
# Define the LP
lp = highspy.HighsLp()
lp.num_col_ = 3
lp.num_row_ = 2
lp.col_cost_ = np.array([0, 0, 0], dtype=np.double)
lp.col_lower_ = np.array([0, 0, 0], dtype=np.double)
lp.col_upper_ = np.array([1, 1, -1], dtype=np.double)
lp.row_lower_ = np.array([1, 0], dtype=np.double)
lp.row_upper_ = np.array([0, 1], dtype=np.double)
lp.a_matrix_.format_ = highspy.HighsMatrixFormat.kRowwise
lp.a_matrix_.start_ = np.array([0, 2, 4])
lp.a_matrix_.index_ = np.array([1, 2, 0, 2])
lp.a_matrix_.value_ = np.array([1, 1, 1, 1], dtype=np.double)

h.clear()
h.passModel(lp)
h.run()
assert h.getModelStatus() == highspy.HighsModelStatus.kInfeasible

# Set IIS strategy to row priority and get IIS
h.setOptionValue("iis_strategy", highspy.HighsIisStrategy.kFromLpRowPriority)

iis = highspy.HighsIis()
assert h.getIis(iis) == highspy.HighsStatus.kOk
assert len(iis.col_index) == 0
assert len(iis.row_index) == 1
assert iis.row_index[0] == 0
assert iis.row_bound[0] == highspy.HighsIisBoundStatus.kBoxed

# Set IIS strategy to column priority and get IIS
h.setOptionValue("iis_strategy", highspy.HighsIisStrategy.kFromLpColPriority)
iis.invalidate()
assert h.getIis(iis) == highspy.HighsStatus.kOk
assert len(iis.col_index) == 1
assert len(iis.row_index) == 0
assert iis.col_index[0] == 2
assert iis.col_bound[0] == highspy.HighsIisBoundStatus.kBoxed

print("IIS computation completed successfully")

print("computing feasibility relaxation")
h.clear()
inf = h.getInfinity()

num_col = 2
num_row = 3
num_nz = 6
a_format = highspy.HighsMatrixFormat.kColwise
sense = highspy.HighsObjSense.kMinimize
offset = 0
col_cost = np.array([1, -2], dtype=np.double)
col_lower = np.array([5, -inf], dtype=np.double)
col_upper = np.array([inf, inf], dtype=np.double)
row_lower = np.array([2, -inf, -inf], dtype=np.double)
row_upper = np.array([inf, 1, 20], dtype=np.double)
a_start = np.array([0, 3])
a_index = np.array([0, 1, 2, 0, 1, 2])
a_value = np.array([-1, -3, 20, 21, 2, 1], dtype=np.double)
integrality = np.array([highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger])

h.passModel(
    num_col, num_row, num_nz, a_format, sense, offset,
    col_cost, col_lower, col_upper,
    row_lower, row_upper,
    a_start, a_index, a_value,
    integrality
)

assert h.feasibilityRelaxation(1, 1, 1) == highspy.HighsStatus.kOk

info = h.getInfo()
objective_function_value = info.objective_function_value

solution = h.getSolution()

assert abs(objective_function_value - 5) < 1e-6, f"Expected objective value 5, got {objective_function_value}"
assert abs(solution.col_value[0] - 1) < 1e-6, f"Expected solution[0] = 1, got {solution.col_value[0]}"
assert abs(solution.col_value[1] - 1) < 1e-6, f"Expected solution[1] = 1, got {solution.col_value[1]}"

print("Feasibility Relaxation Test Passed")

# Using infeasible LP from AMPL documentation
h = highspy.Highs()
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 3
lp.col_cost_ = np.array([1, -2], dtype=np.double)
lp.col_lower_ = np.array([5, -h.getInfinity()], dtype=np.double)
lp.col_upper_ = np.array([h.getInfinity(), h.getInfinity()], dtype=np.double)
lp.col_names_ = ["X", "Y"]
lp.row_lower_ = np.array([2, -h.getInfinity(), -h.getInfinity()], dtype=np.double)
lp.row_upper_ = np.array([h.getInfinity(), 1, 20], dtype=np.double)
lp.row_names_ = ["R0", "R1", "R2"]
lp.a_matrix_.start_ = np.array([0, 3, 6])
lp.a_matrix_.index_ = np.array([0, 1, 2, 0, 1, 2])
lp.a_matrix_.value_ = np.array([-1, -3, 20, 21, 2, 1], dtype=np.double)
lp.integrality_ = np.array([highspy.HighsVarType.kInteger, highspy.HighsVarType.kInteger])
h.setOptionValue("output_flag", False)
h.passModel(lp)

# Vanilla feasibility relaxation
print("Vanilla feasibility relaxation")
h.feasibilityRelaxation(1, 1, 1)
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]})")
print(f"Slacks: ({solution.row_value[0] - lp.row_lower_[0]}, "
      f"{lp.row_upper_[1] - solution.row_value[1]}, "
      f"{lp.row_upper_[2] - solution.row_value[2]})")

# Respect all lower bounds
print("\nRespect all lower bounds")
h.feasibilityRelaxation(-1, 1, 1)
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]})")
print(f"Slacks: ({solution.row_value[0] - lp.row_lower_[0]}, "
      f"{lp.row_upper_[1] - solution.row_value[1]}, "
      f"{lp.row_upper_[2] - solution.row_value[2]})")

# Local penalties RHS {1, -1, 10}
print("\nLocal penalties RHS {1, -1, 10}")
local_rhs_penalty = np.array([1, -1, 10], dtype=np.double)
h.feasibilityRelaxation(1, 1, 0, None, None, local_rhs_penalty)
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]})")
print(f"Slacks: ({solution.row_value[0] - lp.row_lower_[0]}, "
      f"{lp.row_upper_[1] - solution.row_value[1]}, "
      f"{lp.row_upper_[2] - solution.row_value[2]})")

# Local penalties RHS {10, 1, 1}
print("\nLocal penalties RHS {10, 1, 1}")
local_rhs_penalty = np.array([10, 1, 1], dtype=np.double)
h.feasibilityRelaxation(1, 1, 0, None, None, local_rhs_penalty)
solution = h.getSolution()
print(f"Solution: ({solution.col_value[0]}, {solution.col_value[1]})")
print(f"Slacks: ({solution.row_value[0] - lp.row_lower_[0]}, "
      f"{lp.row_upper_[1] - solution.row_value[1]}, "
      f"{lp.row_upper_[2] - solution.row_value[2]})")

iis = highspy.HighsIis()
assert h.getIis(iis) == highspy.HighsStatus.kOk

print("\nIIS")
print("row_index:", iis.row_index)
print("row_bound:", iis.row_bound)

print("col_index:", iis.col_index)
print("col_bound:", iis.col_bound)
