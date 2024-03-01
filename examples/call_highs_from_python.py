# NB Must have installed highspy, either by using pip install highspy,
# or by following the instructions on
# https://github.com/ERGO-Code/HiGHS#python
#
# The paths to MPS file instances assumes that this is run in the
# directory of this file (ie highs/examples) or any other subdirectory
# of HiGHS
import numpy as np
import highspy._highs
import highspy._highs.cb as hscb

h = highspy._highs.Highs()
inf = highspy._highs.kHighsInf
alt_inf = h.getInfinity()
print('highspy._highs.kHighsInf = ', inf, '; h.getInfinity() = ', alt_inf)

# NB The callbacks are not available in
# https://pypi.org/project/highspy/ (HiGHS v1.5.3). To use them,
# highspy must be installed locally using (at least) HiGHS v1.6.0
def user_interrupt_callback(
    callback_type, message, data_out, data_in, user_callback_data
):
    dev_run = True  # or any other condition

    # Constants for iteration limits or objective targets, adjust as per requirement
    ADLITTLE_SIMPLEX_ITERATION_LIMIT = 100
    ADLITTLE_IPM_ITERATION_LIMIT = 100
    EGOUT_OBJECTIVE_TARGET = 1.0

    # Callback for MIP Improving Solution
    if callback_type == hscb.HighsCallbackType.kCallbackMipImprovingSolution:
        assert user_callback_data is not None, "User callback data is None!"
        local_callback_data = user_callback_data[0]  # Assuming it is a list or array

        if dev_run:
            print(
                f"userCallback(type {callback_type}; data {local_callback_data:.4g}): {message} with objective {data_out.objective_function_value} and solution[0] = {data_out.mip_solution[0]}"
            )

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
                print(
                    f"userInterruptCallback(type {callback_type}): {message} with iteration count = {data_out.simplex_iteration_count}"
                )
            data_in.user_interrupt = (
                data_out.simplex_iteration_count > ADLITTLE_SIMPLEX_ITERATION_LIMIT
            )

        elif callback_type == hscb.HighsCallbackType.kCallbackIpmInterrupt:
            if dev_run:
                print(
                    f"userInterruptCallback(type {callback_type}): {message} with iteration count = {data_out.ipm_iteration_count}"
                )
            data_in.user_interrupt = (
                data_out.ipm_iteration_count > ADLITTLE_IPM_ITERATION_LIMIT
            )

        elif callback_type == hscb.HighsCallbackType.kCallbackMipInterrupt:
            if dev_run:
                print(
                    f"userInterruptCallback(type {callback_type}): {message} with Bounds ({data_out.mip_dual_bound:.4g}, {data_out.mip_primal_bound:.4g}); Gap = {data_out.mip_gap:.4g}; Objective = {data_out.objective_function_value:.4g}"
                )
            data_in.user_interrupt = (
                data_out.objective_function_value < EGOUT_OBJECTIVE_TARGET
            )


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
h.setOptionValue("log_to_console", True)
h.setCallback(user_interrupt_callback, None)
h.startCallback(hscb.HighsCallbackType.kCallbackLogging)
h.run()
h.stopCallback(hscb.HighsCallbackType.kCallbackLogging)
num_var = h.getNumCol()
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)
print(
    "Primal solution status = ", h.solutionStatusToString(info.primal_solution_status)
)
print("Dual solution status = ", h.solutionStatusToString(info.dual_solution_status))
print("Basis validity = ", h.basisValidityToString(info.basis_validity))
for icol in range(num_var):
    print(icol, solution.col_value[icol], h.basisStatusToString(basis.col_status[icol]))

# Read in and solve avgas
h.readModel("check/instances/avgas.mps")
# h.writeModel("ml.mps")
h.run()
lp = h.getLp()
num_nz = h.getNumNz()
print(
    "LP has ", lp.num_col_, " columns", lp.num_row_, " rows and ", num_nz, " nonzeros"
)

# Clear so that incumbent model is empty
h.clear()

# Now define the blending model as a HighsLp instance
#
lp = highspy._highs.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 2
lp.sense_ = highspy._highs.ObjSense.kMaximize
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
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)
print(
    "Primal solution status = ", h.solutionStatusToString(info.primal_solution_status)
)
print("Dual solution status = ", h.solutionStatusToString(info.dual_solution_status))
print("Basis validity = ", h.basisValidityToString(info.basis_validity))
num_var = h.getNumCol()
num_row = h.getNumRow()
print("Variables")
for icol in range(num_var):
    print(icol, solution.col_value[icol], h.basisStatusToString(basis.col_status[icol]))
print("Constraints")
for irow in range(num_row):
    print(irow, solution.row_value[irow], h.basisStatusToString(basis.row_status[irow]))

# Clear so that incumbent model is empty
h.clear()
# Now define the test-semi-definite0 model (from TestQpSolver.cpp) as a HighsModel instance
#
model = highspy._highs.HighsModel()
model.lp_.model_name_ = "semi-definite"
model.lp_.num_col_ = 3
model.lp_.num_row_ = 1
model.lp_.col_cost_ = np.array([1.0, 1.0, 2.0], dtype=np.double)
model.lp_.col_lower_ = np.array([0, 0, 0], dtype=np.double)
model.lp_.col_upper_ = np.array([inf, inf, inf], dtype=np.double)
model.lp_.row_lower_ = np.array([2], dtype=np.double)
model.lp_.row_upper_ = np.array([inf], dtype=np.double)
model.lp_.a_matrix_.format_ = highspy._highs.MatrixFormat.kColwise
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

# Clear so that incumbent model is empty
h.clear()
num_col = 3
num_row = 1
sense = highspy._highs.ObjSense.kMinimize
offset = 0
col_cost = np.array([1.0, 1.0, 2.0], dtype=np.double)
col_lower = np.array([0, 0, 0], dtype=np.double)
col_upper = np.array([inf, inf, inf], dtype=np.double)
row_lower = np.array([2], dtype=np.double)
row_upper = np.array([inf], dtype=np.double)
a_matrix_format = highspy._highs.MatrixFormat.kColwise
a_matrix_start = np.array([0, 1, 2, 3])
a_matrix_index = np.array([0, 0, 0])
a_matrix_value = np.array([1.0, 1.0, 1.0], dtype=np.double)
a_matrix_num_nz = a_matrix_start[num_col]
hessian_format = highspy._highs.HessianFormat.kTriangular
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

# Clear so that incumbent model is empty
h.clear()
print("25fv47 as HighsModel")

h.readModel("check/instances/25fv47.mps")
h.presolve()
presolved_lp = h.getPresolvedLp()
# Create a HiGHS instance to solve the presolved LP
print('\nCreate Highs instance to solve presolved LP')
h1 = highspy._highs.Highs()
h1.passModel(presolved_lp)
options = h1.getOptions()
options.presolve = "off"
options.solver = "ipm"
h1.passOptions(options)
h1.writeOptions("")
h1.run()
solution = h1.getSolution()
basis = h1.getBasis()
print(
    "\nCrossover then postsolve LP using solution and basis from other Highs instance"
)
h.postsolve(solution, basis)
info = h.getInfo()
model_status = h.getModelStatus()
print("Model status = ", h.modelStatusToString(model_status))
print("Optimal objective = ", info.objective_function_value)
print("Iteration count = ", info.simplex_iteration_count)

run_time = h.getRunTime()
print("Total HiGHS run time is ", run_time)
