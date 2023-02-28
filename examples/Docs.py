import highspy
import numpy as np

# minimize    f  =  x0 +  x1
# subject to              x1 <= 7
#             5 <=  x0 + 2x1 <= 15
#             6 <= 3x0 + 2x1
#             0 <= x0 <= 4; 1 <= x1
# Highs h
h = highspy.Highs()
inf = highspy.kHighsInf

# Load a model from MPS file model.mps
print("\nLoading the model from an MPS file")
filename = 'model.mps'
h.readModel(filename)
h.run()

print("\nBuilding a model using single variables and constraints")
# Build a model with single variables and constraints
h.clear()
# Define two variables, first using identifiers for the bound values,
# and then using constants
lower = 0
upper = 4
h.addVar(lower, upper)
h.addVar(1, inf)
# Define the objective coefficients (costs) of the two variables,
# identifying the variable by index, and defining its cost
cost = 1
h.changeColCost(0, cost)
h.changeColCost(1, 1)
# Define constraints for the model
#
# The first constraint (x_1<=7) has only one nonzero coefficient,
# identified by variable index 1 and value 1
lower = -inf
upper = 7
num_nz = 1
index = 1
value = 1
h.addRow(lower, upper, num_nz, index, value)
# The second constraint (5 <= x_0 + 2x_1 <= 15) has two nonzero
# coefficients, so arrays of indices and values are required
num_nz = 2
index = np.array([0, 1])
value = np.array([1, 2], dtype=np.double)
h.addRow(5, 15, num_nz, index, value)
# The final constraint (6 <= 3x_0 + 2x_1) has the same indices but different values
num_nz = 2
value = np.array([3, 2], dtype=np.double)
h.addRow(6, inf, num_nz, index, value)

# Access LP 
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns', lp.num_row_, ' rows and ', num_nz, ' nonzeros')

#h.writeModel("")

# Build a model with multiple columns and rows
h.clear()
print("\nBuilding a model using multiple variables and constraints")
# The constraint matrix is defined with the rows below, but parameters
# for an empty (column-wise) matrix must be passed
cost = np.array([1, 1], dtype=np.double)
lower = np.array([0, 1], dtype=np.double)
upper = np.array([4, inf], dtype=np.double)
num_nz = 0
start = 0
index = 0
value = 0
h.addCols(2, cost, lower, upper, num_nz, start, index, value)
# Add the rows, with the constraint matrix row-wise
lower = np.array([-inf, 5, 6], dtype=np.double)
upper = np.array([7, 15, inf], dtype=np.double)
num_nz = 5
start = np.array([0, 1, 3])
index = np.array([1, 0, 1, 0, 1])
value = np.array([1, 1, 2, 3, 2], dtype=np.double)
h.addRows(3, lower, upper, num_nz, start, index, value)

h.writeModel("")
h.run()

# Pass the following model from a HighsLp instance
h.clear()
print("Passing the model via HighsLp")
lp = highspy.HighsLp()
lp.num_col_ = 2;
lp.num_row_ = 3;
lp.col_cost_ = np.array([1, 1], dtype=np.double)
lp.col_lower_ = np.array([0, 1], dtype=np.double)
lp.col_upper_ = np.array([4, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, 5, 6], dtype=np.double)
lp.row_upper_ = np.array([7, 15, inf], dtype=np.double)
# In a HighsLp instsance, the number of nonzeros is given by a fictitious final start
lp.a_matrix_.start_ = np.array([0, 2, 5])
lp.a_matrix_.index_ = np.array([1, 2, 0, 1, 2])
lp.a_matrix_.value_ = np.array([1, 3, 1, 2, 2], dtype=np.double)
h.passModel(lp)
h.writeModel("")
h.run()

