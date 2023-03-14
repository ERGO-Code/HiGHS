## Initialize HiGHS

HiGHS must be initialized before making calls to the HiGHS Python
library, and the examples below assume that it has been done

```
import highspy
import numpy as np

# Highs h
h = highspy.Highs()
```

## Load a model

```
# Load a model from MPS file model.mps
filename = 'model.mps'
h.readModel(filename)
```

## Build a model

Build the model

```
minimize    f  =  x0 +  x1
subject to              x1 <= 7
            5 <=  x0 + 2x1 <= 15
            6 <= 3x0 + 2x1
            0 <= x0 <= 4; 1 <= x1
```

Firstly, one variable at a time, via a sequence of calls to __addVar__ and __addRow__.
```
inf = highspy.kHighsInf
# Define two variables, first using identifiers for the bound values,
# and then using constants
lower = 0
upper = 4
h.addVar(lower, upper)
h.addVar(1, inf)

# Define the objective coefficients (costs) of the two variables,
# identifying the variable by index, and changing its cost from the
# default value of zero
cost = 1
h.changeColCost(0, cost)
h.changeColCost(1, 1)

# Define constraints for the model
#
# The first constraint (x1<=7) has only one nonzero coefficient,
# identified by variable index 1 and value 1
num_nz = 1
index = 1
value = 1
h.addRow(-inf, 7, num_nz, index, value)

# The second constraint (5 <= x0 + 2x1 <= 15) has two nonzero
# coefficients, so arrays of indices and values are required
num_nz = 2
index = np.array([0, 1])
value = np.array([1, 2])
h.addRow(5, 15, num_nz, index, value)

# The final constraint (6 <= 3x0 + 2x1) has the same indices but
# different values
num_nz = 2
value = np.array([3, 2])
h.addRow(6, inf, num_nz, index, value)

# Access LP 
lp = h.getLp()
num_nz = h.getNumNz()
print('LP has ', lp.num_col_, ' columns', lp.num_row_, ' rows and ', num_nz, ' nonzeros')

```
Alternatively, via calls to __addCols__ and __addRows__.

```
inf = highspy.kHighsInf
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
```

## Pass a model

Pass a model from a HighsLp instance
```
inf = highspy.kHighsInf
# Define a HighsLp instance
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
```

## Solve a model

The incumbent model in HiGHS is solved by calling
```
h.run()
```

## Print solution information 
```
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()
print('Model status = ', h.modelStatusToString(model_status))
print()
print('Optimal objective = ', info.objective_function_value)
print('Iteration count = ', info.simplex_iteration_count)
print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
print('Basis validity = ', h.basisValidityToString(info.basis_validity))
```