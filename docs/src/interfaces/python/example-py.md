# [Examples](@id example-py)

## Initialize HiGHS

HiGHS must be initialized before making calls to the HiGHS Python
library, and the examples below assume that it has been done

```python
import highspy
import numpy as np
h = highspy.Highs()
```

## Read a model

To read a model into HiGHS from a MPS files and (CPLEX) LP files pass the file name to `readModel`. 

```python
# Read a model from MPS file model.mps
filename = 'model.mps'
status = h.readModel(filename)
print('Reading model file ', filename, ' returns a status of ', status)
filename = 'model.dat'
status = h.readModel(filename)
print('Reading model file ', filename, ' returns a status of ', status)
```

## Build a model

Build the model

```raw
minimize    f  =  x0 +  x1
subject to              x1 <= 7
            5 <=  x0 + 2x1 <= 15
            6 <= 3x0 + 2x1
            0 <= x0 <= 4; 1 <= x1
```

Using the simplified interface, the model can be built as follows:

```python
x0 = h.addVariable(lb = 0, ub = 4)
x1 = h.addVariable(lb = 1, ub = 7)

h.addConstr(5 <=   x0 + 2*x1 <= 15)
h.addConstr(6 <= 3*x0 + 2*x1)

h.minimize(x0 + x1)
```

Alternatively, the model can be built using the more general interface, which allows the user to specify the model in a more flexible way.

Firstly, one variable at a time, via a sequence of calls to `addVar` and `addRow`:s
```python
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

Alternatively, via calls to `addCols` and `addRows`.

```python
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

 * `passColName`
 * `passRowName`

## Pass a model

Pass a model from a HighsLp instance
```python
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

## Solve the model

The incumbent model in HiGHS is solved by calling
```python
h.run()
```

## Print solution information

```python
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
!!! warning
    The following are markers for documentation that has yet to be added

## Extract results

 * `getModelStatus`
 * `getInfo`
 * `getSolution`
 * `getBasis`

## Report results

 * `writeSolution`

## [Option values](@id example-py-option-values)

 * `setOptionValue`
 * `getOptionValue`

## Get model data

 * `getNumCol`
 * `getNumRow`
 * `getNumNz`
 * `getCol`
 * `getRow`
 * `getColEntries`
 * `getRowEntries`
 * `getCols`
 * `getRows`
 * `getColsEntries`
 * `getRowsEntries`
 * `getColName`
 * `getColByName`
 * `getRowName`
 * `getRowByName`
 * `getCoeff`

## Modify model data

 * `changeObjectiveSense`
 * `changeColCost`
 * `changeColBounds`
 * `changeRowBounds`
 * `changeColsCost`
 * `changeColsBounds`
 * `changeRowsBounds`
 * `changeCoeff`

## Set solution

 * `setSolution`

## Set basis

 * `setBasis`

## Presolve/postsolve

 * `presolve`
 * `getPresolvedLp`
 * `getPresolvedModel`
 * `getPresolveLog`
 * `getPresolveOrigColsIndex`
 * `getPresolveOrigRowsIndex`
 * `getModelPresolveStatus`
 * `writePresolvedModel`
 * `presolveStatusToString`
 * `presolveRuleTypeToString`
 * `postsolve`
 
## Multi-objective optimization

* `passLinearObjectives`
* `addLinearObjective`
* `clearLinearObjectives`




