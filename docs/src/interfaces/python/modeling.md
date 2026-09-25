<a id="modeling-examples"></a>

# Modeling Examples

<a id="simplified-modeling-interface"></a>

## Simplified Modeling Interface

The simplified interface creates variable and constraint wrapper objects that
can be combined into linear expressions.

```python
import highspy

h = highspy.Highs()
x0 = h.addVariable(lb=0, ub=4)
x1 = h.addVariable(lb=1, ub=7)

h.addConstr(x1 <= 7)
h.addConstr(5 <= x0 + 2 * x1 <= 15)
h.addConstr(6 <= 3 * x0 + 2 * x1)

h.minimize(x0 + x1)
print(h.getObjectiveValue())
```

<a id="extracting-variable-values"></a>

## Extracting Variable Values

Use [`highspy.Highs.val()`](api.md#highspy.Highs.val) to read back the value of a single variable
(or linear expression) after solving, or [`highspy.Highs.vals()`](api.md#highspy.Highs.vals) for a
collection of them at once:

```python
print(h.val(x0))
print(h.vals([x0, x1]))
```

<a id="matrix-oriented-interface"></a>

## Matrix-Oriented Interface

The lower-level interface accepts explicit column and row data. Use
`highspy.kHighsInf` for infinite bounds.

```python
import highspy
import numpy as np

h = highspy.Highs()
inf = highspy.kHighsInf

cost = np.array([1, 1], dtype=np.double)
lower = np.array([0, 1], dtype=np.double)
upper = np.array([4, inf], dtype=np.double)
h.addCols(2, cost, lower, upper, 0, [0], [0], [0])

row_lower = np.array([-inf, 5, 6], dtype=np.double)
row_upper = np.array([7, 15, inf], dtype=np.double)
start = np.array([0, 1, 3], dtype=np.int32)
index = np.array([1, 0, 1, 0, 1], dtype=np.int32)
value = np.array([1, 1, 2, 3, 2], dtype=np.double)
h.addRows(3, row_lower, row_upper, 5, start, index, value)

h.run()
```

<a id="adding-variables-and-rows-individually"></a>

## Adding Variables and Rows Individually

Variables and rows can also be added one at a time with
[`highspy.Highs.addVar()`](api.md#highspy.Highs.addVar) and [`highspy.Highs.addRow()`](api.md#highspy.Highs.addRow), with
costs updated afterwards via [`highspy.Highs.changeColCost()`](api.md#highspy.Highs.changeColCost):

```python
import highspy

h = highspy.Highs()
inf = highspy.kHighsInf

# Add two variables, then set their objective coefficients (costs) by
# index; a newly added variable has cost 0 by default.
h.addVar(0, 4)
h.addVar(1, inf)
h.changeColCost(0, 1)
h.changeColCost(1, 1)

# x1 <= 7
h.addRow(-inf, 7, 1, [1], [1])
# 5 <= x0 + 2 * x1 <= 15
h.addRow(5, 15, 2, [0, 1], [1, 2])
# 6 <= 3 * x0 + 2 * x1
h.addRow(6, inf, 2, [0, 1], [3, 2])

lp = h.getLp()
print("LP has", lp.num_col_, "columns,", lp.num_row_, "rows, and", h.getNumNz(), "nonzeros")
```

<a id="passing-a-highslp-object"></a>

## Passing a `HighsLp` Object

For full control, populate a [`highspy.HighsLp`](api.md#highspy.HighsLp) object and pass it
to HiGHS.

```python
import highspy
import numpy as np

inf = highspy.kHighsInf
lp = highspy.HighsLp()
lp.num_col_ = 2
lp.num_row_ = 3
lp.col_cost_ = np.array([1, 1], dtype=np.double)
lp.col_lower_ = np.array([0, 1], dtype=np.double)
lp.col_upper_ = np.array([4, inf], dtype=np.double)
lp.row_lower_ = np.array([-inf, 5, 6], dtype=np.double)
lp.row_upper_ = np.array([7, 15, inf], dtype=np.double)
lp.a_matrix_.start_ = np.array([0, 2, 5], dtype=np.int32)
lp.a_matrix_.index_ = np.array([1, 2, 0, 1, 2], dtype=np.int32)
lp.a_matrix_.value_ = np.array([1, 3, 1, 2, 2], dtype=np.double)

h = highspy.Highs()
h.passModel(lp)
h.run()
```

<a id="inspecting-results"></a>

## Inspecting Results

```python
solution = h.getSolution()
basis = h.getBasis()
info = h.getInfo()
model_status = h.getModelStatus()

print("Model status:", h.modelStatusToString(model_status))
print("Objective:", info.objective_function_value)
print("Column values:", list(solution.col_value))
print("Primal solution status:", h.solutionStatusToString(info.primal_solution_status))
print("Dual solution status:", h.solutionStatusToString(info.dual_solution_status))
print("Basis validity:", h.basisValidityToString(info.basis_validity))
```
