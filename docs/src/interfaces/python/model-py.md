# [Modelling](@id model-py)

HiGHS has a rudimentry modelling language that allows models to be built and run using `highspy`. 

Below is an example of building a mathematical LP. The functions used are documented in detail below
```
# model and solve the LP
#
# maximize 10 x1 + 25x2
# s. t.       x1 +  2x2 <=  80
#             x1 +  4x2 <= 120
#             x1 >= 0; x2 >= 0
import highspy

h = highspy.Highs()

x1 = h.addVar()
x2 = h.addVar()

h.addConstr(x1 + 2*x2 <=  80)
h.addConstr(x1 + 4*x2 <= 120)

h.maximize(10*x1 + 25*x2)

print("x1 = ", h.val(x1))
print("x2 = ", h.val(x2))
```

## addVar

Adds a variable to the model. By default it is continuous,
non-negative, with zero objective coefficient, and has no name
associated with it.

```
addVar(lb = 0, ub = kHighsInf, obj = 0, type=HighsVarType.kContinuous, name = None)
```

## addConstr

Adds a constraint to the model. It must be defined in terms of a
linear function, with `*` used when there are non-unit
coefficients. By default it has a lower bound of -infinity, an upper
bound of +infinity, and no name associated with it.

```
addConstr(cons, name = None)
```

## maximize

Calls HiGHS to maximize the objective. By default it uses the
objective coefficients defined when the variables were added to the
model. However, a linear function can be passed as an argument.

```
maximize(obj=None)
```

## minimize

Calls HiGHS to minimize the objective. By default it uses the
objective coefficients defined when the variables were added to the
model. However, a linear function can be passed as an argument.

```
minimize(obj=None)
```

## val

Extracts the current value of a particular variable

```
val(var)
```

## vals

Extracts the current values of a particular set of variables

```
vals(vars)
```

## MIP Example




