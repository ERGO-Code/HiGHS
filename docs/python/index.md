# highspy Documentation

`highspy` is the Python interface to the HiGHS optimization solver. It
supports reading model files, building linear and mixed-integer models from
Python, solving them with HiGHS, and inspecting solutions, bases, ranging, and
solver information.

Install `highspy` from PyPI or conda-forge, then create a solver instance:

```python
import highspy

h = highspy.Highs()
```

The simplified modeling interface lets you build expressions directly from
variables:

```python
x0 = h.addVariable(lb=0, ub=4)
x1 = h.addVariable(lb=1, ub=7)

h.addConstr(5 <= x0 + 2 * x1 <= 15)
h.addConstr(6 <= 3 * x0 + 2 * x1)
h.minimize(x0 + x1)
```

```{toctree}
:maxdepth: 2
:caption: User Guide

getting-started
modeling
callbacks
hipo
```

```{toctree}
:maxdepth: 2
:caption: API Reference

api
```
