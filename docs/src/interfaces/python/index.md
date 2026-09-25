<a id="highspy-documentation"></a>

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

# User Guide

* [Getting Started](getting-started.md)
  * [Installation](getting-started.md#installation)
  * [Import and Initialize](getting-started.md#import-and-initialize)
  * [Read and Solve a Model](getting-started.md#read-and-solve-a-model)
  * [Solver Logging](getting-started.md#solver-logging)
  * [Return Status Values](getting-started.md#return-status-values)
  * [Efficient Solution Access](getting-started.md#efficient-solution-access)
* [Modeling Examples](modeling.md)
  * [Simplified Modeling Interface](modeling.md#simplified-modeling-interface)
  * [Extracting Variable Values](modeling.md#extracting-variable-values)
  * [Matrix-Oriented Interface](modeling.md#matrix-oriented-interface)
  * [Adding Variables and Rows Individually](modeling.md#adding-variables-and-rows-individually)
  * [Passing a `HighsLp` Object](modeling.md#passing-a-highslp-object)
  * [Inspecting Results](modeling.md#inspecting-results)
* [Callbacks](callbacks.md)
  * [Callback API](callbacks.md#callback-api)
* [HiPO in Python](hipo.md)
  * [Installation](hipo.md#installation)
  * [Usage](hipo.md#usage)
  * [Local Installation Requirements](hipo.md#local-installation-requirements)
  * [Uninstall](hipo.md#uninstall)
  * [License](hipo.md#license)

# API Reference

* [API Reference](api.md)
  * [Top-Level Package](api.md#module-highspy)
  * [High-Level Python Interface](api.md#high-level-python-interface)
  * [Core Solver Bindings](api.md#core-solver-bindings)
  * [Callback Bindings](api.md#module-highspy._core.cb)
  * [Simplex Constants](api.md#module-highspy._core.simplex_constants)
