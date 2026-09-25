<a id="getting-started"></a>

# Getting Started

<a id="installation"></a>

## Installation

Install `highspy` from PyPI:

```shell
python -m pip install highspy
```

The package is also available from conda-forge.

For enhanced performance on large linear and quadratic programs via the HiPO
interior point solver, see [HiPO in Python](hipo.md) for installing the optional
`highspy-extras` package.

<a id="import-and-initialize"></a>

## Import and Initialize

Create a [`highspy.Highs`](api.md#highspy.Highs) instance before calling the solver API:

```python
import highspy

h = highspy.Highs()
```

<a id="read-and-solve-a-model"></a>

## Read and Solve a Model

HiGHS can read common optimization model formats such as MPS and LP files.

```python
import highspy

h = highspy.Highs()
h.readModel("model.mps")
h.run()
print(h.modelStatusToString(h.getModelStatus()))
```

<a id="solver-logging"></a>

## Solver Logging

To disable console output, call [`highspy.Highs.silent()`](api.md#highspy.Highs.silent). To send logs
to a file, set the `log_file` option explicitly:

```python
h.silent()
h.setOptionValue("log_file", "highs.log")
```

<a id="return-status-values"></a>

## Return Status Values

Most methods that can fail return a [`highspy.HighsStatus`](api.md#highspy.HighsStatus) value.
Methods that only retrieve already-available data usually return that data
directly.

<a id="efficient-solution-access"></a>

## Efficient Solution Access

Arrays returned by `highspy` may be backed by C++ storage. Convert them to a
Python list or NumPy array before repeated element-by-element access.

```python
solution = h.getSolution()
col_value = list(solution.col_value)
values = [col_value[i] for i in range(len(col_value))]
```
