# Getting Started

## Installation

Install `highspy` from PyPI:

```shell
python -m pip install highspy
```

The package is also available from conda-forge.

For enhanced performance on large linear and quadratic programs via the HiPO
interior point solver, see {doc}`hipo` for installing the optional
`highspy-extras` package.

## Import and Initialize

Create a {py:class}`highspy.Highs` instance before calling the solver API:

```python
import highspy

h = highspy.Highs()
```

## Read and Solve a Model

HiGHS can read common optimization model formats such as MPS and LP files.

```python
import highspy

h = highspy.Highs()
h.readModel("model.mps")
h.run()
print(h.modelStatusToString(h.getModelStatus()))
```

## Solver Logging

To disable console output, call {py:meth}`highspy.Highs.silent`. To send logs
to a file, set the `log_file` option explicitly:

```python
h.silent()
h.setOptionValue("log_file", "highs.log")
```

## Return Status Values

Most methods that can fail return a {py:class}`highspy.HighsStatus` value.
Methods that only retrieve already-available data usually return that data
directly.

## Efficient Solution Access

Arrays returned by `highspy` may be backed by C++ storage. Convert them to a
Python list or NumPy array before repeated element-by-element access.

```python
solution = h.getSolution()
col_value = list(solution.col_value)
values = [col_value[i] for i in range(len(col_value))]
```
