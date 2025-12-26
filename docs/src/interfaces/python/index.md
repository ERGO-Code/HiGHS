# [Getting started](@id python-getting-started)

## Install

HiGHS is available as `highspy` on [PyPi](https://pypi.org/project/highspy/).

If `highspy` is not already installed, run:

```shell
$ pip install highspy
```

## Import

To use `highspy` within a Python program, it must be imported

```python
import highspy
```

When using `highspy`, it is likely that `numpy` structures will be needed,
so must also be imported

```python
import numpy as np
```

## Initialize

HiGHS must be initialized before making calls to the HiGHS Python library:

```python
h = highspy.Highs()
```

## Logging

When called from C++, or via the C API, console logging is duplicated
to a file that, by default, is `Highs.log`. However, to channel
logging to a file from `highspy`, the name of the file needs to be
specified explicitly via a call to `setOptionValue('log_file',
'foo.bar')`.

## Methods

Detailed documentation of the methods and structures is given in the
[examples section](@ref example-py).

## Return status

Unless a method just returns data from HiGHS, so is guaranteed to run
successfully, each method returns a status to indicate whether it has run
successfully. This value is an instance of the enum [HighsStatus](@ref), and in
the [examples section](@ref example-py), it is referred to as `status`.

## First example

The following Python code reads a model from the file `model.mps`, and then
solves it.

```python
import highspy

h = highspy.Highs()
filename = 'model.mps'
h.readModel(filename)
h.run()
print('Model ', filename, ' has status ', h.getModelStatus())
```

## Extracting values efficiently

When arrays of values are returned by `highspy`, accessing them
entry-by-entry can be very slow. Such arrays should first be converted
into lists. The following example illustrates how the method
`getSolution()` is used to obtain the solution of a model.

```python
import highspy

h = highspy.Highs()
h.readModel('model.mps')
h.run()

solution = h.getSolution()
num_vars = len(solution.col_value)

value = [solution.col_value[icol]
         for icol in range(num_vars)]

col_value = list(solution.col_value)
value = [col_value[icol]
         for icol in range(num_vars)]
```

For an example LP that is solved in 0.025s, accessing the values
directly from `solution.col_value` takes 0.04s. Forming the list
`col_value` and accessing the values directly from it takes 0.0001s.
