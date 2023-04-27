# [Getting started](@id python-getting-started)

## Install

HiGHS is available as `highspy` on [PyPi](https://pypi.org/project/highspy/).

If `highspy` is not already installed, run:

```bash
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
