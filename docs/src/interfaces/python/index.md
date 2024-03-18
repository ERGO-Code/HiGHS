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
