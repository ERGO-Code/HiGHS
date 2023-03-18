### Install
HiGHS is available as __highspy__ on [PyPi](https://pypi.org/project/highspy/).

If __highspy__ is not already installed, run

```
pip install highspy
```

### Import

To use __highspy__ within a Python program, it must be imported

```
import highspy
```

When using __highspy__, it is likely that __numpy__ structures will be needed, so must also be imported

```
import numpy as np
```

### Initialize

HiGHS must be initialized before making calls to the HiGHS Python
library

```
h = highspy.Highs()
```

### Methods

Detailed documentation of the methods and structures is given in the [examples section](http://ergo-code.github.io/HiGHS/python/example-py.html).

### Return status

Unless a method just returns data from HiGHS, so is guaranteed to run
successfully, each method returns a status to indicate whether it has
run successfully. This value is an instance of the enum
[HighsStatus](http://ergo-code.github.io/HiGHS/python/enums.html#HighsStatus),
and in the [examples
section](http://ergo-code.github.io/HiGHS/python/example-py.html), it
is referred to as `status`.

### First example

The following Python code reads a problem from the file `model.mps`, and then solves it.

```
import highspy
import numpy as np

# Highs h
h = highspy.Highs()

# Load a model from MPS file model.mps
filename = 'model.mps'
status = h.readModel(filename)
status = h.run()
print('Model ', filename, ' has return status ', h.modelStatusToString(h.getModelStatus))
```
