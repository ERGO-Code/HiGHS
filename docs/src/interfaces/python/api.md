<a id="api-reference"></a>

# API Reference

The public `highspy` package re-exports the main solver class, model data
structures, enums, constants, and callback helpers.

<a id="module-highspy"></a>

<a id="top-level-package"></a>

## Top-Level Package

<a id="highspy.BasisValidity"></a>

### *class* highspy.BasisValidity

Bases: `pybind11_object`

Members:

kBasisValidityInvalid

kBasisValidityValid

<a id="highspy.BasisValidity.kBasisValidityInvalid"></a>

#### kBasisValidityInvalid *= <BasisValidity.kBasisValidityInvalid: 0>*

<a id="highspy.BasisValidity.kBasisValidityValid"></a>

#### kBasisValidityValid *= <BasisValidity.kBasisValidityValid: 1>*

### BasisValidity.name -> str

<a id="highspy.BasisValidity.value"></a>

#### *property* value

<a id="highspy.HessianFormat"></a>

### *class* highspy.HessianFormat

Bases: `pybind11_object`

Members:

kTriangular

kSquare

<a id="highspy.HessianFormat.kSquare"></a>

#### kSquare *= <HessianFormat.kSquare: 2>*

<a id="highspy.HessianFormat.kTriangular"></a>

#### kTriangular *= <HessianFormat.kTriangular: 1>*

### HessianFormat.name -> str

<a id="highspy.HessianFormat.value"></a>

#### *property* value

<a id="highspy.Highs"></a>

### *class* highspy.Highs

Bases: `_Highs`

HiGHS solver interface

<a id="highspy.Highs.silent"></a>

#### silent(turn_off_output=True)

Disables solver output to the console.

* **Parameters:**
  **turn_off_output** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))

<a id="highspy.Highs.solve"></a>

#### solve()

Runs the solver on the current problem.

* **Returns:**
  A HighsStatus object containing the solve status.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.Highs.startSolve"></a>

#### startSolve()

Starts the solver in a separate thread.  Useful for handling KeyboardInterrupts.
Do not attempt to modify the model while the solver is running.

* **Returns:**
  A Thread object representing the solver thread.
* **Return type:**
  [*Thread*](https://docs.python.org/3/library/threading.html#threading.Thread)

<a id="highspy.Highs.is_solver_running"></a>

#### is_solver_running()

* **Return type:**
  [bool](https://docs.python.org/3/builtins/functions.html#bool)

<a id="highspy.Highs.joinSolve"></a>

#### joinSolve(solver_thread=None, interrupt_limit=5)

Waits for the solver to finish. If solver_thread is provided, it will handle KeyboardInterrupts.

* **Parameters:**
  * **solver_thread** ([*Thread*](https://docs.python.org/3/library/threading.html#threading.Thread) *|* *None*) – A Thread object representing the solver thread (optional).
  * **interrupt_limit** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – The number of times to allow KeyboardInterrupt before forcing termination (optional).
* **Returns:**
  A HighsStatus object containing the solve status.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.Highs.wait"></a>

#### wait(timeout=-1.0)

* **Parameters:**
  **timeout** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[bool](https://docs.python.org/3/builtins/functions.html#bool), [*HighsStatus*](#highspy._core.HighsStatus) | None]

<a id="highspy.Highs.optimize"></a>

#### optimize()

Alias for the solve method.

* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.Highs.getObjective"></a>

#### getObjective()

Retrieves the current objective function (as a linear expression) and sense.

* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[*highs_linear_expression*](#highspy.highs.highs_linear_expression), [*ObjSense*](#highspy._core.ObjSense)]

<a id="highspy.Highs.setObjective"></a>

#### setObjective(obj=None, sense=None)

Updates the costs.

* **Parameters:**
  * **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
  * **sense** ([*ObjSense*](#highspy._core.ObjSense) *|* *None*) – An optional ObjSense value representing the new objective sense.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.

<a id="highspy.Highs.minimize"></a>

#### minimize(obj=None)

Solves a minimization of the objective and optionally updates the costs.

* **Parameters:**
  **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.
* **Returns:**
  A HighsStatus object containing the solve status after minimization.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.Highs.maximize"></a>

#### maximize(obj=None)

Solves a maximization of the objective and optionally updates the costs.

* **Parameters:**
  **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.
* **Returns:**
  A HighsStatus object containing the solve status after maximization.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.Highs.internal_get_value"></a>

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Internal method to get the value of an index from an array of values. Could be value or dual, variable or constraint.

* **Parameters:**
  * **array_values** ([*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*float*](https://docs.python.org/3/builtins/functions.html#float) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[**float64* *]* *]*)
  * **index_collection** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray))
* **Return type:**
  [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool) | [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*Any*](https://docs.python.org/3/library/typing.html#typing.Any)] | [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*float64*]]

<a id="highspy.Highs.val"></a>

#### val(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### val(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### val(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### val(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### val(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Gets the value of a variable/index or expression in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var/index or highs_linear_expression object representing the variable.
* **Returns:**
  The value of the variable in the solution.

<a id="highspy.Highs.vals"></a>

#### vals(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### vals(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### vals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### vals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### vals(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Gets the values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.

<a id="highspy.Highs.variableName"></a>

#### variableName(var)

Retrieves the name of a specific variable.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var)) – A highs_var object representing the variable.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If the variable name cannot be found.
* **Returns:**
  The name of the specified variable.

<a id="highspy.Highs.variableNames"></a>

#### variableNames(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var) | [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

#### variableNames(idxs: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[highs_var](#highspy.highs.highs_var) | [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral)]) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

Retrieves the names of multiple variables.

* **Parameters:**
  **idxs** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *]*) – An iterable of highs_var objects or a mapping where keys are identifiers and values are highs_var objects.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If any variable name cannot be found.
* **Returns:**
  If idxs is a mapping, returns a dict where keys are the same keys from the input idxs and values are the names of the corresponding variables.
  If idxs is an iterable, returns a list of names for the specified variables.

<a id="highspy.Highs.allVariableNames"></a>

#### allVariableNames()

Retrieves the names of all variables in the model.

* **Returns:**
  A list of strings representing the names of all variables.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.Highs.variableValue"></a>

#### variableValue(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableValue(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableValue(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableValue(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableValue(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the value of a specific variable in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var object representing the variable.
* **Returns:**
  The value of the specified variable in the solution.

<a id="highspy.Highs.variableValues"></a>

#### variableValues(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableValues(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableValues(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableValues(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableValues(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.

<a id="highspy.Highs.allVariableValues"></a>

#### allVariableValues()

Retrieves the values of all variables in the solution.

* **Returns:**
  A list of values for all variables in the solution.

<a id="highspy.Highs.variableDual"></a>

#### variableDual(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableDual(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableDual(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableDual(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableDual(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual value of a specific variable/index or expression in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var object representing the variable.
* **Returns:**
  The dual value of the specified variable in the solution.

<a id="highspy.Highs.variableDuals"></a>

#### variableDuals(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableDuals(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableDuals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableDuals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableDuals(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the dual values of the corresponding variables. If idxs is an iterable, returns a list of dual values for the variables.

<a id="highspy.Highs.allVariableDuals"></a>

#### allVariableDuals()

Retrieves the dual values of all variables in the solution.

* **Returns:**
  A list of dual values for all variables in the solution.

<a id="highspy.Highs.constrValue"></a>

#### constrValue(con: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrValue(con: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrValue(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrValue(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrValue(con: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the value of a specific constraint in the solution.

* **Parameters:**
  **con** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_con object representing the constraint.
* **Returns:**
  The value of the specified constraint in the solution.

<a id="highspy.Highs.constrValues"></a>

#### constrValues(cons: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrValues(cons: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrValues(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrValues(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrValues(cons: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the values of multiple constraints in the solution.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.
* **Returns:**
  If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.

<a id="highspy.Highs.allConstrValues"></a>

#### allConstrValues()

Retrieves the values of all constraints in the solution.

* **Returns:**
  A list of values for all constraints in the solution.

<a id="highspy.Highs.constrDual"></a>

#### constrDual(con: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrDual(con: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrDual(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrDual(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrDual(con: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual value of a specific constraint in the solution.

* **Parameters:**
  **con** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_con object representing the constraint.
* **Returns:**
  The dual value of the specified constraint in the solution.

<a id="highspy.Highs.constrDuals"></a>

#### constrDuals(cons: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrDuals(cons: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrDuals(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrDuals(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrDuals(cons: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual values of multiple constraints in the solution.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.
* **Returns:**
  If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.

<a id="highspy.Highs.allConstrDuals"></a>

#### allConstrDuals()

Retrieves the dual values of all constraints in the solution.

* **Returns:**
  A list of dual values for all constraints in the solution.

<a id="highspy.Highs.addVariable"></a>

#### addVariable(lb=0, ub=inf, obj=0.0, type=<HighsVarType.kContinuous: 0>, name=None)

Adds a variable to the model.

* **Parameters:**
  * **lb** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Lower bound of the variable (default is 0).
  * **ub** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Upper bound of the variable (default is infinity).
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Objective coefficient of the variable (default is 0).
  * **type** ([*HighsVarType*](#highspy._core.HighsVarType)) – Type of the variable (continuous, integer; default is continuous).
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – Optional name for the variable.
* **Returns:**
  A highs_var object representing the added variable.

<a id="highspy.Highs.addVariables"></a>

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addVariables(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addVariables(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Adds multiple variables to the model.

* **Parameters:**
  * **\*args** – A sequence of variables to be added. Can be a collection of scalars or indices (or mix).
  * **\*\*kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*) – Optional keyword arguments.  Can be scalars, arrays, or mappings.
    lb: Lower bound of the variables (default is 0).
    ub: Upper bound of the variables (default is infinity).
    obj: Objective coefficient of the variables (default is 0).
    type: Type of the variables (continuous, integer; default is continuous).
    name: A collection of names for the variables (list or mapping).
    name_prefix: Prefix for the variable names.  Constructed name will be name_prefix + index.
    out_array: Return an array of highs_var objects instead of a dictionary.
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **\*\*kwargs**
* **Returns:**
  A highs_var collection (array or dictionary) representing the added variables.
* **Return type:**
  [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*highs_var*](#highspy.highs.highs_var)] | [*HighspyArray*](#highspy.highs.HighspyArray) | None

<a id="highspy.Highs.addIntegrals"></a>

#### addIntegrals(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addIntegrals(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addIntegrals(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addIntegrals(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Alias for the addVariables method, for integer variables.

* **Parameters:**
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)

<a id="highspy.Highs.addBinaries"></a>

#### addBinaries(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addBinaries(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addBinaries(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addBinaries(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Alias for the addVariables method, for binary variables.

* **Parameters:**
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)

<a id="highspy.Highs.addIntegral"></a>

#### addIntegral(lb=0.0, ub=inf, obj=0.0, name=None)

Alias for the addVariable method, for integer variables.

* **Parameters:**
  * **lb** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **ub** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*)

<a id="highspy.Highs.addBinary"></a>

#### addBinary(obj=0.0, name=None)

Alias for the addVariable method, for binary variables.

* **Parameters:**
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*)

<a id="highspy.Highs.deleteVariable"></a>

#### deleteVariable(var_or_index, \*args)

Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

* **Parameters:**
  * **var_or_index** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – A highs_var object or an index representing the variable to be deleted.
  * **\*args** ([*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*highs_var*](#highspy.highs.highs_var) *]*  *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*highs_var*](#highspy.highs.highs_var) *|* [*HighspyArray*](#highspy.highs.HighspyArray) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – Optional collections (lists, dicts, etc.) of highs_var objects whose indices need to be updated.

<a id="highspy.Highs.getVariables"></a>

#### getVariables()

Retrieves all variables in the model.

* **Returns:**
  A list of highs_var objects, each representing a variable in the model.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[*highs_var*](#highspy.highs.highs_var)]

<a id="highspy.Highs.inf"></a>

#### *property* inf *: [float](https://docs.python.org/3/builtins/functions.html#float)*

Represents infinity in the context of the solver.

* **Returns:**
  The value used to represent infinity.

<a id="highspy.Highs.numVariables"></a>

#### *property* numVariables *: [int](https://docs.python.org/3/builtins/functions.html#int)*

Gets the number of variables in the model.

* **Returns:**
  The number of variables.

<a id="highspy.Highs.numConstrs"></a>

#### *property* numConstrs *: [int](https://docs.python.org/3/builtins/functions.html#int)*

Gets the number of constraints in the model.

* **Returns:**
  The number of constraints.

<a id="highspy.Highs.addConstr"></a>

#### addConstr(expr, name=None)

Adds a constraint to the model.

* **Parameters:**
  * **expr** ([*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – A highs_linear_expression to be added.
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – Optional name of the constraint.
* **Returns:**
  A highs_cons object representing the added constraint.
* **Return type:**
  [*highs_cons*](#highspy.highs.highs_cons)

<a id="highspy.Highs.addConstrs"></a>

#### addConstrs(\*args: [highs_linear_expression](#highspy.highs.highs_linear_expression), \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_linear_expression](#highspy.highs.highs_linear_expression)], \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[highs_linear_expression](#highspy.highs.highs_linear_expression)], \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [HighspyArray](#highspy.highs.HighspyArray), \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

Adds multiple constraints to the model.

* **Parameters:**
  * **\*args** ([*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A sequence of highs_linear_expression to be added.
  * **\*\*kwargs** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *]*  *|* *None*) – Optional keyword arguments.
    name_prefix: Prefix for the constraint names.  Constructed name will be name_prefix + index.
    name: A collection of names for the constraints (list or mapping).
* **Returns:**
  A highs_con collection array representing the added constraints.

<a id="highspy.Highs.expr"></a>

#### expr(optional=None)

Creates a new highs_linear_expression object.

* **Returns:**
  A highs_linear_expression object.
* **Parameters:**
  **optional** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.Highs.getExpr"></a>

#### getExpr(cons)

Retrieves the highs_linear_expression of a constraint.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_cons*](#highspy.highs.highs_cons)) – A highs_con object or index representing the constraint.
* **Returns:**
  A highs_linear_expression object representing the expression of the constraint.
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.Highs.chgCoeff"></a>

#### chgCoeff(cons, var, val)

Changes the coefficient of a variable in a constraint.

* **Parameters:**
  * **cons** ([*highs_cons*](#highspy.highs.highs_cons) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_con object representing the constraint.
  * **var** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_var object representing the variable.
  * **val** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – The new coefficient value for the variable in the constraint.

<a id="highspy.Highs.getConstrs"></a>

#### getConstrs()

Retrieves all constraints in the model.

* **Returns:**
  A list of highs_cons objects, each representing a constraint in the model.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[*highs_cons*](#highspy.highs.highs_cons)]

<a id="highspy.Highs.removeConstr"></a>

#### removeConstr(cons_or_index, \*args)

Removes a constraint from the model and updates the indices of subsequent constraints in provided collections.

* **Parameters:**
  * **cons_or_index** ([*highs_cons*](#highspy.highs.highs_cons) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_cons object or an index representing the constraint to be removed.
  * **\*args** ([*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*highs_cons*](#highspy.highs.highs_cons) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*highs_cons*](#highspy.highs.highs_cons) *]*  *|* [*highs_cons*](#highspy.highs.highs_cons)) – Optional collections (lists, dicts, etc.) of highs_cons objects whose indices need to be updated after the removal.

<a id="highspy.Highs.setMinimize"></a>

#### setMinimize()

Sets the objective sense of the model to minimization.

<a id="highspy.Highs.setMaximize"></a>

#### setMaximize()

Sets the objective sense of the model to maximization.

<a id="highspy.Highs.setInteger"></a>

#### setInteger(var_or_collection)

Sets a variable/collection to integer.

* **Parameters:**
  **var_or_collection** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A highs_var object/collection representing the variable to be set as integer.

<a id="highspy.Highs.setContinuous"></a>

#### setContinuous(var_or_collection)

Sets a variable/collection to continuous.

* **Parameters:**
  **var_or_collection** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A highs_var object/collection representing the variable to be set as continuous.

<a id="highspy.Highs.idx"></a>

#### *static* idx(\*args)

Convert highs_var/highs_cons to a flat int32 index array.

Can be called as:
: - `h.idx(array)` with a HighspyArray, numpy array, list, or tuple
  - `h.idx(a, b, c)` with individual highs_var or highs_cons objects

* **Returns:**
  A flat int32 numpy array of the underlying indices.
* **Return type:**
  [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*int32*]]

<a id="highspy.Highs.qsum"></a>

#### *static* qsum(items, initial=None)

Performs a faster sum for highs_linear_expressions.

* **Parameters:**
  * **items** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[*[*object_*](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_) *]* *]*) – A collection of highs_linear_expressions or highs_vars to be summed.
  * **initial** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.Highs.enableCallbacks"></a>

#### enableCallbacks()

Enables callbacks, restarting them if they were previously enabled.

<a id="highspy.Highs.clearCallbacks"></a>

#### clearCallbacks()

Clears all callbacks.

<a id="highspy.Highs.disableCallbacks"></a>

#### disableCallbacks()

Disables all callbacks, but does not clear them.

<a id="highspy.Highs.cancelSolve"></a>

#### cancelSolve()

If HandleUserInterrupt is enabled, this method will signal the solver to stop.

<a id="highspy.Highs.HandleKeyboardInterrupt"></a>

#### *property* HandleKeyboardInterrupt *: [bool](https://docs.python.org/3/builtins/functions.html#bool)*

Get/Set whether the solver should handle KeyboardInterrupt (i.e., cancel solve on Ctrl+C). Also enables/disables HandleUserInterrupt.

<a id="highspy.Highs.HandleUserInterrupt"></a>

#### *property* HandleUserInterrupt *: [bool](https://docs.python.org/3/builtins/functions.html#bool)*

Get/Set whether the solver should handle user interrupts (i.e., cancel solve on user request)

<a id="highspy.Highs.cbLogging"></a>

#### cbLogging *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbSimplexInterrupt"></a>

#### cbSimplexInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbIpmInterrupt"></a>

#### cbIpmInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipSolution"></a>

#### cbMipSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipImprovingSolution"></a>

#### cbMipImprovingSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipLogging"></a>

#### cbMipLogging *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipInterrupt"></a>

#### cbMipInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipGetCutPool"></a>

#### cbMipGetCutPool *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipDefineLazyConstraints"></a>

#### cbMipDefineLazyConstraints *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.cbMipUserSolution"></a>

#### cbMipUserSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.Highs.addCol"></a>

#### addCol(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg3: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addCols"></a>

#### addCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg4: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg6: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg7: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addLinearObjective"></a>

#### addLinearObjective(self: highspy._core._Highs, arg0: [HighsLinearObjective](#highspy.HighsLinearObjective)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addRow"></a>

#### addRow(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addRows"></a>

#### addRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg6: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addVar"></a>

#### addVar(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.addVars"></a>

#### addVars(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.basisStatusToString"></a>

#### basisStatusToString(self: highspy._core._Highs, arg0: [highspy._core.HighsBasisStatus](#highspy._core.HighsBasisStatus)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.basisValidityToString"></a>

#### basisValidityToString(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.changeCoeff"></a>

#### changeCoeff(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColBounds"></a>

#### changeColBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColCost"></a>

#### changeColCost(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColIntegrality"></a>

#### changeColIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [highspy._core.HighsVarType](#highspy._core.HighsVarType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColsBounds"></a>

#### changeColsBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColsCost"></a>

#### changeColsCost(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeColsIntegrality"></a>

#### changeColsIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.uint8]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeObjectiveOffset"></a>

#### changeObjectiveOffset(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeObjectiveSense"></a>

#### changeObjectiveSense(self: highspy._core._Highs, arg0: [highspy._core.ObjSense](#highspy._core.ObjSense)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeRowBounds"></a>

#### changeRowBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.changeRowsBounds"></a>

#### changeRowsBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.clear"></a>

#### clear(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.clearLinearObjectives"></a>

#### clearLinearObjectives(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.clearModel"></a>

#### clearModel(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.clearSolver"></a>

#### clearSolver(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.crossover"></a>

#### crossover(self: highspy._core._Highs, arg0: [HighsSolution](#highspy.HighsSolution)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.deleteCols"></a>

#### deleteCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.deleteRows"></a>

#### deleteRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.deleteVars"></a>

#### deleteVars(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.ensureColwise"></a>

#### ensureColwise(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.ensureRowwise"></a>

#### ensureRowwise(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.feasibilityRelaxation"></a>

#### feasibilityRelaxation(self: highspy._core._Highs, global_lower_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), global_upper_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), global_rhs_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), local_lower_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None, local_upper_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None, local_rhs_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.getBasicVariables"></a>

#### getBasicVariables(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getBasis"></a>

#### getBasis(self: highspy._core._Highs) → [HighsBasis](#highspy.HighsBasis)

<a id="highspy.Highs.getBasisInverseCol"></a>

#### getBasisInverseCol(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getBasisInverseColSparse"></a>

#### getBasisInverseColSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getBasisInverseRow"></a>

#### getBasisInverseRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getBasisInverseRowSparse"></a>

#### getBasisInverseRowSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getBasisSolve"></a>

#### getBasisSolve(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getBasisSolveSparse"></a>

#### getBasisSolveSparse(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getBasisTransposeSolve"></a>

#### getBasisTransposeSolve(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getBasisTransposeSolveSparse"></a>

#### getBasisTransposeSolveSparse(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getCol"></a>

#### getCol(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getColByName"></a>

#### getColByName(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getColEntries"></a>

#### getColEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getColIntegrality"></a>

#### getColIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsVarType](#highspy._core.HighsVarType)]

<a id="highspy.Highs.getColName"></a>

#### getColName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.Highs.getCols"></a>

#### getCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getColsEntries"></a>

#### getColsEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getDualObjectiveValue"></a>

#### getDualObjectiveValue(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.getDualRay"></a>

#### getDualRay(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getDualRayExist"></a>

#### getDualRayExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.Highs.getDualUnboundednessDirection"></a>

#### getDualUnboundednessDirection(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getDualUnboundednessDirectionExist"></a>

#### getDualUnboundednessDirectionExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.Highs.getFixedLp"></a>

#### getFixedLp(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsLp](#highspy._core.HighsLp)]

<a id="highspy.Highs.getHessianNumNz"></a>

#### getHessianNumNz(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.getIis"></a>

#### getIis(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [HighsIis](#highspy.HighsIis)]

<a id="highspy.Highs.getInfinity"></a>

#### getInfinity(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.Highs.getInfo"></a>

#### getInfo(self: highspy._core._Highs) → [highspy._core.HighsInfo](#highspy._core.HighsInfo)

<a id="highspy.Highs.getInfoType"></a>

#### getInfoType(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsInfoType](#highspy._core.HighsInfoType)]

<a id="highspy.Highs.getInfoValue"></a>

#### getInfoValue(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [object](https://docs.python.org/3/builtins/functions.html#object)]

<a id="highspy.Highs.getLinearObjective"></a>

#### getLinearObjective(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [HighsLinearObjective](#highspy.HighsLinearObjective)

<a id="highspy.Highs.getLp"></a>

#### getLp(self: highspy._core._Highs) → [highspy._core.HighsLp](#highspy._core.HighsLp)

<a id="highspy.Highs.getModel"></a>

#### getModel(self: highspy._core._Highs) → [highspy._core.HighsModel](#highspy._core.HighsModel)

<a id="highspy.Highs.getModelPresolveStatus"></a>

#### getModelPresolveStatus(self: highspy._core._Highs) → [highspy._core.HighsPresolveStatus](#highspy._core.HighsPresolveStatus)

<a id="highspy.Highs.getModelStatus"></a>

#### getModelStatus(self: highspy._core._Highs) → [highspy._core.HighsModelStatus](#highspy._core.HighsModelStatus)

<a id="highspy.Highs.getNumCol"></a>

#### getNumCol(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.getNumLinearObjectives"></a>

#### getNumLinearObjectives(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.getNumNz"></a>

#### getNumNz(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.getNumRow"></a>

#### getNumRow(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.getObjectiveOffset"></a>

#### getObjectiveOffset(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float)]

<a id="highspy.Highs.getObjectiveSense"></a>

#### getObjectiveSense(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.ObjSense](#highspy._core.ObjSense)]

<a id="highspy.Highs.getObjectiveValue"></a>

#### getObjectiveValue(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.Highs.getOptionType"></a>

#### getOptionType(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsOptionType](#highspy._core.HighsOptionType)]

<a id="highspy.Highs.getOptionValue"></a>

#### getOptionValue(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [object](https://docs.python.org/3/builtins/functions.html#object)]

<a id="highspy.Highs.getOptions"></a>

#### getOptions(self: highspy._core._Highs) → [highspy._core.HighsOptions](#highspy._core.HighsOptions)

<a id="highspy.Highs.getPresolvedLp"></a>

#### getPresolvedLp(self: highspy._core._Highs) → [highspy._core.HighsLp](#highspy._core.HighsLp)

<a id="highspy.Highs.getPrimalRay"></a>

#### getPrimalRay(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getPrimalRayExist"></a>

#### getPrimalRayExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.Highs.getRanging"></a>

#### getRanging(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [HighsRanging](#highspy.HighsRanging)]

<a id="highspy.Highs.getReducedColumn"></a>

#### getReducedColumn(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getReducedColumnSparse"></a>

#### getReducedColumnSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getReducedRow"></a>

#### getReducedRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getReducedRowSparse"></a>

#### getReducedRowSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.Highs.getRow"></a>

#### getRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getRowByName"></a>

#### getRowByName(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getRowEntries"></a>

#### getRowEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getRowName"></a>

#### getRowName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.Highs.getRows"></a>

#### getRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.Highs.getRowsEntries"></a>

#### getRowsEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.Highs.getRunTime"></a>

#### getRunTime(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.Highs.getSavedMipSolutions"></a>

#### getSavedMipSolutions(self: highspy._core._Highs) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[HighsObjectiveSolution](#highspy.HighsObjectiveSolution)]

<a id="highspy.Highs.getSolution"></a>

#### getSolution(self: highspy._core._Highs) → [HighsSolution](#highspy.HighsSolution)

<a id="highspy.Highs.getThirdPartyNotice"></a>

#### getThirdPartyNotice(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.githash"></a>

#### githash(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.modelStatusToString"></a>

#### modelStatusToString(self: highspy._core._Highs, arg0: [highspy._core.HighsModelStatus](#highspy._core.HighsModelStatus)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.passColName"></a>

#### passColName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.passHessian"></a>

#### passHessian(\*args, \*\*kwargs)

Overloaded function.

1. passHessian(self: highspy._core._Highs, arg0: highspy._core.HighsHessian) -> highspy._core.HighsStatus
2. passHessian(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg4: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg5: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus

<a id="highspy.Highs.passModel"></a>

#### passModel(\*args, \*\*kwargs)

Overloaded function.

1. passModel(self: highspy._core._Highs, arg0: highspy._core.HighsModel) -> highspy._core.HighsStatus
2. passModel(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.SupportsInt | typing.SupportsIndex, arg4: typing.SupportsInt | typing.SupportsIndex, arg5: typing.SupportsInt | typing.SupportsIndex, arg6: typing.SupportsInt | typing.SupportsIndex, arg7: typing.SupportsFloat | typing.SupportsIndex, arg8: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg9: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg10: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg11: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg12: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg13: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg14: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg15: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg16: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg17: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg18: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg19: typing.Annotated[numpy.typing.ArrayLike, numpy.int32]) -> highspy._core.HighsStatus
3. passModel(self: highspy._core._Highs, arg0: highspy._core.HighsLp) -> highspy._core.HighsStatus
4. passModel(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.SupportsInt | typing.SupportsIndex, arg4: typing.SupportsInt | typing.SupportsIndex, arg5: typing.SupportsFloat | typing.SupportsIndex, arg6: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg7: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg8: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg9: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg10: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg11: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg12: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg13: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg14: typing.Annotated[numpy.typing.ArrayLike, numpy.int32]) -> highspy._core.HighsStatus

<a id="highspy.Highs.passOptions"></a>

#### passOptions(self: highspy._core._Highs, arg0: [highspy._core.HighsOptions](#highspy._core.HighsOptions)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.passRowName"></a>

#### passRowName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.postsolve"></a>

#### postsolve(\*args, \*\*kwargs)

Overloaded function.

1. postsolve(self: highspy._core._Highs, arg0: HighsSolution, arg1: HighsBasis) -> highspy._core.HighsStatus
2. postsolve(self: highspy._core._Highs, arg0: HighsSolution) -> highspy._core.HighsStatus

<a id="highspy.Highs.presolve"></a>

#### presolve(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.readBasis"></a>

#### readBasis(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.readModel"></a>

#### readModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.readOptions"></a>

#### readOptions(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.readSolution"></a>

#### readSolution(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.releaseMemory"></a>

#### releaseMemory(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.resetGlobalScheduler"></a>

#### *static* resetGlobalScheduler(arg0: [bool](https://docs.python.org/3/builtins/functions.html#bool)) → [None](https://docs.python.org/3/builtins/constants.html#None)

<a id="highspy.Highs.resetOptions"></a>

#### resetOptions(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.run"></a>

#### run(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.setBasis"></a>

#### setBasis(\*args, \*\*kwargs)

Overloaded function.

1. setBasis(self: highspy._core._Highs, arg0: HighsBasis) -> highspy._core.HighsStatus
2. setBasis(self: highspy._core._Highs) -> highspy._core.HighsStatus

<a id="highspy.Highs.setCallback"></a>

#### setCallback(self: highspy._core._Highs, arg0: [collections.abc.Callable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable)[[[int](https://docs.python.org/3/builtins/functions.html#int), [str](https://docs.python.org/3/builtins/stdtypes.html#str), [HighsCallbackOutput](#highspy._core.cb.HighsCallbackOutput), [HighsCallbackInput](#highspy._core.cb.HighsCallbackInput), [object](https://docs.python.org/3/builtins/functions.html#object)], [None](https://docs.python.org/3/builtins/constants.html#None)], arg1: [object](https://docs.python.org/3/builtins/functions.html#object)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.setOptionValue"></a>

#### setOptionValue(\*args, \*\*kwargs)

Overloaded function.

1. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: bool) -> highspy._core.HighsStatus
2. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: typing.SupportsInt | typing.SupportsIndex) -> highspy._core.HighsStatus
3. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: typing.SupportsFloat | typing.SupportsIndex) -> highspy._core.HighsStatus
4. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: str) -> highspy._core.HighsStatus

<a id="highspy.Highs.setSolution"></a>

#### setSolution(\*args, \*\*kwargs)

Overloaded function.

1. setSolution(self: highspy._core._Highs, arg0: HighsSolution) -> highspy._core.HighsStatus
2. setSolution(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus

<a id="highspy.Highs.solutionStatusToString"></a>

#### solutionStatusToString(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.startCallback"></a>

#### startCallback(self: highspy._core._Highs, arg0: [HighsCallbackType](#highspy._core.cb.HighsCallbackType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.startCallbackInt"></a>

#### startCallbackInt(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.stopCallback"></a>

#### stopCallback(self: highspy._core._Highs, arg0: [HighsCallbackType](#highspy._core.cb.HighsCallbackType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.stopCallbackInt"></a>

#### stopCallbackInt(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.version"></a>

#### version(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.Highs.versionMajor"></a>

#### versionMajor(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.versionMinor"></a>

#### versionMinor(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.versionPatch"></a>

#### versionPatch(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.Highs.writeBasis"></a>

#### writeBasis(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writeIisModel"></a>

#### writeIisModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writeInfo"></a>

#### writeInfo(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writeModel"></a>

#### writeModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writeOptions"></a>

#### writeOptions(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writePresolvedModel"></a>

#### writePresolvedModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.Highs.writeSolution"></a>

#### writeSolution(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.HighsBasis"></a>

### *class* highspy.HighsBasis

Bases: `pybind11_object`

<a id="highspy.HighsBasis.alien"></a>

#### *property* alien

<a id="highspy.HighsBasis.col_status"></a>

#### *property* col_status

<a id="highspy.HighsBasis.debug_id"></a>

#### *property* debug_id

<a id="highspy.HighsBasis.debug_origin_name"></a>

#### *property* debug_origin_name

<a id="highspy.HighsBasis.debug_update_count"></a>

#### *property* debug_update_count

<a id="highspy.HighsBasis.row_status"></a>

#### *property* row_status

<a id="highspy.HighsBasis.valid"></a>

#### *property* valid

<a id="highspy.HighsBasis.was_alien"></a>

#### *property* was_alien

<a id="highspy.HighsBasisStatus"></a>

### *class* highspy.HighsBasisStatus

Bases: `pybind11_object`

Members:

kLower

kBasic

kUpper

kZero

kNonbasic

<a id="highspy.HighsBasisStatus.kBasic"></a>

#### kBasic *= <HighsBasisStatus.kBasic: 1>*

<a id="highspy.HighsBasisStatus.kLower"></a>

#### kLower *= <HighsBasisStatus.kLower: 0>*

<a id="highspy.HighsBasisStatus.kNonbasic"></a>

#### kNonbasic *= <HighsBasisStatus.kNonbasic: 4>*

<a id="highspy.HighsBasisStatus.kUpper"></a>

#### kUpper *= <HighsBasisStatus.kUpper: 2>*

<a id="highspy.HighsBasisStatus.kZero"></a>

#### kZero *= <HighsBasisStatus.kZero: 3>*

### HighsBasisStatus.name -> str

<a id="highspy.HighsBasisStatus.value"></a>

#### *property* value

<a id="highspy.HighsCallback"></a>

### *class* highspy.HighsCallback(callback_type, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

* **Parameters:**
  * **callback_type** ([*cb.HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **highs** ([*Highs*](#highspy.Highs))

<a id="highspy.HighsCallback.callbacks"></a>

#### callbacks *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[Callable](https://docs.python.org/3/library/typing.html#typing.Callable)[[[HighsCallbackEvent](#highspy.highs.HighsCallbackEvent)], [None](https://docs.python.org/3/builtins/constants.html#None)]]*

<a id="highspy.HighsCallback.user_callback_data"></a>

#### user_callback_data *: [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Callable](https://docs.python.org/3/library/typing.html#typing.Callable)[[[HighsCallbackEvent](#highspy.highs.HighsCallbackEvent)], [None](https://docs.python.org/3/builtins/constants.html#None)], [Any](https://docs.python.org/3/library/typing.html#typing.Any)]*

<a id="highspy.HighsCallback.callback_type"></a>

#### callback_type

<a id="highspy.HighsCallback.highs"></a>

#### highs

<a id="highspy.HighsCallback.subscribe"></a>

#### subscribe(callback, user_data=None)

Subscribes a callback to the event.

* **Parameters:**
  * **callback** ([*Callable*](https://docs.python.org/3/library/typing.html#typing.Callable) *[* *[*[*HighsCallbackEvent*](#highspy.highs.HighsCallbackEvent) *]* *,* *None* *]*) – The callback function to be executed.
  * **user_data** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*) – Optional user data to be passed to the callback.

<a id="highspy.HighsCallback.unsubscribe"></a>

#### unsubscribe(callback)

Unsubscribes a callback from the event.

* **Parameters:**
  **callback** ([*Callable*](https://docs.python.org/3/library/typing.html#typing.Callable) *[* *[*[*HighsCallbackEvent*](#highspy.highs.HighsCallbackEvent) *]* *,* *None* *]*) – The callback function to be removed.

<a id="highspy.HighsCallback.unsubscribe_by_data"></a>

#### unsubscribe_by_data(user_data)

Unsubscribes a callback by user data.

* **Parameters:**
  **user_data** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*) – The user data corresponding to the callback(s) to be removed.

<a id="highspy.HighsCallback.clear"></a>

#### clear()

Unsubscribes all callbacks from the event.

<a id="highspy.HighsCallback.fire"></a>

#### fire(callback_type, message, data_out, data_in)

Fires the event, executing all subscribed callbacks.

* **Parameters:**
  * **callback_type** ([*HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **message** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **data_out** ([*HighsCallbackOutput*](#highspy._core.cb.HighsCallbackOutput))
  * **data_in** ([*HighsCallbackInput*](#highspy._core.cb.HighsCallbackInput))

<a id="highspy.HighsCallbackEvent"></a>

### *class* highspy.HighsCallbackEvent(callback_type, message, data_out, data_in, user_data)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

* **Parameters:**
  * **callback_type** ([*cb.HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **message** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **data_out** ([*cb.HighsCallbackOutput*](#highspy._core.cb.HighsCallbackOutput))
  * **data_in** ([*cb.HighsCallbackInput*](#highspy._core.cb.HighsCallbackInput) *|* *None*)
  * **user_data** (*Any* *|* *None*)

<a id="highspy.HighsCallbackEvent.callback_type"></a>

#### callback_type

<a id="highspy.HighsCallbackEvent.message"></a>

#### message

<a id="highspy.HighsCallbackEvent.data_out"></a>

#### data_out

<a id="highspy.HighsCallbackEvent.data_in"></a>

#### data_in

<a id="highspy.HighsCallbackEvent.user_data"></a>

#### user_data

<a id="highspy.HighsCallbackEvent.interrupt"></a>

#### interrupt(interrupt_value=True)

Sets the user interrupt flag in the callback data.

* **Parameters:**
  **interrupt_value** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))

<a id="highspy.HighsCallbackEvent.val"></a>

#### val(var_expr)

Gets the value(s) of a variable/index or expression in the callback solution.

* **Parameters:**
  **var_expr** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray))

<a id="highspy.HighsCallbackEvent.cut"></a>

#### cut(index)

Gets the cut pool for the given index.

* **Parameters:**
  **index** ([*int*](https://docs.python.org/3/builtins/functions.html#int))

<a id="highspy.HighsCallbackEvent.cuts"></a>

#### *property* cuts

Gets all cuts in the cut pool.

<a id="highspy.HighsError"></a>

### *exception* highspy.HighsError

Bases: [`Exception`](https://docs.python.org/3/builtins/exceptions.html#Exception)

Base exception for highspy errors.

<a id="highspy.HighsHessian"></a>

### *class* highspy.HighsHessian

Bases: `pybind11_object`

<a id="highspy.HighsHessian.dim_"></a>

#### *property* dim_

<a id="highspy.HighsHessian.format_"></a>

#### *property* format_

<a id="highspy.HighsHessian.index_"></a>

#### *property* index_

<a id="highspy.HighsHessian.start_"></a>

#### *property* start_

<a id="highspy.HighsHessian.value_"></a>

#### *property* value_

<a id="highspy.HighsIis"></a>

### *class* highspy.HighsIis

Bases: `pybind11_object`

<a id="highspy.HighsIis.clear"></a>

#### clear(self: [highspy._core.HighsIis](#highspy._core.HighsIis)) → [None](https://docs.python.org/3/builtins/constants.html#None)

<a id="highspy.HighsIis.col_bound_"></a>

#### *property* col_bound_

<a id="highspy.HighsIis.col_index_"></a>

#### *property* col_index_

<a id="highspy.HighsIis.col_status_"></a>

#### *property* col_status_

<a id="highspy.HighsIis.info_"></a>

#### *property* info_

<a id="highspy.HighsIis.model_"></a>

#### *property* model_

<a id="highspy.HighsIis.row_bound_"></a>

#### *property* row_bound_

<a id="highspy.HighsIis.row_index_"></a>

#### *property* row_index_

<a id="highspy.HighsIis.row_status_"></a>

#### *property* row_status_

<a id="highspy.HighsIis.status_"></a>

#### *property* status_

<a id="highspy.HighsIis.strategy_"></a>

#### *property* strategy_

<a id="highspy.HighsIis.valid_"></a>

#### *property* valid_

<a id="highspy.HighsInfo"></a>

### *class* highspy.HighsInfo

Bases: `pybind11_object`

<a id="highspy.HighsInfo.basis_validity"></a>

#### *property* basis_validity

<a id="highspy.HighsInfo.crossover_iteration_count"></a>

#### *property* crossover_iteration_count

<a id="highspy.HighsInfo.dual_solution_status"></a>

#### *property* dual_solution_status

<a id="highspy.HighsInfo.ipm_iteration_count"></a>

#### *property* ipm_iteration_count

<a id="highspy.HighsInfo.max_complementarity_violation"></a>

#### *property* max_complementarity_violation

<a id="highspy.HighsInfo.max_dual_infeasibility"></a>

#### *property* max_dual_infeasibility

<a id="highspy.HighsInfo.max_dual_residual_error"></a>

#### *property* max_dual_residual_error

<a id="highspy.HighsInfo.max_integrality_violation"></a>

#### *property* max_integrality_violation

<a id="highspy.HighsInfo.max_primal_infeasibility"></a>

#### *property* max_primal_infeasibility

<a id="highspy.HighsInfo.max_primal_residual_error"></a>

#### *property* max_primal_residual_error

<a id="highspy.HighsInfo.max_relative_dual_infeasibility"></a>

#### *property* max_relative_dual_infeasibility

<a id="highspy.HighsInfo.max_relative_dual_residual_error"></a>

#### *property* max_relative_dual_residual_error

<a id="highspy.HighsInfo.max_relative_primal_infeasibility"></a>

#### *property* max_relative_primal_infeasibility

<a id="highspy.HighsInfo.max_relative_primal_residual_error"></a>

#### *property* max_relative_primal_residual_error

<a id="highspy.HighsInfo.mip_dual_bound"></a>

#### *property* mip_dual_bound

<a id="highspy.HighsInfo.mip_gap"></a>

#### *property* mip_gap

<a id="highspy.HighsInfo.mip_node_count"></a>

#### *property* mip_node_count

<a id="highspy.HighsInfo.num_complementarity_violations"></a>

#### *property* num_complementarity_violations

<a id="highspy.HighsInfo.num_dual_infeasibilities"></a>

#### *property* num_dual_infeasibilities

<a id="highspy.HighsInfo.num_dual_residual_errors"></a>

#### *property* num_dual_residual_errors

<a id="highspy.HighsInfo.num_primal_infeasibilities"></a>

#### *property* num_primal_infeasibilities

<a id="highspy.HighsInfo.num_primal_residual_errors"></a>

#### *property* num_primal_residual_errors

<a id="highspy.HighsInfo.num_relative_dual_infeasibilities"></a>

#### *property* num_relative_dual_infeasibilities

<a id="highspy.HighsInfo.num_relative_dual_residual_errors"></a>

#### *property* num_relative_dual_residual_errors

<a id="highspy.HighsInfo.num_relative_primal_infeasibilities"></a>

#### *property* num_relative_primal_infeasibilities

<a id="highspy.HighsInfo.num_relative_primal_residual_errors"></a>

#### *property* num_relative_primal_residual_errors

<a id="highspy.HighsInfo.objective_function_value"></a>

#### *property* objective_function_value

<a id="highspy.HighsInfo.pdlp_iteration_count"></a>

#### *property* pdlp_iteration_count

<a id="highspy.HighsInfo.primal_dual_integral"></a>

#### *property* primal_dual_integral

<a id="highspy.HighsInfo.primal_dual_objective_error"></a>

#### *property* primal_dual_objective_error

<a id="highspy.HighsInfo.primal_solution_status"></a>

#### *property* primal_solution_status

<a id="highspy.HighsInfo.qp_iteration_count"></a>

#### *property* qp_iteration_count

<a id="highspy.HighsInfo.simplex_iteration_count"></a>

#### *property* simplex_iteration_count

<a id="highspy.HighsInfo.sum_dual_infeasibilities"></a>

#### *property* sum_dual_infeasibilities

<a id="highspy.HighsInfo.sum_primal_infeasibilities"></a>

#### *property* sum_primal_infeasibilities

<a id="highspy.HighsInfo.valid"></a>

#### *property* valid

<a id="highspy.HighsInfoType"></a>

### *class* highspy.HighsInfoType

Bases: `pybind11_object`

Members:

kInt64

kInt

kDouble

<a id="highspy.HighsInfoType.kDouble"></a>

#### kDouble *= <HighsInfoType.kDouble: 2>*

<a id="highspy.HighsInfoType.kInt"></a>

#### kInt *= <HighsInfoType.kInt: 1>*

<a id="highspy.HighsInfoType.kInt64"></a>

#### kInt64 *= <HighsInfoType.kInt64: -1>*

### HighsInfoType.name -> str

<a id="highspy.HighsInfoType.value"></a>

#### *property* value

<a id="highspy.HighsLinearObjective"></a>

### *class* highspy.HighsLinearObjective

Bases: `pybind11_object`

<a id="highspy.HighsLinearObjective.abs_tolerance"></a>

#### *property* abs_tolerance

<a id="highspy.HighsLinearObjective.coefficients"></a>

#### *property* coefficients

<a id="highspy.HighsLinearObjective.offset"></a>

#### *property* offset

<a id="highspy.HighsLinearObjective.priority"></a>

#### *property* priority

<a id="highspy.HighsLinearObjective.rel_tolerance"></a>

#### *property* rel_tolerance

<a id="highspy.HighsLinearObjective.weight"></a>

#### *property* weight

<a id="highspy.HighsLogType"></a>

### *class* highspy.HighsLogType

Bases: `pybind11_object`

Members:

kInfo

kDetailed

kVerbose

kWarning

kError

<a id="highspy.HighsLogType.kDetailed"></a>

#### kDetailed *= <HighsLogType.kDetailed: 2>*

<a id="highspy.HighsLogType.kError"></a>

#### kError *= <HighsLogType.kError: 5>*

<a id="highspy.HighsLogType.kInfo"></a>

#### kInfo *= <HighsLogType.kInfo: 1>*

<a id="highspy.HighsLogType.kVerbose"></a>

#### kVerbose *= <HighsLogType.kVerbose: 3>*

<a id="highspy.HighsLogType.kWarning"></a>

#### kWarning *= <HighsLogType.kWarning: 4>*

### HighsLogType.name -> str

<a id="highspy.HighsLogType.value"></a>

#### *property* value

<a id="highspy.HighsLp"></a>

### *class* highspy.HighsLp

Bases: `pybind11_object`

<a id="highspy.HighsLp.a_matrix_"></a>

#### *property* a_matrix_

<a id="highspy.HighsLp.col_cost_"></a>

#### *property* col_cost_

<a id="highspy.HighsLp.col_lower_"></a>

#### *property* col_lower_

<a id="highspy.HighsLp.col_names_"></a>

#### *property* col_names_

<a id="highspy.HighsLp.col_upper_"></a>

#### *property* col_upper_

<a id="highspy.HighsLp.integrality_"></a>

#### *property* integrality_

<a id="highspy.HighsLp.is_moved_"></a>

#### *property* is_moved_

<a id="highspy.HighsLp.is_scaled_"></a>

#### *property* is_scaled_

<a id="highspy.HighsLp.model_name_"></a>

#### *property* model_name_

<a id="highspy.HighsLp.mods_"></a>

#### *property* mods_

<a id="highspy.HighsLp.num_col_"></a>

#### *property* num_col_

<a id="highspy.HighsLp.num_row_"></a>

#### *property* num_row_

<a id="highspy.HighsLp.offset_"></a>

#### *property* offset_

<a id="highspy.HighsLp.row_lower_"></a>

#### *property* row_lower_

<a id="highspy.HighsLp.row_names_"></a>

#### *property* row_names_

<a id="highspy.HighsLp.row_upper_"></a>

#### *property* row_upper_

<a id="highspy.HighsLp.scale_"></a>

#### *property* scale_

<a id="highspy.HighsLp.sense_"></a>

#### *property* sense_

<a id="highspy.HighsModel"></a>

### *class* highspy.HighsModel

Bases: `pybind11_object`

<a id="highspy.HighsModel.hessian_"></a>

#### *property* hessian_

<a id="highspy.HighsModel.lp_"></a>

#### *property* lp_

<a id="highspy.HighsModelStatus"></a>

### *class* highspy.HighsModelStatus

Bases: `pybind11_object`

Members:

kNotset

kLoadError

kModelError

kPresolveError

kSolveError

kPostsolveError

kModelEmpty

kOptimal

kInfeasible

kUnboundedOrInfeasible

kUnbounded

kObjectiveBound

kObjectiveTarget

kTimeLimit

kIterationLimit

kUnknown

kSolutionLimit

kInterrupt

kMemoryLimit

kHighsInterrupt

<a id="highspy.HighsModelStatus.kHighsInterrupt"></a>

#### kHighsInterrupt *= <HighsModelStatus.kHighsInterrupt: 19>*

<a id="highspy.HighsModelStatus.kInfeasible"></a>

#### kInfeasible *= <HighsModelStatus.kInfeasible: 8>*

<a id="highspy.HighsModelStatus.kInterrupt"></a>

#### kInterrupt *= <HighsModelStatus.kInterrupt: 17>*

<a id="highspy.HighsModelStatus.kIterationLimit"></a>

#### kIterationLimit *= <HighsModelStatus.kIterationLimit: 14>*

<a id="highspy.HighsModelStatus.kLoadError"></a>

#### kLoadError *= <HighsModelStatus.kLoadError: 1>*

<a id="highspy.HighsModelStatus.kMemoryLimit"></a>

#### kMemoryLimit *= <HighsModelStatus.kMemoryLimit: 18>*

<a id="highspy.HighsModelStatus.kModelEmpty"></a>

#### kModelEmpty *= <HighsModelStatus.kModelEmpty: 6>*

<a id="highspy.HighsModelStatus.kModelError"></a>

#### kModelError *= <HighsModelStatus.kModelError: 2>*

<a id="highspy.HighsModelStatus.kNotset"></a>

#### kNotset *= <HighsModelStatus.kNotset: 0>*

<a id="highspy.HighsModelStatus.kObjectiveBound"></a>

#### kObjectiveBound *= <HighsModelStatus.kObjectiveBound: 11>*

<a id="highspy.HighsModelStatus.kObjectiveTarget"></a>

#### kObjectiveTarget *= <HighsModelStatus.kObjectiveTarget: 12>*

<a id="highspy.HighsModelStatus.kOptimal"></a>

#### kOptimal *= <HighsModelStatus.kOptimal: 7>*

<a id="highspy.HighsModelStatus.kPostsolveError"></a>

#### kPostsolveError *= <HighsModelStatus.kPostsolveError: 5>*

<a id="highspy.HighsModelStatus.kPresolveError"></a>

#### kPresolveError *= <HighsModelStatus.kPresolveError: 3>*

<a id="highspy.HighsModelStatus.kSolutionLimit"></a>

#### kSolutionLimit *= <HighsModelStatus.kSolutionLimit: 16>*

<a id="highspy.HighsModelStatus.kSolveError"></a>

#### kSolveError *= <HighsModelStatus.kSolveError: 4>*

<a id="highspy.HighsModelStatus.kTimeLimit"></a>

#### kTimeLimit *= <HighsModelStatus.kTimeLimit: 13>*

<a id="highspy.HighsModelStatus.kUnbounded"></a>

#### kUnbounded *= <HighsModelStatus.kUnbounded: 10>*

<a id="highspy.HighsModelStatus.kUnboundedOrInfeasible"></a>

#### kUnboundedOrInfeasible *= <HighsModelStatus.kUnboundedOrInfeasible: 9>*

<a id="highspy.HighsModelStatus.kUnknown"></a>

#### kUnknown *= <HighsModelStatus.kUnknown: 15>*

### HighsModelStatus.name -> str

<a id="highspy.HighsModelStatus.value"></a>

#### *property* value

<a id="highspy.HighsObjectiveSolution"></a>

### *class* highspy.HighsObjectiveSolution

Bases: `pybind11_object`

<a id="highspy.HighsObjectiveSolution.col_value"></a>

#### *property* col_value

<a id="highspy.HighsObjectiveSolution.objective"></a>

#### *property* objective

<a id="highspy.HighsOptionType"></a>

### *class* highspy.HighsOptionType

Bases: `pybind11_object`

Members:

kBool

kInt

kDouble

kString

<a id="highspy.HighsOptionType.kBool"></a>

#### kBool *= <HighsOptionType.kBool: 0>*

<a id="highspy.HighsOptionType.kDouble"></a>

#### kDouble *= <HighsOptionType.kDouble: 2>*

<a id="highspy.HighsOptionType.kInt"></a>

#### kInt *= <HighsOptionType.kInt: 1>*

<a id="highspy.HighsOptionType.kString"></a>

#### kString *= <HighsOptionType.kString: 3>*

### HighsOptionType.name -> str

<a id="highspy.HighsOptionType.value"></a>

#### *property* value

<a id="highspy.HighsOptions"></a>

### *class* highspy.HighsOptions

Bases: `pybind11_object`

<a id="highspy.HighsOptions.allow_unbounded_or_infeasible"></a>

#### *property* allow_unbounded_or_infeasible

<a id="highspy.HighsOptions.allowed_matrix_scale_factor"></a>

#### *property* allowed_matrix_scale_factor

<a id="highspy.HighsOptions.blend_multi_objectives"></a>

#### *property* blend_multi_objectives

<a id="highspy.HighsOptions.dual_feasibility_tolerance"></a>

#### *property* dual_feasibility_tolerance

<a id="highspy.HighsOptions.dual_residual_tolerance"></a>

#### *property* dual_residual_tolerance

<a id="highspy.HighsOptions.glpsol_cost_row_location"></a>

#### *property* glpsol_cost_row_location

<a id="highspy.HighsOptions.highs_analysis_level"></a>

#### *property* highs_analysis_level

<a id="highspy.HighsOptions.highs_debug_level"></a>

#### *property* highs_debug_level

<a id="highspy.HighsOptions.infinite_bound"></a>

#### *property* infinite_bound

<a id="highspy.HighsOptions.infinite_cost"></a>

#### *property* infinite_cost

<a id="highspy.HighsOptions.ipm_iteration_limit"></a>

#### *property* ipm_iteration_limit

<a id="highspy.HighsOptions.ipm_optimality_tolerance"></a>

#### *property* ipm_optimality_tolerance

<a id="highspy.HighsOptions.ipx_dualize_strategy"></a>

#### *property* ipx_dualize_strategy

<a id="highspy.HighsOptions.kkt_tolerance"></a>

#### *property* kkt_tolerance

<a id="highspy.HighsOptions.large_matrix_value"></a>

#### *property* large_matrix_value

<a id="highspy.HighsOptions.log_dev_level"></a>

#### *property* log_dev_level

<a id="highspy.HighsOptions.log_file"></a>

#### *property* log_file

<a id="highspy.HighsOptions.log_githash"></a>

#### *property* log_githash

<a id="highspy.HighsOptions.log_to_console"></a>

#### *property* log_to_console

<a id="highspy.HighsOptions.mip_abs_gap"></a>

#### *property* mip_abs_gap

<a id="highspy.HighsOptions.mip_detect_symmetry"></a>

#### *property* mip_detect_symmetry

<a id="highspy.HighsOptions.mip_feasibility_tolerance"></a>

#### *property* mip_feasibility_tolerance

<a id="highspy.HighsOptions.mip_heuristic_effort"></a>

#### *property* mip_heuristic_effort

<a id="highspy.HighsOptions.mip_heuristic_run_feasibility_jump"></a>

#### *property* mip_heuristic_run_feasibility_jump

<a id="highspy.HighsOptions.mip_heuristic_run_rens"></a>

#### *property* mip_heuristic_run_rens

<a id="highspy.HighsOptions.mip_heuristic_run_rins"></a>

#### *property* mip_heuristic_run_rins

<a id="highspy.HighsOptions.mip_heuristic_run_root_reduced_cost"></a>

#### *property* mip_heuristic_run_root_reduced_cost

<a id="highspy.HighsOptions.mip_heuristic_run_shifting"></a>

#### *property* mip_heuristic_run_shifting

<a id="highspy.HighsOptions.mip_heuristic_run_zi_round"></a>

#### *property* mip_heuristic_run_zi_round

<a id="highspy.HighsOptions.mip_lp_age_limit"></a>

#### *property* mip_lp_age_limit

<a id="highspy.HighsOptions.mip_max_improving_sols"></a>

#### *property* mip_max_improving_sols

<a id="highspy.HighsOptions.mip_max_leaves"></a>

#### *property* mip_max_leaves

<a id="highspy.HighsOptions.mip_max_nodes"></a>

#### *property* mip_max_nodes

<a id="highspy.HighsOptions.mip_max_stall_nodes"></a>

#### *property* mip_max_stall_nodes

<a id="highspy.HighsOptions.mip_min_cliquetable_entries_for_parallelism"></a>

#### *property* mip_min_cliquetable_entries_for_parallelism

<a id="highspy.HighsOptions.mip_min_logging_interval"></a>

#### *property* mip_min_logging_interval

<a id="highspy.HighsOptions.mip_pool_age_limit"></a>

#### *property* mip_pool_age_limit

<a id="highspy.HighsOptions.mip_pool_soft_limit"></a>

#### *property* mip_pool_soft_limit

<a id="highspy.HighsOptions.mip_pscost_minreliable"></a>

#### *property* mip_pscost_minreliable

<a id="highspy.HighsOptions.mip_rel_gap"></a>

#### *property* mip_rel_gap

<a id="highspy.HighsOptions.mip_report_level"></a>

#### *property* mip_report_level

<a id="highspy.HighsOptions.objective_bound"></a>

#### *property* objective_bound

<a id="highspy.HighsOptions.objective_target"></a>

#### *property* objective_target

<a id="highspy.HighsOptions.optimality_tolerance"></a>

#### *property* optimality_tolerance

<a id="highspy.HighsOptions.output_flag"></a>

#### *property* output_flag

<a id="highspy.HighsOptions.parallel"></a>

#### *property* parallel

<a id="highspy.HighsOptions.pdlp_cupdlpc_restart_method"></a>

#### *property* pdlp_cupdlpc_restart_method

<a id="highspy.HighsOptions.pdlp_iteration_limit"></a>

#### *property* pdlp_iteration_limit

<a id="highspy.HighsOptions.pdlp_optimality_tolerance"></a>

#### *property* pdlp_optimality_tolerance

<a id="highspy.HighsOptions.pdlp_scaling_mode"></a>

#### *property* pdlp_scaling_mode

<a id="highspy.HighsOptions.presolve"></a>

#### *property* presolve

<a id="highspy.HighsOptions.primal_feasibility_tolerance"></a>

#### *property* primal_feasibility_tolerance

<a id="highspy.HighsOptions.primal_residual_tolerance"></a>

#### *property* primal_residual_tolerance

<a id="highspy.HighsOptions.qp_iteration_limit"></a>

#### *property* qp_iteration_limit

<a id="highspy.HighsOptions.qp_nullspace_limit"></a>

#### *property* qp_nullspace_limit

<a id="highspy.HighsOptions.qp_regularization_value"></a>

#### *property* qp_regularization_value

<a id="highspy.HighsOptions.random_seed"></a>

#### *property* random_seed

<a id="highspy.HighsOptions.ranging"></a>

#### *property* ranging

<a id="highspy.HighsOptions.read_basis_file"></a>

#### *property* read_basis_file

<a id="highspy.HighsOptions.read_solution_file"></a>

#### *property* read_solution_file

<a id="highspy.HighsOptions.run_crossover"></a>

#### *property* run_crossover

<a id="highspy.HighsOptions.simplex_crash_strategy"></a>

#### *property* simplex_crash_strategy

<a id="highspy.HighsOptions.simplex_dual_edge_weight_strategy"></a>

#### *property* simplex_dual_edge_weight_strategy

<a id="highspy.HighsOptions.simplex_dualize_strategy"></a>

#### *property* simplex_dualize_strategy

<a id="highspy.HighsOptions.simplex_iteration_limit"></a>

#### *property* simplex_iteration_limit

<a id="highspy.HighsOptions.simplex_max_concurrency"></a>

#### *property* simplex_max_concurrency

<a id="highspy.HighsOptions.simplex_min_concurrency"></a>

#### *property* simplex_min_concurrency

<a id="highspy.HighsOptions.simplex_permute_strategy"></a>

#### *property* simplex_permute_strategy

<a id="highspy.HighsOptions.simplex_price_strategy"></a>

#### *property* simplex_price_strategy

<a id="highspy.HighsOptions.simplex_primal_edge_weight_strategy"></a>

#### *property* simplex_primal_edge_weight_strategy

<a id="highspy.HighsOptions.simplex_scale_strategy"></a>

#### *property* simplex_scale_strategy

<a id="highspy.HighsOptions.simplex_strategy"></a>

#### *property* simplex_strategy

<a id="highspy.HighsOptions.simplex_update_limit"></a>

#### *property* simplex_update_limit

<a id="highspy.HighsOptions.small_matrix_value"></a>

#### *property* small_matrix_value

<a id="highspy.HighsOptions.solution_file"></a>

#### *property* solution_file

<a id="highspy.HighsOptions.solve_relaxation"></a>

#### *property* solve_relaxation

<a id="highspy.HighsOptions.solver"></a>

#### *property* solver

<a id="highspy.HighsOptions.threads"></a>

#### *property* threads

<a id="highspy.HighsOptions.time_limit"></a>

#### *property* time_limit

<a id="highspy.HighsOptions.timeless_log"></a>

#### *property* timeless_log

<a id="highspy.HighsOptions.user_bound_scale"></a>

#### *property* user_bound_scale

<a id="highspy.HighsOptions.user_objective_scale"></a>

#### *property* user_objective_scale

<a id="highspy.HighsOptions.write_basis_file"></a>

#### *property* write_basis_file

<a id="highspy.HighsOptions.write_model_file"></a>

#### *property* write_model_file

<a id="highspy.HighsOptions.write_model_to_file"></a>

#### *property* write_model_to_file

<a id="highspy.HighsOptions.write_presolved_model_file"></a>

#### *property* write_presolved_model_file

<a id="highspy.HighsOptions.write_solution_style"></a>

#### *property* write_solution_style

<a id="highspy.HighsOptions.write_solution_to_file"></a>

#### *property* write_solution_to_file

<a id="highspy.HighsPresolveStatus"></a>

### *class* highspy.HighsPresolveStatus

Bases: `pybind11_object`

Members:

kNotPresolved

kNotReduced

kInfeasible

kUnboundedOrInfeasible

kReduced

kReducedToEmpty

kTimeout

kNullError

kOptionsError

<a id="highspy.HighsPresolveStatus.kInfeasible"></a>

#### kInfeasible *= <HighsPresolveStatus.kInfeasible: 1>*

<a id="highspy.HighsPresolveStatus.kNotPresolved"></a>

#### kNotPresolved *= <HighsPresolveStatus.kNotPresolved: -1>*

<a id="highspy.HighsPresolveStatus.kNotReduced"></a>

#### kNotReduced *= <HighsPresolveStatus.kNotReduced: 0>*

<a id="highspy.HighsPresolveStatus.kNullError"></a>

#### kNullError *= <HighsPresolveStatus.kNullError: 6>*

<a id="highspy.HighsPresolveStatus.kOptionsError"></a>

#### kOptionsError *= <HighsPresolveStatus.kOptionsError: 7>*

<a id="highspy.HighsPresolveStatus.kReduced"></a>

#### kReduced *= <HighsPresolveStatus.kReduced: 3>*

<a id="highspy.HighsPresolveStatus.kReducedToEmpty"></a>

#### kReducedToEmpty *= <HighsPresolveStatus.kReducedToEmpty: 4>*

<a id="highspy.HighsPresolveStatus.kTimeout"></a>

#### kTimeout *= <HighsPresolveStatus.kTimeout: 5>*

<a id="highspy.HighsPresolveStatus.kUnboundedOrInfeasible"></a>

#### kUnboundedOrInfeasible *= <HighsPresolveStatus.kUnboundedOrInfeasible: 2>*

### HighsPresolveStatus.name -> str

<a id="highspy.HighsPresolveStatus.value"></a>

#### *property* value

<a id="highspy.HighsRanging"></a>

### *class* highspy.HighsRanging

Bases: `pybind11_object`

<a id="highspy.HighsRanging.col_bound_dn"></a>

#### *property* col_bound_dn

<a id="highspy.HighsRanging.col_bound_up"></a>

#### *property* col_bound_up

<a id="highspy.HighsRanging.col_cost_dn"></a>

#### *property* col_cost_dn

<a id="highspy.HighsRanging.col_cost_up"></a>

#### *property* col_cost_up

<a id="highspy.HighsRanging.row_bound_dn"></a>

#### *property* row_bound_dn

<a id="highspy.HighsRanging.row_bound_up"></a>

#### *property* row_bound_up

<a id="highspy.HighsRanging.valid"></a>

#### *property* valid

<a id="highspy.HighsRangingRecord"></a>

### *class* highspy.HighsRangingRecord

Bases: `pybind11_object`

<a id="highspy.HighsRangingRecord.in_var_"></a>

#### *property* in_var_

<a id="highspy.HighsRangingRecord.objective_"></a>

#### *property* objective_

<a id="highspy.HighsRangingRecord.ou_var_"></a>

#### *property* ou_var_

<a id="highspy.HighsRangingRecord.value_"></a>

#### *property* value_

<a id="highspy.HighsSolution"></a>

### *class* highspy.HighsSolution

Bases: `pybind11_object`

<a id="highspy.HighsSolution.col_dual"></a>

#### *property* col_dual

<a id="highspy.HighsSolution.col_value"></a>

#### *property* col_value

<a id="highspy.HighsSolution.dual_valid"></a>

#### *property* dual_valid

<a id="highspy.HighsSolution.row_dual"></a>

#### *property* row_dual

<a id="highspy.HighsSolution.row_value"></a>

#### *property* row_value

<a id="highspy.HighsSolution.value_valid"></a>

#### *property* value_valid

<a id="highspy.HighsSparseMatrix"></a>

### *class* highspy.HighsSparseMatrix

Bases: `pybind11_object`

<a id="highspy.HighsSparseMatrix.format_"></a>

#### *property* format_

<a id="highspy.HighsSparseMatrix.index_"></a>

#### *property* index_

<a id="highspy.HighsSparseMatrix.num_col_"></a>

#### *property* num_col_

<a id="highspy.HighsSparseMatrix.num_row_"></a>

#### *property* num_row_

<a id="highspy.HighsSparseMatrix.p_end_"></a>

#### *property* p_end_

<a id="highspy.HighsSparseMatrix.start_"></a>

#### *property* start_

<a id="highspy.HighsSparseMatrix.value_"></a>

#### *property* value_

<a id="highspy.HighsStatus"></a>

### *class* highspy.HighsStatus

Bases: `pybind11_object`

Members:

kError

kOk

kWarning

<a id="highspy.HighsStatus.kError"></a>

#### kError *= <HighsStatus.kError: -1>*

<a id="highspy.HighsStatus.kOk"></a>

#### kOk *= <HighsStatus.kOk: 0>*

<a id="highspy.HighsStatus.kWarning"></a>

#### kWarning *= <HighsStatus.kWarning: 1>*

### HighsStatus.name -> str

<a id="highspy.HighsStatus.value"></a>

#### *property* value

<a id="highspy.HighsStatusError"></a>

### *exception* highspy.HighsStatusError(operation, status)

Bases: [`HighsError`](#highspy.highs.HighsError)

Raised when a HiGHS operation returns an unsuccessful status.

* **Parameters:**
  * **operation** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **status** ([*HighsStatus*](#highspy.HighsStatus))

<a id="highspy.HighsVarType"></a>

### *class* highspy.HighsVarType

Bases: `pybind11_object`

Members:

kContinuous

kInteger

kSemiContinuous

kSemiInteger

kImplicitInteger

<a id="highspy.HighsVarType.kContinuous"></a>

#### kContinuous *= <HighsVarType.kContinuous: 0>*

<a id="highspy.HighsVarType.kImplicitInteger"></a>

#### kImplicitInteger *= <HighsVarType.kImplicitInteger: 4>*

<a id="highspy.HighsVarType.kInteger"></a>

#### kInteger *= <HighsVarType.kInteger: 1>*

<a id="highspy.HighsVarType.kSemiContinuous"></a>

#### kSemiContinuous *= <HighsVarType.kSemiContinuous: 2>*

<a id="highspy.HighsVarType.kSemiInteger"></a>

#### kSemiInteger *= <HighsVarType.kSemiInteger: 3>*

### HighsVarType.name -> str

<a id="highspy.HighsVarType.value"></a>

#### *property* value

<a id="highspy.HighspyArray"></a>

### *class* highspy.HighspyArray(input_array, highs)

Bases: [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[`Any`](https://docs.python.org/3/library/typing.html#typing.Any), [`dtype`](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[`object_`](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_)]]

A numpy array wrapper for highs_var/highs_linear_expression objects.

This provides additional type information for static analysis, and also allows faster sum operations.

* **Parameters:**
  * **input_array** (*np.ndarray* *[**Any* *,* *np.dtype* *[**np.object_* *]* *]*)
  * **highs** ([*Highs*](#highspy.highs.Highs) *|* *None*)
* **Return type:**
  Self

<a id="highspy.HighspyArray.highs"></a>

#### highs *: [Highs](#highspy.highs.Highs) | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.HighspyArray.sum"></a>

#### sum(axis: [None](https://docs.python.org/3/builtins/constants.html#None) = None, dtype: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, out: [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [highs_linear_expression](#highspy.highs.highs_linear_expression)

#### sum(axis: [Any](https://docs.python.org/3/library/typing.html#typing.Any), dtype: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, out: [HighspyArray](#highspy.highs.HighspyArray) = None) → [HighspyArray](#highspy.highs.HighspyArray)

Return the sum of the array elements over the given axis.

Refer to numpy.sum for full documentation.

#### SEE ALSO
[`numpy.sum`](https://numpy.org/doc/stable/reference/generated/numpy.sum.html#numpy.sum)
: equivalent function

* **Parameters:**
  * **axis** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
  * **dtype** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*)
  * **out** ([*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[*[*object_*](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_) *]* *]*  *|* *None*)
  * **unused_kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any))
* **Return type:**
  [*HighspyArray*](#highspy.highs.HighspyArray) | [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.HighspyArray.idx"></a>

#### idx()

Convert to a flat int32 index array for passing to the HiGHS C++ API.

Each element’s `__index__` method is called to extract its integer
index (e.g., `highs_var.index` or `highs_cons.index`).
The result is always 1-D, regardless of the array’s shape.

* **Returns:**
  A new flat int32 numpy array of the underlying indices.
* **Return type:**
  [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*int32*]]

<a id="highspy.IisBoundStatus"></a>

### *class* highspy.IisBoundStatus

Bases: `pybind11_object`

Members:

kIisBoundStatusDropped

kIisBoundStatusNull

kIisBoundStatusFree

kIisBoundStatusLower

kIisBoundStatusUpper

kIisBoundStatusBoxed

<a id="highspy.IisBoundStatus.kIisBoundStatusBoxed"></a>

#### kIisBoundStatusBoxed *= <IisBoundStatus.kIisBoundStatusBoxed: 4>*

<a id="highspy.IisBoundStatus.kIisBoundStatusDropped"></a>

#### kIisBoundStatusDropped *= <IisBoundStatus.kIisBoundStatusDropped: -1>*

<a id="highspy.IisBoundStatus.kIisBoundStatusFree"></a>

#### kIisBoundStatusFree *= <IisBoundStatus.kIisBoundStatusFree: 1>*

<a id="highspy.IisBoundStatus.kIisBoundStatusLower"></a>

#### kIisBoundStatusLower *= <IisBoundStatus.kIisBoundStatusLower: 2>*

<a id="highspy.IisBoundStatus.kIisBoundStatusNull"></a>

#### kIisBoundStatusNull *= <IisBoundStatus.kIisBoundStatusNull: 0>*

<a id="highspy.IisBoundStatus.kIisBoundStatusUpper"></a>

#### kIisBoundStatusUpper *= <IisBoundStatus.kIisBoundStatusUpper: 3>*

### IisBoundStatus.name -> str

<a id="highspy.IisBoundStatus.value"></a>

#### *property* value

<a id="highspy.IisStatus"></a>

### *class* highspy.IisStatus

Bases: `pybind11_object`

Members:

kIisStatusNotInConflict

kIisStatusMaybeInConflict

kIisStatusInConflict

<a id="highspy.IisStatus.kIisStatusInConflict"></a>

#### kIisStatusInConflict *= <IisStatus.kIisStatusInConflict: 1>*

<a id="highspy.IisStatus.kIisStatusMaybeInConflict"></a>

#### kIisStatusMaybeInConflict *= <IisStatus.kIisStatusMaybeInConflict: 0>*

<a id="highspy.IisStatus.kIisStatusNotInConflict"></a>

#### kIisStatusNotInConflict *= <IisStatus.kIisStatusNotInConflict: -1>*

### IisStatus.name -> str

<a id="highspy.IisStatus.value"></a>

#### *property* value

<a id="highspy.IisStrategy"></a>

### *class* highspy.IisStrategy

Bases: `pybind11_object`

Members:

kIisStrategyMin

kIisStrategyLight

kIisStrategyFromRay

kIisStrategyFromLp

kIisStrategyIrreducible

kIisStrategyColPriority

kIisStrategyRelaxation

kIisStrategyMax

<a id="highspy.IisStrategy.kIisStrategyColPriority"></a>

#### kIisStrategyColPriority *= <IisStrategy.kIisStrategyColPriority: 8>*

<a id="highspy.IisStrategy.kIisStrategyFromLp"></a>

#### kIisStrategyFromLp *= <IisStrategy.kIisStrategyFromLp: 2>*

<a id="highspy.IisStrategy.kIisStrategyFromRay"></a>

#### kIisStrategyFromRay *= <IisStrategy.kIisStrategyFromRay: 1>*

<a id="highspy.IisStrategy.kIisStrategyIrreducible"></a>

#### kIisStrategyIrreducible *= <IisStrategy.kIisStrategyIrreducible: 4>*

<a id="highspy.IisStrategy.kIisStrategyLight"></a>

#### kIisStrategyLight *= <IisStrategy.kIisStrategyMin: 0>*

<a id="highspy.IisStrategy.kIisStrategyMax"></a>

#### kIisStrategyMax *= <IisStrategy.kIisStrategyMax: 31>*

<a id="highspy.IisStrategy.kIisStrategyMin"></a>

#### kIisStrategyMin *= <IisStrategy.kIisStrategyMin: 0>*

<a id="highspy.IisStrategy.kIisStrategyRelaxation"></a>

#### kIisStrategyRelaxation *= <IisStrategy.kIisStrategyRelaxation: 16>*

### IisStrategy.name -> str

<a id="highspy.IisStrategy.value"></a>

#### *property* value

<a id="highspy.MatrixFormat"></a>

### *class* highspy.MatrixFormat

Bases: `pybind11_object`

Members:

kColwise

kRowwise

kRowwisePartitioned

<a id="highspy.MatrixFormat.kColwise"></a>

#### kColwise *= <MatrixFormat.kColwise: 1>*

<a id="highspy.MatrixFormat.kRowwise"></a>

#### kRowwise *= <MatrixFormat.kRowwise: 2>*

<a id="highspy.MatrixFormat.kRowwisePartitioned"></a>

#### kRowwisePartitioned *= <MatrixFormat.kRowwisePartitioned: 3>*

### MatrixFormat.name -> str

<a id="highspy.MatrixFormat.value"></a>

#### *property* value

<a id="highspy.ObjSense"></a>

### *class* highspy.ObjSense

Bases: `pybind11_object`

Members:

kMinimize

kMaximize

<a id="highspy.ObjSense.kMaximize"></a>

#### kMaximize *= <ObjSense.kMaximize: -1>*

<a id="highspy.ObjSense.kMinimize"></a>

#### kMinimize *= <ObjSense.kMinimize: 1>*

### ObjSense.name -> str

<a id="highspy.ObjSense.value"></a>

#### *property* value

<a id="highspy.SolutionStatus"></a>

### *class* highspy.SolutionStatus

Bases: `pybind11_object`

Members:

kSolutionStatusNone

kSolutionStatusInfeasible

kSolutionStatusFeasible

<a id="highspy.SolutionStatus.kSolutionStatusFeasible"></a>

#### kSolutionStatusFeasible *= <SolutionStatus.kSolutionStatusFeasible: 2>*

<a id="highspy.SolutionStatus.kSolutionStatusInfeasible"></a>

#### kSolutionStatusInfeasible *= <SolutionStatus.kSolutionStatusInfeasible: 1>*

<a id="highspy.SolutionStatus.kSolutionStatusNone"></a>

#### kSolutionStatusNone *= <SolutionStatus.kSolutionStatusNone: 0>*

### SolutionStatus.name -> str

<a id="highspy.SolutionStatus.value"></a>

#### *property* value

<a id="highspy.highs_cons"></a>

### *class* highspy.highs_cons(i, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Constraint index wrapper for HiGHS

* **Parameters:**
  * **i** ([*int*](https://docs.python.org/3/builtins/functions.html#int))
  * **highs** ([*Highs*](#highspy.Highs))

<a id="highspy.highs_cons.index"></a>

#### index

<a id="highspy.highs_cons.highs"></a>

#### highs

<a id="highspy.highs_cons.expr"></a>

#### expr()

Retrieves the expression of the constraint.

* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs_cons.name"></a>

#### *property* name *: [str](https://docs.python.org/3/builtins/stdtypes.html#str)*

<a id="highspy.highs_linear_expression"></a>

### *class* highspy.highs_linear_expression(other: [None](https://docs.python.org/3/builtins/constants.html#None) = None)

### *class* highspy.highs_linear_expression(other: [float](https://docs.python.org/3/builtins/functions.html#float))

### *class* highspy.highs_linear_expression(other: [highs_var](#highspy.highs.highs_var))

### *class* highspy.highs_linear_expression(other: [highs_linear_expression](#highspy.highs.highs_linear_expression))

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Linear constraint builder for HiGHS

* **Parameters:**
  **other** ([*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*highs_var*](#highspy.highs_var) *|* [*highs_linear_expression*](#highspy.highs_linear_expression) *|* *None*)

<a id="highspy.highs_linear_expression.bounds"></a>

#### bounds *: [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float)] | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.highs_linear_expression.idxs"></a>

#### idxs *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]*

<a id="highspy.highs_linear_expression.vals"></a>

#### vals *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[float](https://docs.python.org/3/builtins/functions.html#float)]*

<a id="highspy.highs_linear_expression.constant"></a>

#### constant *: [float](https://docs.python.org/3/builtins/functions.html#float) | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.highs_linear_expression.simplify"></a>

#### simplify()

Simplifies the linear expression by combining duplicate variables.

<a id="highspy.highs_linear_expression.copy"></a>

#### copy()

Creates a copy of the linear expression.

<a id="highspy.highs_linear_expression.evaluate"></a>

#### evaluate(values)

Evaluates the linear expression given a solution array (values).

* **Parameters:**
  **values** ([*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*float*](https://docs.python.org/3/builtins/functions.html#float) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[**float64* *]* *]*)
* **Return type:**
  [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

<a id="highspy.highs_linear_expression.unique_elements"></a>

#### unique_elements()

Collects unique variables and sums their corresponding values.  Keeps all values (including zeros).

<a id="highspy.highs_linear_expression.reduced_elements"></a>

#### reduced_elements()

Similar to unique_elements, except keeps only non-zero values

<a id="highspy.highs_var"></a>

### *class* highspy.highs_var(i, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Variable index wrapper for HiGHS

* **Parameters:**
  * **i** ([*int*](https://docs.python.org/3/builtins/functions.html#int))
  * **highs** ([*Highs*](#highspy.Highs))

<a id="highspy.highs_var.index"></a>

#### index

<a id="highspy.highs_var.highs"></a>

#### highs

<a id="highspy.highs_var.name"></a>

#### *property* name *: [str](https://docs.python.org/3/builtins/stdtypes.html#str)*

<a id="highspy.qsum"></a>

### highspy.qsum(items, initial=None)

Performs a faster sum for highs_linear_expressions.

* **Parameters:**
  * **items** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*) – A collection of highs_linear_expressions or highs_vars to be summed.
  * **initial** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)

<a id="high-level-python-interface"></a>

## High-Level Python Interface

`highspy.Highs` inherits most of its solver API (`readModel`, `run`,
`getSolution`, `setOptionValue`, and so on) from the pybind11-bound `_Highs`
base class, so those inherited members are included below too.

<a id="module-highspy.highs"></a>

<a id="highspy.highs.HighspyLeafTypes"></a>

### highspy.highs.HighspyLeafTypes *: [TypeAlias](https://docs.python.org/3/library/typing.html#typing.TypeAlias)*

<a id="highspy.highs.HighspyNestedIndex"></a>

### highspy.highs.HighspyNestedIndex *: [TypeAlias](https://docs.python.org/3/library/typing.html#typing.TypeAlias)*

<a id="highspy.highs.HighspyNestedResult"></a>

### highspy.highs.HighspyNestedResult *: [TypeAlias](https://docs.python.org/3/library/typing.html#typing.TypeAlias)*

<a id="highspy.highs.HighsError"></a>

### *exception* highspy.highs.HighsError

Bases: [`Exception`](https://docs.python.org/3/builtins/exceptions.html#Exception)

Base exception for highspy errors.

<a id="highspy.highs.HighsStatusError"></a>

### *exception* highspy.highs.HighsStatusError(operation, status)

Bases: [`HighsError`](#highspy.highs.HighsError)

Raised when a HiGHS operation returns an unsuccessful status.

* **Parameters:**
  * **operation** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **status** ([*HighsStatus*](#highspy.HighsStatus))

<a id="highspy.highs.Highs"></a>

### *class* highspy.highs.Highs

Bases: `_Highs`

HiGHS solver interface

<a id="highspy.highs.Highs.silent"></a>

#### silent(turn_off_output=True)

Disables solver output to the console.

* **Parameters:**
  **turn_off_output** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))

<a id="highspy.highs.Highs.solve"></a>

#### solve()

Runs the solver on the current problem.

* **Returns:**
  A HighsStatus object containing the solve status.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.highs.Highs.startSolve"></a>

#### startSolve()

Starts the solver in a separate thread.  Useful for handling KeyboardInterrupts.
Do not attempt to modify the model while the solver is running.

* **Returns:**
  A Thread object representing the solver thread.
* **Return type:**
  [*Thread*](https://docs.python.org/3/library/threading.html#threading.Thread)

<a id="highspy.highs.Highs.is_solver_running"></a>

#### is_solver_running()

* **Return type:**
  [bool](https://docs.python.org/3/builtins/functions.html#bool)

<a id="highspy.highs.Highs.joinSolve"></a>

#### joinSolve(solver_thread=None, interrupt_limit=5)

Waits for the solver to finish. If solver_thread is provided, it will handle KeyboardInterrupts.

* **Parameters:**
  * **solver_thread** ([*Thread*](https://docs.python.org/3/library/threading.html#threading.Thread) *|* *None*) – A Thread object representing the solver thread (optional).
  * **interrupt_limit** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – The number of times to allow KeyboardInterrupt before forcing termination (optional).
* **Returns:**
  A HighsStatus object containing the solve status.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.highs.Highs.wait"></a>

#### wait(timeout=-1.0)

* **Parameters:**
  **timeout** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[bool](https://docs.python.org/3/builtins/functions.html#bool), [*HighsStatus*](#highspy._core.HighsStatus) | None]

<a id="highspy.highs.Highs.optimize"></a>

#### optimize()

Alias for the solve method.

* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.highs.Highs.getObjective"></a>

#### getObjective()

Retrieves the current objective function (as a linear expression) and sense.

* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[*highs_linear_expression*](#highspy.highs.highs_linear_expression), [*ObjSense*](#highspy._core.ObjSense)]

<a id="highspy.highs.Highs.setObjective"></a>

#### setObjective(obj=None, sense=None)

Updates the costs.

* **Parameters:**
  * **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
  * **sense** ([*ObjSense*](#highspy._core.ObjSense) *|* *None*) – An optional ObjSense value representing the new objective sense.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.

<a id="highspy.highs.Highs.minimize"></a>

#### minimize(obj=None)

Solves a minimization of the objective and optionally updates the costs.

* **Parameters:**
  **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.
* **Returns:**
  A HighsStatus object containing the solve status after minimization.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.highs.Highs.maximize"></a>

#### maximize(obj=None)

Solves a maximization of the objective and optionally updates the costs.

* **Parameters:**
  **obj** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*) – An optional highs_linear_expression representing the new objective function.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If obj is an inequality or not a highs_linear_expression.
* **Returns:**
  A HighsStatus object containing the solve status after maximization.
* **Return type:**
  [*HighsStatus*](#highspy._core.HighsStatus) | None

<a id="highspy.highs.Highs.internal_get_value"></a>

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### *static* internal_get_value(array_values: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[float](https://docs.python.org/3/builtins/functions.html#float)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]], index_collection: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Internal method to get the value of an index from an array of values. Could be value or dual, variable or constraint.

* **Parameters:**
  * **array_values** ([*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*float*](https://docs.python.org/3/builtins/functions.html#float) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[**float64* *]* *]*)
  * **index_collection** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray))
* **Return type:**
  [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool) | [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*Any*](https://docs.python.org/3/library/typing.html#typing.Any)] | [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*float64*]]

<a id="highspy.highs.Highs.val"></a>

#### val(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### val(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### val(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### val(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### val(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Gets the value of a variable/index or expression in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var/index or highs_linear_expression object representing the variable.
* **Returns:**
  The value of the variable in the solution.

<a id="highspy.highs.Highs.vals"></a>

#### vals(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### vals(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### vals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### vals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### vals(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Gets the values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.

<a id="highspy.highs.Highs.variableName"></a>

#### variableName(var)

Retrieves the name of a specific variable.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var)) – A highs_var object representing the variable.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If the variable name cannot be found.
* **Returns:**
  The name of the specified variable.

<a id="highspy.highs.Highs.variableNames"></a>

#### variableNames(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var) | [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

#### variableNames(idxs: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[highs_var](#highspy.highs.highs_var) | [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral)]) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

Retrieves the names of multiple variables.

* **Parameters:**
  **idxs** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *]*) – An iterable of highs_var objects or a mapping where keys are identifiers and values are highs_var objects.
* **Raises:**
  [**Exception**](https://docs.python.org/3/builtins/exceptions.html#Exception) – If any variable name cannot be found.
* **Returns:**
  If idxs is a mapping, returns a dict where keys are the same keys from the input idxs and values are the names of the corresponding variables.
  If idxs is an iterable, returns a list of names for the specified variables.

<a id="highspy.highs.Highs.allVariableNames"></a>

#### allVariableNames()

Retrieves the names of all variables in the model.

* **Returns:**
  A list of strings representing the names of all variables.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.highs.Highs.variableValue"></a>

#### variableValue(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableValue(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableValue(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableValue(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableValue(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the value of a specific variable in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var object representing the variable.
* **Returns:**
  The value of the specified variable in the solution.

<a id="highspy.highs.Highs.variableValues"></a>

#### variableValues(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableValues(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableValues(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableValues(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableValues(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.

<a id="highspy.highs.Highs.allVariableValues"></a>

#### allVariableValues()

Retrieves the values of all variables in the solution.

* **Returns:**
  A list of values for all variables in the solution.

<a id="highspy.highs.Highs.variableDual"></a>

#### variableDual(var: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableDual(var: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableDual(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableDual(var: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableDual(var: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual value of a specific variable/index or expression in the solution.

* **Parameters:**
  **var** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_var object representing the variable.
* **Returns:**
  The dual value of the specified variable in the solution.

<a id="highspy.highs.Highs.variableDuals"></a>

#### variableDuals(idxs: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### variableDuals(idxs: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### variableDuals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_var](#highspy.highs.highs_var) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### variableDuals(idxs: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### variableDuals(idxs: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual values of multiple variables in the solution.

* **Parameters:**
  **idxs** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.
* **Returns:**
  If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the dual values of the corresponding variables. If idxs is an iterable, returns a list of dual values for the variables.

<a id="highspy.highs.Highs.allVariableDuals"></a>

#### allVariableDuals()

Retrieves the dual values of all variables in the solution.

* **Returns:**
  A list of dual values for all variables in the solution.

<a id="highspy.highs.Highs.constrValue"></a>

#### constrValue(con: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrValue(con: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrValue(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrValue(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrValue(con: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the value of a specific constraint in the solution.

* **Parameters:**
  **con** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_con object representing the constraint.
* **Returns:**
  The value of the specified constraint in the solution.

<a id="highspy.highs.Highs.constrValues"></a>

#### constrValues(cons: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrValues(cons: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrValues(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrValues(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrValues(cons: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the values of multiple constraints in the solution.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.
* **Returns:**
  If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.

<a id="highspy.highs.Highs.allConstrValues"></a>

#### allConstrValues()

Retrieves the values of all constraints in the solution.

* **Returns:**
  A list of values for all constraints in the solution.

<a id="highspy.highs.Highs.constrDual"></a>

#### constrDual(con: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrDual(con: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrDual(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrDual(con: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrDual(con: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual value of a specific constraint in the solution.

* **Parameters:**
  **con** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A highs_con object representing the constraint.
* **Returns:**
  The dual value of the specified constraint in the solution.

<a id="highspy.highs.Highs.constrDuals"></a>

#### constrDuals(cons: [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### constrDuals(cons: [highs_linear_expression](#highspy.highs.highs_linear_expression)) → [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### constrDuals(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [int](https://docs.python.org/3/builtins/functions.html#int) | [Integral](https://docs.python.org/3/library/numbers.html#numbers.Integral) | [highs_cons](#highspy.highs.highs_cons)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [float](https://docs.python.org/3/builtins/functions.html#float)]

#### constrDuals(cons: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]

#### constrDuals(cons: [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]]) → [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [dtype](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[float64]]

Retrieves the dual values of multiple constraints in the solution.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.
* **Returns:**
  If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.

<a id="highspy.highs.Highs.allConstrDuals"></a>

#### allConstrDuals()

Retrieves the dual values of all constraints in the solution.

* **Returns:**
  A list of dual values for all constraints in the solution.

<a id="highspy.highs.Highs.addVariable"></a>

#### addVariable(lb=0, ub=inf, obj=0.0, type=<HighsVarType.kContinuous: 0>, name=None)

Adds a variable to the model.

* **Parameters:**
  * **lb** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Lower bound of the variable (default is 0).
  * **ub** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Upper bound of the variable (default is infinity).
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – Objective coefficient of the variable (default is 0).
  * **type** ([*HighsVarType*](#highspy._core.HighsVarType)) – Type of the variable (continuous, integer; default is continuous).
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – Optional name for the variable.
* **Returns:**
  A highs_var object representing the added variable.

<a id="highspy.highs.Highs.addVariables"></a>

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addVariables(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addVariables(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addVariables(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Adds multiple variables to the model.

* **Parameters:**
  * **\*args** – A sequence of variables to be added. Can be a collection of scalars or indices (or mix).
  * **\*\*kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*) – Optional keyword arguments.  Can be scalars, arrays, or mappings.
    lb: Lower bound of the variables (default is 0).
    ub: Upper bound of the variables (default is infinity).
    obj: Objective coefficient of the variables (default is 0).
    type: Type of the variables (continuous, integer; default is continuous).
    name: A collection of names for the variables (list or mapping).
    name_prefix: Prefix for the variable names.  Constructed name will be name_prefix + index.
    out_array: Return an array of highs_var objects instead of a dictionary.
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **\*\*kwargs**
* **Returns:**
  A highs_var collection (array or dictionary) representing the added variables.
* **Return type:**
  [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*highs_var*](#highspy.highs.highs_var)] | [*HighspyArray*](#highspy.highs.HighspyArray) | None

<a id="highspy.highs.Highs.addIntegrals"></a>

#### addIntegrals(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addIntegrals(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addIntegrals(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addIntegrals(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Alias for the addVariables method, for integer variables.

* **Parameters:**
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)

<a id="highspy.highs.Highs.addBinaries"></a>

#### addBinaries(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int), out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addBinaries(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[False] = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)]

#### addBinaries(\*nvars: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [Literal](https://docs.python.org/3/library/typing.html#typing.Literal)[True], \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [HighspyArray](#highspy.highs.HighspyArray)

#### addBinaries(\*nvars: [int](https://docs.python.org/3/builtins/functions.html#int) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)], out_array: [bool](https://docs.python.org/3/builtins/functions.html#bool) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, \*\*kwargs: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [HighsVarType](#highspy._core.HighsVarType) | [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [Any](https://docs.python.org/3/library/typing.html#typing.Any)] | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_var](#highspy.highs.highs_var)] | [HighspyArray](#highspy.highs.HighspyArray) | [None](https://docs.python.org/3/builtins/constants.html#None)

Alias for the addVariables method, for binary variables.

* **Parameters:**
  * **nvars** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)
  * **out_array** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*)
  * **kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* [*HighsVarType*](#highspy._core.HighsVarType) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*)

<a id="highspy.highs.Highs.addIntegral"></a>

#### addIntegral(lb=0.0, ub=inf, obj=0.0, name=None)

Alias for the addVariable method, for integer variables.

* **Parameters:**
  * **lb** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **ub** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*)

<a id="highspy.highs.Highs.addBinary"></a>

#### addBinary(obj=0.0, name=None)

Alias for the addVariable method, for binary variables.

* **Parameters:**
  * **obj** ([*float*](https://docs.python.org/3/builtins/functions.html#float))
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*)

<a id="highspy.highs.Highs.deleteVariable"></a>

#### deleteVariable(var_or_index, \*args)

Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

* **Parameters:**
  * **var_or_index** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – A highs_var object or an index representing the variable to be deleted.
  * **\*args** ([*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*highs_var*](#highspy.highs.highs_var) *]*  *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*highs_var*](#highspy.highs.highs_var) *|* [*HighspyArray*](#highspy.highs.HighspyArray) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – Optional collections (lists, dicts, etc.) of highs_var objects whose indices need to be updated.

<a id="highspy.highs.Highs.getVariables"></a>

#### getVariables()

Retrieves all variables in the model.

* **Returns:**
  A list of highs_var objects, each representing a variable in the model.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[*highs_var*](#highspy.highs.highs_var)]

<a id="highspy.highs.Highs.inf"></a>

#### *property* inf *: [float](https://docs.python.org/3/builtins/functions.html#float)*

Represents infinity in the context of the solver.

* **Returns:**
  The value used to represent infinity.

<a id="highspy.highs.Highs.numVariables"></a>

#### *property* numVariables *: [int](https://docs.python.org/3/builtins/functions.html#int)*

Gets the number of variables in the model.

* **Returns:**
  The number of variables.

<a id="highspy.highs.Highs.numConstrs"></a>

#### *property* numConstrs *: [int](https://docs.python.org/3/builtins/functions.html#int)*

Gets the number of constraints in the model.

* **Returns:**
  The number of constraints.

<a id="highspy.highs.Highs.addConstr"></a>

#### addConstr(expr, name=None)

Adds a constraint to the model.

* **Parameters:**
  * **expr** ([*highs_linear_expression*](#highspy.highs.highs_linear_expression)) – A highs_linear_expression to be added.
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – Optional name of the constraint.
* **Returns:**
  A highs_cons object representing the added constraint.
* **Return type:**
  [*highs_cons*](#highspy.highs.highs_cons)

<a id="highspy.highs.Highs.addConstrs"></a>

#### addConstrs(\*args: [highs_linear_expression](#highspy.highs.highs_linear_expression), \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_linear_expression](#highspy.highs.highs_linear_expression)], \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Any](https://docs.python.org/3/library/typing.html#typing.Any), [highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[highs_linear_expression](#highspy.highs.highs_linear_expression)], \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

#### addConstrs(\*args: [HighspyArray](#highspy.highs.HighspyArray), \*\*kwargs: [str](https://docs.python.org/3/builtins/stdtypes.html#str) | [Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)] | [None](https://docs.python.org/3/builtins/constants.html#None)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[highs_cons](#highspy.highs.highs_cons)]

Adds multiple constraints to the model.

* **Parameters:**
  * **\*args** ([*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A sequence of highs_linear_expression to be added.
  * **\*\*kwargs** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *]*  *|* *None*) – Optional keyword arguments.
    name_prefix: Prefix for the constraint names.  Constructed name will be name_prefix + index.
    name: A collection of names for the constraints (list or mapping).
* **Returns:**
  A highs_con collection array representing the added constraints.

<a id="highspy.highs.Highs.expr"></a>

#### expr(optional=None)

Creates a new highs_linear_expression object.

* **Returns:**
  A highs_linear_expression object.
* **Parameters:**
  **optional** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs.Highs.getExpr"></a>

#### getExpr(cons)

Retrieves the highs_linear_expression of a constraint.

* **Parameters:**
  **cons** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_cons*](#highspy.highs.highs_cons)) – A highs_con object or index representing the constraint.
* **Returns:**
  A highs_linear_expression object representing the expression of the constraint.
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs.Highs.chgCoeff"></a>

#### chgCoeff(cons, var, val)

Changes the coefficient of a variable in a constraint.

* **Parameters:**
  * **cons** ([*highs_cons*](#highspy.highs.highs_cons) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_con object representing the constraint.
  * **var** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_var object representing the variable.
  * **val** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – The new coefficient value for the variable in the constraint.

<a id="highspy.highs.Highs.getConstrs"></a>

#### getConstrs()

Retrieves all constraints in the model.

* **Returns:**
  A list of highs_cons objects, each representing a constraint in the model.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[*highs_cons*](#highspy.highs.highs_cons)]

<a id="highspy.highs.Highs.removeConstr"></a>

#### removeConstr(cons_or_index, \*args)

Removes a constraint from the model and updates the indices of subsequent constraints in provided collections.

* **Parameters:**
  * **cons_or_index** ([*highs_cons*](#highspy.highs.highs_cons) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral)) – A highs_cons object or an index representing the constraint to be removed.
  * **\*args** ([*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*highs_cons*](#highspy.highs.highs_cons) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*highs_cons*](#highspy.highs.highs_cons) *]*  *|* [*highs_cons*](#highspy.highs.highs_cons)) – Optional collections (lists, dicts, etc.) of highs_cons objects whose indices need to be updated after the removal.

<a id="highspy.highs.Highs.setMinimize"></a>

#### setMinimize()

Sets the objective sense of the model to minimization.

<a id="highspy.highs.Highs.setMaximize"></a>

#### setMaximize()

Sets the objective sense of the model to maximization.

<a id="highspy.highs.Highs.setInteger"></a>

#### setInteger(var_or_collection)

Sets a variable/collection to integer.

* **Parameters:**
  **var_or_collection** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A highs_var object/collection representing the variable to be set as integer.

<a id="highspy.highs.Highs.setContinuous"></a>

#### setContinuous(var_or_collection)

Sets a variable/collection to continuous.

* **Parameters:**
  **var_or_collection** ([*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*HighspyArray*](#highspy.highs.HighspyArray)) – A highs_var object/collection representing the variable to be set as continuous.

<a id="highspy.highs.Highs.idx"></a>

#### *static* idx(\*args)

Convert highs_var/highs_cons to a flat int32 index array.

Can be called as:
: - `h.idx(array)` with a HighspyArray, numpy array, list, or tuple
  - `h.idx(a, b, c)` with individual highs_var or highs_cons objects

* **Returns:**
  A flat int32 numpy array of the underlying indices.
* **Return type:**
  [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*int32*]]

<a id="highspy.highs.Highs.qsum"></a>

#### *static* qsum(items, initial=None)

Performs a faster sum for highs_linear_expressions.

* **Parameters:**
  * **items** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[*[*object_*](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_) *]* *]*) – A collection of highs_linear_expressions or highs_vars to be summed.
  * **initial** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs.Highs.enableCallbacks"></a>

#### enableCallbacks()

Enables callbacks, restarting them if they were previously enabled.

<a id="highspy.highs.Highs.clearCallbacks"></a>

#### clearCallbacks()

Clears all callbacks.

<a id="highspy.highs.Highs.disableCallbacks"></a>

#### disableCallbacks()

Disables all callbacks, but does not clear them.

<a id="highspy.highs.Highs.cancelSolve"></a>

#### cancelSolve()

If HandleUserInterrupt is enabled, this method will signal the solver to stop.

<a id="highspy.highs.Highs.HandleKeyboardInterrupt"></a>

#### *property* HandleKeyboardInterrupt *: [bool](https://docs.python.org/3/builtins/functions.html#bool)*

Get/Set whether the solver should handle KeyboardInterrupt (i.e., cancel solve on Ctrl+C). Also enables/disables HandleUserInterrupt.

<a id="highspy.highs.Highs.HandleUserInterrupt"></a>

#### *property* HandleUserInterrupt *: [bool](https://docs.python.org/3/builtins/functions.html#bool)*

Get/Set whether the solver should handle user interrupts (i.e., cancel solve on user request)

<a id="highspy.highs.Highs.cbLogging"></a>

#### cbLogging *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbSimplexInterrupt"></a>

#### cbSimplexInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbIpmInterrupt"></a>

#### cbIpmInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipSolution"></a>

#### cbMipSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipImprovingSolution"></a>

#### cbMipImprovingSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipLogging"></a>

#### cbMipLogging *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipInterrupt"></a>

#### cbMipInterrupt *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipGetCutPool"></a>

#### cbMipGetCutPool *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipDefineLazyConstraints"></a>

#### cbMipDefineLazyConstraints *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.cbMipUserSolution"></a>

#### cbMipUserSolution *: [HighsCallback](#highspy.highs.HighsCallback)*

<a id="highspy.highs.Highs.addCol"></a>

#### addCol(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg3: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addCols"></a>

#### addCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg4: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg6: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg7: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addLinearObjective"></a>

#### addLinearObjective(self: highspy._core._Highs, arg0: [HighsLinearObjective](#highspy.HighsLinearObjective)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addRow"></a>

#### addRow(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addRows"></a>

#### addRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg4: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg5: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg6: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addVar"></a>

#### addVar(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.addVars"></a>

#### addVars(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.basisStatusToString"></a>

#### basisStatusToString(self: highspy._core._Highs, arg0: [highspy._core.HighsBasisStatus](#highspy._core.HighsBasisStatus)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.basisValidityToString"></a>

#### basisValidityToString(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.changeCoeff"></a>

#### changeCoeff(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColBounds"></a>

#### changeColBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColCost"></a>

#### changeColCost(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColIntegrality"></a>

#### changeColIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [highspy._core.HighsVarType](#highspy._core.HighsVarType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColsBounds"></a>

#### changeColsBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColsCost"></a>

#### changeColsCost(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeColsIntegrality"></a>

#### changeColsIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.uint8]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeObjectiveOffset"></a>

#### changeObjectiveOffset(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeObjectiveSense"></a>

#### changeObjectiveSense(self: highspy._core._Highs, arg0: [highspy._core.ObjSense](#highspy._core.ObjSense)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeRowBounds"></a>

#### changeRowBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.changeRowsBounds"></a>

#### changeRowsBounds(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64], arg3: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.clear"></a>

#### clear(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.clearLinearObjectives"></a>

#### clearLinearObjectives(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.clearModel"></a>

#### clearModel(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.clearSolver"></a>

#### clearSolver(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.crossover"></a>

#### crossover(self: highspy._core._Highs, arg0: [HighsSolution](#highspy.HighsSolution)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.deleteCols"></a>

#### deleteCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.deleteRows"></a>

#### deleteRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.deleteVars"></a>

#### deleteVars(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.ensureColwise"></a>

#### ensureColwise(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.ensureRowwise"></a>

#### ensureRowwise(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.feasibilityRelaxation"></a>

#### feasibilityRelaxation(self: highspy._core._Highs, global_lower_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), global_upper_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), global_rhs_penalty: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), local_lower_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None, local_upper_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None, local_rhs_penalty: [object](https://docs.python.org/3/builtins/functions.html#object) = None) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.getBasicVariables"></a>

#### getBasicVariables(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getBasis"></a>

#### getBasis(self: highspy._core._Highs) → [HighsBasis](#highspy.HighsBasis)

<a id="highspy.highs.Highs.getBasisInverseCol"></a>

#### getBasisInverseCol(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getBasisInverseColSparse"></a>

#### getBasisInverseColSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getBasisInverseRow"></a>

#### getBasisInverseRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getBasisInverseRowSparse"></a>

#### getBasisInverseRowSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getBasisSolve"></a>

#### getBasisSolve(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getBasisSolveSparse"></a>

#### getBasisSolveSparse(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getBasisTransposeSolve"></a>

#### getBasisTransposeSolve(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getBasisTransposeSolveSparse"></a>

#### getBasisTransposeSolveSparse(self: highspy._core._Highs, arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getCol"></a>

#### getCol(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getColByName"></a>

#### getColByName(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getColEntries"></a>

#### getColEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getColIntegrality"></a>

#### getColIntegrality(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsVarType](#highspy._core.HighsVarType)]

<a id="highspy.highs.Highs.getColName"></a>

#### getColName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.highs.Highs.getCols"></a>

#### getCols(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getColsEntries"></a>

#### getColsEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getDualObjectiveValue"></a>

#### getDualObjectiveValue(self: highspy._core._Highs, arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.getDualRay"></a>

#### getDualRay(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getDualRayExist"></a>

#### getDualRayExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.highs.Highs.getDualUnboundednessDirection"></a>

#### getDualUnboundednessDirection(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getDualUnboundednessDirectionExist"></a>

#### getDualUnboundednessDirectionExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.highs.Highs.getFixedLp"></a>

#### getFixedLp(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsLp](#highspy._core.HighsLp)]

<a id="highspy.highs.Highs.getHessianNumNz"></a>

#### getHessianNumNz(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.getIis"></a>

#### getIis(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [HighsIis](#highspy.HighsIis)]

<a id="highspy.highs.Highs.getInfinity"></a>

#### getInfinity(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.highs.Highs.getInfo"></a>

#### getInfo(self: highspy._core._Highs) → [highspy._core.HighsInfo](#highspy._core.HighsInfo)

<a id="highspy.highs.Highs.getInfoType"></a>

#### getInfoType(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsInfoType](#highspy._core.HighsInfoType)]

<a id="highspy.highs.Highs.getInfoValue"></a>

#### getInfoValue(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [object](https://docs.python.org/3/builtins/functions.html#object)]

<a id="highspy.highs.Highs.getLinearObjective"></a>

#### getLinearObjective(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [HighsLinearObjective](#highspy.HighsLinearObjective)

<a id="highspy.highs.Highs.getLp"></a>

#### getLp(self: highspy._core._Highs) → [highspy._core.HighsLp](#highspy._core.HighsLp)

<a id="highspy.highs.Highs.getModel"></a>

#### getModel(self: highspy._core._Highs) → [highspy._core.HighsModel](#highspy._core.HighsModel)

<a id="highspy.highs.Highs.getModelPresolveStatus"></a>

#### getModelPresolveStatus(self: highspy._core._Highs) → [highspy._core.HighsPresolveStatus](#highspy._core.HighsPresolveStatus)

<a id="highspy.highs.Highs.getModelStatus"></a>

#### getModelStatus(self: highspy._core._Highs) → [highspy._core.HighsModelStatus](#highspy._core.HighsModelStatus)

<a id="highspy.highs.Highs.getNumCol"></a>

#### getNumCol(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.getNumLinearObjectives"></a>

#### getNumLinearObjectives(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.getNumNz"></a>

#### getNumNz(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.getNumRow"></a>

#### getNumRow(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.getObjectiveOffset"></a>

#### getObjectiveOffset(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float)]

<a id="highspy.highs.Highs.getObjectiveSense"></a>

#### getObjectiveSense(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.ObjSense](#highspy._core.ObjSense)]

<a id="highspy.highs.Highs.getObjectiveValue"></a>

#### getObjectiveValue(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.highs.Highs.getOptionType"></a>

#### getOptionType(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [highspy._core.HighsOptionType](#highspy._core.HighsOptionType)]

<a id="highspy.highs.Highs.getOptionValue"></a>

#### getOptionValue(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [object](https://docs.python.org/3/builtins/functions.html#object)]

<a id="highspy.highs.Highs.getOptions"></a>

#### getOptions(self: highspy._core._Highs) → [highspy._core.HighsOptions](#highspy._core.HighsOptions)

<a id="highspy.highs.Highs.getPresolvedLp"></a>

#### getPresolvedLp(self: highspy._core._Highs) → [highspy._core.HighsLp](#highspy._core.HighsLp)

<a id="highspy.highs.Highs.getPrimalRay"></a>

#### getPrimalRay(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getPrimalRayExist"></a>

#### getPrimalRayExist(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [bool](https://docs.python.org/3/builtins/functions.html#bool)]

<a id="highspy.highs.Highs.getRanging"></a>

#### getRanging(self: highspy._core._Highs) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [HighsRanging](#highspy.HighsRanging)]

<a id="highspy.highs.Highs.getReducedColumn"></a>

#### getReducedColumn(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getReducedColumnSparse"></a>

#### getReducedColumnSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getReducedRow"></a>

#### getReducedRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getReducedRowSparse"></a>

#### getReducedRowSparse(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.int32]]

<a id="highspy.highs.Highs.getRow"></a>

#### getRow(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getRowByName"></a>

#### getRowByName(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getRowEntries"></a>

#### getRowEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getRowName"></a>

#### getRowName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [str](https://docs.python.org/3/builtins/stdtypes.html#str)]

<a id="highspy.highs.Highs.getRows"></a>

#### getRows(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), [int](https://docs.python.org/3/builtins/functions.html#int), numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64], [int](https://docs.python.org/3/builtins/functions.html#int)]

<a id="highspy.highs.Highs.getRowsEntries"></a>

#### getRowsEntries(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.int32]) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[highspy._core.HighsStatus](#highspy._core.HighsStatus), numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.int32], numpy.typing.NDArray[numpy.float64]]

<a id="highspy.highs.Highs.getRunTime"></a>

#### getRunTime(self: highspy._core._Highs) → [float](https://docs.python.org/3/builtins/functions.html#float)

<a id="highspy.highs.Highs.getSavedMipSolutions"></a>

#### getSavedMipSolutions(self: highspy._core._Highs) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[HighsObjectiveSolution](#highspy.HighsObjectiveSolution)]

<a id="highspy.highs.Highs.getSolution"></a>

#### getSolution(self: highspy._core._Highs) → [HighsSolution](#highspy.HighsSolution)

<a id="highspy.highs.Highs.getThirdPartyNotice"></a>

#### getThirdPartyNotice(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.githash"></a>

#### githash(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.modelStatusToString"></a>

#### modelStatusToString(self: highspy._core._Highs, arg0: [highspy._core.HighsModelStatus](#highspy._core.HighsModelStatus)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.passColName"></a>

#### passColName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.passHessian"></a>

#### passHessian(\*args, \*\*kwargs)

Overloaded function.

1. passHessian(self: highspy._core._Highs, arg0: highspy._core.HighsHessian) -> highspy._core.HighsStatus
2. passHessian(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg4: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg5: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.passModel"></a>

#### passModel(\*args, \*\*kwargs)

Overloaded function.

1. passModel(self: highspy._core._Highs, arg0: highspy._core.HighsModel) -> highspy._core.HighsStatus
2. passModel(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.SupportsInt | typing.SupportsIndex, arg4: typing.SupportsInt | typing.SupportsIndex, arg5: typing.SupportsInt | typing.SupportsIndex, arg6: typing.SupportsInt | typing.SupportsIndex, arg7: typing.SupportsFloat | typing.SupportsIndex, arg8: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg9: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg10: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg11: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg12: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg13: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg14: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg15: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg16: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg17: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg18: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg19: typing.Annotated[numpy.typing.ArrayLike, numpy.int32]) -> highspy._core.HighsStatus
3. passModel(self: highspy._core._Highs, arg0: highspy._core.HighsLp) -> highspy._core.HighsStatus
4. passModel(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.SupportsInt | typing.SupportsIndex, arg4: typing.SupportsInt | typing.SupportsIndex, arg5: typing.SupportsFloat | typing.SupportsIndex, arg6: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg7: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg8: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg9: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg10: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg11: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg12: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg13: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg14: typing.Annotated[numpy.typing.ArrayLike, numpy.int32]) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.passOptions"></a>

#### passOptions(self: highspy._core._Highs, arg0: [highspy._core.HighsOptions](#highspy._core.HighsOptions)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.passRowName"></a>

#### passRowName(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arg1: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.postsolve"></a>

#### postsolve(\*args, \*\*kwargs)

Overloaded function.

1. postsolve(self: highspy._core._Highs, arg0: HighsSolution, arg1: HighsBasis) -> highspy._core.HighsStatus
2. postsolve(self: highspy._core._Highs, arg0: HighsSolution) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.presolve"></a>

#### presolve(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.readBasis"></a>

#### readBasis(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.readModel"></a>

#### readModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.readOptions"></a>

#### readOptions(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.readSolution"></a>

#### readSolution(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.releaseMemory"></a>

#### releaseMemory(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.resetGlobalScheduler"></a>

#### *static* resetGlobalScheduler(arg0: [bool](https://docs.python.org/3/builtins/functions.html#bool)) → [None](https://docs.python.org/3/builtins/constants.html#None)

<a id="highspy.highs.Highs.resetOptions"></a>

#### resetOptions(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.run"></a>

#### run(self: highspy._core._Highs) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.setBasis"></a>

#### setBasis(\*args, \*\*kwargs)

Overloaded function.

1. setBasis(self: highspy._core._Highs, arg0: HighsBasis) -> highspy._core.HighsStatus
2. setBasis(self: highspy._core._Highs) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.setCallback"></a>

#### setCallback(self: highspy._core._Highs, arg0: [collections.abc.Callable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable)[[[int](https://docs.python.org/3/builtins/functions.html#int), [str](https://docs.python.org/3/builtins/stdtypes.html#str), [HighsCallbackOutput](#highspy._core.cb.HighsCallbackOutput), [HighsCallbackInput](#highspy._core.cb.HighsCallbackInput), [object](https://docs.python.org/3/builtins/functions.html#object)], [None](https://docs.python.org/3/builtins/constants.html#None)], arg1: [object](https://docs.python.org/3/builtins/functions.html#object)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.setOptionValue"></a>

#### setOptionValue(\*args, \*\*kwargs)

Overloaded function.

1. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: bool) -> highspy._core.HighsStatus
2. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: typing.SupportsInt | typing.SupportsIndex) -> highspy._core.HighsStatus
3. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: typing.SupportsFloat | typing.SupportsIndex) -> highspy._core.HighsStatus
4. setOptionValue(self: highspy._core._Highs, arg0: str, arg1: str) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.setSolution"></a>

#### setSolution(\*args, \*\*kwargs)

Overloaded function.

1. setSolution(self: highspy._core._Highs, arg0: HighsSolution) -> highspy._core.HighsStatus
2. setSolution(self: highspy._core._Highs, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus

<a id="highspy.highs.Highs.solutionStatusToString"></a>

#### solutionStatusToString(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.startCallback"></a>

#### startCallback(self: highspy._core._Highs, arg0: [HighsCallbackType](#highspy._core.cb.HighsCallbackType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.startCallbackInt"></a>

#### startCallbackInt(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.stopCallback"></a>

#### stopCallback(self: highspy._core._Highs, arg0: [HighsCallbackType](#highspy._core.cb.HighsCallbackType)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.stopCallbackInt"></a>

#### stopCallbackInt(self: highspy._core._Highs, arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.version"></a>

#### version(self: highspy._core._Highs) → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="highspy.highs.Highs.versionMajor"></a>

#### versionMajor(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.versionMinor"></a>

#### versionMinor(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.versionPatch"></a>

#### versionPatch(self: highspy._core._Highs) → [int](https://docs.python.org/3/builtins/functions.html#int)

<a id="highspy.highs.Highs.writeBasis"></a>

#### writeBasis(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writeIisModel"></a>

#### writeIisModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writeInfo"></a>

#### writeInfo(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writeModel"></a>

#### writeModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writeOptions"></a>

#### writeOptions(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writePresolvedModel"></a>

#### writePresolvedModel(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.Highs.writeSolution"></a>

#### writeSolution(self: highspy._core._Highs, arg0: [str](https://docs.python.org/3/builtins/stdtypes.html#str), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy.highs.HighsCallbackEvent"></a>

### *class* highspy.highs.HighsCallbackEvent(callback_type, message, data_out, data_in, user_data)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

* **Parameters:**
  * **callback_type** ([*cb.HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **message** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **data_out** ([*cb.HighsCallbackOutput*](#highspy._core.cb.HighsCallbackOutput))
  * **data_in** ([*cb.HighsCallbackInput*](#highspy._core.cb.HighsCallbackInput) *|* *None*)
  * **user_data** (*Any* *|* *None*)

<a id="highspy.highs.HighsCallbackEvent.callback_type"></a>

#### callback_type

<a id="highspy.highs.HighsCallbackEvent.message"></a>

#### message

<a id="highspy.highs.HighsCallbackEvent.data_out"></a>

#### data_out

<a id="highspy.highs.HighsCallbackEvent.data_in"></a>

#### data_in

<a id="highspy.highs.HighsCallbackEvent.user_data"></a>

#### user_data

<a id="highspy.highs.HighsCallbackEvent.interrupt"></a>

#### interrupt(interrupt_value=True)

Sets the user interrupt flag in the callback data.

* **Parameters:**
  **interrupt_value** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))

<a id="highspy.highs.HighsCallbackEvent.val"></a>

#### val(var_expr)

Gets the value(s) of a variable/index or expression in the callback solution.

* **Parameters:**
  **var_expr** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*Integral*](https://docs.python.org/3/library/numbers.html#numbers.Integral) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_cons*](#highspy.highs.highs_cons) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*Mapping*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray))

<a id="highspy.highs.HighsCallbackEvent.cut"></a>

#### cut(index)

Gets the cut pool for the given index.

* **Parameters:**
  **index** ([*int*](https://docs.python.org/3/builtins/functions.html#int))

<a id="highspy.highs.HighsCallbackEvent.cuts"></a>

#### *property* cuts

Gets all cuts in the cut pool.

<a id="highspy.highs.HighsCallback"></a>

### *class* highspy.highs.HighsCallback(callback_type, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

* **Parameters:**
  * **callback_type** ([*cb.HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **highs** ([*Highs*](#highspy.highs.Highs))

<a id="highspy.highs.HighsCallback.callbacks"></a>

#### callbacks *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[Callable](https://docs.python.org/3/library/typing.html#typing.Callable)[[[HighsCallbackEvent](#highspy.highs.HighsCallbackEvent)], [None](https://docs.python.org/3/builtins/constants.html#None)]]*

<a id="highspy.highs.HighsCallback.user_callback_data"></a>

#### user_callback_data *: [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[Callable](https://docs.python.org/3/library/typing.html#typing.Callable)[[[HighsCallbackEvent](#highspy.highs.HighsCallbackEvent)], [None](https://docs.python.org/3/builtins/constants.html#None)], [Any](https://docs.python.org/3/library/typing.html#typing.Any)]*

<a id="highspy.highs.HighsCallback.callback_type"></a>

#### callback_type

<a id="highspy.highs.HighsCallback.highs"></a>

#### highs

<a id="highspy.highs.HighsCallback.subscribe"></a>

#### subscribe(callback, user_data=None)

Subscribes a callback to the event.

* **Parameters:**
  * **callback** ([*Callable*](https://docs.python.org/3/library/typing.html#typing.Callable) *[* *[*[*HighsCallbackEvent*](#highspy.highs.HighsCallbackEvent) *]* *,* *None* *]*) – The callback function to be executed.
  * **user_data** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*) – Optional user data to be passed to the callback.

<a id="highspy.highs.HighsCallback.unsubscribe"></a>

#### unsubscribe(callback)

Unsubscribes a callback from the event.

* **Parameters:**
  **callback** ([*Callable*](https://docs.python.org/3/library/typing.html#typing.Callable) *[* *[*[*HighsCallbackEvent*](#highspy.highs.HighsCallbackEvent) *]* *,* *None* *]*) – The callback function to be removed.

<a id="highspy.highs.HighsCallback.unsubscribe_by_data"></a>

#### unsubscribe_by_data(user_data)

Unsubscribes a callback by user data.

* **Parameters:**
  **user_data** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*) – The user data corresponding to the callback(s) to be removed.

<a id="highspy.highs.HighsCallback.clear"></a>

#### clear()

Unsubscribes all callbacks from the event.

<a id="highspy.highs.HighsCallback.fire"></a>

#### fire(callback_type, message, data_out, data_in)

Fires the event, executing all subscribed callbacks.

* **Parameters:**
  * **callback_type** ([*HighsCallbackType*](#highspy._core.cb.HighsCallbackType))
  * **message** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str))
  * **data_out** ([*HighsCallbackOutput*](#highspy._core.cb.HighsCallbackOutput))
  * **data_in** ([*HighsCallbackInput*](#highspy._core.cb.HighsCallbackInput))

<a id="highspy.highs.HighspyArray"></a>

### *class* highspy.highs.HighspyArray(input_array, highs)

Bases: [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[`Any`](https://docs.python.org/3/library/typing.html#typing.Any), [`dtype`](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[[`object_`](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_)]]

A numpy array wrapper for highs_var/highs_linear_expression objects.

This provides additional type information for static analysis, and also allows faster sum operations.

* **Parameters:**
  * **input_array** (*np.ndarray* *[**Any* *,* *np.dtype* *[**np.object_* *]* *]*)
  * **highs** ([*Highs*](#highspy.highs.Highs) *|* *None*)
* **Return type:**
  Self

<a id="highspy.highs.HighspyArray.highs"></a>

#### highs *: [Highs](#highspy.highs.Highs) | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.highs.HighspyArray.sum"></a>

#### sum(axis: [None](https://docs.python.org/3/builtins/constants.html#None) = None, dtype: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, out: [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [highs_linear_expression](#highspy.highs.highs_linear_expression)

#### sum(axis: [Any](https://docs.python.org/3/library/typing.html#typing.Any), dtype: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, out: [HighspyArray](#highspy.highs.HighspyArray) = None) → [HighspyArray](#highspy.highs.HighspyArray)

Return the sum of the array elements over the given axis.

Refer to numpy.sum for full documentation.

#### SEE ALSO
[`numpy.sum`](https://numpy.org/doc/stable/reference/generated/numpy.sum.html#numpy.sum)
: equivalent function

* **Parameters:**
  * **axis** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)
  * **dtype** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *|* *None*)
  * **out** ([*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[*[*object_*](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.object_) *]* *]*  *|* *None*)
  * **unused_kwargs** ([*Any*](https://docs.python.org/3/library/typing.html#typing.Any))
* **Return type:**
  [*HighspyArray*](#highspy.highs.HighspyArray) | [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs.HighspyArray.idx"></a>

#### idx()

Convert to a flat int32 index array for passing to the HiGHS C++ API.

Each element’s `__index__` method is called to extract its integer
index (e.g., `highs_var.index` or `highs_cons.index`).
The result is always 1-D, regardless of the array’s shape.

* **Returns:**
  A new flat int32 numpy array of the underlying indices.
* **Return type:**
  [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)[[*Any*](https://docs.python.org/3/library/typing.html#typing.Any), [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype)[*int32*]]

<a id="highspy.highs.highs_var"></a>

### *class* highspy.highs.highs_var(i, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Variable index wrapper for HiGHS

* **Parameters:**
  * **i** ([*int*](https://docs.python.org/3/builtins/functions.html#int))
  * **highs** ([*Highs*](#highspy.highs.Highs))

<a id="highspy.highs.highs_var.index"></a>

#### index

<a id="highspy.highs.highs_var.highs"></a>

#### highs

<a id="highspy.highs.highs_var.name"></a>

#### *property* name *: [str](https://docs.python.org/3/builtins/stdtypes.html#str)*

<a id="highspy.highs.highs_cons"></a>

### *class* highspy.highs.highs_cons(i, highs)

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Constraint index wrapper for HiGHS

* **Parameters:**
  * **i** ([*int*](https://docs.python.org/3/builtins/functions.html#int))
  * **highs** ([*Highs*](#highspy.highs.Highs))

<a id="highspy.highs.highs_cons.index"></a>

#### index

<a id="highspy.highs.highs_cons.highs"></a>

#### highs

<a id="highspy.highs.highs_cons.expr"></a>

#### expr()

Retrieves the expression of the constraint.

* **Return type:**
  [*highs_linear_expression*](#highspy.highs.highs_linear_expression)

<a id="highspy.highs.highs_cons.name"></a>

#### *property* name *: [str](https://docs.python.org/3/builtins/stdtypes.html#str)*

<a id="highspy.highs.highs_linear_expression"></a>

### *class* highspy.highs.highs_linear_expression(other: [None](https://docs.python.org/3/builtins/constants.html#None) = None)

### *class* highspy.highs.highs_linear_expression(other: [float](https://docs.python.org/3/builtins/functions.html#float))

### *class* highspy.highs.highs_linear_expression(other: [highs_var](#highspy.highs.highs_var))

### *class* highspy.highs.highs_linear_expression(other: [highs_linear_expression](#highspy.highs.highs_linear_expression))

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Linear constraint builder for HiGHS

* **Parameters:**
  **other** ([*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* *None*)

<a id="highspy.highs.highs_linear_expression.bounds"></a>

#### bounds *: [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[float](https://docs.python.org/3/builtins/functions.html#float), [float](https://docs.python.org/3/builtins/functions.html#float)] | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.highs.highs_linear_expression.idxs"></a>

#### idxs *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]*

<a id="highspy.highs.highs_linear_expression.vals"></a>

#### vals *: [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[float](https://docs.python.org/3/builtins/functions.html#float)]*

<a id="highspy.highs.highs_linear_expression.constant"></a>

#### constant *: [float](https://docs.python.org/3/builtins/functions.html#float) | [None](https://docs.python.org/3/builtins/constants.html#None)*

<a id="highspy.highs.highs_linear_expression.simplify"></a>

#### simplify()

Simplifies the linear expression by combining duplicate variables.

<a id="highspy.highs.highs_linear_expression.copy"></a>

#### copy()

Creates a copy of the linear expression.

<a id="highspy.highs.highs_linear_expression.evaluate"></a>

#### evaluate(values)

Evaluates the linear expression given a solution array (values).

* **Parameters:**
  **values** ([*Sequence*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence) *[*[*float*](https://docs.python.org/3/builtins/functions.html#float) *]*  *|* [*ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray) *[*[*Any*](https://docs.python.org/3/library/typing.html#typing.Any) *,* [*dtype*](https://numpy.org/doc/stable/reference/generated/numpy.dtype.html#numpy.dtype) *[**float64* *]* *]*)
* **Return type:**
  [float](https://docs.python.org/3/builtins/functions.html#float) | [bool](https://docs.python.org/3/builtins/functions.html#bool)

<a id="highspy.highs.highs_linear_expression.unique_elements"></a>

#### unique_elements()

Collects unique variables and sums their corresponding values.  Keeps all values (including zeros).

<a id="highspy.highs.highs_linear_expression.reduced_elements"></a>

#### reduced_elements()

Similar to unique_elements, except keeps only non-zero values

<a id="highspy.highs.qsum"></a>

### highspy.highs.qsum(items, initial=None)

Performs a faster sum for highs_linear_expressions.

* **Parameters:**
  * **items** ([*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable) *[*[*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *]*) – A collection of highs_linear_expressions or highs_vars to be summed.
  * **initial** ([*highs_var*](#highspy.highs.highs_var) *|* [*highs_linear_expression*](#highspy.highs.highs_linear_expression) *|* [*float*](https://docs.python.org/3/builtins/functions.html#float) *|* [*int*](https://docs.python.org/3/builtins/functions.html#int) *|* *None*)

<a id="core-solver-bindings"></a>

## Core Solver Bindings

The [`highspy._core`](#module-highspy._core) module contains the pybind11 bindings to the
HiGHS C++ API. In an installed build these entries include the signatures
exposed by the compiled extension.

<a id="module-highspy._core"></a>

<a id="highspy._core.BasisValidity"></a>

### *class* highspy._core.BasisValidity

Bases: `pybind11_object`

Members:

kBasisValidityInvalid

kBasisValidityValid

<a id="highspy._core.BasisValidity.kBasisValidityInvalid"></a>

#### kBasisValidityInvalid *= <BasisValidity.kBasisValidityInvalid: 0>*

<a id="highspy._core.BasisValidity.kBasisValidityValid"></a>

#### kBasisValidityValid *= <BasisValidity.kBasisValidityValid: 1>*

### BasisValidity.name -> str

<a id="highspy._core.BasisValidity.value"></a>

#### *property* value

<a id="highspy._core.HessianFormat"></a>

### *class* highspy._core.HessianFormat

Bases: `pybind11_object`

Members:

kTriangular

kSquare

<a id="highspy._core.HessianFormat.kSquare"></a>

#### kSquare *= <HessianFormat.kSquare: 2>*

<a id="highspy._core.HessianFormat.kTriangular"></a>

#### kTriangular *= <HessianFormat.kTriangular: 1>*

### HessianFormat.name -> str

<a id="highspy._core.HessianFormat.value"></a>

#### *property* value

<a id="highspy._core.HighsBasis"></a>

### *class* highspy._core.HighsBasis

Bases: `pybind11_object`

<a id="highspy._core.HighsBasis.alien"></a>

#### *property* alien

<a id="highspy._core.HighsBasis.col_status"></a>

#### *property* col_status

<a id="highspy._core.HighsBasis.debug_id"></a>

#### *property* debug_id

<a id="highspy._core.HighsBasis.debug_origin_name"></a>

#### *property* debug_origin_name

<a id="highspy._core.HighsBasis.debug_update_count"></a>

#### *property* debug_update_count

<a id="highspy._core.HighsBasis.row_status"></a>

#### *property* row_status

<a id="highspy._core.HighsBasis.valid"></a>

#### *property* valid

<a id="highspy._core.HighsBasis.was_alien"></a>

#### *property* was_alien

<a id="highspy._core.HighsBasisStatus"></a>

### *class* highspy._core.HighsBasisStatus

Bases: `pybind11_object`

Members:

kLower

kBasic

kUpper

kZero

kNonbasic

<a id="highspy._core.HighsBasisStatus.kBasic"></a>

#### kBasic *= <HighsBasisStatus.kBasic: 1>*

<a id="highspy._core.HighsBasisStatus.kLower"></a>

#### kLower *= <HighsBasisStatus.kLower: 0>*

<a id="highspy._core.HighsBasisStatus.kNonbasic"></a>

#### kNonbasic *= <HighsBasisStatus.kNonbasic: 4>*

<a id="highspy._core.HighsBasisStatus.kUpper"></a>

#### kUpper *= <HighsBasisStatus.kUpper: 2>*

<a id="highspy._core.HighsBasisStatus.kZero"></a>

#### kZero *= <HighsBasisStatus.kZero: 3>*

### HighsBasisStatus.name -> str

<a id="highspy._core.HighsBasisStatus.value"></a>

#### *property* value

<a id="highspy._core.HighsDebugLevel"></a>

### *class* highspy._core.HighsDebugLevel

Bases: `pybind11_object`

Members:

kHighsDebugLevelNone

kHighsDebugLevelCheap

kHighsDebugLevelCostly

kHighsDebugLevelExpensive

kHighsDebugLevelMin

kHighsDebugLevelMax

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelCheap"></a>

#### kHighsDebugLevelCheap *= <HighsDebugLevel.kHighsDebugLevelCheap: 1>*

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelCostly"></a>

#### kHighsDebugLevelCostly *= <HighsDebugLevel.kHighsDebugLevelCostly: 2>*

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelExpensive"></a>

#### kHighsDebugLevelExpensive *= <HighsDebugLevel.kHighsDebugLevelExpensive: 3>*

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelMax"></a>

#### kHighsDebugLevelMax *= <HighsDebugLevel.kHighsDebugLevelExpensive: 3>*

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelMin"></a>

#### kHighsDebugLevelMin *= <HighsDebugLevel.kHighsDebugLevelNone: 0>*

<a id="highspy._core.HighsDebugLevel.kHighsDebugLevelNone"></a>

#### kHighsDebugLevelNone *= <HighsDebugLevel.kHighsDebugLevelNone: 0>*

### HighsDebugLevel.name -> str

<a id="highspy._core.HighsDebugLevel.value"></a>

#### *property* value

<a id="highspy._core.HighsHessian"></a>

### *class* highspy._core.HighsHessian

Bases: `pybind11_object`

<a id="highspy._core.HighsHessian.dim_"></a>

#### *property* dim_

<a id="highspy._core.HighsHessian.format_"></a>

#### *property* format_

<a id="highspy._core.HighsHessian.index_"></a>

#### *property* index_

<a id="highspy._core.HighsHessian.start_"></a>

#### *property* start_

<a id="highspy._core.HighsHessian.value_"></a>

#### *property* value_

<a id="highspy._core.HighsIis"></a>

### *class* highspy._core.HighsIis

Bases: `pybind11_object`

<a id="highspy._core.HighsIis.clear"></a>

#### clear(self: [highspy._core.HighsIis](#highspy._core.HighsIis)) → [None](https://docs.python.org/3/builtins/constants.html#None)

<a id="highspy._core.HighsIis.col_bound_"></a>

#### *property* col_bound_

<a id="highspy._core.HighsIis.col_index_"></a>

#### *property* col_index_

<a id="highspy._core.HighsIis.col_status_"></a>

#### *property* col_status_

<a id="highspy._core.HighsIis.info_"></a>

#### *property* info_

<a id="highspy._core.HighsIis.model_"></a>

#### *property* model_

<a id="highspy._core.HighsIis.row_bound_"></a>

#### *property* row_bound_

<a id="highspy._core.HighsIis.row_index_"></a>

#### *property* row_index_

<a id="highspy._core.HighsIis.row_status_"></a>

#### *property* row_status_

<a id="highspy._core.HighsIis.status_"></a>

#### *property* status_

<a id="highspy._core.HighsIis.strategy_"></a>

#### *property* strategy_

<a id="highspy._core.HighsIis.valid_"></a>

#### *property* valid_

<a id="highspy._core.HighsIisInfo"></a>

### *class* highspy._core.HighsIisInfo

Bases: `pybind11_object`

<a id="highspy._core.HighsIisInfo.max_simplex_iteration_count"></a>

#### *property* max_simplex_iteration_count

<a id="highspy._core.HighsIisInfo.max_simplex_time"></a>

#### *property* max_simplex_time

<a id="highspy._core.HighsIisInfo.min_simplex_iteration_count"></a>

#### *property* min_simplex_iteration_count

<a id="highspy._core.HighsIisInfo.min_simplex_time"></a>

#### *property* min_simplex_time

<a id="highspy._core.HighsIisInfo.num_lp_solved"></a>

#### *property* num_lp_solved

<a id="highspy._core.HighsIisInfo.sum_simplex_iteration_counts"></a>

#### *property* sum_simplex_iteration_counts

<a id="highspy._core.HighsIisInfo.sum_simplex_times"></a>

#### *property* sum_simplex_times

<a id="highspy._core.HighsInfo"></a>

### *class* highspy._core.HighsInfo

Bases: `pybind11_object`

<a id="highspy._core.HighsInfo.basis_validity"></a>

#### *property* basis_validity

<a id="highspy._core.HighsInfo.crossover_iteration_count"></a>

#### *property* crossover_iteration_count

<a id="highspy._core.HighsInfo.dual_solution_status"></a>

#### *property* dual_solution_status

<a id="highspy._core.HighsInfo.ipm_iteration_count"></a>

#### *property* ipm_iteration_count

<a id="highspy._core.HighsInfo.max_complementarity_violation"></a>

#### *property* max_complementarity_violation

<a id="highspy._core.HighsInfo.max_dual_infeasibility"></a>

#### *property* max_dual_infeasibility

<a id="highspy._core.HighsInfo.max_dual_residual_error"></a>

#### *property* max_dual_residual_error

<a id="highspy._core.HighsInfo.max_integrality_violation"></a>

#### *property* max_integrality_violation

<a id="highspy._core.HighsInfo.max_primal_infeasibility"></a>

#### *property* max_primal_infeasibility

<a id="highspy._core.HighsInfo.max_primal_residual_error"></a>

#### *property* max_primal_residual_error

<a id="highspy._core.HighsInfo.max_relative_dual_infeasibility"></a>

#### *property* max_relative_dual_infeasibility

<a id="highspy._core.HighsInfo.max_relative_dual_residual_error"></a>

#### *property* max_relative_dual_residual_error

<a id="highspy._core.HighsInfo.max_relative_primal_infeasibility"></a>

#### *property* max_relative_primal_infeasibility

<a id="highspy._core.HighsInfo.max_relative_primal_residual_error"></a>

#### *property* max_relative_primal_residual_error

<a id="highspy._core.HighsInfo.mip_dual_bound"></a>

#### *property* mip_dual_bound

<a id="highspy._core.HighsInfo.mip_gap"></a>

#### *property* mip_gap

<a id="highspy._core.HighsInfo.mip_node_count"></a>

#### *property* mip_node_count

<a id="highspy._core.HighsInfo.num_complementarity_violations"></a>

#### *property* num_complementarity_violations

<a id="highspy._core.HighsInfo.num_dual_infeasibilities"></a>

#### *property* num_dual_infeasibilities

<a id="highspy._core.HighsInfo.num_dual_residual_errors"></a>

#### *property* num_dual_residual_errors

<a id="highspy._core.HighsInfo.num_primal_infeasibilities"></a>

#### *property* num_primal_infeasibilities

<a id="highspy._core.HighsInfo.num_primal_residual_errors"></a>

#### *property* num_primal_residual_errors

<a id="highspy._core.HighsInfo.num_relative_dual_infeasibilities"></a>

#### *property* num_relative_dual_infeasibilities

<a id="highspy._core.HighsInfo.num_relative_dual_residual_errors"></a>

#### *property* num_relative_dual_residual_errors

<a id="highspy._core.HighsInfo.num_relative_primal_infeasibilities"></a>

#### *property* num_relative_primal_infeasibilities

<a id="highspy._core.HighsInfo.num_relative_primal_residual_errors"></a>

#### *property* num_relative_primal_residual_errors

<a id="highspy._core.HighsInfo.objective_function_value"></a>

#### *property* objective_function_value

<a id="highspy._core.HighsInfo.pdlp_iteration_count"></a>

#### *property* pdlp_iteration_count

<a id="highspy._core.HighsInfo.primal_dual_integral"></a>

#### *property* primal_dual_integral

<a id="highspy._core.HighsInfo.primal_dual_objective_error"></a>

#### *property* primal_dual_objective_error

<a id="highspy._core.HighsInfo.primal_solution_status"></a>

#### *property* primal_solution_status

<a id="highspy._core.HighsInfo.qp_iteration_count"></a>

#### *property* qp_iteration_count

<a id="highspy._core.HighsInfo.simplex_iteration_count"></a>

#### *property* simplex_iteration_count

<a id="highspy._core.HighsInfo.sum_dual_infeasibilities"></a>

#### *property* sum_dual_infeasibilities

<a id="highspy._core.HighsInfo.sum_primal_infeasibilities"></a>

#### *property* sum_primal_infeasibilities

<a id="highspy._core.HighsInfo.valid"></a>

#### *property* valid

<a id="highspy._core.HighsInfoType"></a>

### *class* highspy._core.HighsInfoType

Bases: `pybind11_object`

Members:

kInt64

kInt

kDouble

<a id="highspy._core.HighsInfoType.kDouble"></a>

#### kDouble *= <HighsInfoType.kDouble: 2>*

<a id="highspy._core.HighsInfoType.kInt"></a>

#### kInt *= <HighsInfoType.kInt: 1>*

<a id="highspy._core.HighsInfoType.kInt64"></a>

#### kInt64 *= <HighsInfoType.kInt64: -1>*

### HighsInfoType.name -> str

<a id="highspy._core.HighsInfoType.value"></a>

#### *property* value

<a id="highspy._core.HighsLinearObjective"></a>

### *class* highspy._core.HighsLinearObjective

Bases: `pybind11_object`

<a id="highspy._core.HighsLinearObjective.abs_tolerance"></a>

#### *property* abs_tolerance

<a id="highspy._core.HighsLinearObjective.coefficients"></a>

#### *property* coefficients

<a id="highspy._core.HighsLinearObjective.offset"></a>

#### *property* offset

<a id="highspy._core.HighsLinearObjective.priority"></a>

#### *property* priority

<a id="highspy._core.HighsLinearObjective.rel_tolerance"></a>

#### *property* rel_tolerance

<a id="highspy._core.HighsLinearObjective.weight"></a>

#### *property* weight

<a id="highspy._core.HighsLogType"></a>

### *class* highspy._core.HighsLogType

Bases: `pybind11_object`

Members:

kInfo

kDetailed

kVerbose

kWarning

kError

<a id="highspy._core.HighsLogType.kDetailed"></a>

#### kDetailed *= <HighsLogType.kDetailed: 2>*

<a id="highspy._core.HighsLogType.kError"></a>

#### kError *= <HighsLogType.kError: 5>*

<a id="highspy._core.HighsLogType.kInfo"></a>

#### kInfo *= <HighsLogType.kInfo: 1>*

<a id="highspy._core.HighsLogType.kVerbose"></a>

#### kVerbose *= <HighsLogType.kVerbose: 3>*

<a id="highspy._core.HighsLogType.kWarning"></a>

#### kWarning *= <HighsLogType.kWarning: 4>*

### HighsLogType.name -> str

<a id="highspy._core.HighsLogType.value"></a>

#### *property* value

<a id="highspy._core.HighsLp"></a>

### *class* highspy._core.HighsLp

Bases: `pybind11_object`

<a id="highspy._core.HighsLp.a_matrix_"></a>

#### *property* a_matrix_

<a id="highspy._core.HighsLp.col_cost_"></a>

#### *property* col_cost_

<a id="highspy._core.HighsLp.col_lower_"></a>

#### *property* col_lower_

<a id="highspy._core.HighsLp.col_names_"></a>

#### *property* col_names_

<a id="highspy._core.HighsLp.col_upper_"></a>

#### *property* col_upper_

<a id="highspy._core.HighsLp.integrality_"></a>

#### *property* integrality_

<a id="highspy._core.HighsLp.is_moved_"></a>

#### *property* is_moved_

<a id="highspy._core.HighsLp.is_scaled_"></a>

#### *property* is_scaled_

<a id="highspy._core.HighsLp.model_name_"></a>

#### *property* model_name_

<a id="highspy._core.HighsLp.mods_"></a>

#### *property* mods_

<a id="highspy._core.HighsLp.num_col_"></a>

#### *property* num_col_

<a id="highspy._core.HighsLp.num_row_"></a>

#### *property* num_row_

<a id="highspy._core.HighsLp.offset_"></a>

#### *property* offset_

<a id="highspy._core.HighsLp.row_lower_"></a>

#### *property* row_lower_

<a id="highspy._core.HighsLp.row_names_"></a>

#### *property* row_names_

<a id="highspy._core.HighsLp.row_upper_"></a>

#### *property* row_upper_

<a id="highspy._core.HighsLp.scale_"></a>

#### *property* scale_

<a id="highspy._core.HighsLp.sense_"></a>

#### *property* sense_

<a id="highspy._core.HighsLpMods"></a>

### *class* highspy._core.HighsLpMods

Bases: `pybind11_object`

<a id="highspy._core.HighsModel"></a>

### *class* highspy._core.HighsModel

Bases: `pybind11_object`

<a id="highspy._core.HighsModel.hessian_"></a>

#### *property* hessian_

<a id="highspy._core.HighsModel.lp_"></a>

#### *property* lp_

<a id="highspy._core.HighsModelStatus"></a>

### *class* highspy._core.HighsModelStatus

Bases: `pybind11_object`

Members:

kNotset

kLoadError

kModelError

kPresolveError

kSolveError

kPostsolveError

kModelEmpty

kOptimal

kInfeasible

kUnboundedOrInfeasible

kUnbounded

kObjectiveBound

kObjectiveTarget

kTimeLimit

kIterationLimit

kUnknown

kSolutionLimit

kInterrupt

kMemoryLimit

kHighsInterrupt

<a id="highspy._core.HighsModelStatus.kHighsInterrupt"></a>

#### kHighsInterrupt *= <HighsModelStatus.kHighsInterrupt: 19>*

<a id="highspy._core.HighsModelStatus.kInfeasible"></a>

#### kInfeasible *= <HighsModelStatus.kInfeasible: 8>*

<a id="highspy._core.HighsModelStatus.kInterrupt"></a>

#### kInterrupt *= <HighsModelStatus.kInterrupt: 17>*

<a id="highspy._core.HighsModelStatus.kIterationLimit"></a>

#### kIterationLimit *= <HighsModelStatus.kIterationLimit: 14>*

<a id="highspy._core.HighsModelStatus.kLoadError"></a>

#### kLoadError *= <HighsModelStatus.kLoadError: 1>*

<a id="highspy._core.HighsModelStatus.kMemoryLimit"></a>

#### kMemoryLimit *= <HighsModelStatus.kMemoryLimit: 18>*

<a id="highspy._core.HighsModelStatus.kModelEmpty"></a>

#### kModelEmpty *= <HighsModelStatus.kModelEmpty: 6>*

<a id="highspy._core.HighsModelStatus.kModelError"></a>

#### kModelError *= <HighsModelStatus.kModelError: 2>*

<a id="highspy._core.HighsModelStatus.kNotset"></a>

#### kNotset *= <HighsModelStatus.kNotset: 0>*

<a id="highspy._core.HighsModelStatus.kObjectiveBound"></a>

#### kObjectiveBound *= <HighsModelStatus.kObjectiveBound: 11>*

<a id="highspy._core.HighsModelStatus.kObjectiveTarget"></a>

#### kObjectiveTarget *= <HighsModelStatus.kObjectiveTarget: 12>*

<a id="highspy._core.HighsModelStatus.kOptimal"></a>

#### kOptimal *= <HighsModelStatus.kOptimal: 7>*

<a id="highspy._core.HighsModelStatus.kPostsolveError"></a>

#### kPostsolveError *= <HighsModelStatus.kPostsolveError: 5>*

<a id="highspy._core.HighsModelStatus.kPresolveError"></a>

#### kPresolveError *= <HighsModelStatus.kPresolveError: 3>*

<a id="highspy._core.HighsModelStatus.kSolutionLimit"></a>

#### kSolutionLimit *= <HighsModelStatus.kSolutionLimit: 16>*

<a id="highspy._core.HighsModelStatus.kSolveError"></a>

#### kSolveError *= <HighsModelStatus.kSolveError: 4>*

<a id="highspy._core.HighsModelStatus.kTimeLimit"></a>

#### kTimeLimit *= <HighsModelStatus.kTimeLimit: 13>*

<a id="highspy._core.HighsModelStatus.kUnbounded"></a>

#### kUnbounded *= <HighsModelStatus.kUnbounded: 10>*

<a id="highspy._core.HighsModelStatus.kUnboundedOrInfeasible"></a>

#### kUnboundedOrInfeasible *= <HighsModelStatus.kUnboundedOrInfeasible: 9>*

<a id="highspy._core.HighsModelStatus.kUnknown"></a>

#### kUnknown *= <HighsModelStatus.kUnknown: 15>*

### HighsModelStatus.name -> str

<a id="highspy._core.HighsModelStatus.value"></a>

#### *property* value

<a id="highspy._core.HighsObjectiveSolution"></a>

### *class* highspy._core.HighsObjectiveSolution

Bases: `pybind11_object`

<a id="highspy._core.HighsObjectiveSolution.col_value"></a>

#### *property* col_value

<a id="highspy._core.HighsObjectiveSolution.objective"></a>

#### *property* objective

<a id="highspy._core.HighsOptionType"></a>

### *class* highspy._core.HighsOptionType

Bases: `pybind11_object`

Members:

kBool

kInt

kDouble

kString

<a id="highspy._core.HighsOptionType.kBool"></a>

#### kBool *= <HighsOptionType.kBool: 0>*

<a id="highspy._core.HighsOptionType.kDouble"></a>

#### kDouble *= <HighsOptionType.kDouble: 2>*

<a id="highspy._core.HighsOptionType.kInt"></a>

#### kInt *= <HighsOptionType.kInt: 1>*

<a id="highspy._core.HighsOptionType.kString"></a>

#### kString *= <HighsOptionType.kString: 3>*

### HighsOptionType.name -> str

<a id="highspy._core.HighsOptionType.value"></a>

#### *property* value

<a id="highspy._core.HighsOptions"></a>

### *class* highspy._core.HighsOptions

Bases: `pybind11_object`

<a id="highspy._core.HighsOptions.allow_unbounded_or_infeasible"></a>

#### *property* allow_unbounded_or_infeasible

<a id="highspy._core.HighsOptions.allowed_matrix_scale_factor"></a>

#### *property* allowed_matrix_scale_factor

<a id="highspy._core.HighsOptions.blend_multi_objectives"></a>

#### *property* blend_multi_objectives

<a id="highspy._core.HighsOptions.dual_feasibility_tolerance"></a>

#### *property* dual_feasibility_tolerance

<a id="highspy._core.HighsOptions.dual_residual_tolerance"></a>

#### *property* dual_residual_tolerance

<a id="highspy._core.HighsOptions.glpsol_cost_row_location"></a>

#### *property* glpsol_cost_row_location

<a id="highspy._core.HighsOptions.highs_analysis_level"></a>

#### *property* highs_analysis_level

<a id="highspy._core.HighsOptions.highs_debug_level"></a>

#### *property* highs_debug_level

<a id="highspy._core.HighsOptions.infinite_bound"></a>

#### *property* infinite_bound

<a id="highspy._core.HighsOptions.infinite_cost"></a>

#### *property* infinite_cost

<a id="highspy._core.HighsOptions.ipm_iteration_limit"></a>

#### *property* ipm_iteration_limit

<a id="highspy._core.HighsOptions.ipm_optimality_tolerance"></a>

#### *property* ipm_optimality_tolerance

<a id="highspy._core.HighsOptions.ipx_dualize_strategy"></a>

#### *property* ipx_dualize_strategy

<a id="highspy._core.HighsOptions.kkt_tolerance"></a>

#### *property* kkt_tolerance

<a id="highspy._core.HighsOptions.large_matrix_value"></a>

#### *property* large_matrix_value

<a id="highspy._core.HighsOptions.log_dev_level"></a>

#### *property* log_dev_level

<a id="highspy._core.HighsOptions.log_file"></a>

#### *property* log_file

<a id="highspy._core.HighsOptions.log_githash"></a>

#### *property* log_githash

<a id="highspy._core.HighsOptions.log_to_console"></a>

#### *property* log_to_console

<a id="highspy._core.HighsOptions.mip_abs_gap"></a>

#### *property* mip_abs_gap

<a id="highspy._core.HighsOptions.mip_detect_symmetry"></a>

#### *property* mip_detect_symmetry

<a id="highspy._core.HighsOptions.mip_feasibility_tolerance"></a>

#### *property* mip_feasibility_tolerance

<a id="highspy._core.HighsOptions.mip_heuristic_effort"></a>

#### *property* mip_heuristic_effort

<a id="highspy._core.HighsOptions.mip_heuristic_run_feasibility_jump"></a>

#### *property* mip_heuristic_run_feasibility_jump

<a id="highspy._core.HighsOptions.mip_heuristic_run_rens"></a>

#### *property* mip_heuristic_run_rens

<a id="highspy._core.HighsOptions.mip_heuristic_run_rins"></a>

#### *property* mip_heuristic_run_rins

<a id="highspy._core.HighsOptions.mip_heuristic_run_root_reduced_cost"></a>

#### *property* mip_heuristic_run_root_reduced_cost

<a id="highspy._core.HighsOptions.mip_heuristic_run_shifting"></a>

#### *property* mip_heuristic_run_shifting

<a id="highspy._core.HighsOptions.mip_heuristic_run_zi_round"></a>

#### *property* mip_heuristic_run_zi_round

<a id="highspy._core.HighsOptions.mip_lp_age_limit"></a>

#### *property* mip_lp_age_limit

<a id="highspy._core.HighsOptions.mip_max_improving_sols"></a>

#### *property* mip_max_improving_sols

<a id="highspy._core.HighsOptions.mip_max_leaves"></a>

#### *property* mip_max_leaves

<a id="highspy._core.HighsOptions.mip_max_nodes"></a>

#### *property* mip_max_nodes

<a id="highspy._core.HighsOptions.mip_max_stall_nodes"></a>

#### *property* mip_max_stall_nodes

<a id="highspy._core.HighsOptions.mip_min_cliquetable_entries_for_parallelism"></a>

#### *property* mip_min_cliquetable_entries_for_parallelism

<a id="highspy._core.HighsOptions.mip_min_logging_interval"></a>

#### *property* mip_min_logging_interval

<a id="highspy._core.HighsOptions.mip_pool_age_limit"></a>

#### *property* mip_pool_age_limit

<a id="highspy._core.HighsOptions.mip_pool_soft_limit"></a>

#### *property* mip_pool_soft_limit

<a id="highspy._core.HighsOptions.mip_pscost_minreliable"></a>

#### *property* mip_pscost_minreliable

<a id="highspy._core.HighsOptions.mip_rel_gap"></a>

#### *property* mip_rel_gap

<a id="highspy._core.HighsOptions.mip_report_level"></a>

#### *property* mip_report_level

<a id="highspy._core.HighsOptions.objective_bound"></a>

#### *property* objective_bound

<a id="highspy._core.HighsOptions.objective_target"></a>

#### *property* objective_target

<a id="highspy._core.HighsOptions.optimality_tolerance"></a>

#### *property* optimality_tolerance

<a id="highspy._core.HighsOptions.output_flag"></a>

#### *property* output_flag

<a id="highspy._core.HighsOptions.parallel"></a>

#### *property* parallel

<a id="highspy._core.HighsOptions.pdlp_cupdlpc_restart_method"></a>

#### *property* pdlp_cupdlpc_restart_method

<a id="highspy._core.HighsOptions.pdlp_iteration_limit"></a>

#### *property* pdlp_iteration_limit

<a id="highspy._core.HighsOptions.pdlp_optimality_tolerance"></a>

#### *property* pdlp_optimality_tolerance

<a id="highspy._core.HighsOptions.pdlp_scaling_mode"></a>

#### *property* pdlp_scaling_mode

<a id="highspy._core.HighsOptions.presolve"></a>

#### *property* presolve

<a id="highspy._core.HighsOptions.primal_feasibility_tolerance"></a>

#### *property* primal_feasibility_tolerance

<a id="highspy._core.HighsOptions.primal_residual_tolerance"></a>

#### *property* primal_residual_tolerance

<a id="highspy._core.HighsOptions.qp_iteration_limit"></a>

#### *property* qp_iteration_limit

<a id="highspy._core.HighsOptions.qp_nullspace_limit"></a>

#### *property* qp_nullspace_limit

<a id="highspy._core.HighsOptions.qp_regularization_value"></a>

#### *property* qp_regularization_value

<a id="highspy._core.HighsOptions.random_seed"></a>

#### *property* random_seed

<a id="highspy._core.HighsOptions.ranging"></a>

#### *property* ranging

<a id="highspy._core.HighsOptions.read_basis_file"></a>

#### *property* read_basis_file

<a id="highspy._core.HighsOptions.read_solution_file"></a>

#### *property* read_solution_file

<a id="highspy._core.HighsOptions.run_crossover"></a>

#### *property* run_crossover

<a id="highspy._core.HighsOptions.simplex_crash_strategy"></a>

#### *property* simplex_crash_strategy

<a id="highspy._core.HighsOptions.simplex_dual_edge_weight_strategy"></a>

#### *property* simplex_dual_edge_weight_strategy

<a id="highspy._core.HighsOptions.simplex_dualize_strategy"></a>

#### *property* simplex_dualize_strategy

<a id="highspy._core.HighsOptions.simplex_iteration_limit"></a>

#### *property* simplex_iteration_limit

<a id="highspy._core.HighsOptions.simplex_max_concurrency"></a>

#### *property* simplex_max_concurrency

<a id="highspy._core.HighsOptions.simplex_min_concurrency"></a>

#### *property* simplex_min_concurrency

<a id="highspy._core.HighsOptions.simplex_permute_strategy"></a>

#### *property* simplex_permute_strategy

<a id="highspy._core.HighsOptions.simplex_price_strategy"></a>

#### *property* simplex_price_strategy

<a id="highspy._core.HighsOptions.simplex_primal_edge_weight_strategy"></a>

#### *property* simplex_primal_edge_weight_strategy

<a id="highspy._core.HighsOptions.simplex_scale_strategy"></a>

#### *property* simplex_scale_strategy

<a id="highspy._core.HighsOptions.simplex_strategy"></a>

#### *property* simplex_strategy

<a id="highspy._core.HighsOptions.simplex_update_limit"></a>

#### *property* simplex_update_limit

<a id="highspy._core.HighsOptions.small_matrix_value"></a>

#### *property* small_matrix_value

<a id="highspy._core.HighsOptions.solution_file"></a>

#### *property* solution_file

<a id="highspy._core.HighsOptions.solve_relaxation"></a>

#### *property* solve_relaxation

<a id="highspy._core.HighsOptions.solver"></a>

#### *property* solver

<a id="highspy._core.HighsOptions.threads"></a>

#### *property* threads

<a id="highspy._core.HighsOptions.time_limit"></a>

#### *property* time_limit

<a id="highspy._core.HighsOptions.timeless_log"></a>

#### *property* timeless_log

<a id="highspy._core.HighsOptions.user_bound_scale"></a>

#### *property* user_bound_scale

<a id="highspy._core.HighsOptions.user_objective_scale"></a>

#### *property* user_objective_scale

<a id="highspy._core.HighsOptions.write_basis_file"></a>

#### *property* write_basis_file

<a id="highspy._core.HighsOptions.write_model_file"></a>

#### *property* write_model_file

<a id="highspy._core.HighsOptions.write_model_to_file"></a>

#### *property* write_model_to_file

<a id="highspy._core.HighsOptions.write_presolved_model_file"></a>

#### *property* write_presolved_model_file

<a id="highspy._core.HighsOptions.write_solution_style"></a>

#### *property* write_solution_style

<a id="highspy._core.HighsOptions.write_solution_to_file"></a>

#### *property* write_solution_to_file

<a id="highspy._core.HighsPresolveStatus"></a>

### *class* highspy._core.HighsPresolveStatus

Bases: `pybind11_object`

Members:

kNotPresolved

kNotReduced

kInfeasible

kUnboundedOrInfeasible

kReduced

kReducedToEmpty

kTimeout

kNullError

kOptionsError

<a id="highspy._core.HighsPresolveStatus.kInfeasible"></a>

#### kInfeasible *= <HighsPresolveStatus.kInfeasible: 1>*

<a id="highspy._core.HighsPresolveStatus.kNotPresolved"></a>

#### kNotPresolved *= <HighsPresolveStatus.kNotPresolved: -1>*

<a id="highspy._core.HighsPresolveStatus.kNotReduced"></a>

#### kNotReduced *= <HighsPresolveStatus.kNotReduced: 0>*

<a id="highspy._core.HighsPresolveStatus.kNullError"></a>

#### kNullError *= <HighsPresolveStatus.kNullError: 6>*

<a id="highspy._core.HighsPresolveStatus.kOptionsError"></a>

#### kOptionsError *= <HighsPresolveStatus.kOptionsError: 7>*

<a id="highspy._core.HighsPresolveStatus.kReduced"></a>

#### kReduced *= <HighsPresolveStatus.kReduced: 3>*

<a id="highspy._core.HighsPresolveStatus.kReducedToEmpty"></a>

#### kReducedToEmpty *= <HighsPresolveStatus.kReducedToEmpty: 4>*

<a id="highspy._core.HighsPresolveStatus.kTimeout"></a>

#### kTimeout *= <HighsPresolveStatus.kTimeout: 5>*

<a id="highspy._core.HighsPresolveStatus.kUnboundedOrInfeasible"></a>

#### kUnboundedOrInfeasible *= <HighsPresolveStatus.kUnboundedOrInfeasible: 2>*

### HighsPresolveStatus.name -> str

<a id="highspy._core.HighsPresolveStatus.value"></a>

#### *property* value

<a id="highspy._core.HighsRanging"></a>

### *class* highspy._core.HighsRanging

Bases: `pybind11_object`

<a id="highspy._core.HighsRanging.col_bound_dn"></a>

#### *property* col_bound_dn

<a id="highspy._core.HighsRanging.col_bound_up"></a>

#### *property* col_bound_up

<a id="highspy._core.HighsRanging.col_cost_dn"></a>

#### *property* col_cost_dn

<a id="highspy._core.HighsRanging.col_cost_up"></a>

#### *property* col_cost_up

<a id="highspy._core.HighsRanging.row_bound_dn"></a>

#### *property* row_bound_dn

<a id="highspy._core.HighsRanging.row_bound_up"></a>

#### *property* row_bound_up

<a id="highspy._core.HighsRanging.valid"></a>

#### *property* valid

<a id="highspy._core.HighsRangingRecord"></a>

### *class* highspy._core.HighsRangingRecord

Bases: `pybind11_object`

<a id="highspy._core.HighsRangingRecord.in_var_"></a>

#### *property* in_var_

<a id="highspy._core.HighsRangingRecord.objective_"></a>

#### *property* objective_

<a id="highspy._core.HighsRangingRecord.ou_var_"></a>

#### *property* ou_var_

<a id="highspy._core.HighsRangingRecord.value_"></a>

#### *property* value_

<a id="highspy._core.HighsScale"></a>

### *class* highspy._core.HighsScale

Bases: `pybind11_object`

<a id="highspy._core.HighsSolution"></a>

### *class* highspy._core.HighsSolution

Bases: `pybind11_object`

<a id="highspy._core.HighsSolution.col_dual"></a>

#### *property* col_dual

<a id="highspy._core.HighsSolution.col_value"></a>

#### *property* col_value

<a id="highspy._core.HighsSolution.dual_valid"></a>

#### *property* dual_valid

<a id="highspy._core.HighsSolution.row_dual"></a>

#### *property* row_dual

<a id="highspy._core.HighsSolution.row_value"></a>

#### *property* row_value

<a id="highspy._core.HighsSolution.value_valid"></a>

#### *property* value_valid

<a id="highspy._core.HighsSparseMatrix"></a>

### *class* highspy._core.HighsSparseMatrix

Bases: `pybind11_object`

<a id="highspy._core.HighsSparseMatrix.format_"></a>

#### *property* format_

<a id="highspy._core.HighsSparseMatrix.index_"></a>

#### *property* index_

<a id="highspy._core.HighsSparseMatrix.num_col_"></a>

#### *property* num_col_

<a id="highspy._core.HighsSparseMatrix.num_row_"></a>

#### *property* num_row_

<a id="highspy._core.HighsSparseMatrix.p_end_"></a>

#### *property* p_end_

<a id="highspy._core.HighsSparseMatrix.start_"></a>

#### *property* start_

<a id="highspy._core.HighsSparseMatrix.value_"></a>

#### *property* value_

<a id="highspy._core.HighsStatus"></a>

### *class* highspy._core.HighsStatus

Bases: `pybind11_object`

Members:

kError

kOk

kWarning

<a id="highspy._core.HighsStatus.kError"></a>

#### kError *= <HighsStatus.kError: -1>*

<a id="highspy._core.HighsStatus.kOk"></a>

#### kOk *= <HighsStatus.kOk: 0>*

<a id="highspy._core.HighsStatus.kWarning"></a>

#### kWarning *= <HighsStatus.kWarning: 1>*

### HighsStatus.name -> str

<a id="highspy._core.HighsStatus.value"></a>

#### *property* value

<a id="highspy._core.HighsVarType"></a>

### *class* highspy._core.HighsVarType

Bases: `pybind11_object`

Members:

kContinuous

kInteger

kSemiContinuous

kSemiInteger

kImplicitInteger

<a id="highspy._core.HighsVarType.kContinuous"></a>

#### kContinuous *= <HighsVarType.kContinuous: 0>*

<a id="highspy._core.HighsVarType.kImplicitInteger"></a>

#### kImplicitInteger *= <HighsVarType.kImplicitInteger: 4>*

<a id="highspy._core.HighsVarType.kInteger"></a>

#### kInteger *= <HighsVarType.kInteger: 1>*

<a id="highspy._core.HighsVarType.kSemiContinuous"></a>

#### kSemiContinuous *= <HighsVarType.kSemiContinuous: 2>*

<a id="highspy._core.HighsVarType.kSemiInteger"></a>

#### kSemiInteger *= <HighsVarType.kSemiInteger: 3>*

### HighsVarType.name -> str

<a id="highspy._core.HighsVarType.value"></a>

#### *property* value

<a id="highspy._core.IisBoundStatus"></a>

### *class* highspy._core.IisBoundStatus

Bases: `pybind11_object`

Members:

kIisBoundStatusDropped

kIisBoundStatusNull

kIisBoundStatusFree

kIisBoundStatusLower

kIisBoundStatusUpper

kIisBoundStatusBoxed

<a id="highspy._core.IisBoundStatus.kIisBoundStatusBoxed"></a>

#### kIisBoundStatusBoxed *= <IisBoundStatus.kIisBoundStatusBoxed: 4>*

<a id="highspy._core.IisBoundStatus.kIisBoundStatusDropped"></a>

#### kIisBoundStatusDropped *= <IisBoundStatus.kIisBoundStatusDropped: -1>*

<a id="highspy._core.IisBoundStatus.kIisBoundStatusFree"></a>

#### kIisBoundStatusFree *= <IisBoundStatus.kIisBoundStatusFree: 1>*

<a id="highspy._core.IisBoundStatus.kIisBoundStatusLower"></a>

#### kIisBoundStatusLower *= <IisBoundStatus.kIisBoundStatusLower: 2>*

<a id="highspy._core.IisBoundStatus.kIisBoundStatusNull"></a>

#### kIisBoundStatusNull *= <IisBoundStatus.kIisBoundStatusNull: 0>*

<a id="highspy._core.IisBoundStatus.kIisBoundStatusUpper"></a>

#### kIisBoundStatusUpper *= <IisBoundStatus.kIisBoundStatusUpper: 3>*

### IisBoundStatus.name -> str

<a id="highspy._core.IisBoundStatus.value"></a>

#### *property* value

<a id="highspy._core.IisStatus"></a>

### *class* highspy._core.IisStatus

Bases: `pybind11_object`

Members:

kIisStatusNotInConflict

kIisStatusMaybeInConflict

kIisStatusInConflict

<a id="highspy._core.IisStatus.kIisStatusInConflict"></a>

#### kIisStatusInConflict *= <IisStatus.kIisStatusInConflict: 1>*

<a id="highspy._core.IisStatus.kIisStatusMaybeInConflict"></a>

#### kIisStatusMaybeInConflict *= <IisStatus.kIisStatusMaybeInConflict: 0>*

<a id="highspy._core.IisStatus.kIisStatusNotInConflict"></a>

#### kIisStatusNotInConflict *= <IisStatus.kIisStatusNotInConflict: -1>*

### IisStatus.name -> str

<a id="highspy._core.IisStatus.value"></a>

#### *property* value

<a id="highspy._core.IisStrategy"></a>

### *class* highspy._core.IisStrategy

Bases: `pybind11_object`

Members:

kIisStrategyMin

kIisStrategyLight

kIisStrategyFromRay

kIisStrategyFromLp

kIisStrategyIrreducible

kIisStrategyColPriority

kIisStrategyRelaxation

kIisStrategyMax

<a id="highspy._core.IisStrategy.kIisStrategyColPriority"></a>

#### kIisStrategyColPriority *= <IisStrategy.kIisStrategyColPriority: 8>*

<a id="highspy._core.IisStrategy.kIisStrategyFromLp"></a>

#### kIisStrategyFromLp *= <IisStrategy.kIisStrategyFromLp: 2>*

<a id="highspy._core.IisStrategy.kIisStrategyFromRay"></a>

#### kIisStrategyFromRay *= <IisStrategy.kIisStrategyFromRay: 1>*

<a id="highspy._core.IisStrategy.kIisStrategyIrreducible"></a>

#### kIisStrategyIrreducible *= <IisStrategy.kIisStrategyIrreducible: 4>*

<a id="highspy._core.IisStrategy.kIisStrategyLight"></a>

#### kIisStrategyLight *= <IisStrategy.kIisStrategyMin: 0>*

<a id="highspy._core.IisStrategy.kIisStrategyMax"></a>

#### kIisStrategyMax *= <IisStrategy.kIisStrategyMax: 31>*

<a id="highspy._core.IisStrategy.kIisStrategyMin"></a>

#### kIisStrategyMin *= <IisStrategy.kIisStrategyMin: 0>*

<a id="highspy._core.IisStrategy.kIisStrategyRelaxation"></a>

#### kIisStrategyRelaxation *= <IisStrategy.kIisStrategyRelaxation: 16>*

### IisStrategy.name -> str

<a id="highspy._core.IisStrategy.value"></a>

#### *property* value

<a id="highspy._core.MatrixFormat"></a>

### *class* highspy._core.MatrixFormat

Bases: `pybind11_object`

Members:

kColwise

kRowwise

kRowwisePartitioned

<a id="highspy._core.MatrixFormat.kColwise"></a>

#### kColwise *= <MatrixFormat.kColwise: 1>*

<a id="highspy._core.MatrixFormat.kRowwise"></a>

#### kRowwise *= <MatrixFormat.kRowwise: 2>*

<a id="highspy._core.MatrixFormat.kRowwisePartitioned"></a>

#### kRowwisePartitioned *= <MatrixFormat.kRowwisePartitioned: 3>*

### MatrixFormat.name -> str

<a id="highspy._core.MatrixFormat.value"></a>

#### *property* value

<a id="highspy._core.ObjSense"></a>

### *class* highspy._core.ObjSense

Bases: `pybind11_object`

Members:

kMinimize

kMaximize

<a id="highspy._core.ObjSense.kMaximize"></a>

#### kMaximize *= <ObjSense.kMaximize: -1>*

<a id="highspy._core.ObjSense.kMinimize"></a>

#### kMinimize *= <ObjSense.kMinimize: 1>*

### ObjSense.name -> str

<a id="highspy._core.ObjSense.value"></a>

#### *property* value

<a id="highspy._core.SolutionStatus"></a>

### *class* highspy._core.SolutionStatus

Bases: `pybind11_object`

Members:

kSolutionStatusNone

kSolutionStatusInfeasible

kSolutionStatusFeasible

<a id="highspy._core.SolutionStatus.kSolutionStatusFeasible"></a>

#### kSolutionStatusFeasible *= <SolutionStatus.kSolutionStatusFeasible: 2>*

<a id="highspy._core.SolutionStatus.kSolutionStatusInfeasible"></a>

#### kSolutionStatusInfeasible *= <SolutionStatus.kSolutionStatusInfeasible: 1>*

<a id="highspy._core.SolutionStatus.kSolutionStatusNone"></a>

#### kSolutionStatusNone *= <SolutionStatus.kSolutionStatusNone: 0>*

### SolutionStatus.name -> str

<a id="highspy._core.SolutionStatus.value"></a>

#### *property* value

<a id="highspy._core.getExtrasLoadStatus"></a>

### highspy._core.getExtrasLoadStatus() → [str](https://docs.python.org/3/builtins/stdtypes.html#str)

<a id="module-highspy._core.cb"></a>

<a id="callback-bindings"></a>

## Callback Bindings

Callback interface submodule

<a id="highspy._core.cb.HighsCallbackInput"></a>

### *class* highspy._core.cb.HighsCallbackInput

Bases: `pybind11_object`

<a id="highspy._core.cb.HighsCallbackInput.repairSolution"></a>

#### repairSolution(self: [highspy._core.cb.HighsCallbackInput](#highspy._core.cb.HighsCallbackInput)) → [highspy._core.HighsStatus](#highspy._core.HighsStatus)

<a id="highspy._core.cb.HighsCallbackInput.setSolution"></a>

#### setSolution(\*args, \*\*kwargs)

Overloaded function.

1. setSolution(self: highspy._core.cb.HighsCallbackInput, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus
2. setSolution(self: highspy._core.cb.HighsCallbackInput, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> highspy._core.HighsStatus

<a id="highspy._core.cb.HighsCallbackInput.user_has_solution"></a>

#### *property* user_has_solution

<a id="highspy._core.cb.HighsCallbackInput.user_interrupt"></a>

#### *property* user_interrupt

<a id="highspy._core.cb.HighsCallbackInput.user_solution"></a>

#### *property* user_solution

<a id="highspy._core.cb.HighsCallbackOutput"></a>

### *class* highspy._core.cb.HighsCallbackOutput

Bases: `pybind11_object`

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_index"></a>

#### *property* cutpool_index

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_lower"></a>

#### *property* cutpool_lower

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_num_col"></a>

#### *property* cutpool_num_col

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_num_cut"></a>

#### *property* cutpool_num_cut

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_start"></a>

#### *property* cutpool_start

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_upper"></a>

#### *property* cutpool_upper

<a id="highspy._core.cb.HighsCallbackOutput.cutpool_value"></a>

#### *property* cutpool_value

<a id="highspy._core.cb.HighsCallbackOutput.ipm_iteration_count"></a>

#### *property* ipm_iteration_count

<a id="highspy._core.cb.HighsCallbackOutput.log_type"></a>

#### *property* log_type

<a id="highspy._core.cb.HighsCallbackOutput.mip_dual_bound"></a>

#### *property* mip_dual_bound

<a id="highspy._core.cb.HighsCallbackOutput.mip_gap"></a>

#### *property* mip_gap

<a id="highspy._core.cb.HighsCallbackOutput.mip_node_count"></a>

#### *property* mip_node_count

<a id="highspy._core.cb.HighsCallbackOutput.mip_primal_bound"></a>

#### *property* mip_primal_bound

<a id="highspy._core.cb.HighsCallbackOutput.mip_solution"></a>

#### *property* mip_solution

<a id="highspy._core.cb.HighsCallbackOutput.objective_function_value"></a>

#### *property* objective_function_value

<a id="highspy._core.cb.HighsCallbackOutput.pdlp_iteration_count"></a>

#### *property* pdlp_iteration_count

<a id="highspy._core.cb.HighsCallbackOutput.running_time"></a>

#### *property* running_time

<a id="highspy._core.cb.HighsCallbackOutput.simplex_iteration_count"></a>

#### *property* simplex_iteration_count

<a id="highspy._core.cb.HighsCallbackType"></a>

### *class* highspy._core.cb.HighsCallbackType

Bases: `pybind11_object`

Members:

kCallbackMin

kCallbackLogging

kCallbackSimplexInterrupt

kCallbackIpmInterrupt

kCallbackMipSolution

kCallbackMipImprovingSolution

kCallbackMipLogging

kCallbackMipInterrupt

kCallbackMipGetCutPool

kCallbackMipDefineLazyConstraints

kCallbackMipUserSolution

kHighsCallbackQpFirstFeasiblePoint

kHighsCallbackQpInterrupt

kCallbackMax

kNumCallbackType

<a id="highspy._core.cb.HighsCallbackType.kCallbackIpmInterrupt"></a>

#### kCallbackIpmInterrupt *= <HighsCallbackType.kCallbackIpmInterrupt: 2>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackLogging"></a>

#### kCallbackLogging *= <HighsCallbackType.kCallbackMin: 0>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMax"></a>

#### kCallbackMax *= <HighsCallbackType.kHighsCallbackQpInterrupt: 11>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMin"></a>

#### kCallbackMin *= <HighsCallbackType.kCallbackMin: 0>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipDefineLazyConstraints"></a>

#### kCallbackMipDefineLazyConstraints *= <HighsCallbackType.kCallbackMipDefineLazyConstraints: 8>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipGetCutPool"></a>

#### kCallbackMipGetCutPool *= <HighsCallbackType.kCallbackMipGetCutPool: 7>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipImprovingSolution"></a>

#### kCallbackMipImprovingSolution *= <HighsCallbackType.kCallbackMipImprovingSolution: 4>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipInterrupt"></a>

#### kCallbackMipInterrupt *= <HighsCallbackType.kCallbackMipInterrupt: 6>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipLogging"></a>

#### kCallbackMipLogging *= <HighsCallbackType.kCallbackMipLogging: 5>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipSolution"></a>

#### kCallbackMipSolution *= <HighsCallbackType.kCallbackMipSolution: 3>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackMipUserSolution"></a>

#### kCallbackMipUserSolution *= <HighsCallbackType.kCallbackMipUserSolution: 9>*

<a id="highspy._core.cb.HighsCallbackType.kCallbackSimplexInterrupt"></a>

#### kCallbackSimplexInterrupt *= <HighsCallbackType.kCallbackSimplexInterrupt: 1>*

<a id="highspy._core.cb.HighsCallbackType.kHighsCallbackQpFirstFeasiblePoint"></a>

#### kHighsCallbackQpFirstFeasiblePoint *= <HighsCallbackType.kHighsCallbackQpFirstFeasiblePoint: 10>*

<a id="highspy._core.cb.HighsCallbackType.kHighsCallbackQpInterrupt"></a>

#### kHighsCallbackQpInterrupt *= <HighsCallbackType.kHighsCallbackQpInterrupt: 11>*

<a id="highspy._core.cb.HighsCallbackType.kNumCallbackType"></a>

#### kNumCallbackType *= <HighsCallbackType.kNumCallbackType: 12>*

### HighsCallbackType.name -> str

<a id="highspy._core.cb.HighsCallbackType.value"></a>

#### *property* value

<a id="module-highspy._core.simplex_constants"></a>

<a id="simplex-constants"></a>

## Simplex Constants

Submodule for simplex constants

<a id="highspy._core.simplex_constants.EdgeWeightMode"></a>

### *class* highspy._core.simplex_constants.EdgeWeightMode

Bases: `pybind11_object`

Members:

kDantzig

kDevex

kSteepestEdge

kCount

<a id="highspy._core.simplex_constants.EdgeWeightMode.kCount"></a>

#### kCount *= <EdgeWeightMode.kCount: 3>*

<a id="highspy._core.simplex_constants.EdgeWeightMode.kDantzig"></a>

#### kDantzig *= <EdgeWeightMode.kDantzig: 0>*

<a id="highspy._core.simplex_constants.EdgeWeightMode.kDevex"></a>

#### kDevex *= <EdgeWeightMode.kDevex: 1>*

<a id="highspy._core.simplex_constants.EdgeWeightMode.kSteepestEdge"></a>

#### kSteepestEdge *= <EdgeWeightMode.kSteepestEdge: 2>*

### EdgeWeightMode.name -> str

<a id="highspy._core.simplex_constants.EdgeWeightMode.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexEdgeWeightStrategy

Bases: `pybind11_object`

Members:

kSimplexEdgeWeightStrategyMin

kSimplexEdgeWeightStrategyChoose

kSimplexEdgeWeightStrategyDantzig

kSimplexEdgeWeightStrategyDevex

kSimplexEdgeWeightStrategySteepestEdge

kSimplexEdgeWeightStrategyMax

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyChoose"></a>

#### kSimplexEdgeWeightStrategyChoose *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyMin: -1>*

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyDantzig"></a>

#### kSimplexEdgeWeightStrategyDantzig *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyDantzig: 0>*

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyDevex"></a>

#### kSimplexEdgeWeightStrategyDevex *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyDevex: 1>*

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyMax"></a>

#### kSimplexEdgeWeightStrategyMax *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategySteepestEdge: 2>*

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyMin"></a>

#### kSimplexEdgeWeightStrategyMin *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategyMin: -1>*

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategySteepestEdge"></a>

#### kSimplexEdgeWeightStrategySteepestEdge *= <SimplexEdgeWeightStrategy.kSimplexEdgeWeightStrategySteepestEdge: 2>*

### SimplexEdgeWeightStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexEdgeWeightStrategy.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexNlaOperation"></a>

### *class* highspy._core.simplex_constants.SimplexNlaOperation

Bases: `pybind11_object`

Members:

kSimplexNlaNull

kSimplexNlaBtranFull

kSimplexNlaPriceFull

kSimplexNlaBtranBasicFeasibilityChange

kSimplexNlaBtranEp

kSimplexNlaPriceAp

kSimplexNlaFtran

kSimplexNlaFtranBfrt

kSimplexNlaFtranDse

kSimplexNlaBtranPse

kNumSimplexNlaOperation

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kNumSimplexNlaOperation"></a>

#### kNumSimplexNlaOperation *= <SimplexNlaOperation.kNumSimplexNlaOperation: 10>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaBtranBasicFeasibilityChange"></a>

#### kSimplexNlaBtranBasicFeasibilityChange *= <SimplexNlaOperation.kSimplexNlaBtranBasicFeasibilityChange: 2>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaBtranEp"></a>

#### kSimplexNlaBtranEp *= <SimplexNlaOperation.kSimplexNlaBtranEp: 4>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaBtranFull"></a>

#### kSimplexNlaBtranFull *= <SimplexNlaOperation.kSimplexNlaBtranFull: 0>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaBtranPse"></a>

#### kSimplexNlaBtranPse *= <SimplexNlaOperation.kSimplexNlaBtranPse: 9>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaFtran"></a>

#### kSimplexNlaFtran *= <SimplexNlaOperation.kSimplexNlaFtran: 6>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaFtranBfrt"></a>

#### kSimplexNlaFtranBfrt *= <SimplexNlaOperation.kSimplexNlaFtranBfrt: 7>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaFtranDse"></a>

#### kSimplexNlaFtranDse *= <SimplexNlaOperation.kSimplexNlaFtranDse: 8>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaNull"></a>

#### kSimplexNlaNull *= <SimplexNlaOperation.kSimplexNlaNull: -1>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaPriceAp"></a>

#### kSimplexNlaPriceAp *= <SimplexNlaOperation.kSimplexNlaPriceAp: 5>*

<a id="highspy._core.simplex_constants.SimplexNlaOperation.kSimplexNlaPriceFull"></a>

#### kSimplexNlaPriceFull *= <SimplexNlaOperation.kSimplexNlaPriceFull: 1>*

### SimplexNlaOperation.name -> str

<a id="highspy._core.simplex_constants.SimplexNlaOperation.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy

Bases: `pybind11_object`

Members:

kSimplexInfeasibilityProofRefinementMin

kSimplexInfeasibilityProofRefinementNo

kSimplexInfeasibilityProofRefinementUnscaledLp

kSimplexInfeasibilityProofRefinementAlsoScaledLp

kSimplexInfeasibilityProofRefinementMax

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementAlsoScaledLp"></a>

#### kSimplexInfeasibilityProofRefinementAlsoScaledLp *= <SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementAlsoScaledLp: 2>*

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementMax"></a>

#### kSimplexInfeasibilityProofRefinementMax *= <SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementAlsoScaledLp: 2>*

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementMin"></a>

#### kSimplexInfeasibilityProofRefinementMin *= <SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementMin: 0>*

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementNo"></a>

#### kSimplexInfeasibilityProofRefinementNo *= <SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementMin: 0>*

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementUnscaledLp"></a>

#### kSimplexInfeasibilityProofRefinementUnscaledLp *= <SimplexPivotalRowRefinementStrategy.kSimplexInfeasibilityProofRefinementUnscaledLp: 1>*

### SimplexPivotalRowRefinementStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexPivotalRowRefinementStrategy.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexPriceStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexPriceStrategy

Bases: `pybind11_object`

Members:

kSimplexPriceStrategyMin

kSimplexPriceStrategyCol

kSimplexPriceStrategyRow

kSimplexPriceStrategyRowSwitch

kSimplexPriceStrategyRowSwitchColSwitch

kSimplexPriceStrategyMax

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyCol"></a>

#### kSimplexPriceStrategyCol *= <SimplexPriceStrategy.kSimplexPriceStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyMax"></a>

#### kSimplexPriceStrategyMax *= <SimplexPriceStrategy.kSimplexPriceStrategyRowSwitchColSwitch: 3>*

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyMin"></a>

#### kSimplexPriceStrategyMin *= <SimplexPriceStrategy.kSimplexPriceStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyRow"></a>

#### kSimplexPriceStrategyRow *= <SimplexPriceStrategy.kSimplexPriceStrategyRow: 1>*

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyRowSwitch"></a>

#### kSimplexPriceStrategyRowSwitch *= <SimplexPriceStrategy.kSimplexPriceStrategyRowSwitch: 2>*

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.kSimplexPriceStrategyRowSwitchColSwitch"></a>

#### kSimplexPriceStrategyRowSwitchColSwitch *= <SimplexPriceStrategy.kSimplexPriceStrategyRowSwitchColSwitch: 3>*

### SimplexPriceStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexPriceStrategy.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy

Bases: `pybind11_object`

Members:

kSimplexPrimalCorrectionStrategyNone

kSimplexPrimalCorrectionStrategyInRebuild

kSimplexPrimalCorrectionStrategyAlways

<a id="highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyAlways"></a>

#### kSimplexPrimalCorrectionStrategyAlways *= <SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyAlways: 2>*

<a id="highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyInRebuild"></a>

#### kSimplexPrimalCorrectionStrategyInRebuild *= <SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyInRebuild: 1>*

<a id="highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyNone"></a>

#### kSimplexPrimalCorrectionStrategyNone *= <SimplexPrimalCorrectionStrategy.kSimplexPrimalCorrectionStrategyNone: 0>*

### SimplexPrimalCorrectionStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexPrimalCorrectionStrategy.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexSolvePhase"></a>

### *class* highspy._core.simplex_constants.SimplexSolvePhase

Bases: `pybind11_object`

Members:

kSolvePhaseMin

kSolvePhaseError

kSolvePhaseExit

kSolvePhaseUnknown

kSolvePhaseOptimal

kSolvePhase1

kSolvePhase2

kSolvePhasePrimalInfeasibleCleanup

kSolvePhaseOptimalCleanup

kSolvePhaseTabooBasis

kSolvePhaseMax

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhase1"></a>

#### kSolvePhase1 *= <SimplexSolvePhase.kSolvePhase1: 1>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhase2"></a>

#### kSolvePhase2 *= <SimplexSolvePhase.kSolvePhase2: 2>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseError"></a>

#### kSolvePhaseError *= <SimplexSolvePhase.kSolvePhaseMin: -3>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseExit"></a>

#### kSolvePhaseExit *= <SimplexSolvePhase.kSolvePhaseExit: -2>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseMax"></a>

#### kSolvePhaseMax *= <SimplexSolvePhase.kSolvePhaseTabooBasis: 5>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseMin"></a>

#### kSolvePhaseMin *= <SimplexSolvePhase.kSolvePhaseMin: -3>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseOptimal"></a>

#### kSolvePhaseOptimal *= <SimplexSolvePhase.kSolvePhaseOptimal: 0>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseOptimalCleanup"></a>

#### kSolvePhaseOptimalCleanup *= <SimplexSolvePhase.kSolvePhaseOptimalCleanup: 4>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhasePrimalInfeasibleCleanup"></a>

#### kSolvePhasePrimalInfeasibleCleanup *= <SimplexSolvePhase.kSolvePhasePrimalInfeasibleCleanup: 3>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseTabooBasis"></a>

#### kSolvePhaseTabooBasis *= <SimplexSolvePhase.kSolvePhaseTabooBasis: 5>*

<a id="highspy._core.simplex_constants.SimplexSolvePhase.kSolvePhaseUnknown"></a>

#### kSolvePhaseUnknown *= <SimplexSolvePhase.kSolvePhaseUnknown: -1>*

### SimplexSolvePhase.name -> str

<a id="highspy._core.simplex_constants.SimplexSolvePhase.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexStrategy

Bases: `pybind11_object`

Members:

kSimplexStrategyMin

kSimplexStrategyChoose

kSimplexStrategyDual

kSimplexStrategyDualPlain

kSimplexStrategyDualTasks

kSimplexStrategyDualMulti

kSimplexStrategyPrimal

kSimplexStrategyMax

kSimplexStrategyNum

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyChoose"></a>

#### kSimplexStrategyChoose *= <SimplexStrategy.kSimplexStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyDual"></a>

#### kSimplexStrategyDual *= <SimplexStrategy.kSimplexStrategyDual: 1>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyDualMulti"></a>

#### kSimplexStrategyDualMulti *= <SimplexStrategy.kSimplexStrategyDualMulti: 3>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyDualPlain"></a>

#### kSimplexStrategyDualPlain *= <SimplexStrategy.kSimplexStrategyDual: 1>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyDualTasks"></a>

#### kSimplexStrategyDualTasks *= <SimplexStrategy.kSimplexStrategyDualTasks: 2>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyMax"></a>

#### kSimplexStrategyMax *= <SimplexStrategy.kSimplexStrategyPrimal: 4>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyMin"></a>

#### kSimplexStrategyMin *= <SimplexStrategy.kSimplexStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyNum"></a>

#### kSimplexStrategyNum *= <SimplexStrategy.kSimplexStrategyNum: 5>*

<a id="highspy._core.simplex_constants.SimplexStrategy.kSimplexStrategyPrimal"></a>

#### kSimplexStrategyPrimal *= <SimplexStrategy.kSimplexStrategyPrimal: 4>*

### SimplexStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexStrategy.value"></a>

#### *property* value

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy"></a>

### *class* highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy

Bases: `pybind11_object`

Members:

kSimplexUnscaledSolutionStrategyMin

kSimplexUnscaledSolutionStrategyNone

kSimplexUnscaledSolutionStrategyRefine

kSimplexUnscaledSolutionStrategyDirect

kSimplexUnscaledSolutionStrategyMax

kSimplexUnscaledSolutionStrategyNum

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyDirect"></a>

#### kSimplexUnscaledSolutionStrategyDirect *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyDirect: 2>*

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyMax"></a>

#### kSimplexUnscaledSolutionStrategyMax *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyDirect: 2>*

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyMin"></a>

#### kSimplexUnscaledSolutionStrategyMin *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyNone"></a>

#### kSimplexUnscaledSolutionStrategyNone *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyMin: 0>*

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyNum"></a>

#### kSimplexUnscaledSolutionStrategyNum *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyNum: 3>*

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyRefine"></a>

#### kSimplexUnscaledSolutionStrategyRefine *= <SimplexUnscaledSolutionStrategy.kSimplexUnscaledSolutionStrategyRefine: 1>*

### SimplexUnscaledSolutionStrategy.name -> str

<a id="highspy._core.simplex_constants.SimplexUnscaledSolutionStrategy.value"></a>

#### *property* value
