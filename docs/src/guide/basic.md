# [Basic features](@id guide-basic)

The minimal use of HiGHS has the following three stages.

* [Define a model](@ref Defining-a-model)
* [Solve the model](@ref Solving-the-model)
* [Extract results](@ref Extracting-results)

Although its default actions will be sufficient for most users, HiGHS
can be controlled by setting [Option values](@ref Option-values).

_Intro to other basic features_

### HiGHS data structures

There are several specialist data structures that can be used to
interact with HiGHS when using [`C++`](@ref cpp-getting-started) and
`highspy`. These are defined in the sections on [enums](@ref structures-enums)
and [classes](@ref classes-overview), and are referred to below.

#### [Enums](@id guide-basic-enums)

Enums are scalar identifier types that can take only a limited range of values.

#### [Classes](@id guide-basic-classes) The advantage of using the
native `C++` classes in HiGHS is that many fewer parameters are needed
when passing data to and from HiGHS. The binding of the data members
of these classes to `highspy` structures allows them to be used when
calling HiGHS from Python, although they are not necessary for the
basic use of `highspy`. As with the `C` and `Fortran` interfaces,
there are equivalent methods that use simple scalars and vectors of
data.

## Defining a model

HiGHS has comprehensive tools for defining models. This can be done by
either reading a model using a data file created elsewhere, or by
passing model data to HiGHS. Once a model has been defined in HiGHS,
it is referred to as the `incumbent model`.

### Reading a model from a file

The simplest way to define a model in HiGHS is to read it from a file using
the method [`readModel`](@ref Read-a-model). HiGHS infers the file type by the extension. Supported extensions are:

 * `.mps`: for an MPS file
 * `.lp`: for a CPLEX LP file
 
HiGHS can read compressed files that end in the `.gz` extension, but
not (yet) files that end in the `.zip` extension.

### Building a model

The model in HiGHS can be built using a sequence of calls to add
variables and constraints. This is most easily done one-by-one using
the methods [`addCol` and `addRow`](@ref
Build-a-model). Alternatively, calls to [`addVar`](@ref Build-a-model)
can be used to add variables, with calls to [`changeColCost`](@ref
Build-a-model) used to define each objective coefficient.

Addition of multiple variables and constraints can be achieved using
[`addCols` and `addRows`](@ref Build-a-model). Alternatively,
[`addVars`](@ref Build-a-model) can be used to add variables, with
[`changeColsCost`](@ref Build-a-model) used to define objective
coefficients. Note that defining multiple variables and constraints requires
vectors of data and the specification of constraint coefficients as
compressed
[row-wise](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format))
or
[column-wise](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS))
matrices.

### Passing a model

If the entire definition of a model is known, then it can be passed to
HiGHS via individual data arrays using the method [`passModel`](@ref
Pass-a-model). In languages where HiGHS data structures can be used,
an instance of the [`HighsLp`](@ref HighsLp) class can be populated
with data and then passed.

## Solving the model

The incumbent model in HiGHS is solved by a call to the method
[`run`](@ref Solve-the-model). By default, HiGHS minimizes the model's
objective function, although this can be [changed](@ref
Modifying-model-data). Where possible, HiGHS will reduce the solution
time by using data obtained on previous runs, or supplied by the
user. More information on this process of hot starting solvers is
[available](@ref hot-start).

## Extracting results

After solving a model, it is important to know whether it has been
solved to optimality, shown to be infeasible or unbounded, or why the
solver may have terminated early. This model status is given by the
value returned by the method [`getModelStatus`](@ref
Extract-results). This value is of type [`HighsModelStatus`](@ref HighsModelStatus).
Scalar information about a solved model is obtained using the method
[`getInfo`](@ref Extract-results).  The solution and (any) basis are
returned by the methods [`getSolution`](@ref Extract-results) and
[`getBasis`](@ref Extract-results) respectively. HiGHS can also be used
to write the solution to a file using the method [`writeSolution`](@ref
Report-results).

## Option values

The option values that control HiGHS are of type `string`, `bool`,
`int` and `double`. Options are referred to by a `string` identical to
the name of their identifier. A full specification of the options is
given in the [list of options](@ref option-definitions). An option
value is changed by passing its name and value to the method
[`setOptionValue`](@ref example-py-option-values).  The current value
of an option is obtained by passing its name to the method
[`getOptionValue`](@ref example-py-option-values).

