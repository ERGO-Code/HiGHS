# Basic features

This guide describes the basic features of HiGHS that are available
when it is called from [`Python`](@ref python-getting-started),
[`C++`](@ref cpp-getting-started), [`C`](@ref c-api) and
[`Fortran`](@ref fortran-api). Although references to methods link to
`Python` examples, the method names and functionality for other
interfaces are as close as possible.

This basic guide will be sufficient for most users, and covers all the
methods in the `Python` interface `highspy`. A guide to more advanced
methods available via other interfaces is [available](@ref
Advanced-features).

The minimal use of HiGHS has the following three stages.

* [Define a model](@ref Defining-a-model)
* [Solve the model](@ref Solving-the-model)
* [Extract results](@ref Extracting-results)

Although its default actions will be sufficient for most users, HiGHS can be controlled by setting [Option values](@ref Option-values).

_Intro to other basic features_

## HiGHS data structures

There are several specialist data structures that can be used to
interact with HiGHS when using [`C++`](@ref cpp-getting-started) and
`highspy`. These are defined in the sections on [enums](@ref Enums)
and [classes](@ref classes-overview), and are referred to below. However, the
use of classes is not necessary for the basic use of `highspy`. As
with the `C` and `Fortran` interfaces, there are equivalent methods
that use simple scalars and vectors of data.

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
 
HiGHS can read compressed files that end in the `.gz` extension.

### Building a model

The model in HiGHS can be built using a sequence of calls to add
variables and constraints. This is most easily done one-by-one using
the methods [`addCol` and `addRow`](@ref Build-a-model). Alterntively,
[`addVar` and `addRow`](@ref Build-a-model) can be used, with
[`changeColCost`](@ref Build-a-model) used to define each objective
coefficient.

Addition of multiple variables and constraints can be achieved using
[`addVars` and `addRows`](@ref Build-a-model), with
[`changeColsCost`](@ref Build-a-model) used to define objective
coefficients. Note that defining the model in this way requires
vectors of data and the specification of constraint coefficients as
compressed
[row-wise](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format))
or
[column-wise](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS))
matrices.

### Passing a model

If the entire definition of a model is known, then is can be passed to
HiGHS via individual data arrays using the method [`passModel`](@ref
Pass-a-model). In languages where HiGHS data structures can be used,
an instance of the [`HighsLp`](@ref HighsLp) class can be passed.

## Solving the model

The incumbent model in HiGHS is solved by a call to the method [run](@ref Solve-the-model).
By default, HiGHS minimizes the model's objective function. Where possible,
HiGHS will hot start the solver using solution information obtained on previous
runs, or supplied by the user.

## Extracting results

After solving a model, its status is the value returnedby the method
[getModelStatus](@ref Extract-results). This value is of type [HighsModelStatus](@ref).
Scalar information about a solved model is obtained using the method [getInfo](@ref Extract-results).
The solution and (any) basis are returned by the methods [getSolution](@ref Extract-results)
and [getBasis](@ref Extract-results) respectively. HiGHS can also be used to
write the solution to a file using the method [writeSolution](@ref Report-results).

## Option values

The option values that control HiGHS are of type `string`, `bool`, `int` and
`double`. Options are referred to by a `string` identical to the name of their
identifier.

A full specification of the options is given in the [List of options](@ref).

An option value is changed by passing its name and value to the method [setOptionValue](@ref example-py-option-values).
The current value of an option is obtained by passing its name to the method
[getOptionValue](@ref example-py-option-values).

## Model and solution management

HiGHS has comprehensive tools for defining and extracting models. This can be
done either to/from MPS or (CPLEX) format LP files, or via method calls. HiGHS
also has methods that permit the incumbent model to be modified. Solutions can
be supplied and extracted using either files or method calls.

### Extracting model data

The numbers of column, rows and nonzeros in the model are returned by the
methods [getNumCols](@ref Get-model-data), [getNumRows](@ref Get-model-data),
and [getNumEntries](@ref Get-model-data) respectively.

Model data can be extracted for a single column or row by specifying the index
of the column or row and calling the methods [getCol](@ref Get-model-data) and
[getRow](@ref Get-model-data).

As well as returning the value of the cost and bounds, these methods also return
the number of nonzeros in the corresponding column or row of the constraint
matrix. The indices and values of the nonzeros can be obtained using the methods
[getColEntries](@ref Get-model-data) and [getRowEntries](@ref Get-model-data).

For multiple columns and rows defined by a set of indices, the corresponding
data can be extracted using the methods [getCols](@ref Get-model-data),
[getRows](@ref Get-model-data), [getColsEntries](@ref Get-model-data) and
[getRowsEntries](@ref Get-model-data).

Specific matrix coefficients obtained using the method [getCoeff](@ref Get-model-data).

### Modifying model data

The most immediate model modification is to change the sense of the objective.
By default, HiGHS minimizes the model's objective function. The objective sense
can be set to minimize (maximize) using [changeObjectiveSense](@ref Modify-model-data).

Model data for can be changed for one column or row by specifying the index of
the column or row, together with the new scalar value for the cost or bounds,
the specific methods being [changeColCost](@ref Modify-model-data),
[changeColBounds](@ref Modify-model-data). The corresponding method for a row is
[changeRowBounds](@ref Modify-model-data). Changes for multiple columns or rows
are defined by supplying a list of indices, together with arrays of new values,
using the methods [changeColsCost](@ref Modify-model-data),
[changeColsBounds](@ref Modify-model-data). The corresponding method for a row
is [changeRowsBounds](@ref Modify-model-data). An individual matrix coefficient
is changed by passing its row index, column index and new value to
[changeCoeff](@ref Modify-model-data).

### Other operations

To run HiGHS from a user-defined solution or basis, this is passed to HiGHS
using the methods [setSolution](@ref Set-solution) or [setBasis](@ref Set-basis).
