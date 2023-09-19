# [Further features](@id guide-further)

## Model and solution management

HiGHS has comprehensive tools for defining and extracting models. This can be
done either to/from MPS or (CPLEX) format LP files, or via method calls. HiGHS
also has methods that permit the incumbent model to be modified. Solutions can
be supplied and extracted using either files or method calls.

### Extracting model data

The numbers of column, rows and nonzeros in the model are returned by the
methods [`getNumCols`](@ref Get-model-data), [`getNumRows`](@ref Get-model-data),
and [`getNumEntries`](@ref Get-model-data) respectively.

Model data can be extracted for a single column or row by specifying the index
of the column or row and calling the methods [`getCol`](@ref Get-model-data) and
[`getRow`](@ref Get-model-data).

As well as returning the value of the cost and bounds, these methods also return
the number of nonzeros in the corresponding column or row of the constraint
matrix. The indices and values of the nonzeros can be obtained using the methods
[`getColEntries`](@ref Get-model-data) and [`getRowEntries`](@ref Get-model-data).

For multiple columns and rows defined by a set of indices, the corresponding
data can be extracted using the methods [`getCols`](@ref Get-model-data),
[`getRows`](@ref Get-model-data), [`getColsEntries`](@ref Get-model-data) and
[`getRowsEntries`](@ref Get-model-data).

Specific matrix coefficients obtained using the method [`getCoeff`](@ref Get-model-data).

## Modifying model data

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

### [Hot start](@id hot-start)

It may be possible for HiGHS to start solving a model using data
obtained by solving a related model, or supplied by a user. Whether
this is possible depends on the the class of model being solved, the
solver to be used, and the modifications (if any) that have been to
the incumbent model since it was last solved.

#### LP

To run HiGHS from a user-defined solution or basis, this is passed to HiGHS
using the methods [setSolution](@ref Set-solution) or [setBasis](@ref Set-basis). The basis passed to HiGHS need not be complete

* There can be more basic variables then the number of rows in the
  model. HiGHS will identify a set of basic variables of the correct
  dimension by making some basic variables nonbasic

* There can be fewer basic variables then the number of rows in the
  model.  HiGHS will identify a set of basic variables of the correct
  dimension by adding basic variables corresponding to slacks.

* For nonbasic variables, it is unnecessary to specify whether they
  are at their lower or upper bound unless they are "boxed" variables.

#### MIP

If a (partial) feasible assignment of the integer variables is known,
this can be passed to HiGHS via [setSolution](@ref Set-solution). If
integer variables are set to integer values, HiGHS will solve the LP
with these integer variables fixed (or MIP if the assignment of the
integer variables is not complete). If a feasible solution is
obtained, it will be used to provide the MIP solver with an initial
primal bound when it run to solve for all integer variables.


