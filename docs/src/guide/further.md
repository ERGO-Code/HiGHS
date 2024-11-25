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

## [Hot start](@id hot-start)

It may be possible for HiGHS to start solving a model using data
obtained by solving a related model, or supplied by a user. Whether
this is possible depends on the the class of model being solved, the
solver to be used, and the modifications (if any) that have been to
the incumbent model since it was last solved.

### LP

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

### MIP

If a (partial) feasible assignment of the integer variables is known,
this can be passed to HiGHS via [setSolution](@ref Set-solution). If
integer variables are set to integer values, HiGHS will solve the LP
with these integer variables fixed (or MIP if the assignment of the
integer variables is not complete). If a feasible solution is
obtained, it will be used to provide the MIP solver with an initial
primal bound when it run to solve for all integer variables.

## [Presolve](@id guide-presolve)

HiGHS has a sophisticated presolve procedure for LPs and MIPs that
aims to reduce the dimension of the model that must be solved. In most
cases, the time saved by solving the reduced model is very much
greater than the time taken to perform presolve. Once he presolved
model is solved, a postsolve procedure (of minimal computational cost)
deduces the optimal solution to the original model. Hence presolve is
performed by default. The only exception occus when there is a valid
basis for an LP and the simplex solver is used. In this case the
original LP is solved, starting from this basis. In cases where the
use of presolve is found not to be advantageous, its use can be
switched off by setting the [presolve](@ref option-presolve) option to
"off".

HiGHS has a method [presolve](@ref Presolve/postsolve) that performs presolve on
the incumbent model, allowing the presolved model to be extracted or
written to a file. This is intended for users who have their own
solution technique that they wish to test using models presolved by
HiGHS. Note that this does not affect how the incumbent model is
solved. There are two corresponding [postsolve](@ref Presolve/postsolve)
methods, according to whether there are just solution values, or also
a basis.

## [Multi-objective optimization](@id guide-multi-objective-optimization)

Users can specify multiple linear objectives with respect to which
HiGHS will optimize by either blending them, or by performing
lexicographic optimization according to the truth of the
[blend\_multi\_objectives](@ref blend_multi_objectives) option. Each
linear objective is represented by the following data, held in the
[HighsLinearObjective](@ref HighsLinearObjective) structure

- weight: Scalar of type double - The weight of this objective when blending 
- offset: Scalar of type double - The offset of this objective
- coefficients: Vector of type double - The coefficients of this objective
- abs\_tolerance: Scalar of type double - The absolute tolerance on this objective when performing lexicographic optimization 
- rel\_tolerance: Scalar of type double - The relative tolerance on this objective when performing lexicographic optimization 
- priority: Scalar of type HighsInt - The priority of this objective when performing lexicographic optimization

### Methods

Multi-objective optimization in HiGHS is defined by the following methods

- [addLinearObjective](@ref Multi-objective-optimization] - Add a single `HighsLinearObjective` instance to any already stored in HiGHS
- [clearLinearObjectives](@ref Multi-objective-optimization] - Clears any linear objectives stored in HiGHS

When there is at least one `HighsLinearObjective` instance in HiGHS,
the `col_cost_` data in the incumbent model is ignored.

### Blending multiple linear objectives

When [blend\_multi\_objectives](@ref blend_multi_objectives) is `true`,
as it is by default, any `HighsLinearObjective` instances will be
combined according to the `weight` values, and the resulting objective
will be minimized. Hence, any objectives that should be maximized
within the combination must have a negative `weight` value.

### Lexicographic optimization of multiple linear objectives

When [blend\_multi\_objectives](@ref blend_multi_objectives) is `false`,
HiGHS will optimize lexicographically with respect to any
`HighsLinearObjective` instances. This is carried out as follows, according to the
`priority` values in `HighsLinearObjective` instances. Note that _all
priority values must be distinct_.

* Minimize/maximize with respect to the linear objective of highest priority value, according to whether its `weight` is positive/negative

* Add a constraint to the model so that the value of the linear objective of highest priority satsifies a bound given by the values of `abs_tolerance` and/or `rel_tolerance`.
    + If the objective was minimized to a value ``f^*\ge0``, then the constraint ensures that the this objective value is no greater than
```math
\min(f^*+`abs\_tolerance`,~f^*\times[1+`rel\_tolerance`]).
```
    + If the objective was minimized to a value ``f^*\lt0``, then the constraint ensures that the this objective value is no greater than
```math
\min(f^*+`abs\_tolerance`,~f^*\times[1-`rel\_tolerance`]).
```
    + If the objective was maximized to a value ``f^*\ge0``, then the constraint ensures that the this objective value is no less than
```math
\max(f^*-`abs\_tolerance`,~f^*\times[1-`rel\_tolerance`]).
```
    + If the objective was maximized to a value ``f^*\lt0``, then the constraint ensures that the this objective value is no less than
```math
\max(f^*-`abs\_tolerance`,~f^*\times[1+`rel\_tolerance`]).
```
* Minimize/maximize with respect to the linear objective of next highest priority, and then add a corresponding objective constraint to the model, repeating until optimization with respect to the linear objective of lowest priority has taken place.

Note

* Negative values of `abs_tolerance` and `rel_tolerance` will be ignored. This is a convenient way of "switching off" a bounding technique that is not of interest.
* When the model is continuous, no dual information will be returned if there is more than one linear objective.


