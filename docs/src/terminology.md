# Terminology

Any linear optimization model will have __decision variables__, a
linear or quadratic __objective function__, and linear __constraints__
and __bounds__ on the values of the decision variables. A
__mixed-integer__ optimization model will require some or all of the
decision variables to take integer values. The model may require the
objective function to be maximized or minimized whilst satisfying the
constraints and bounds. By default, HiGHS minimizes the objective
function.

## Bounds and the objective function

The bounds on a decision variable are the least and greatest values
that it may take, and infinite bounds can be specified. A linear
objective function is given by a set of coefficients, one for each
decision variable, and its value is the sum of products of
coefficients and values of decision variables. The objective
coefficients are often referred to as __costs__, and some may be
zero. When a model has been solved, the optimal values of the
decision variables are referred to as the __(primal) solution__.

## Constraints and the feasible region

Linear constraints require linear functions of decision variables to
lie between bounds, and infinite bounds can be specified. If the
bounds are equal, then the constraint is an __equation__. If the
bounds are both finite, then the constraint is said to be __boxed__ or
__two-sided__. The set of points satisfying linear constraints and
bounds is known as the __feasible region__. Geometrically, this is a
multi-dimensional convex polyhedron, whose extreme points are referred
to as __vertices__.

## The constraint matrix

The coefficients of the linear constraints are naturally viewed as
rows of a __matrix__. The constraint coefficients associated with a
particular decision variable form a column of the constraint
matrix. Hence constraints are sometimes referred to as __rows__, and
decision variables as __columns__. Constraint matrix coefficients may
be zero. Indeed, for large practical models it is typical for most
of the coefficients to be zero. When this property can be exploited to
computational advantage, the matrix is said to be __sparse__. When the
constraint matrix is not sparse, the solution of large models is
normally intractable computationally.

## Optimization outcomes

It is possible to define a set of constraints and bounds that cannot
be satisfied, in which case the model is said to be
__infeasible__. Conversely, it is possible that the value of the
objective function can be improved without bound whilst satisfying the
constraints and bounds, in which case the model is said to be
__unbounded__. If a model is neither infeasible, nor unbounded, it
has an __optimal solution__. The optimal objective function value for
a linear optimization model may be achieved at more than point, in
which case the optimal solution is said to be __non-unique__.

## Primal values

The values of the decision variables are referred to as __primal__ values to distingush them from __dual__ values.

## Dual values

When none of the decision variables is required to take integer
values, the model is said to be __continuous__. For
continuous models, each variable and constraint has an
associated __dual variable__. The values of the dual
variables constitute the __dual solution__, and it is for
this reason that the term __primal solution__ is used to
distinguish the optimal values of the decision variables. At the
optimal solution of a continuous model, some of the decision
variables and values of constraint functions will be equal to their
lower or upper bounds. Such a bound is said to
be __active__. If a variable or constraint is at a bound,
its corresponding dual solution value will generally be non-zero: when
at a lower bound the dual value will be non-negative; when at an upper
bound the dual value will be non-positive. When maximizing the
objective the required signs of the dual values are reversed. Due to
their economic interpretation, the dual values associated with
constraints are often referred to as __shadow prices__
or __fair prices__. Mathematically, the dual values
associated with variables are often referred to as __reduced
costs__, and the dual values associated with constraints are
often referred to as __Lagrange multipliers__.

## Basic solution

An LP model that is neither infeasible, nor unbounded, has an
optimal solution at a vertex. At a vertex, the decision variables can
be partitioned into as many __basic variables__ as there are
constraints, and __nonbasic variables__. Such a solution is known as a
__basic solution__, and the partition referred to as a __basis__.

## Sensitivity

Analysis of the change in optimal objective value of a continuous
linear optimization model as the cost coefficients and bounds are
changed is referred to in HiGHS as __ranging__. For an
active bound, the corresponding dual value gives the change in the
objective if that bound is increased or decreased. This level of
analysis is often referred to as __sensitivity__. In
general, the change in the objective is only known for a limited range
of values for the active bound. HiGHS will return the limits of
these __bound ranges__ ranges, the objective value at
both limits and the index of a variable or constraint that will
acquire an active bound at both limits. For each variable with an
active bound, the solution will remain optimal for a range of values
of its cost coefficient. HiGHS will return the values of
these __cost ranges__. For a variable or constraint whose
value is not at a bound, HiGHS will return the range of values that
the variable or constraint can take, the objective values at the
limits of the range, and the index of a variable or constraint with a
bound that will become in active at both limits.

## [MIP](@id terminology-mip)

When solving a MIP, some or all the variables must take discrete values. In HiGHS there are three types of discrete variables.

- Integer: those that must take integer values between their bounds
- Semi-continuous: those that must be zero or take continuous values between their bounds
- Semi-integer: those that must be zero or take integer values between their bounds

In the following discussion, for ease of reference to relative
objective values, it is assumed that the objective is being minimized

Any point for which the discrete variables satisfy their requirements,
is said to be __integer feasible__. The objective value at such a
point is an upper bound on the optimal objective value. The least such
bound is known as the __primal bound__. The MIP solver generates a
sequence of LPs, each of which has bounds on the variables that are
tighter than those of the original model. When unsolved, there is a
bound on the optimal objective value for each such LP and, the
greatest such bound is known as the __dual bound__. The optimal
objective value of the MIP cannot be less than the dual bound. Hence
the gap between the primal and dual bounds is a measure of progress of
the MIP solver. Although the absolute gap is of some interest, the gap
relative to the primal bound is a better measure. When the gap reaches
zero then the MIP is solved to optimality. However, it is often
preferable to stop the MIP solver when the relative gap is below a
specified tolerance.
