# [Feasibilty and optimality](@id kkt)

Mathematically, continuous optimization problems have exact feasibilty
and optimality conditions. However, since solvers cannot always
satisfy these conditions exactly when using floating-point arithmetic,
they do so to within tolerances. As explored below, some solvers aim
to satisfy those tolerances absolutely, and others aim to satisfy
tolerances relative to problem data. When tolerances are satisfied
relatively, they are generally not satisfied absolutely. The use of
tolerances relative to problem data is not consistent across solvers,
and can give a misleading claim of optimality. To achieve consistency,
HiGHS reassesses the optimal solution claimed by such a solver in a
reasonable and uniform manner.

### Feasibilty and optimality conditions

To discuss tolerances and their use in different solvers, consider the
standard form linear programming (LP) problem with ``n`` variables and
``m`` equations (``n\ge m``).

```math
\begin{aligned}
\textrm{minimize}   \quad & c^T\! x        \\
\textrm{subject to} \quad & Ax = b  \\
                          & x \ge 0,
\end{aligned}
```

The feasibilty and optimality conditions (KKT conditions) are that, at
a point ``x``, there exist (row) dual values ``y`` and reduced costs
(column dual values) ``s`` such that

```math
\begin{aligned}
Ax=b&\qquad\textrm{Primal~equations}\\
A^Ty+s=c&\qquad\textrm{Dual~equations}\\
x\ge0&\qquad\textrm{Primal~feasibility}\\
s\ge0&\qquad\textrm{Dual~feasibility}\\
c^Tx-b^Ty=0&\qquad\textrm{Optimality}
\end{aligned}
```

The optimality condition is equivalent to the complementarity
condition that ``x^Ts=0``. Since any LP problem can be transformed
into standard form, the following discussion loses no generality. This
discussion also largely applies to quadratic programming (QP)
problems, with the differences explored below.

### The HiGHS feasibility and optimality tolerances

HiGHS has separate tolerances for the following, listed with convenient mathematical notation

- [Primal feasibility](@ref option-primal-feasibility-tolerance) (``\epsilon_P``)
- [Dual feasibility](@ref option-dual-feasibility-tolerance) (``\epsilon_D``)
- Residual errors in the [primal equations](@ref option-primal-residual-tolerance) (``\epsilon_R``)
- Residual errors in the [dual equations](@ref option-dual-residual-tolerance) (``\epsilon_C``)
- [Optimality](@ref option-optimality-tolerance) (``\epsilon_{O}``)

All are set to the same default value of ``10^{-7}``. Although each
can be set to different values by the user, if the user wishes to
solve LPs to a general lower or higher tolerance, the value of the
[KKT tolerance](@ref option-kkt-tolerance) can be changed from this
default value.

### When HiGHS yields an optimal solution

When HiGHS returns a model status of optimal, the solution will
satisfy feasibility and optimality tolerances absolutely or relatively
according to whether the solver yields a basic solution.

### Solutions with a corresponding basis

The HiGHS simplex solvers and the interior point solver after
crossover yield an optimal basic solution of the LP, consisting of
``m`` basic variables and ``n-m`` nonbasic variables. At any basis,
the nonbasic variables are zero, and values for the basic variables
are given by solving a linear system of equations. Values for the row
dual values (``y``) can be computed by solving a linear system of equations,
and the column dual values are then given by ``s=c-A^Ty``. With exact
arithmetic, the basic dual values are zero by construction.

When primal and dual values are computed using floating-point
arithmetic, the basic dual values are set to zero so the optimality
condition holds by construction. However, the primal and dual
equations may not be satisfied exactly, so have nonzero
residuals. Fortunately, when solving a linear system of equations
using a stable technique, any residuals are small relative to the RHS
of the equations, whatever the condition of the matrix of
coefficients. Hence HiGHS does not assess the primal residuals, or the
dual residuals for basic variables. Thus optimality for a basic
solution is assessed by HiGHS according to whether the following
conditions hold

```math
\begin{aligned}
x_i\ge-\epsilon_P&\qquad\forall i=1,\ldots,n\\
s_i\ge-\epsilon_D&\qquad\forall i=1,\ldots,n.
\end{aligned}
```

The HiGHS active set QP solver has an objective function ``(1/2)x^TQx + c^Tx``,
and maintains the QP equivalent of a basis in which a subset
of (up to ``n``) variables are zero. However, there are variables
that are off their bounds whose reduced costs are not zero by
construction. At an optimal solution they will only be less than a
dual feasibility tolerance in magnitude, so the optimality condition
will not be satisfied by construction. The primal and dual equations
(where the latter is ``A^Ty+s=Qx+c``) will be satisfied with small
residuals. Optimality is assessed by HiGHS according to whether primal
and dual feasibility is satisfied to within the corresponding
tolerance.


### Solutions without a corresponding basis

The HiGHS PDLP solver and the interior point solver without crossover
(IPX) yield "optimal" primal and dual values that satisfy internal
conditions for termination of the underlying algorithm. These
conditions are discussed below, and are used for good reason. However
they can lead to a misleading claim of optimality.

#### Interior point solutions

The interior point algorithm uses a single feasibility tolerance
``\epsilon=\min(\epsilon_P, \epsilon_D)``, and an independent
[optimality tolerance](@ref option-ipm-optimality-tolerance)
(``\epsilon_{IPM}``) that, by default, is (currently) ten times
smaller than the other feasibility and optimality tolerances used by
HiGHS. It terminates when

```math
\begin{aligned}
\|Ax-b\|_\infty&\le(1+\|b\|_\infty)\epsilon_R\\
\|c-A^Ty+s\|_\infty&\le(1+\|c\|_\infty)\epsilon_C\\
-x_i&\le\epsilon\qquad\forall i=1,\ldots,n\\
-s_i&\le		      \epsilon\qquad\forall i=1,\ldots,n\\
|c^Tx-b^Ty|&\le(1+|c^Tx+b^Ty|/2)\epsilon_{IPM}.
\end{aligned}
```

#### PDLP solutions

The PDLP algorithm uses an independent [optimality tolerance](@ref
option-pdlp-optimality-tolerance) (``\epsilon_{PDLP}``) that is equal
to the other feasibility and optimality tolerances used by HiGHS. It
determines values of ``x\ge0`` and ``y``, and chooses ``s`` to be the
non-negative values of ``c-A^Ty``. Hence it guarantees primal and dual
feasibility by construction. It terminates when

```math
\begin{aligned}
\|Ax-b\|_2&\le (1+\|b\|_2)\epsilon_P\\
\|c-A^Ty-s\|_2&\le (1+\|c\|_2)\epsilon_D\\
|c^Tx-b^Ty|&\le (1+|c^Tx|+|b^Ty|)\epsilon_{PDLP}.
\end{aligned}
```

#### HiGHS solutions

The relative measures used by PDLP and IPX assume that all components
of the cost and RHS vectors are relevant. When an LP problem is in
standard form this is true for ``b``, but not necessarily for the cost
vector ``c``. Consider a large component of ``c`` for which the
corresponding reduced cost value in ``s`` is also large, in which case
the LP solution is insensitive to the cost. This component will
contribute significantly to ``\|c\|`` and, hence, the RHS of the dual
residual condition, allowing large values of ``\|c-A^Ty-s\|`` to be
accepted. However, this can lead to unacceptably large absolute
residual errors and non-optimal solutions being deemed "optimal". When
equations in ``Ax=b`` correspond to inequality constraints with large
RHS values and a slack variable (so the constraint is redundant) the
same issue occurs in the case of primal residual errors. The solution
of the LP is not sensitive to this large RHS value, but its
contribution to ``||b||`` can allow large absolute primal residual
errors to be overlooked.

To make an informed assessment of whether an "optimal" solution
obtained by IPX or PDLP is acceptable, HiGHS computes infinity norm
measures of ``b`` and ``c`` corresponding to the components that
define the optimal solution. For ``c`` these are the components
corresponding to positive values of ``x`` and reduced costs that are
close to zero. For ``b``, these are the components corresponding to
constraints that are (close to being) satisfied exactly. The resulting
measures are smaller than ``\|b\|`` or ``\|c\|``, and may lead to
relative measures of primal/dual residual errors or infeasibilities
not being satisfied, so the status of the solver's "optimal" solution
may be reduced to "unknown". When this happens - and possibly if
tolerances on relative measures _have_ been satisfied - users can
consult the absolute and relative measures available via
[HighsInfo](@ref info-num-primal-infeasibilities).

### Discrete optimization problems

Discrete optimization problems, such as the mixed-integer programming
(MIP) problems solved by HiGHS, have no local optimality
conditions. Variables required to take integer values will do so to
within the `mip_feasibility_tolerance`. Since MIP sub-problems are
solved with the simplex solver, the values of the variables and
constraints will satisfy absolute feasibility tolerances. Within the
MIP solver, the value of `mip_feasibility_tolerance` is used for
`primal_feasibility_tolerance` when solving LP sub-problems, and one
tenth of this value is used for `dual_feasibility_tolerance`. Hence
any value of `primal_feasibility_tolerance` (or
`dual_feasibility_tolerance`) set by the user has no effect of the MIP
solver.

