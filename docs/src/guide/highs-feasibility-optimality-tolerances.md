# [HiGHS feasibilty and optimality tolerances](@id highs-feasibility-optimality-tolerances)

Mathematically, optimization problems have exact feasibilty conditions and, for continuous optimization problems, exact optimality conditions. However, since solvers cannot always satisfy these conditions exactly when using floating-point arithmetic, they do so to within tolerances. To the user, HiGHS interprets feasibilty and optimality tolerances in the same way for all classes of problem. However, internally, the HiGHS solvers use tolerances to determine feasibilty and optimality in different ways. This, together with internal problem scaling, can lead to solutions that are deemed optimal by a solver not satisfying the HiGHS "quality control" criteria and flagged up as non-optimal.

### Optimality conditions

To discuss tolerances, their use in different solvers, and how to assess the situation where solutions deemed optimal by a solver are flagged up as non-optimal by HiGHS, consider the standard form LP problem with ``n`` variables and ``m`` equations (``n\ge m``) that is assumed to have  an optimal solution.
```math
\begin{aligned}
\min                \quad & c^T\! x        \\
\textrm{subject to} \quad & Ax = b  \\
                          & x \ge 0,
\end{aligned}
```
The optimality conditions are that, at a point ``x``, there exist (row) dual values ``y`` and reduced costs (column dual values) ``s`` such that
```math
\begin{aligned}
Ax=b&\textrm{Primal~equations}\\
A^Ty+s=c&\textrm{Dual~equations}\\
x\ge0&\textrm{Primal~feasibility}\\
s\ge0&\textrm{Dual~feasibility}\\
c^Tx-b^Ty=0&\textrm{Primal-dual~gap}
\end{aligned}
```
The primal-dual gap is equivalent to the complementarity condition that `x^Ts=0`. Since any LP problem can be transformed into standard form, the following discussion loses no generality.

### The HiGHS feasibility and optimality tolerances

Within


### When HiGHS yields an optimal solution

When HiGHS returns a model status of optimal, it can be assumed that

```math
\begin{aligned}
\|Ax-b\|_\infty\le\epsilon_P&\\
\|c-A^Ty+s\|_\infty\le\epsilon_D\\
x_i\ge-\epsilon_P&\forall i=1,\ldots,n\\
s_i\ge-\epsilon_D&\forall i=1,\ldots,n\\
|c^Tx-b^Ty|\le\epsilon_{PD}
\end{aligned}
```
### Solutions with a corresponding basis

The HiGHS simplex solver and the interior point solver after crossover yield an optimal basic solution of the LP. There are ``m`` basic variables and ``n-m`` nonbasic variables. At this basis, the nonbasic variables are zero, and values for the basic variables are given by solving a linear system of equations. Values for the row dual values can be computed by solving a linear system of equations, and the column dual values are then given by ``s=c-A^Ty``. By construction, the dual values for basic variables are zero.

When primal and dual values are computed using floating-point arithmetic, the primal equations may not hold exactly so have nonzero residuals. However, numerical linear algebra theory is such that it is reasonable to assume that any residuals are small, so HiGHS does not measure them. By construction, the dual equations are satisfied exactly and the complementarity condition holds exactly. Hence optimality for a basic solution is assessed by HiGHS according to whether the following conditions hold
```math
\begin{aligned}
x_i\ge-\epsilon_P&\forall i=1,\ldots,n\\
s_i\ge-\epsilon_D&\forall i=1,\ldots,n
```
When HiGHS returns a model status of optimal, it can be assumed that these conditions hold.

### Solutions without a corresponding basis

The HiGHS PDLP solver and the interior point solver without crossover yield "optimal" primal and dual values that satisfy internal conditions for termination of the underlying algorithm. These conditions are discussed below, and are used for good reason. However they can lead to 


#### Interior point solutions

#### PDLP solutions

The PDLP algorithm determines values of ``x\ge0`` and ``y``, and chooses ``s`` to be the non-negative values of ``c-A^Ty``. Hence it guarantees primal and dual feasibility by construction. PDLP terminates when
```math
\begin{aligned}
\|Ax-b\|_2&\le \epsilon_P(1+\|b\|_2)\\
\|c-A^Ty-s\|_2&\le \epsilon_D(1+\|c\|_2)\\
|c^Tx-b^Ty|&\le \epsilon_{PD}(|c^Tx|+|b^Ty|)
\end{aligned}
```
where the value of ``\epsilon_P`` used is the HiGHS option [`primal_feasibility_tolerance`](@ref dprimal_feasibility_tolerance), the value of ``\epsilon_D`` used is the HiGHS option [`dual_feasibility_tolerance`](@ref ddual_feasibility_tolerance), and the value of ``\epsilon_{PD}`` used is the HiGHS option [`pdlp_d_gap_tol`](@ref dp_gap_tol). 

### Discrete optimization problems

Discrete optimization problems have no local optimality conditions. 

