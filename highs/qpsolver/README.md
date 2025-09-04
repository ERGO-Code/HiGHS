# HiGHS Active set QP solver (QUASS)

This document is not a full accout of how active set methods solver
QPs (see, for example, Fletcher's _Practical Methods of Optimization_,
Chapter 10), rather some background information to help understand the
essence of the HiGHS solver. For a

## Active set vs simplex for LP

To document the active set QP solver, it is convenient to start with
observations on the active set method (ASM) for LP. Although the ASM
is equivalent to the primal simplex algorithm, this is obscured by
differences in terminology, description and underlying linear algebra.

Whilst the simplex algorithm considers an LP in standard form with
$n+m$ variables ($n$ original and $m$ slack), the ASM considers an LP
with $m+n$ constraints ($m$ original and $n$ corresponding to variable
bounds). Hence, whilst the simplex algorithm definition relates to
"variables" (that are original or slack) the ASM relates to
"constraints" (that are original or variable bounds).

Whereas the simplex method is based on the factorization of a matrix
of dimension $m$ (whose columns correspond to the basic variables) the
ASM is based on the factorization of a matrix of dimension $n$. The
rows of this matrix correspond to a set of active constraints
(original constraints at a bound and variables at a bound). As such,
the set of active constraints in the ASM is identical to the set of
nonbasic variables in the simplex algorithm. Confusingly, and due to
the role of this matrix, the set of $n$ active constraints in the ASM
is referred to as the "basis". The remaining $m$ constraints in the
ASM are analogous to the set of basic variables in the simplex
algorithm, in that their values are solved for.

For the simplex method, below is an illustration of the partition of
the original and slack variables into basic and nonbasic variables

| Basic original | Basic slack | Nonbasic slack | Nonbasic original |
| --- | --- | --- | --- |
| $B_0$ | $0$ | $I_N$ | $N_0$ |
| $B_1$ | $I_B$ | $0$ | $N_1$ |

Hence the simplex method factorizes the matrix
$\left[\begin{matrix}B_0 & 0\\ B_1 &
I_B\end{matrix}\right]\in\R^{m\times m}$.  For the ASM, below is an
illustration of the partition of the original constraints and variable
bounds into active and inactive constraints

ASM | Simplex basic | Simplex nonbasic |
| --- | --- | --- |
Active original constraints | $B_0$ | $N_0$ |
Active variable bounds| $0$ | $I_N$ |
Inactive variable bounds | $I_B$ | $0$ |
Inactive original constraints | $B_1$ | $N_1$ |

Hence the ASM factorizes the matrix $\left[\begin{matrix}B_0 & N_0\\ 0
& I_N\end{matrix}\right]\in\R^{n\times n}$.  Note that the
factorization requirement for simplex and ASM is essentially the same,
since it is restricted to the "kernel" matrix $B_0$. This is due to
the matrices $I_B$ for simplex and $I_N$ for ASM yielding singleton
pivots in Gaussian elimination.

## An active set method for QP

The QP problem solved by HiGHS is

$$
\min \quad q(x) = \dfrac{1}{2}x^TQx + c^Tx \qquad \textrm{s.t.}~ \quad L \leq Ax \leq U; \quad l \leq x \leq u,\qquad(1)
$$

where $Q$ must be positive semi-definite. Unlike LP where (if an
optimal solution exists) there is an optimal solution at a vertex of
the feasible region, in general the optimal solution of a QP is not a
vertex. This is trivially so in the case of unconstrained convex QP
problems. When a QP has linear constraints, so long as the no solution
of the unconstrained QP is feasible, there will be at least one active
constraint at an optimal solution. The null space orthogonal to these
active constraints is fundamental to the properties of an optimal
solution, and the characterisation of constraint null spaces is
critical to finding an optimal solution.

### Karush-Kuhn-Tucker (KKT) conditions for convex QP

A point $x$ is an optimal solution of (1) iff $x$ is feasible and
there is a set $\cal{A}$ of indices of active constraints (and
bounds), and Lagrange multipliers $\lambda$, such that
$Qx+c=A_{\cal{A}}^T\lambda$ (where $A_{\cal{A}}$ is the matrix
corresponding to the active constraints) and the signs of the
components of $\lambda$ are consistent with the bounds at which the
constraint is active. Note that $\cal{A}$ need only include the
indices of active constraints with nonzero Lagrange multiplier, so it
can be assumed that $A_{\cal{A}}$ is of full rank.

The aim of a QP active set method is to find $\cal{A}$ such that
$A_{\cal{A}}$ is of full rank and the KKT conditions hold. 

### Some linear algebra and calculus

For a set $\cal{A}$ and matrix $A_{\cal{A}}$, there is a "null space"
normal to the row space of $A_{\cal{A}}$. This null space plays a
major role in defining an active set method for QP. Suppose that
$r=|\cal{A}|$, so $A_{\cal{A}}\in\R^{r\times n}$, and the matrix
$Y^T\in\R^{r\times n}$ spans the row space of $A_{\cal{A}}$ such that
$A_{\cal{A}}Y=I$ (so $Y$ is the right inverse of $A_{\cal{A}}$). Let
$Z^T\in\R^{(n-r)\times n}$ span the null space of $A_{\cal{A}}$, so
$A_{\cal{A}}Z=0$.

With these definitions, note that the point $x=Yb_{\cal{A}} + Zy$,
where $b_{\cal{A}}$ are the active bounds for the constraints indexed
by $\cal{A}$, satisfies
$$
A_{\cal{A}}x=A_{\cal{A}}Yb_{\cal{A}} + A_{\cal{A}}Zy = b_{\cal{A}}.
$$
Hence $y$ parameterises $x$ in the null space.

Active set methods form and operate with $Y$ and $Z$ by identifying a
matrix $V\in\R^{(n-r)\times n}$ such that
$B_{\cal{A}}\displaystyle\left[\begin{matrix}A_{\cal{A}}\\V\end{matrix}\right]$
is well conditioned, so
$\left[\begin{matrix}Y&Z\end{matrix}\right]=B_{\cal{A}}^{-1}$. Note
that $Y$ and $Z$ are not formed explicitly since, for example,
$Yb_{\cal{A}}$ can be formed by factorizing $B$ and solving
$$
B_{\cal{A}}y=\left[\begin{matrix}b_{\cal{A}}\\0\end{matrix}\right]\qquad(2)
$$

Note that substituting $x=Yb_{\cal{A}} + Zy$ into the objective yields
$$
q(y) = \frac{1}{2} y^TZ^TQZy +(c+QYb_{\cal{A}})^TZy +\left(c+\frac{1}{2}QYb_{\cal{A}}\right)^TYb_{\cal{A}}
$$
Hence the minimizer of $q(y)$ is a solution of
$$
Z^TQZy = -Z^T(c+QYb_{\cal{A}})\qquad(3)
$$
The matrix $Z^TQZ$ is known as the **reduced Hessian**.



and its properties are critical to the performance of 