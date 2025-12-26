# [Numerical considerations](@id numerics)

Optimization solvers cannot be expected to find the exact solution of
a problem, since it may not be possible to represent that solution
using floating-point arithemtic. However, solvers will typically run
faster and find more accurate solutions if the problem has good
numerical properties. Ideally the optimal value of all primal
variables (and dual variables when relevant) will be of order
unity. This typically occurs if all objective and constraint matrix
coefficients, as well as finite variable and constraint bounds, are of
order unity. Whilst there may be some reasons why this ideal cannot be
achieved in all models, there are many pitfalls to avoid. For an
insight into reasons why a model may have bad numerical properties and
how to avoid them, users are recommended to study this [JuMP
tutorial](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/tolerances/). Improving
the numerical properties of a model will typically lead to it being
solved faster and more accurately/reliably, so the investment should
pay off!

Internally, the HiGHS continuous optimization solvers scale the
constraint matrix to improve the numerical properties of the problem,
[feasibility and optimality tolerances](@ref kkt) are determined with
respect to the original, unscaled problem. However, faced with a model
with bad numerical properties, there is only so much that HiGHS can do
to solve it efficiently and accurately.

If the optimal values of many variables in a model are typically very
large, this can correspond to very large values of the objective
coefficients and finite variable and constraint bounds. Since most of
the HiGHS solvers terminate according to small absolute [feasibility
tolerances](@ref kkt), large objective coefficients and bounds force
the solvers to achieve an accuracy that may be unrealsitic in the
context of a model. As well as having an impact on efficiency, the
solver may ultimately be unable to achieve the required accuracy and
fail. Objective coefficients and bounds that are less than the
feasibility and optimality tolerances can also be problematic,
although this is less common and less serious.

HiGHS offers a facility to enable users to assess the consequences of
better problem scaling, in cases where some objective coefficients or
bounds are large, or if all objective coefficients or bounds are
small. By setting the options [__user\_objective\_scale__](@ref
option-user_objective_scale) and/or [__user\_bound\_scale__](@ref
option-user_bound_scale), HiGHS will solve the given model with
uniform scaling of the objective coefficients or bounds. Note that
these options define the exponent in power-of-two scaling factors so
that model accuracy is not compromised. After solving the problem,
feasibility and optimality will be assessed for the original model,
with a warning given if the tolerances are not satisfied. Note that
uniform scaling of bounds on discrete variables is not possible, and
is achieved implicitly by scaling their cost and matrix
coefficients. Also, when bounds on variables in a quadratic
programming problem are scaled up (down), the values in the Hessian
matrix must be scaled down (up) so that the overall scaling of the
objective is uniform.