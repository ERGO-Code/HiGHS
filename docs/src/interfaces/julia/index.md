# HiGHS.jl

[HiGHS.jl](https://github.com/jump-dev/HiGHS.jl) is a Julia package that
interfaces with HiGHS.

HiGHS.jl has two components:

 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `HiGHS.Highs_xxx` functions, where the names and
arguments are identical to the C API.

## Installation

Install HiGHS as follows:
```julia
import Pkg
Pkg.add("HiGHS")
```

In addition to installing the HiGHS.jl package, this command will also download
and install the HiGHS binaries. (You do not need to install or compile HiGHS
separately.)

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

Pass `HiGHS.Optimizer` to `JuMP.Model` to create a JuMP model with HiGHS as the
optimizer. Set options using `set_optimizer_attribute`.

```julia
using JuMP
import HiGHS
model = Model(HiGHS.Optimizer)
set_optimizer_attribute(model, "presolve", "on")
set_optimizer_attribute(model, "time_limit", 60.0)
```

For more details, including a range of tutorials and examples using HiGHS, see
the [JuMP documentation](https://jump.dev/JuMP.jl/stable/).

## Issues and feedback

HiGHS.jl is maintained by the JuMP community and is not officially maintained
or supported by the HiGHS developers.

To report a problem (e.g., incorrect results, or a crash of the solver),
or make a suggestion for how to improve HiGHS.jl, please
[file a GitHub issue at HiGHS.jl](https://github.com/jump-dev/HiGHS.jl).

If you use HiGHS from JuMP, use `JuMP.write_to_file(model, "filename.mps")`
to write your model an MPS file, then upload the MPS file to [https://gist.github.com](https://gist.github.com)
and provide a link to the gist in the GitHub issue.

## C API

HiGHS.jl is a thin wrapper around the complete [HiGHS C API](@ref c-api).

As a basic example, we solve the model:

```math
\begin{aligned}
\min                \quad & x + y \\
\textrm{subject to} \quad & 5 \le x + 2y \le 15   \\
                          & 6 \le 3x + 2y         \\
                          & 0 \le x \le 4         \\
                          & 1 \le y               \\
                          & y \in \mathbb{Z}.

\end{aligned}
```

Here is the corresponding Julia code:

```julia
julia> using HiGHS

julia> highs = Highs_create()
Ptr{Nothing} @0x00007fc4557d3200

julia> ret = Highs_setBoolOptionValue(highs, "log_to_console", false)
0

julia> @assert ret == 0  # If ret != 0, something went wrong

julia> Highs_addCol(highs, 1.0, 0.0, 4.0, 0, C_NULL, C_NULL)   # x is column 0
0

julia> Highs_addCol(highs, 1.0, 1.0, Inf, 0, C_NULL, C_NULL)   # y is column 1
0

julia> Highs_changeColIntegrality(highs, 1, kHighsVarTypeInteger)
0

julia> Highs_changeObjectiveSense(highs, kHighsObjSenseMinimize)
0

julia> senseP = Ref{Cint}(0)  # Instead of passing `&sense`, pass a Julia `Ref`
Base.RefValue{Int32}(0)

julia> Highs_getObjectiveSense(model, senseP)
0

julia> senseP[] == kHighsObjSenseMinimize  # Dereference with senseP[]
true

julia> Highs_addRow(highs, 5.0, 15.0, 2, Cint[0, 1], [1.0, 2.0])
0

julia> Highs_addRow(highs, 6.0, Inf, 2, Cint[0, 1], [3.0, 2.0])
0

julia> Highs_run(highs)
0

julia> col_value = zeros(Cdouble, 2);

julia> Highs_getSolution(highs, col_value, C_NULL, C_NULL, C_NULL)
0

julia> col_value
2-element Vector{Float64}:
 1.0
 2.0

julia> Highs_destroy(highs)
```
