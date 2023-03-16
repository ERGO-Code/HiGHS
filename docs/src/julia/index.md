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

## Issues and feedback

HiGHS.jl is maintained by the JuMP community and is not officially maintained
or supported by the HiGHS developers.

To report a problem (e.g., incorrect results, or a crash of the solver),
or make a suggestion for how to improve HiGHS.jl, please
[file a GitHub issue at HiGHS.jl](https://github.com/jump-dev/HiGHS.jl).

If you use HiGHS from JuMP, use `JuMP.write_to_file(model, "filename.mps")`
to write your model an MPS file, then upload the MPS file to [https://gist.github.com](https://gist.github.com)
and provide a link to the gist in the GitHub issue.
