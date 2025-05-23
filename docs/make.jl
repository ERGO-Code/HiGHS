#=* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                       *
*    This file is part of the HiGHS linear optimization suite           *
*                                                                       *
*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    *
*    Leona Gottwald and Michael Feldmeier                               *
*                                                                       *
*    Available as open-source under the MIT License                     *
*                                                                       *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *=#

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

import Documenter

# ==============================================================================
#  Parse and build docstrings from the C API
# ==============================================================================

const libhighs = ""
include(joinpath(@__DIR__, "c_api_gen", "build.jl"))
include(joinpath(@__DIR__, "c_api_gen", "libhighs.jl"))

# ==============================================================================
#  Make the documentation
# ==============================================================================

Documenter.makedocs(
    sitename = "HiGHS Documentation",
    authors = "Julian Hall and Ivet Galabova",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        # prettyurls = !("local" in ARGS),
        prettyurls = get(ENV, "CI", nothing) == "true",
        highlights = ["yaml"],
    ),
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
    repo = "https://github.com/ERGO-Code/HiGHS/tree/latest{path}",
    linkcheck = false,
    linkcheck_ignore = [
        "https://crates.io/crates/highs",
        "https://crates.io/crates/good_lp",
        "https://link.springer.com/article/10.1007/s12532-017-0130-5",
	"https://link.springer.com/article/10.1007/s12532-020-00181-8",
	"https://github.com/ERGO-Code/HiGHS/blob/master/highs/Highs.h"
    ],
    pages = [
        "About" => "index.md",
        "installation.md",
        "Executable" => "executable.md",
        "Guide" => Any[
            "guide/index.md",
            "guide/basic.md",
            "guide/further.md",
            "guide/advanced.md",
            "guide/gpu.md",
            "guide/kkt.md"
        ],
	"Data structures" => Any[
	    "structures/index.md",
	    "structures/enums.md",
	    "Classes" => Any[
                "structures/classes/index.md",
                "structures/classes/HighsSparseMatrix.md",
                "structures/classes/HighsLp.md",
                "structures/classes/HighsHessian.md",
                "structures/classes/HighsModel.md"
            ],
	    "Structures" => Any[
                "structures/structs/index.md",
                "structures/structs/HighsSolution.md",
                "structures/structs/HighsBasis.md",
                "structures/structs/HighsInfo.md",
                "structures/structs/HighsLinearObjective.md"
            ],
	],
        "Callbacks" => "callbacks.md",
        "Interfaces" => Any[
            "C++" => Any[
                "interfaces/cpp/index.md",
                "The HiGHS library" => "interfaces/cpp/library.md",
                "Examples" => "interfaces/cpp/examples.md",
            ],
            "C" => "interfaces/c_api.md",
            "Fortran" => "interfaces/fortran.md",
            "Python" => Any[
                "interfaces/python/index.md",
                "interfaces/python/example-py.md",
            ],
            "CSharp" => "interfaces/csharp.md",
            "Julia" => "interfaces/julia/index.md",
            "Other" => "interfaces/other.md",
        ],
        "Options" => Any[
            "options/intro.md",
            "options/definitions.md"
        ],
        "Parallel" => "parallel.md",
        "Solvers" => "solvers.md",
        "Terminology" => "terminology.md",
    ],
)

# ==============================================================================
#  Deploy everything in `build`
# ==============================================================================

Documenter.deploydocs(;
    repo = "github.com/ERGO-Code/HiGHS.git",
    push_preview = true,
    devbranch = "latest",
)
