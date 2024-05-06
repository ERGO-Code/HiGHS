#=* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                       *
*    This file is part of the HiGHS linear optimization suite           *
*                                                                       *
*    Written and engineered 2008-2023 by Julian Hall, Ivet Galabova,    *
*    Leona Gottwald and Michael Feldmeier                               *
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
    linkcheck = true,
    linkcheck_ignore = [
        "https://crates.io/crates/highs",
        "https://crates.io/crates/good_lp",
        "https://link.springer.com/article/10.1007/s12532-017-0130-5",
    ],
    pages = [
        "About" => "index.md",
        "installation.md",
        "Executable" => "executable.md",
        "Guide" => Any[
            "guide/index.md",
            "guide/basic.md",
            "guide/further.md",
            "guide/advanced.md"
        ],
	"Data structures" => Any[
	    "structures/index.md",
	    "structures/enums.md",
	    "Classes" => Any[
                "structures/classes/index.md",
                "structures/classes/HighsSparseMatrix.md",
                "structures/classes/HighsLp.md",
                "structures/classes/HighsSolution.md",
                "structures/classes/HighsBasis.md",
                "structures/classes/HighsInfo.md",
            ],
	],
        "Callbacks" => "callbacks.md",
        "Interfaces" => Any[
            "Python" => Any[
                "interfaces/python/index.md",
                "interfaces/python/example-py.md",
            ],
            "C++" => Any[
                "interfaces/cpp/index.md",
                "The HiGHS library" => "interfaces/cpp/library.md",
                "Linking" => "interfaces/cpp/link.md",
                "Examples" => "interfaces/cpp/examples.md",
            ],
            "C" => "interfaces/c/index.md",
            "Julia" => "interfaces/julia/index.md",
            "Other" => "interfaces/other.md",
        ],
        "Options" => Any[
            "options/intro.md",
            "options/definitions.md"
        ],
        "Parallel" => "parallel.md",
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
