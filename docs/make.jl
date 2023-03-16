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
    pages = [
        "About" => "index.md",
        "Guide" => "guide.md",
        "HiGHS in Python" => Any[
            "Get started in Python" => "python/pip.md",
            "Enums" => "python/enums.md",
            "Classes" => Any[
                "Introduction" => "python/classes/Index.md",
                "HighsSparseMatrix" => "python/classes/HighsSparseMatrix.md",
                "HighsLp" => "python/classes/HighsLp.md",
                "HighsSolution" => "python/classes/HighsSolution.md",
                "HighsBasis" => "python/classes/HighsBasis.md",
                "HighsInfo" => "python/classes/HighsInfo.md",
                "Other" => "python/classes/Other.md",
            ],
            "Examples" => "python/example-py.md",
            "Notebooks" => "python/notebooks.md",
        ],
        "HiGHS in C++" => Any[
            "Get started in C++" => "cpp/get-started.md",
            "The HiGHS library" => "cpp/library.md",
            "Linking" => "cpp/link.md",
            "Examples" => "cpp/examples.md",
        ],
        "HiGHS in Julia" => "julia/index.md",
        "Binaries" => "binaries.md",
        "Executable" => "executable.md",
        "Options" => Any[
            "Introduction" => "options/intro.md",
            "Definitions" => "options/definitions.md"
        ],
        "Parallel" => "parallel.md",
        "Interfaces" => "interfaces.md",
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
