push!(LOAD_PATH,"src/")

Pkg.instantiate()

using Documenter, HighsDocs

# makedocs(sitename="HiGHS Documentation",format = Documenter.HTML(
# ))


makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        # prettyurls = !("local" in ARGS),
        prettyurls = get(ENV, "CI", nothing) == "true",
        highlights = ["yaml"],
    ),
    clean = false,
    sitename = "HiGHS Documentation",
    authors = "Julian Hall, Ivet Galabova, Leona Gottwald and Michael Feldmeier.",
    pages = [
        "Home" => "index.md",
        "About" => "about.md", 
        "Get Started" => "get-started.md",
        "Running HiGHS" => "run-executable.md",
        "Examples" =>"examples.md",
        "HiGHS Library" => Any[
            "Guide" => "man/guide.md",
            "Library" => "man/library.md",
            "Model" => "man/model-definition.md",
            "Options" => "man/options.md",
            "Linking" => "man/link.md",
        ],
    ],
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
)