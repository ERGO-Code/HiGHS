# HighsDocs.jl

Documentation for [HiGHS](https://github.com/ERGO-Code/HiGHS)

## Editing the Docs

To edit the documentation, simply edit the Markdown files in [docs/src](https://github.com/galabovaa/HighsDocs.jl)
The new documentation files will be pushed to the [`gh-pages`](https://github.com/galabovaa/HighsDocs.jl/tree/gh-pages) branch in this repository (useful for debugging the `HighsDocs.jl` documentation build through CI).

## Building the docs

With proper dependencies installed, run `GKSwstype=nul julia --project=docs/ docs/make.jl`.