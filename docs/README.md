# Documentation

This directory contains the source files for the [documentation](https://ergo-code.github.io/HiGHS).

## Editing the documentation

To edit the documentation, checkout a branch and edit the Markdown files in the
`src` directory.

## Building the documentation

To build locally, [install Julia](https://julialang.org/downloads/), then, from
the repository root, run:

``` bash
$ julia --project=docs -e 'using Pkg; Pkg.instantiate()'
$ julia --project=docs docs/make.jl
```

The first time you run this command, Julia will download and install the
necessary packages. This may take a couple of minutes.

The website is generated in the `build/` folder. To check it out, load
`build/index.html` in your browser.

## Deploying the documentation

The documentation is automatically built and deployed by a GitHub action. You
should not check the `build/` directory into git.

## Building the Python (highspy) API documentation

The `highspy` API reference under `python/` is authored in Markdown (using
the [MyST](https://myst-parser.readthedocs.io/) Sphinx extension) and
generated separately with Sphinx, rendered directly to Markdown. The
generated Markdown lives directly under `src/interfaces/python/` and **is
checked into git** — it is not rebuilt automatically by CI. Whenever you
change a file under `docs/python/` (or the `highspy` docstrings/API), you
must rebuild and commit the regenerated Markdown yourself. From the
repository root, run:

``` bash
$ pip install sphinx sphinx_markdown_builder myst-parser
$ pip install ./highspy  # or otherwise ensure `import highspy` resolves to a real, built copy
$ sphinx-build -b markdown -d docs/python/_doctrees docs/python docs/src/interfaces/python
```

To remove the doctree build cache (not the generated Markdown, which is
committed), run:

``` bash
$ rm -rf docs/python/_doctrees
```
