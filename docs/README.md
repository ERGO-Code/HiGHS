# Documentation

This directory contains the source files for the [documentation](https://ergo-code.github.io/HiGHS).

## Editing the documentation

To edit the documentation, checkout a branch and edit the Markdown files in the
`src` directory.

## Building the documentation

To build locally, [install Julia](https://julialang.org/downloads/), then (from the `docs` directory) run:

```
$ julia make.jl
```

The first time you run this command, Julia will download and install the
necessary packages. This may take a couple of minutes.

The website is generated in the `build/` folder. To check it out, load
`build/index.html` in your browser.

## Deploying the documentation

The documentation is automatically built and deployed by a GitHub action. You
should not check the `build/` directory into git.
