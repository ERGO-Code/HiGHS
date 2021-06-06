---
title: Running the HiGHS executable
permalink: /run-executable/
---
#### Input file formats

HiGHS can parse .mps and .lp files. Models can also be loaded at runtime from another program using the library interface.

# Running the executable

Assuming the executable was created following the [build]() instructions, and the lp model is specified in `model.mps` (see more on LP [Input file formats])In the following discussion, the name of the executable file created in build/bin when building HiGHS is assumed to be highs. HiGHS can read plain text MPS files and LP files (but not compressed files), and the following command solves the model in model.mps

```bash
  /build/bin/highs /path/to/model.mps
```

### Command line options

When HiGHS is run from the command line, some fundamental option values may be specified directly. Many more may be specified via a file. Formally, the usage is

```bash
./build/bin/highs --help
Running HiGHS 1.0.0 [date: 2020-03-23, git hash: b83a87d9]
Copyright (c) 2020 ERGO-Code under MIT licence terms

HiGHS options
Usage:
  ./build/bin/highs [OPTION...] [file]

      --model_file arg    File of model to solve.
      --presolve arg      Presolve: "choose" by default - "on"/"off" are
                          alternatives.
      --solver arg        Solver: "choose" by default - "simplex"/"ipm" are
                          alternatives.
      --parallel arg      Parallel solve: "choose" by default - "on"/"off"
                          are alternatives.
      --time_limit arg    Run time limit (double).
      --options_file arg  File containing HiGHS options.
  -h, --help              Print help.
```

