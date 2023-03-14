For convenience, the executable is assumed to be `bin/highs`

### Running the executable

The model given by the MPS file `model.mps` is solved by the command

```bash
  bin/highs model.mps
```

If the model file is not in the folder from which the command was issued, then a path name can be given

### Command line options

When HiGHS is run from the command line, some fundamental option values may be specified directly. Many more may be specified via a file. Formally, the usage is

```bash
  bin/highs --help
HiGHS options
Usage:
  bin/highs [OPTION...] [file]

      --model_file arg          File of model to solve.
      --read_solution_file arg  File of solution to read.
      --options_file arg        File containing HiGHS options.
      --presolve arg            Presolve: "choose" by default - "on"/"off"
                                are alternatives.
      --solver arg              Solver: "choose" by default - "simplex"/"ipm"
                                are alternatives.
      --parallel arg            Parallel solve: "choose" by default -
                                "on"/"off" are alternatives.
      --run_crossover arg       Run crossover: "on" by default -
                                "choose"/"off" are alternatives.
      --time_limit arg          Run time limit (seconds - double).
      --solution_file arg       File for writing out model solution.
      --write_model_file arg    File for writing out model.
      --random_seed arg         Seed to initialize random number generation.
      --ranging arg             Compute cost, bound, RHS and basic solution
                                ranging.
      --version                 Print version.
  -h, --help                    Print help.
```

The [Options](https://ergo-code.github.io/HiGHS/options/definitions.html) section gives a full list of options, and the format in which they are specified.

