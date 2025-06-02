# [Executable](@id executable)

HiGHS can run as a standalone program with a command-line
interface. It solves an optimization problem provided by either an
[MPS](https://docs.gurobi.com/projects/optimizer/en/current/reference/fileformats/modelformats.html#mps-format)
file, or
[LP](https://docs.gurobi.com/projects/optimizer/en/current/reference/fileformats/modelformats.html#lp-format)
file. Note that HiGHS cannot read the [lpsolve LP file
format](https://lpsolve.sourceforge.net/5.5/lp-format.htm).

### Running the executable

For convenience, the executable is assumed to be `bin/highs`.
The model given by the MPS file `model.mps` is solved by the command:

```shell
$ bin/highs model.mps
```

If the model file is not in the folder from which the command was issued, then a
path name can be given.

### Command line options

When HiGHS is run from the command line, some fundamental option
values may be specified directly. Many more may be specified via a
file containing HiGHS options settings. Formally, the usage is:

```shell
$ bin/highs --help
usage:
      ./bin/highs [options] [file]

options:
      --model_file file          File of model to solve.
      --options_file file        File containing HiGHS options.
      --read_solution_file file  File of solution to read.
      --read_basis_file text     File of initial basis to read. 
      --write_model_file text    File for writing out model.
      --solution_file text       File for writing out solution.
      --write_basis_file text    File for writing out final basis.
      --presolve text            Set presolve option to:
                                   "choose" * default 
                                   "on"
                                   "off"
      --solver text              Set solver option to: 
                                   "choose" * default 
                                   "simplex"
                                   "ipm" 
      --parallel text            Set parallel option to: 
                                   "choose" * default 
                                   "on" 
                                   "off" 
      --run_crossover text       Set run_crossover option to: 
                                   "choose" 
                                   "on" * default 
                                   "off" 
      --time_limit float         Run time limit (seconds - double).
      --random_seed int          Seed to initialize random number 
                                 generation.
      --ranging text             Compute cost, bound, RHS and basic 
                                 solution ranging:
                                   "on" 
                                   "off" * default 
  -v, --version                  Print version.
  -h, --help                     Print help.
```

The [list of options](@ref option-definitions) section gives a full
list of options, and their default values. In a file containing HiGHS
options they are specified as `name = value`, with one per line, and
any line beginning with `#` treated as a comment.

### Return code values

Consistent with the callable methods in HiGHS, there are three possible return codes

 * -1: An error has occurred in HiGHS
 * 0: HiGHS has run successfully
 * 1: HiGHS has recovered from an unusual event, or has terminated due to reaching a time or iteration limit
