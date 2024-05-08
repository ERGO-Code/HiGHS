# [List of options](@id option-definitions)

## presolve
- Presolve option: "off", "choose" or "on"
- Type: string
- Default: "choose"

## solver
- Solver option: "simplex", "choose", "ipm" or "pdlp". If "simplex"/"ipm"/"pdlp" is chosen then, for a MIP (QP) the integrality constraint (quadratic term) will be ignored
- Type: string
- Default: "choose"

## parallel
- Parallel option: "off", "choose" or "on"
- Type: string
- Default: "choose"

## run\_crossover
- Run IPM crossover: "off", "choose" or "on"
- Type: string
- Default: "on"

## time\_limit
- Time limit (seconds)
- Type: double
- Range: [0, inf]
- Default: inf

## ranging
- Compute cost, bound, RHS and basic solution ranging: "off" or "on"
- Type: string
- Default: "off"

## infinite\_cost
- Limit on |cost coefficient|: values greater than or equal to this will be treated as infinite
- Type: double
- Range: [1e+15, inf]
- Default: 1e+20

## infinite\_bound
- Limit on |constraint bound|: values greater than or equal to this will be treated as infinite
- Type: double
- Range: [1e+15, inf]
- Default: 1e+20

## small\_matrix\_value
- Lower limit on |matrix entries|: values less than or equal to this will be treated as zero
- Type: double
- Range: [1e-12, inf]
- Default: 1e-09

## large\_matrix\_value
- Upper limit on |matrix entries|: values greater than or equal to this will be treated as infinite
- Type: double
- Range: [1, inf]
- Default: 1e+15

## primal\_feasibility\_tolerance
- Primal feasibility tolerance
- Type: double
- Range: [1e-10, inf]
- Default: 1e-07

## dual\_feasibility\_tolerance
- Dual feasibility tolerance
- Type: double
- Range: [1e-10, inf]
- Default: 1e-07

## ipm\_optimality\_tolerance
- IPM optimality tolerance
- Type: double
- Range: [1e-12, inf]
- Default: 1e-08

## objective\_bound
- Objective bound for termination of the dual simplex solver
- Type: double
- Range: [-inf, inf]
- Default: inf

## objective\_target
- Objective target for termination of the MIP solver
- Type: double
- Range: [-inf, inf]
- Default: -inf

## random\_seed
- Random seed used in HiGHS
- Type: integer
- Range: {0, 2147483647}
- Default: 0

## threads
- Number of threads used by HiGHS (0: automatic)
- Type: integer
- Range: {0, 2147483647}
- Default: 0

## user\_bound\_scale
- Exponent of power-of-two bound scaling for model
- Type: integer
- Range: {-2147483647, 2147483647}
- Default: 0

## user\_cost\_scale
- Exponent of power-of-two cost scaling for model
- Type: integer
- Range: {-2147483647, 2147483647}
- Default: 0

## simplex\_strategy
- Strategy for simplex solver 0 => Choose; 1 => Dual (serial); 2 => Dual (PAMI); 3 => Dual (SIP); 4 => Primal
- Type: integer
- Range: {0, 4}
- Default: 1

## simplex\_scale\_strategy
- Simplex scaling strategy: off / choose / equilibration / forced equilibration / max value 0 / max value 1 (0/1/2/3/4/5)
- Type: integer
- Range: {0, 5}
- Default: 1

## simplex\_dual\_edge\_weight\_strategy
- Strategy for simplex dual edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)
- Type: integer
- Range: {-1, 2}
- Default: -1

## simplex\_primal\_edge\_weight\_strategy
- Strategy for simplex primal edge weights: Choose / Dantzig / Devex / Steepest Edge (-1/0/1/2)
- Type: integer
- Range: {-1, 2}
- Default: -1

## simplex\_iteration\_limit
- Iteration limit for simplex solver when solving LPs, but not subproblems in the MIP solver
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## simplex\_update\_limit
- Limit on the number of simplex UPDATE operations
- Type: integer
- Range: {0, 2147483647}
- Default: 5000

## simplex\_max\_concurrency
- Maximum level of concurrency in parallel simplex
- Type: integer
- Range: {1, 8}
- Default: 8

## output\_flag
- Enables or disables solver output
- Type: boolean
- Default: "true"

## log\_to\_console
- Enables or disables console logging
- Type: boolean
- Default: "true"

## solution\_file
- Solution file
- Type: string
- Default: ""

## log\_file
- Log file
- Type: string
- Default: ""

## write\_solution\_to\_file
- Write the primal and dual solution to a file
- Type: boolean
- Default: "false"

## write\_solution\_style
- Style of solution file (raw = computer-readable, pretty = human-readable): -1 => HiGHS old raw (deprecated); 0 => HiGHS raw; 1 => HiGHS pretty; 2 => Glpsol raw; 3 => Glpsol pretty; 4 => HiGHS sparse raw
- Type: integer
- Range: {-1, 4}
- Default: 0

## glpsol\_cost\_row\_location
- Location of cost row for Glpsol file: -2 => Last; -1 => None; 0 => None if empty, otherwise data file location; 1 <= n <= num\_row => Location n; n > num\_row => Last
- Type: integer
- Range: {-2, 2147483647}
- Default: 0

## write\_model\_file
- Write model file
- Type: string
- Default: ""

## write\_model\_to\_file
- Write the model to a file
- Type: boolean
- Default: "false"

## mip\_detect\_symmetry
- Whether MIP symmetry should be detected
- Type: boolean
- Default: "true"

## mip\_allow\_restart
- Whether MIP restart is permitted
- Type: boolean
- Default: "true"

## mip\_max\_nodes
- MIP solver max number of nodes
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## mip\_max\_stall\_nodes
- MIP solver max number of nodes where estimate is above cutoff bound
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## mip\_improving\_solution\_save
- Whether improving MIP solutions should be saved
- Type: boolean
- Default: "false"

## mip\_improving\_solution\_report\_sparse
- Whether improving MIP solutions should be reported in sparse format
- Type: boolean
- Default: "false"

## mip\_improving\_solution\_file
- File for reporting improving MIP solutions: not reported for an empty string \"\"
- Type: string
- Default: ""

## mip\_max\_leaves
- MIP solver max number of leave nodes
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## mip\_max\_improving\_sols
- Limit on the number of improving solutions found to stop the MIP solver prematurely
- Type: integer
- Range: {1, 2147483647}
- Default: 2147483647

## mip\_lp\_age\_limit
- Maximal age of dynamic LP rows before they are removed from the LP relaxation in the MIP solver
- Type: integer
- Range: {0, 32767}
- Default: 10

## mip\_pool\_age\_limit
- Maximal age of rows in the MIP solver cutpool before they are deleted
- Type: integer
- Range: {0, 1000}
- Default: 30

## mip\_pool\_soft\_limit
- Soft limit on the number of rows in the MIP solver cutpool for dynamic age adjustment
- Type: integer
- Range: {1, 2147483647}
- Default: 10000

## mip\_pscost\_minreliable
- Minimal number of observations before MIP solver pseudo costs are considered reliable
- Type: integer
- Range: {0, 2147483647}
- Default: 8

## mip\_min\_cliquetable\_entries\_for\_parallelism
- Minimal number of entries in the MIP solver cliquetable before neighbourhood queries of the conflict graph use parallel processing
- Type: integer
- Range: {0, 2147483647}
- Default: 100000

## mip\_feasibility\_tolerance
- MIP feasibility tolerance
- Type: double
- Range: [1e-10, inf]
- Default: 1e-06

## mip\_heuristic\_effort
- Effort spent for MIP heuristics
- Type: double
- Range: [0, 1]
- Default: 0.05

## mip\_rel\_gap
- Tolerance on relative gap, |ub-lb|/|ub|, to determine whether optimality has been reached for a MIP instance
- Type: double
- Range: [0, inf]
- Default: 0.0001

## mip\_abs\_gap
- Tolerance on absolute gap of MIP, |ub-lb|, to determine whether optimality has been reached for a MIP instance
- Type: double
- Range: [0, inf]
- Default: 1e-06

## mip\_min\_logging\_interval
- MIP minimum logging interval
- Type: double
- Range: [0, inf]
- Default: 5

## ipm\_iteration\_limit
- Iteration limit for IPM solver
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## pdlp\_native\_termination
- Use native termination for PDLP solver: Default = false
- Type: boolean
- Default: "false"

## pdlp\_scaling
- Scaling option for PDLP solver: Default = true
- Type: boolean
- Default: "true"

## pdlp\_iteration\_limit
- Iteration limit for PDLP solver
- Type: integer
- Range: {0, 2147483647}
- Default: 2147483647

## pdlp\_e\_restart\_method
- Restart mode for PDLP solver: 0 => none; 1 => GPU (default); 2 => CPU 
- Type: integer
- Range: {0, 2}
- Default: 1

## pdlp\_d\_gap\_tol
- Duality gap tolerance for PDLP solver: Default = 1e-4
- Type: double
- Range: [1e-12, inf]
- Default: 0.0001

