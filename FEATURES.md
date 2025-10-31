## Code changes

The HiPO release exposed various issues flagged up via GitHub and email.
- Fix some overflows when computing statistics of analyse phase.
- Free memory used for normal equations, if augmented system is preferred.
- Fix bug in supernode amalgamation.
- Printing of BLAS library moved to HiGHS header, so it is printed when using HiPO without logging.
- Recommend to install Metis branch `521-ts` due to better oredring quality on many problems. Update workflows and documentation accordingly.
- Add option `hipo_metis_no2hop` to control option `no2hop` of Metis and add warning if the fill-in is large.

Added singleton column stuffing to MIP presolve - see Gamrath et al., Progress in presolving for mixed integer programming. Math. Prog. Comp. 7, 367–398 (2015).

## Build changes
