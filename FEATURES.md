## Build changes

The python wrapper highspy is now available for aarch64 on manylinux

This allows highs to be run through Python on AWS arm64 

## Code changes

The accessor function Highs_getCallbackDataOutItem in the C API means
that `pdlp_iteration_count` can be moved back to where it was inserted
into the `HighsCallbackDataOut` struct in v1.7.0, which broke the C
API. This fixes #1812

Some duplicate code has been eliminated from the MIP solver, and
modifications made to eliminate compiler warnings

Declaration of the (deprecated) method `char* highsCompilationDate()`
has been corrected

Fixed bug when describing integrality status during the human-readable solution write

