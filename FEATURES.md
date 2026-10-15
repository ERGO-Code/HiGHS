The highlights of v1.16 are 

## Code changes

Rather than provide an explicit Hessian matrix, users can provide an
"oracle" callback that provides the necessary data associated with the
Hessian.

The LP file reader has been rewritten as a parser for a formal grammar
of the LP file format. Invalid files are now rejected, rather than
being read silently as a different model, and syntax errors are
reported with the line and column of the error. The reader is also
faster and uses less memory.

## Build changes

