* Optimal objective is -4
NAME          Comment
* Comment
ROWS
 N  COST $Comment
 G  R0
 G  R1 *Comment
 G  R2 $Comment
COLUMNS
    C0        R0            1.0        COST       -1.0 $Comment
    C0        R1           -1.0 $Comment
    C0        R2           -1.0 $Comment
    C1        R0            1.0        COST       -2.0 *Comment
    C1        R1           -1.0 *Comment
    C1        R2           -1.0 *Comment
RHS
    DEMANDS   R0            1.0        R1         -2.0 $Comment
    DEMANDS   R2           -2.0 $Comment
RANGES
    test      R0        1  $Comment
BOUNDS
 UP SERVINGS  C0            1.0 $Comment
ENDATA
