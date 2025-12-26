NAME          WARNINGS
ROWS
 G  R0
 G  R1
 G  R2
 G  R3
 G  R4
 G  R5
 N  COST
COLUMNS
    C0        R0           -1.0        COST       -1.0
    C0        R11           1.0        R12         1.0
    C0        R13           1.0        R14         1.0
    C0        R15           1.0        COST        1.0
    C0        COST          2.0        COST        3.0
    C0        COST          4.0        COST        5.0
    C0        R0            1.0        R0          2.0
    C0        R0            3.0        R0          4.0
    C0        R0            5.0        R0          6.0
    MARK0000  'MARKER'                 'INTORG'
    C1        R0           -1.0        COST       -1.0
    C2        R0           -1.0        COST       -1.0
    C3        R0           -1.0        COST       -1.0
    C4        R0           -1.0        COST       -1.0
    C5        R0           -1.0        COST       -1.0
    MARK0001  'MARKER'                 'INTEND'
    C6        R0           -1.0        COST       -1.0
    C7        R0           -1.0        COST       -1.0
    C8        R0           -1.0        COST       -1.0
RHS
    DEMANDS   R0           -1.0        R11        -1.0
    DEMANDS   R12          -1.0        R13        -1.0
    DEMANDS   R14          -1.0        R15        -1.0
    DEMANDS   R0            1.0        R0          2.0
    DEMANDS   R0            3.0        R0          4.0
    DEMANDS   R0            5.0
BOUNDS
 UP SERVINGS  C0           -1.0
 UP SERVINGS  C0            1.0
 UP SERVINGS  C0            2.0
 UP SERVINGS  C0            3.0
 UP SERVINGS  C0            4.0
 UP SERVINGS  C0            5.0
 UP SERVINGS  C10           1.0
 LI SERVINGS  C1            1.1
 UI SERVINGS  C2            2.1
 SI SERVINGS  C3            3.1
 LI SERVINGS  C4            4.1
 UI SERVINGS  C5            5.1
 LO SERVINGS  C4            1.0
 LO SERVINGS  C6            1.0
 LO SERVINGS  C7            1.0
 LO SERVINGS  C8            1.0
 UP SERVINGS  C4           -1.0
 UP SERVINGS  C6           -1.0
 UP SERVINGS  C7           -1.0
 UP SERVINGS  C8           -1.0
RANGES
    RNG       R0            0.1       R0            1.0
    RNG       R0            2.0       R0            3.0
    RNG       R0            4.0       R0            5.0
    RNG       R11           0.1       R12           1.0
    RNG       R13           2.0       R14           3.0
    RNG       R15           4.0       COST          5.0
    RNG       COST          6.0       COST          7.0
    RNG       COST          8.0       COST          9.0
    
ENDATA
