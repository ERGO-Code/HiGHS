NAME          FOLDMAX
OBJSENSE
    MAX
ROWS
 N  COST
 L  R1
 L  R2
 L  R3
 L  S
COLUMNS
    X1        COST      1.0   R1        1.0
    X1        R3        1.0   S         1.0
    X2        COST      1.0   R1        1.0
    X2        R2        1.0   S         1.0
    X3        COST      1.0   R2        1.0
    X3        R3        1.0   S         1.0
    Z         COST      0.5   S         1.0
    W1        COST      0.1   S         1.0
    W2        COST      0.1   S         1.0
RHS
    RHS       R1        4.0   R2        4.0
    RHS       R3        4.0   S         10.0
ENDATA