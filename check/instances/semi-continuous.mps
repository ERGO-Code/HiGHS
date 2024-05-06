NAME lp_solve_sc
OBJSENSE
  MAX
ROWS
 N  r_0
 L  r_1
 G  r_2
 G  r_3
 G  r_4
COLUMNS
    x1        r_0                  1   r_1                  1
    x1        r_2                  2   r_3                 -1
    x2        r_0                  2   r_1                  1
    x2        r_2                 -1   r_3                  3
    MARK0000  'MARKER'                 'INTORG'
    x3        r_0               -0.1   r_4                  1
    MARK0001  'MARKER'                 'INTEND'
    x4        r_0                 -3   r_4                  1
RHS
    RHS       r_1                  5   r_4                0.5
BOUNDS
 SC BND       x3                  10
 LO BND       x3                 1.1
ENDATA
