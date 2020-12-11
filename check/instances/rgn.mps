*NAME:         rgn
*ROWS:         24
*COLUMNS:      180
*INTEGER:      100
*NONZERO:      460
*BEST SOLN:    82.1999 (opt)
*LP SOLN:      48.7999
*SOURCE:       Linus E. Schrage (U. Chicago) 
*              Laurence A. Wolsey (University of Louvain)
*              Martin W. P. Savelsbergh (Eindhoven Univ. of Technology)
*APPLICATION:  unknown
*COMMENTS:     all integer variables are binary
*       
*       
NAME          RGN
ROWS
 N  1       
 L  2       
 L  3       
 L  4       
 L  5       
 E  6       
 E  7       
 E  8       
 E  9       
 E  10      
 E  11      
 E  12      
 E  13      
 E  14      
 E  15      
 E  16      
 E  17      
 E  18      
 E  19      
 E  20      
 E  21      
 E  22      
 E  23      
 E  24      
 E  25      
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    A1        2                    1   6          -4.60000002
    B1        2                    1   7          -4.60000002
    C1        2                    1   8          -4.60000002
    D1        2                    1   9          -4.60000002
    E1        2                    1   10         -4.60000002
    AB1       2                    1   6          -4.60000002
    AB1       7          -4.60000002
    AC1       2                    1   6          -4.60000002
    AC1       8          -4.60000002
    AD1       2                    1   6          -4.60000002
    AD1       9          -4.60000002
    AE1       2                    1   6          -4.60000002
    AE1       10         -4.60000002
    BC1       2                    1   7          -4.60000002
    BC1       8          -4.60000002
    BD1       2                    1   7          -4.60000002
    BD1       9          -4.60000002
    BE1       2                    1   7          -4.60000002
    BE1       10         -4.60000002
    CD1       2                    1   8          -4.60000002
    CD1       9          -4.60000002
    CE1       2                    1   8          -4.60000002
    CE1       10         -4.60000002
    DE1       2                    1   9          -4.60000002
    DE1       10         -4.60000002
    ABC1      2                    1   6          -4.60000002
    ABC1      7          -4.60000002   8          -4.60000002
    ABD1      2                    1   6          -4.60000002
    ABD1      7          -4.60000002   9          -4.60000002
    ABE1      2                    1   6          -4.60000002
    ABE1      7          -4.60000002   10         -4.60000002
    ACD1      2                    1   6          -4.60000002
    ACD1      8          -4.60000002   9          -4.60000002
    ACE1      2                    1   6          -4.60000002
    ACE1      8          -4.60000002   10         -4.60000002
    ADE1      2                    1   6          -4.60000002
    ADE1      9          -4.60000002   10         -4.60000002
    BCD1      2                    1   7          -4.60000002
    BCD1      8          -4.60000002   9          -4.60000002
    BCE1      2                    1   7          -4.60000002
    BCE1      8          -4.60000002   10         -4.60000002
    BDE1      2                    1   7          -4.60000002
    BDE1      9          -4.60000002   10         -4.60000002
    CDE1      2                    1   8          -4.60000002
    CDE1      9          -4.60000002   10         -4.60000002
    A2        3                    1   11         -4.60000002
    B2        3                    1   12         -4.60000002
    C2        3                    1   13         -4.60000002
    D2        3                    1   14         -4.60000002
    E2        3                    1   15         -4.60000002
    AB2       3                    1   11         -4.60000002
    AB2       12         -4.60000002
    AC2       3                    1   11         -4.60000002
    AC2       13         -4.60000002
    AD2       3                    1   11         -4.60000002
    AD2       14         -4.60000002
    AE2       3                    1   11         -4.60000002
    AE2       15         -4.60000002
    BC2       3                    1   12         -4.60000002
    BC2       13         -4.60000002
    BD2       3                    1   12         -4.60000002
    BD2       14         -4.60000002
    BE2       3                    1   12         -4.60000002
    BE2       15         -4.60000002
    CD2       3                    1   13         -4.60000002
    CD2       14         -4.60000002
    CE2       3                    1   13         -4.60000002
    CE2       15         -4.60000002
    DE2       3                    1   14         -4.60000002
    DE2       15         -4.60000002
    ABC2      3                    1   11         -4.60000002
    ABC2      12         -4.60000002   13         -4.60000002
    ABD2      3                    1   11         -4.60000002
    ABD2      12         -4.60000002   14         -4.60000002
    ABE2      3                    1   11         -4.60000002
    ABE2      12         -4.60000002   15         -4.60000002
    ACD2      3                    1   11         -4.60000002
    ACD2      13         -4.60000002   14         -4.60000002
    ACE2      3                    1   11         -4.60000002
    ACE2      13         -4.60000002   15         -4.60000002
    ADE2      3                    1   11         -4.60000002
    ADE2      14         -4.60000002   15         -4.60000002
    BCD2      3                    1   12         -4.60000002
    BCD2      13         -4.60000002   14         -4.60000002
    BCE2      3                    1   12         -4.60000002
    BCE2      13         -4.60000002   15         -4.60000002
    BDE2      3                    1   12         -4.60000002
    BDE2      14         -4.60000002   15         -4.60000002
    CDE2      3                    1   13         -4.60000002
    CDE2      14         -4.60000002   15         -4.60000002
    A3        4                    1   16         -4.60000002
    B3        4                    1   17         -4.60000002
    C3        4                    1   18         -4.60000002
    D3        4                    1   19         -4.60000002
    E3        4                    1   20         -4.60000002
    AB3       4                    1   16         -4.60000002
    AB3       17         -4.60000002
    AC3       4                    1   16         -4.60000002
    AC3       18         -4.60000002
    AD3       4                    1   16         -4.60000002
    AD3       19         -4.60000002
    AE3       4                    1   16         -4.60000002
    AE3       20         -4.60000002
    BC3       4                    1   17         -4.60000002
    BC3       18         -4.60000002
    BD3       4                    1   17         -4.60000002
    BD3       19         -4.60000002
    BE3       4                    1   17         -4.60000002
    BE3       20         -4.60000002
    CD3       4                    1   18         -4.60000002
    CD3       19         -4.60000002
    CE3       4                    1   18         -4.60000002
    CE3       20         -4.60000002
    DE3       4                    1   19         -4.60000002
    DE3       20         -4.60000002
    ABC3      4                    1   16         -4.60000002
    ABC3      17         -4.60000002   18         -4.60000002
    ABD3      4                    1   16         -4.60000002
    ABD3      17         -4.60000002   19         -4.60000002
    ABE3      4                    1   16         -4.60000002
    ABE3      17         -4.60000002   20         -4.60000002
    ACD3      4                    1   16         -4.60000002
    ACD3      18         -4.60000002   19         -4.60000002
    ACE3      4                    1   16         -4.60000002
    ACE3      18         -4.60000002   20         -4.60000002
    ADE3      4                    1   16         -4.60000002
    ADE3      19         -4.60000002   20         -4.60000002
    BCD3      4                    1   17         -4.60000002
    BCD3      18         -4.60000002   19         -4.60000002
    BCE3      4                    1   17         -4.60000002
    BCE3      18         -4.60000002   20         -4.60000002
    BDE3      4                    1   17         -4.60000002
    BDE3      19         -4.60000002   20         -4.60000002
    CDE3      4                    1   18         -4.60000002
    CDE3      19         -4.60000002   20         -4.60000002
    A4        5                    1   21         -4.60000002
    B4        5                    1   22         -4.60000002
    C4        5                    1   23         -4.60000002
    D4        5                    1   24         -4.60000002
    E4        5                    1   25         -4.60000002
    AB4       5                    1   21         -4.60000002
    AB4       22         -4.60000002
    AC4       5                    1   21         -4.60000002
    AC4       23         -4.60000002
    AD4       5                    1   21         -4.60000002
    AD4       24         -4.60000002
    AE4       5                    1   21         -4.60000002
    AE4       25         -4.60000002
    BC4       5                    1   22         -4.60000002
    BC4       23         -4.60000002
    BD4       5                    1   22         -4.60000002
    BD4       24         -4.60000002
    BE4       5                    1   22         -4.60000002
    BE4       25         -4.60000002
    CD4       5                    1   23         -4.60000002
    CD4       24         -4.60000002
    CE4       5                    1   23         -4.60000002
    CE4       25         -4.60000002
    DE4       5                    1   24         -4.60000002
    DE4       25         -4.60000002
    ABC4      5                    1   21         -4.60000002
    ABC4      22         -4.60000002   23         -4.60000002
    ABD4      5                    1   21         -4.60000002
    ABD4      22         -4.60000002   24         -4.60000002
    ABE4      5                    1   21         -4.60000002
    ABE4      22         -4.60000002   25         -4.60000002
    ACD4      5                    1   21         -4.60000002
    ACD4      23         -4.60000002   24         -4.60000002
    ACE4      5                    1   21         -4.60000002
    ACE4      23         -4.60000002   25         -4.60000002
    ADE4      5                    1   21         -4.60000002
    ADE4      24         -4.60000002   25         -4.60000002
    BCD4      5                    1   22         -4.60000002
    BCD4      23         -4.60000002   24         -4.60000002
    BCE4      5                    1   22         -4.60000002
    BCE4      23         -4.60000002   25         -4.60000002
    BDE4      5                    1   22         -4.60000002
    BDE4      24         -4.60000002   25         -4.60000002
    CDE4      5                    1   23         -4.60000002
    CDE4      24         -4.60000002   25         -4.60000002
    MARK0001  'MARKER'                 'INTEND'
    TA1       1                    3   6                    1
    TA1       11                  -1
    TA2       1                    3   11                   1
    TA2       16                  -1
    TA3       1                    3   16                   1
    TA3       21                  -1
    TA4       1                    3   21                   1
    TB1       1                    3   7                    1
    TB1       12                  -1
    TB2       1                    3   12                   1
    TB2       17                  -1
    TB3       1                    3   17                   1
    TB3       22                  -1
    TB4       1                    3   22                   1
    TC1       1                    3   8                    1
    TC1       13                  -1
    TC2       1                    3   13                   1
    TC2       18                  -1
    TC3       1                    3   18                   1
    TC3       23                  -1
    TC4       1                    3   23                   1
    TD1       1                    3   9                    1
    TD1       14                  -1
    TD2       1                    3   14                   1
    TD2       19                  -1
    TD3       1                    3   19                   1
    TD3       24                  -1
    TD4       1                    3   24                   1
    TE1       1                    3   10                   1
    TE1       15                  -1
    TE2       1                    3   15                   1
    TE2       20                  -1
    TE3       1                    3   20                   1
    TE3       25                  -1
    TE4       1                    3   25                   1
    UA1       1                    1   6                    1
    UA1       11                  -1
    UA2       1                    1   11                   1
    UA2       16                  -1
    UA3       1                    1   16                   1
    UA3       21                  -1
    UA4       1                    1   21                   1
    UB1       1                    1   7                    1
    UB1       12                  -1
    UB2       1                    1   12                   1
    UB2       17                  -1
    UB3       1                    1   17                   1
    UB3       22                  -1
    UB4       1                    1   22                   1
    UC1       1                    1   8                    1
    UC1       13                  -1
    UC2       1                    1   13                   1
    UC2       18                  -1
    UC3       1                    1   18                   1
    UC3       23                  -1
    UC4       1                    1   23                   1
    UD1       1                    1   9                    1
    UD1       14                  -1
    UD2       1                    1   14                   1
    UD2       19                  -1
    UD3       1                    1   19                   1
    UD3       24                  -1
    UD4       1                    1   24                   1
    UE1       1                    1   10                   1
    UE1       15                  -1
    UE2       1                    1   15                   1
    UE2       20                  -1
    UE3       1                    1   20                   1
    UE3       25                  -1
    UE4       1                    1   25                   1
    VA1       1                    1   6                   -1
    VA1       11                   1
    VA2       1                    1   11                  -1
    VA2       16                   1
    VA3       1                    1   16                  -1
    VA3       21                   1
    VA4       1                    1   21                  -1
    VB1       1                    1   7                   -1
    VB1       12                   1
    VB2       1                    1   12                  -1
    VB2       17                   1
    VB3       1                    1   17                  -1
    VB3       22                   1
    VB4       1                    1   22                  -1
    VC1       1                    1   8                   -1
    VC1       13                   1
    VC2       1                    1   13                  -1
    VC2       18                   1
    VC3       1                    1   18                  -1
    VC3       23                   1
    VC4       1                    1   23                  -1
    VD1       1                    1   9                   -1
    VD1       14                   1
    VD2       1                    1   14                  -1
    VD2       19                   1
    VD3       1                    1   19                  -1
    VD3       24                   1
    VD4       1                    1   24                  -1
    VE1       1                    1   10                  -1
    VE1       15                   1
    VE2       1                    1   15                  -1
    VE2       20                   1
    VE3       1                    1   20                  -1
    VE3       25                   1
    VE4       1                    1   25                  -1
    WA1       1                    3   6                   -1
    WA1       11                   1
    WA2       1                    3   11                  -1
    WA2       16                   1
    WA3       1                    3   16                  -1
    WA3       21                   1
    WA4       1                    3   21                  -1
    WB1       1                    3   7                   -1
    WB1       12                   1
    WB2       1                    3   12                  -1
    WB2       17                   1
    WB3       1                    3   17                  -1
    WB3       22                   1
    WB4       1                    3   22                  -1
    WC1       1                    3   8                   -1
    WC1       13                   1
    WC2       1                    3   13                  -1
    WC2       18                   1
    WC3       1                    3   18                  -1
    WC3       23                   1
    WC4       1                    3   23                  -1
    WD1       1                    3   9                   -1
    WD1       14                   1
    WD2       1                    3   14                  -1
    WD2       19                   1
    WD3       1                    3   19                  -1
    WD3       24                   1
    WD4       1                    3   24                  -1
    WE1       1                    3   10                  -1
    WE1       15                   1
    WE2       1                    3   15                  -1
    WE2       20                   1
    WE3       1                    3   20                  -1
    WE3       25                   1
    WE4       1                    3   25                  -1
RHS
    RHS       2                    1   3                    1
    RHS       4                    1   5                    1
    RHS       6                 -3.5   7                 -3.5
    RHS       8                 -3.5   9                 -3.5
    RHS       10                -3.5   11                -3.5
    RHS       12                -3.5   13                -3.5
    RHS       14                -3.5   15                -3.5
    RHS       16                -3.5   17                -3.5
    RHS       18                -3.5   19                -3.5
    RHS       20                -3.5   21                -3.5
    RHS       22                -3.5   23                -3.5
    RHS       24                -3.5   25                -3.5
BOUNDS
 UP LINDOBND  A1                   1
 UP LINDOBND  B1                   1
 UP LINDOBND  C1                   1
 UP LINDOBND  D1                   1
 UP LINDOBND  E1                   1
 UP LINDOBND  AB1                  1
 UP LINDOBND  AC1                  1
 UP LINDOBND  AD1                  1
 UP LINDOBND  AE1                  1
 UP LINDOBND  BC1                  1
 UP LINDOBND  BD1                  1
 UP LINDOBND  BE1                  1
 UP LINDOBND  CD1                  1
 UP LINDOBND  CE1                  1
 UP LINDOBND  DE1                  1
 UP LINDOBND  ABC1                 1
 UP LINDOBND  ABD1                 1
 UP LINDOBND  ABE1                 1
 UP LINDOBND  ACD1                 1
 UP LINDOBND  ACE1                 1
 UP LINDOBND  ADE1                 1
 UP LINDOBND  BCD1                 1
 UP LINDOBND  BCE1                 1
 UP LINDOBND  BDE1                 1
 UP LINDOBND  CDE1                 1
 UP LINDOBND  A2                   1
 UP LINDOBND  B2                   1
 UP LINDOBND  C2                   1
 UP LINDOBND  D2                   1
 UP LINDOBND  E2                   1
 UP LINDOBND  AB2                  1
 UP LINDOBND  AC2                  1
 UP LINDOBND  AD2                  1
 UP LINDOBND  AE2                  1
 UP LINDOBND  BC2                  1
 UP LINDOBND  BD2                  1
 UP LINDOBND  BE2                  1
 UP LINDOBND  CD2                  1
 UP LINDOBND  CE2                  1
 UP LINDOBND  DE2                  1
 UP LINDOBND  ABC2                 1
 UP LINDOBND  ABD2                 1
 UP LINDOBND  ABE2                 1
 UP LINDOBND  ACD2                 1
 UP LINDOBND  ACE2                 1
 UP LINDOBND  ADE2                 1
 UP LINDOBND  BCD2                 1
 UP LINDOBND  BCE2                 1
 UP LINDOBND  BDE2                 1
 UP LINDOBND  CDE2                 1
 UP LINDOBND  A3                   1
 UP LINDOBND  B3                   1
 UP LINDOBND  C3                   1
 UP LINDOBND  D3                   1
 UP LINDOBND  E3                   1
 UP LINDOBND  AB3                  1
 UP LINDOBND  AC3                  1
 UP LINDOBND  AD3                  1
 UP LINDOBND  AE3                  1
 UP LINDOBND  BC3                  1
 UP LINDOBND  BD3                  1
 UP LINDOBND  BE3                  1
 UP LINDOBND  CD3                  1
 UP LINDOBND  CE3                  1
 UP LINDOBND  DE3                  1
 UP LINDOBND  ABC3                 1
 UP LINDOBND  ABD3                 1
 UP LINDOBND  ABE3                 1
 UP LINDOBND  ACD3                 1
 UP LINDOBND  ACE3                 1
 UP LINDOBND  ADE3                 1
 UP LINDOBND  BCD3                 1
 UP LINDOBND  BCE3                 1
 UP LINDOBND  BDE3                 1
 UP LINDOBND  CDE3                 1
 UP LINDOBND  A4                   1
 UP LINDOBND  B4                   1
 UP LINDOBND  C4                   1
 UP LINDOBND  D4                   1
 UP LINDOBND  E4                   1
 UP LINDOBND  AB4                  1
 UP LINDOBND  AC4                  1
 UP LINDOBND  AD4                  1
 UP LINDOBND  AE4                  1
 UP LINDOBND  BC4                  1
 UP LINDOBND  BD4                  1
 UP LINDOBND  BE4                  1
 UP LINDOBND  CD4                  1
 UP LINDOBND  CE4                  1
 UP LINDOBND  DE4                  1
 UP LINDOBND  ABC4                 1
 UP LINDOBND  ABD4                 1
 UP LINDOBND  ABE4                 1
 UP LINDOBND  ACD4                 1
 UP LINDOBND  ACE4                 1
 UP LINDOBND  ADE4                 1
 UP LINDOBND  BCD4                 1
 UP LINDOBND  BCE4                 1
 UP LINDOBND  BDE4                 1
 UP LINDOBND  CDE4                 1
 UP LINDOBND  TA1                100
 UP LINDOBND  TA2                100
 UP LINDOBND  TA3                100
 UP LINDOBND  TA4                100
 UP LINDOBND  TB1                100
 UP LINDOBND  TB2                100
 UP LINDOBND  TB3                100
 UP LINDOBND  TB4                100
 UP LINDOBND  TC1                100
 UP LINDOBND  TC2                100
 UP LINDOBND  TC3                100
 UP LINDOBND  TC4                100
 UP LINDOBND  TD1                100
 UP LINDOBND  TD2                100
 UP LINDOBND  TD3                100
 UP LINDOBND  TD4                100
 UP LINDOBND  TE1                100
 UP LINDOBND  TE2                100
 UP LINDOBND  TE3                100
 UP LINDOBND  TE4                100
 UP LINDOBND  UA1                  2
 UP LINDOBND  UA2                  2
 UP LINDOBND  UA3                  2
 UP LINDOBND  UA4                  2
 UP LINDOBND  UB1                  2
 UP LINDOBND  UB2                  2
 UP LINDOBND  UB3                  2
 UP LINDOBND  UB4                  2
 UP LINDOBND  UC1                  2
 UP LINDOBND  UC2                  2
 UP LINDOBND  UC3                  2
 UP LINDOBND  UC4                  2
 UP LINDOBND  UD1                  2
 UP LINDOBND  UD2                  2
 UP LINDOBND  UD3                  2
 UP LINDOBND  UD4                  2
 UP LINDOBND  UE1                  2
 UP LINDOBND  UE2                  2
 UP LINDOBND  UE3                  2
 UP LINDOBND  UE4                  2
 UP LINDOBND  VA1                  2
 UP LINDOBND  VA2                  2
 UP LINDOBND  VA3                  2
 UP LINDOBND  VA4                  2
 UP LINDOBND  VB1                  2
 UP LINDOBND  VB2                  2
 UP LINDOBND  VB3                  2
 UP LINDOBND  VB4                  2
 UP LINDOBND  VC1                  2
 UP LINDOBND  VC2                  2
 UP LINDOBND  VC3                  2
 UP LINDOBND  VC4                  2
 UP LINDOBND  VD1                  2
 UP LINDOBND  VD2                  2
 UP LINDOBND  VD3                  2
 UP LINDOBND  VD4                  2
 UP LINDOBND  VE1                  2
 UP LINDOBND  VE2                  2
 UP LINDOBND  VE3                  2
 UP LINDOBND  VE4                  2
 UP LINDOBND  WA1                100
 UP LINDOBND  WA2                100
 UP LINDOBND  WA3                100
 UP LINDOBND  WA4                100
 UP LINDOBND  WB1                100
 UP LINDOBND  WB2                100
 UP LINDOBND  WB3                100
 UP LINDOBND  WB4                100
 UP LINDOBND  WC1                100
 UP LINDOBND  WC2                100
 UP LINDOBND  WC3                100
 UP LINDOBND  WC4                100
 UP LINDOBND  WD1                100
 UP LINDOBND  WD2                100
 UP LINDOBND  WD3                100
 UP LINDOBND  WD4                100
 UP LINDOBND  WE1                100
 UP LINDOBND  WE2                100
 UP LINDOBND  WE3                100
 UP LINDOBND  WE4                100
ENDATA
