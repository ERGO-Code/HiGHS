NAME          indicator1
ROWS
 N  obj
 G  con1
COLUMNS
    INT1      'MARKER'                 'INTORG'
    z         obj            0.0
    INT1END   'MARKER'                 'INTEND'
    x         obj            1.0        con1           1.0
RHS
    rhs       con1           5.0
BOUNDS
 UP bnd       x              10.0
 UP bnd       z              1.0
INDICATORS
 IF con1 z 1
ENDATA
