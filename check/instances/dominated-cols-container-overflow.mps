NAME        dominated_columns_container_overflow
ROWS
 N  Obj
 L  total_disk
 L  total_gcu
 L  table_disk
 L  table_gcu
 E  req2_served
 E  req1_served
 L  avail_2_3
 L  avail_1_3
 L  avail_2_1
 L  avail_1_2
 L  avail_1_1
 L  avail_2_5
 L  avail_1_5
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    t1        Obj                      1
    t1        total_disk               200
    t1        total_gcu                200
    t1        table_disk               200
    t1        table_gcu                200
    t1        avail_2_1                -1
    t1        avail_1_1                -1
    t2        Obj                      1
    t2        total_disk               100
    t2        total_gcu                100
    t2        table_disk               100
    t2        table_gcu                100
    t2        avail_1_2                -1
    t3        Obj                      -1
    t3        total_disk               200
    t3        total_gcu                200
    t3        avail_2_3                -1
    t3        avail_1_3                -1
    t4        Obj                      -1
    t4        total_disk               100
    t4        total_gcu                100
    t5        avail_2_5                -1
    t5        avail_1_5                -1
    d1_1      Obj                      2000
    d1_1      req1_served              1
    d1_1      avail_1_1                1
    d1_2      req1_served              1
    d1_2      avail_1_2                1
    d1_3      Obj                      2000
    d1_3      req1_served              1
    d1_3      avail_1_3                1
    d1_5      Obj                      20000
    d1_5      req1_served              1
    d1_5      avail_1_5                1
    d2_1      req2_served              1
    d2_1      avail_2_1                1
    d2_3      req2_served              1
    d2_3      avail_2_3                1
    d2_5      Obj                      40000
    d2_5      req2_served              1
    d2_5      avail_2_5                1
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS_V     Obj                      -2
    RHS_V     total_disk               300.000001
    RHS_V     total_gcu                300.000001
    RHS_V     table_disk               100.000001
    RHS_V     table_gcu                100.000001
    RHS_V     req2_served              1
    RHS_V     req1_served              1
BOUNDS
 BV BOUND     t1
 BV BOUND     t2
 BV BOUND     t3
 BV BOUND     t4
 BV BOUND     t5
 BV BOUND     d1_1
 BV BOUND     d1_2
 BV BOUND     d1_3
 BV BOUND     d1_5
 BV BOUND     d2_1
 BV BOUND     d2_3
 BV BOUND     d2_5
ENDATA
