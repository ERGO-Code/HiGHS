* Synthetic worst-case execution time (WCET) model in the IPET formulation
* (implicit path enumeration: longest path through a control-flow graph as an ILP).
* Generated from a made-up structured control-flow graph; no third-party data.
*
* Graph: basic blocks b0..b34, a DAG; control enters at b0 and leaves at b25.
*   x_bV   = how often block V executes                (integer >= 0)
*   e_A_B  = how often control passes from A to B      (integer >= 0)
*   in_bV  : x_bV - (sum of e_*_V)  = 0   (= 1 for the entry block)
*   out_bV : x_bV - (sum of e_V_*)  = 0   (not stated for the exit block)
* Objective: maximise the cost of the executed blocks = the most expensive
* entry-to-exit path. Expected: Optimal, objective 25.
*
* Row and column ORDER matter for reproducing the presolve failure: keep them.
NAME equality_row_addition_repro
OBJSENSE
    MAX
ROWS
 N  obj
 E  in_b0
 E  in_b1
 E  in_b2
 E  in_b3
 E  in_b4
 E  in_b5
 E  in_b6
 E  in_b7
 E  in_b8
 E  in_b9
 E  in_b10
 E  in_b11
 E  in_b12
 E  in_b13
 E  in_b14
 E  in_b15
 E  in_b16
 E  in_b17
 E  in_b18
 E  in_b19
 E  in_b20
 E  in_b21
 E  in_b22
 E  in_b23
 E  in_b24
 E  in_b25
 E  in_b26
 E  in_b27
 E  in_b28
 E  in_b29
 E  in_b30
 E  in_b31
 E  in_b32
 E  in_b33
 E  in_b34
 E  out_b0
 E  out_b1
 E  out_b2
 E  out_b3
 E  out_b4
 E  out_b5
 E  out_b6
 E  out_b7
 E  out_b8
 E  out_b9
 E  out_b10
 E  out_b11
 E  out_b12
 E  out_b13
 E  out_b14
 E  out_b15
 E  out_b16
 E  out_b17
 E  out_b18
 E  out_b19
 E  out_b20
 E  out_b21
 E  out_b22
 E  out_b23
 E  out_b24
 E  out_b26
 E  out_b27
 E  out_b28
 E  out_b29
 E  out_b30
 E  out_b31
 E  out_b32
 E  out_b33
 E  out_b34
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x_b0  in_b0  1
    x_b0  out_b0  1
    x_b1  in_b1  1
    x_b1  out_b1  1
    x_b2  in_b2  1
    x_b2  out_b2  1
    x_b3  in_b3  1
    x_b3  out_b3  1
    x_b4  in_b4  1
    x_b4  out_b4  1
    x_b5  obj  11
    x_b5  in_b5  1
    x_b5  out_b5  1
    x_b6  obj  11
    x_b6  in_b6  1
    x_b6  out_b6  1
    x_b7  in_b7  1
    x_b7  out_b7  1
    x_b8  in_b8  1
    x_b8  out_b8  1
    x_b9  in_b9  1
    x_b9  out_b9  1
    x_b10  in_b10  1
    x_b10  out_b10  1
    x_b11  in_b11  1
    x_b11  out_b11  1
    x_b12  in_b12  1
    x_b12  out_b12  1
    x_b13  in_b13  1
    x_b13  out_b13  1
    x_b14  in_b14  1
    x_b14  out_b14  1
    x_b15  in_b15  1
    x_b15  out_b15  1
    x_b16  in_b16  1
    x_b16  out_b16  1
    x_b17  in_b17  1
    x_b17  out_b17  1
    x_b18  in_b18  1
    x_b18  out_b18  1
    x_b19  in_b19  1
    x_b19  out_b19  1
    x_b20  in_b20  1
    x_b20  out_b20  1
    x_b21  obj  1
    x_b21  in_b21  1
    x_b21  out_b21  1
    x_b22  in_b22  1
    x_b22  out_b22  1
    x_b23  obj  1
    x_b23  in_b23  1
    x_b23  out_b23  1
    x_b24  in_b24  1
    x_b24  out_b24  1
    x_b25  in_b25  1
    x_b26  in_b26  1
    x_b26  out_b26  1
    x_b27  in_b27  1
    x_b27  out_b27  1
    x_b28  obj  1
    x_b28  in_b28  1
    x_b28  out_b28  1
    x_b29  in_b29  1
    x_b29  out_b29  1
    x_b30  in_b30  1
    x_b30  out_b30  1
    x_b31  obj  1
    x_b31  in_b31  1
    x_b31  out_b31  1
    x_b32  in_b32  1
    x_b32  out_b32  1
    x_b33  obj  1
    x_b33  in_b33  1
    x_b33  out_b33  1
    x_b34  in_b34  1
    x_b34  out_b34  1
    e_2_3  in_b3  -1
    e_2_3  out_b2  -1
    e_2_4  in_b4  -1
    e_2_4  out_b2  -1
    e_4_3  in_b3  -1
    e_4_3  out_b4  -1
    e_5_6  in_b6  -1
    e_5_6  out_b5  -1
    e_5_7  in_b7  -1
    e_5_7  out_b5  -1
    e_7_6  in_b6  -1
    e_7_6  out_b7  -1
    e_2_5  in_b5  -1
    e_2_5  out_b2  -1
    e_6_3  in_b3  -1
    e_6_3  out_b6  -1
    e_8_9  in_b9  -1
    e_8_9  out_b8  -1
    e_8_10  in_b10  -1
    e_8_10  out_b8  -1
    e_10_9  in_b9  -1
    e_10_9  out_b10  -1
    e_8_11  in_b11  -1
    e_8_11  out_b8  -1
    e_11_9  in_b9  -1
    e_11_9  out_b11  -1
    e_8_12  in_b12  -1
    e_8_12  out_b8  -1
    e_12_9  in_b9  -1
    e_12_9  out_b12  -1
    e_2_8  in_b8  -1
    e_2_8  out_b2  -1
    e_9_3  in_b3  -1
    e_9_3  out_b9  -1
    e_0_2  in_b2  -1
    e_0_2  out_b0  -1
    e_3_1  in_b1  -1
    e_3_1  out_b3  -1
    e_0_1  in_b1  -1
    e_0_1  out_b0  -1
    e_13_15  in_b15  -1
    e_13_15  out_b13  -1
    e_15_14  in_b14  -1
    e_15_14  out_b15  -1
    e_13_14  in_b14  -1
    e_13_14  out_b13  -1
    e_16_18  in_b18  -1
    e_16_18  out_b16  -1
    e_18_17  in_b17  -1
    e_18_17  out_b18  -1
    e_19_21  in_b21  -1
    e_19_21  out_b19  -1
    e_21_20  in_b20  -1
    e_21_20  out_b21  -1
    e_19_22  in_b22  -1
    e_19_22  out_b19  -1
    e_22_20  in_b20  -1
    e_22_20  out_b22  -1
    e_19_23  in_b23  -1
    e_19_23  out_b19  -1
    e_23_20  in_b20  -1
    e_23_20  out_b23  -1
    e_19_20  in_b20  -1
    e_19_20  out_b19  -1
    e_16_19  in_b19  -1
    e_16_19  out_b16  -1
    e_20_17  in_b17  -1
    e_20_17  out_b20  -1
    e_14_16  in_b16  -1
    e_14_16  out_b14  -1
    e_1_13  in_b13  -1
    e_1_13  out_b1  -1
    e_24_25  in_b25  -1
    e_24_25  out_b24  -1
    e_26_28  in_b28  -1
    e_26_28  out_b26  -1
    e_28_27  in_b27  -1
    e_28_27  out_b28  -1
    e_26_27  in_b27  -1
    e_26_27  out_b26  -1
    e_30_31  in_b31  -1
    e_30_31  out_b30  -1
    e_30_32  in_b32  -1
    e_30_32  out_b30  -1
    e_32_31  in_b31  -1
    e_32_31  out_b32  -1
    e_30_33  in_b33  -1
    e_30_33  out_b30  -1
    e_33_31  in_b31  -1
    e_33_31  out_b33  -1
    e_30_34  in_b34  -1
    e_30_34  out_b30  -1
    e_34_31  in_b31  -1
    e_34_31  out_b34  -1
    e_29_30  in_b30  -1
    e_29_30  out_b29  -1
    e_26_29  in_b29  -1
    e_26_29  out_b26  -1
    e_31_27  in_b27  -1
    e_31_27  out_b31  -1
    e_24_26  in_b26  -1
    e_24_26  out_b24  -1
    e_27_25  in_b25  -1
    e_27_25  out_b27  -1
    e_17_24  in_b24  -1
    e_17_24  out_b17  -1
    MARKER  'MARKER'  'INTEND'
RHS
    RHS  in_b0  1
BOUNDS
 PL BND  x_b0
 PL BND  x_b1
 PL BND  x_b2
 PL BND  x_b3
 PL BND  x_b4
 PL BND  x_b5
 PL BND  x_b6
 PL BND  x_b7
 PL BND  x_b8
 PL BND  x_b9
 PL BND  x_b10
 PL BND  x_b11
 PL BND  x_b12
 PL BND  x_b13
 PL BND  x_b14
 PL BND  x_b15
 PL BND  x_b16
 PL BND  x_b17
 PL BND  x_b18
 PL BND  x_b19
 PL BND  x_b20
 PL BND  x_b21
 PL BND  x_b22
 PL BND  x_b23
 PL BND  x_b24
 PL BND  x_b25
 PL BND  x_b26
 PL BND  x_b27
 PL BND  x_b28
 PL BND  x_b29
 PL BND  x_b30
 PL BND  x_b31
 PL BND  x_b32
 PL BND  x_b33
 PL BND  x_b34
 PL BND  e_2_3
 PL BND  e_2_4
 PL BND  e_4_3
 PL BND  e_5_6
 PL BND  e_5_7
 PL BND  e_7_6
 PL BND  e_2_5
 PL BND  e_6_3
 PL BND  e_8_9
 PL BND  e_8_10
 PL BND  e_10_9
 PL BND  e_8_11
 PL BND  e_11_9
 PL BND  e_8_12
 PL BND  e_12_9
 PL BND  e_2_8
 PL BND  e_9_3
 PL BND  e_0_2
 PL BND  e_3_1
 PL BND  e_0_1
 PL BND  e_13_15
 PL BND  e_15_14
 PL BND  e_13_14
 PL BND  e_16_18
 PL BND  e_18_17
 PL BND  e_19_21
 PL BND  e_21_20
 PL BND  e_19_22
 PL BND  e_22_20
 PL BND  e_19_23
 PL BND  e_23_20
 PL BND  e_19_20
 PL BND  e_16_19
 PL BND  e_20_17
 PL BND  e_14_16
 PL BND  e_1_13
 PL BND  e_24_25
 PL BND  e_26_28
 PL BND  e_28_27
 PL BND  e_26_27
 PL BND  e_30_31
 PL BND  e_30_32
 PL BND  e_32_31
 PL BND  e_30_33
 PL BND  e_33_31
 PL BND  e_30_34
 PL BND  e_34_31
 PL BND  e_29_30
 PL BND  e_26_29
 PL BND  e_31_27
 PL BND  e_24_26
 PL BND  e_27_25
 PL BND  e_17_24
ENDATA
