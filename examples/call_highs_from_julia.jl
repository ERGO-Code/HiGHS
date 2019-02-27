include("../src/interfaces/highs_lp_solver.jl")

cc = cc = [1, -2]
cl = [0,0]
cu = [10, 10]
ru = [2,1]
rl = [0,0]
astart = [0,2,4]
aindex = [0,1,0,1]
avalue = [1.0,2.0,1.0,3.0)]
call_highs(cc, cl, cu, rl, ru, astart, aindex, avalue)
