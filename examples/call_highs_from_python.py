file = "../HiGHS/src/interfaces/highs_lp_solver.py"
exec(compile(open(file).read(), file, 'exec'))

# Problem is
#
# min 8 x_1 + 10 x_2
#
# st  2 x_1 +  2 x_2 >= 7
#     3 x_1 +  4 x_2 >= 12
#     2 x_1 +    x_2 >= 6
#
#       x_1 >= 0 x_2 >= 0
#
inf = 1e30
colcost = (8.0, 10.0)
collower = (0.0, 0.0)
colupper = (inf, inf)
rowlower = (7.0, 12.0, 6.0)
rowupper = (inf, inf, inf)
astart = (0, 3)
aindex = (0, 1, 2, 0, 1, 2)
avalue = (2.0, 3.0, 2.0, 2.0, 4.0, 1.0)

# Find the continuous solution of the LP
return_code, col_value, col_dual, row_value, row_dual , col_basis, row_basis = Highs_lpCall(colcost, collower, colupper,
                                                                                            rowlower, rowupper,
                                                                                            astart, aindex, avalue)

print (return_code, col_value, col_dual, row_value, row_dual, col_basis, row_basis)

file = "../HiGHS/src/interfaces/highs_mip_solver.py"
exec(compile(open(file).read(), file, 'exec'))

# Now find the integer solution of the LP
# Add integrality for the two variables and call the MIP solver
integrality = (1, 1)
return_code, col_value, row_value = Highs_mipCall(colcost, collower, colupper,
                                                  rowlower, rowupper,
                                                  astart, aindex, avalue,
                                                  integrality)

print (return_code, col_value, row_value)
