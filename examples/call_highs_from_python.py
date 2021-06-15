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

# Illustrate the solution of a QP
#
# minimize -x_2 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
#
# subject to x_1 + x_2 + x_3 >= 1; x>=0
# Add a spurious second row so len(rowlower) is OK

file = "../HiGHS/src/interfaces/highs_qp_solver.py"
exec(compile(open(file).read(), file, 'exec'))

colcost = (0.0, -1.0, 0.0)
collower = (0.0, 0.0, 0.0)
colupper = (inf, inf, inf)
rowlower = (1.0, -inf)
rowupper = (inf, inf)
astart = (0, 1, 2)
aindex = (0, 0, 0)
avalue = (1.0, 1.0, 1.0)
qstart = (0, 2, 3)
qindex = (0, 2, 1, 0, 2)
qvalue = (2.0, -1.0, 0.2, -1.0, 2.0)
return_code, col_value, col_dual, row_value, row_dual , col_basis, row_basis = Highs_qpCall(colcost, collower, colupper,
                                                 rowlower, rowupper,
                                                 astart, aindex, avalue,
                                                 qstart, qindex, qvalue)

print (return_code, col_value, col_dual, row_value, row_dual, col_basis, row_basis)
