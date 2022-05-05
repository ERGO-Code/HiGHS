file = "../src/interfaces/highs_python_api.py"
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
col_cost = (8.0, 10.0)
col_lower = (0.0, 0.0)
col_upper = (inf, inf)
row_lower = (7.0, 12.0, 6.0)
row_upper = (inf, inf, inf)
a_start = (0, 3)
a_index = (0, 1, 2, 0, 1, 2)
a_value = (2.0, 3.0, 2.0, 2.0, 4.0, 1.0)

# Find the continuous solution of the LP
return_status, model_status, col_value, col_dual, row_value, row_dual , col_basis, row_basis = Highs_lpCall(
    col_cost, col_lower, col_upper,
    row_lower, row_upper,
    a_start, a_index, a_value)

print ("return_status = ", return_status)
print ("model_status = ", model_status)
print (col_value, col_dual, row_value, row_dual, col_basis, row_basis)

# Now find the integer solution of the LP
# Add integrality for the two variables and call the MIP solver
integrality = (1, 1)
return_status, model_status, col_value, row_value = Highs_mipCall(
    col_cost, col_lower, col_upper,
    row_lower, row_upper,
    a_start, a_index, a_value,
    integrality)

print ("return_status = ", return_status)
print ("model_status = ", model_status)
print (col_value, row_value)

# Illustrate the solution of a QP
#
# minimize -x_2 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
#
# subject to x_1 + x_2 + x_3 >= 1; x>=0
# Add a spurious second row so len(row_lower) is OK

col_cost = (0.0, -1.0, 0.0)
col_lower = (0.0, 0.0, 0.0)
col_upper = (inf, inf, inf)
row_lower = (1.0, -inf)
row_upper = (inf, inf)
a_start = (0, 1, 2)
a_index = (0, 0, 0)
a_value = (1.0, 1.0, 1.0)
q_start = (0, 2, 3)
q_index = (0, 2, 1, 0, 2)
q_value = (2.0, -1.0, 0.2, -1.0, 2.0)
return_status, model_status, col_value, col_dual, row_value, row_dual , col_basis, row_basis = Highs_qpCall(
    col_cost, col_lower, col_upper,
    row_lower, row_upper,
    a_start, a_index, a_value,
    q_start, q_index, q_value)

print ("return_status = ", return_status)
print ("model_status = ", model_status)
print (col_value, col_dual, row_value, row_dual, col_basis, row_basis)
