# model and solve the LP
#
# maximize 10 x1 + 25x2
# s. t.       x1 +  2x2 <= 80
#             x1 +  4x2 <= 120
#             x1 >= 0; x2 >= 0
import highspy

h = highspy.Highs()

x1 = h.addVar(obj = 10)
x2 = h.addVar(obj = 25)

h.addConstr(x1 + 2*x2 <= 80)
h.addConstr(x1 + 4*x2 <= 120)

h.writeModel("")
print("x1 = ", h.val(x1))
print("x2 = ", h.val(x2))

h.maximize()

print("x1 = ", h.val(x1))
print("x2 = ", h.val(x2))

print("x = ", h.vals([x1, x2]))



