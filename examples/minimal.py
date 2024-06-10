import highspy

h = highspy.Highs()

x1 = h.addVariable(lb = -h.inf)
x2 = h.addVariable(lb = -h.inf)

h.addConstr(x2 - x1 >= 2)
h.addConstr(x1 + x2 >= 0)

h.minimize(x2)
