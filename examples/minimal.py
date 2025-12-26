import highspy

h = highspy.Highs()

x1 = h.addVariable(lb = -h.inf)
x2 = h.addVariable(lb = -h.inf)

h.addConstrs(x2 - x1 >= 2, 
             x1 + x2 >= 0)

h.minimize(x2)
