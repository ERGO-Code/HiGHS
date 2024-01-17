# This example cannot be run with the version (1.5.3) of highspy
# available from PyPI. It requires a local installation of highspy for
# (at least) HiGHS version 1.6.0
import highspy

h = highspy.Highs()

x1 = h.addVar(lb=-h.inf)
x2 = h.addVar(lb=-h.inf)

h.addConstr(x2 - x1 >= 2)
h.addConstr(x1 + x2 >= 0)

h.minimize(x2)
