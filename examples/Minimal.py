import highspy
import numpy as np
inf = highspy.kHighsInf
h = highspy.Highs()
h.addVars(2, np.array([-inf, -inf]), np.array([inf, inf]))
h.changeColsCost(2, np.array([0, 1]), np.array([0, 1], dtype=np.double))
num_cons = 2
lower = np.array([2, 0], dtype=np.double)
upper = np.array([inf, inf], dtype=np.double)
num_new_nz = 4
starts = np.array([0, 2])
indices = np.array([0, 1, 0, 1])
values = np.array([-1, 1, 1, 1], dtype=np.double)
h.addRows(num_cons, lower, upper, num_new_nz, starts, indices, values)
h.setOptionValue('log_to_console', True)
h.run()
