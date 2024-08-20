# This is an example of the N-Queens problem, which is a classic combinatorial problem.
# The problem is to place N queens on an N x N chessboard so that no two queens attack each other.
#
# We show how to model the problem as a MIP and solve it using highspy.
# Using numpy can simplify the construction of the constraints (i.e., diagonal).

import highspy
import numpy as np

N = 8
h = highspy.Highs()
h.silent()

x = np.reshape(h.addBinaries(N*N), (N, N))

h.addConstrs(sum(x[i,:]) == 1 for i in range(N))    # each row has exactly one queen
h.addConstrs(sum(x[:,j]) == 1 for j in range(N))    # each col has exactly one queen

y = np.fliplr(x)
h.addConstrs(x.diagonal(k).sum() <= 1 for k in range(-N + 1, N))   # each diagonal has at most one queen
h.addConstrs(y.diagonal(k).sum() <= 1 for k in range(-N + 1, N))   # each 'reverse' diagonal has at most one queen

h.solve()
sol = np.array(h.vals(x))

print('Queens:')

for i in range(N):
    print(''.join('Q' if sol[i, j] > 0.5 else '*' for j in range(N)))