# An example of solving the Generalized Assignment Problem (GAP) using highspy
# Also demonstrates how to use callbacks to print cuts as they are found
from highspy import *
import numpy as np

#
# GAP instances can be taken from:
# http://people.brunel.ac.uk/~mastjjb/jeb/orlib/gapinfo.html
#
# Expected format:
# - number machines
# _ number jobs
# - cost of each job on each machine
# - size of each job on each machine
# - capacity of each machine

input_data = ''' 5 15
17 21 22 18 24 15 20 18 19 18 16 22 24 24 16
23 16 21 16 17 16 19 25 18 21 17 15 25 17 24
16 20 16 25 24 16 17 19 19 18 20 16 17 21 24
19 19 22 22 20 16 19 17 21 19 25 23 25 25 25
18 19 15 15 21 25 16 16 23 15 22 17 19 22 24
8 15 14 23 8 16 8 25 9 17 25 15 10 8 24
15 7 23 22 11 11 12 10 17 16 7 16 10 18 22
21 20 6 22 24 10 24 9 21 14 11 14 11 19 16
20 11 8 14 9 5 6 19 19 7 6 6 13 9 18
8 13 13 13 10 20 25 16 16 17 10 10 5 12 23
36 34 38 27 33
'''.split()

# read from file
# with open(r"filename", 'r') as input_file:
#     input_data = input_file.read().split()

# parse input
M = int(input_data[0])
J = int(input_data[1])
idx = np.cumsum([2, M * J, M * J, M])
cost = np.array(input_data[idx[0]:idx[1]], dtype=np.float64).reshape((M, J))
size = np.array(input_data[idx[1]:idx[2]], dtype=np.float64).reshape((M, J))
capacity = np.array(input_data[idx[2]:idx[3]], dtype=np.float64)

# build model
model = Highs()

X = model.addBinaries(M, J)

# assign each job to exactly one machine
model.addConstrs(X[:, j].sum() == 1 for j in range(J))

# each machine can only take jobs that fit
model.addConstrs((size[m, :] * X[m, :]).sum() <= capacity[m] for m in range(M))

# print out the cuts as we solve
def printCuts(e):
    for c in e.cuts:
        print(c)

model.cbMipGetCutPool += printCuts

# minimize total cost
model.minimize((cost * X).sum())

# print out solution (i.e., which jobs are assigned to which machines)
print(model.getObjectiveValue())

np.set_printoptions(formatter={'all': lambda x: str(x)})

for m in range(M):
    jobs_on_machine = np.nonzero(model.vals(X[m, :]) > 0)[0]
    print(jobs_on_machine)
