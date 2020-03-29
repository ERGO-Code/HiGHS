from linprog_mps import linprog_mps
from linprog import linprog

import numpy as np
from scipy.sparse import csc_matrix

if __name__ == '__main__':

    #linprog_mps(model_file='ex1.mps', solver='simplex')
    #linprog_mps(model_file='ex1.mps', solver='ipm', run_quiet=True)
    #linprog_mps(model_file='25fv47', presolve=None, solver='choose', run_quiet=True)

    c = np.array([5, 4, 3]).astype('double')
    A = np.array([
        [2, 3, 1],
        [4, 1, 2],
        [3, 4, 2],
    ]).astype('double')
    b = np.array([5, 11, 8]).astype('double')
    res = linprog(c, A, b, presolve=False, sense=-1)
    print(res)
