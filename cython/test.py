from linprog import highs_wrapper

import numpy as np
from scipy.sparse import csc_matrix

if __name__ == '__main__':

    c = np.array([5, 4, 3]).astype('double')
    A = np.array([
        [2, 3, 1],
        [4, 1, 2],
        [3, 4, 2],
    ]).astype('double')
    b = np.array([5, 11, 8]).astype('double')
    options = {
        'presolve': False,
        'sense': -1,
        'solver': None,
        'parallel': True,
        'time_limit': 1,
        'message_level': 0,
    }
    res = highs_wrapper(c, A, b, options=options)
    print(res)
