'''Use cython wrapper to solve problem described by MPS file.'''

import pathlib

from pyHiGHS.linprog_mps import linprog_mps

if __name__ == '__main__':

    mpsfile = str(pathlib.Path(__file__).parent / '25fv47')
    linprog_mps(mpsfile)
