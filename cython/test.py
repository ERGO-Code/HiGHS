from linprog_mps import linprog_mps

if __name__ == '__main__':

    #linprog_mps(model_file='ex1.mps', solver='simplex')
    #linprog_mps(model_file='ex1.mps', solver='ipm', run_quiet=True)
    linprog_mps(model_file='25fv47', presolve=None, solver='choose', run_quiet=True)
