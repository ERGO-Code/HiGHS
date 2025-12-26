# Multi-objective binary knapsack example
from highspy import Highs, HighsLinearObjective
import numpy as np

# Parameters
capacity = 13

# chosen such that applying each objective lexicographically gives a different solution
profit = np.asarray([
        [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1],
        [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
        [1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0],
        [1,1,1,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,0,0],

    ], dtype=np.float64).reshape(4, 20)

OBJECTIVES = profit.shape[0]
ITEMS = profit.shape[1]

def pretty_print(objective_values):
    return '[ ' + ', '.join(map(lambda x: "{:2.0f}".format(x), objective_values)) + ' ]'

#
# individual optimization of each objective
def individual_objectives(h, X):
    print('## Individual:\n')
    objective_values = []

    for k in range(OBJECTIVES):
        h.maximize(np.dot(X, profit[k,:]))
        print(f'  Obj {k+1}: [{"".join(map(lambda x: "{:1.0f}".format(x), profit[k,:]))}]')
        print(f'  Sol {k+1}: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]\n')
        objective_values.append(h.getObjectiveValue())

    print(f'  OBJ: {pretty_print(objective_values)}\n\n')

#
# manual implementation of lexicographic multi-objective optimization
def manual_lexicographic(h, X):
    print('## Manual lexicographic:\n')
    cons = []
    objs = np.dot(profit, X)

    for k in range(OBJECTIVES):
        h.maximize(objs[k])
        print(f'  Obj {k+1}: [{"".join(map(lambda x: "{:1.0f}".format(x), profit[k,:]))}]')
        print(f'  Sol {k+1}: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]\n')

        # add constraint to ensure next solution are at least as good in this objective
        cons.append(h.addConstr(np.dot(profit[k,:], X) >= h.getObjectiveValue()))

    objective_values = h.vals(objs)
    h.deleteRows(len(cons), cons)  # clean up constraints
    print(f'  SOL: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]')
    print(f'  OBJ: {pretty_print(objective_values)}\n\n')


#
# built-in lexicographic multi-objective optimization
def highs_lexicographic(h, X):
    print('## Built-in lexicographic:\n')
    h.setOptionValue('blend_multi_objectives', False) # use lexicographic

    for k in range(OBJECTIVES):
        obj = HighsLinearObjective()
        obj.coefficients = profit[k,:].tolist()
        obj.weight = -1    # maximize
        obj.priority = -k  # higher priority for lower k
        obj.abs_tolerance = 0.01
        obj.rel_tolerance = 0.001

        h.addLinearObjective(obj)

    h.solve()
    print(f'  SOL: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]')
    print(f'  OBJ: {pretty_print(np.dot(profit, h.vals(X)))}\n\n')

    print(f'  Number of objectives: {h.getNumLinearObjectives()}')
    for k in range(h.getNumLinearObjectives()):
        obj = h.getLinearObjective(k)
        print(f'  Obj {k+1}: weight={obj.weight}, priority={obj.priority}, abs_tol={obj.abs_tolerance}, rel_tol={obj.rel_tolerance}')
    print('\n')

    h.clearLinearObjectives()


#
# manual implementation of weighted multi-objective optimization
def manual_weighted(h, X):
    print('## Manual Weighted:\n')
    weights = np.asarray([1.0/(k+1) for k in range(OBJECTIVES)], dtype=np.float64)
    h.maximize(np.dot(weights[:,None] * profit, X).sum())

    print(f'  SOL: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]')
    print(f'  OBJ: {pretty_print(np.dot(profit, h.vals(X)))}\n')

    print(f'  Obj: {h.getObjective()[0]}\n\n')


#
# built-in weighted multi-objective optimization
def highs_weighted(h, X):
    print('## Built-in Weighted:\n')
    h.setOptionValue('blend_multi_objectives', True) # use weighted

    for k in range(OBJECTIVES):
        obj = HighsLinearObjective()
        obj.coefficients = profit[k,:].tolist()
        obj.weight = -1.0/(k+1)

        h.addLinearObjective(obj)

    h.solve()
    print(f'  SOL: [{"".join(map(lambda x: "{:1.0f}".format(x), abs(h.vals(X))))}]')
    print(f'  OBJ: {pretty_print(np.dot(profit, h.vals(X)))}\n')

    print(f'  Number of objectives: {h.getNumLinearObjectives()}')
    for k in range(h.getNumLinearObjectives()):
        obj = h.getLinearObjective(k)
        print(f'  Obj {k+1}: weight={obj.weight}, priority={obj.priority}, abs_tol={obj.abs_tolerance}, rel_tol={obj.rel_tolerance}')
    print('\n')

    h.clearLinearObjectives()


if __name__ == "__main__":
    h = Highs()
    h.silent(True)

    X = h.addBinaries(ITEMS)
    h.addConstr(X.sum() <= capacity)

    individual_objectives(h, X)

    manual_lexicographic(h, X)
    manual_weighted(h, X)

    highs_lexicographic(h, X)
    highs_weighted(h, X)
