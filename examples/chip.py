# model and solve the LP
#
# maximize 10 Tables +  25 SetsOfChairs
# s. t.       Tables +   2 SetsOfChairs <=  80: Assembly
#             Tables +   4 SetsOfChairs <= 120: Finishing
#             Tables >= 0; SetsOfChairs >= 0
import highspy

h = highspy.Highs()
h.silent()

varNames = list()
varNames.append('Tables')
varNames.append('Sets of chairs')

x1 = h.addVar(obj = 10, name = varNames[0])
x2 = h.addVar(obj = 25, name = varNames[1])

vars = list()
vars.append(x1)
vars.append(x2)

constrNames = list()
constrNames.append('Assembly')
constrNames.append('Finishing')

h.addConstr(x1 + 2*x2 <=  80, name = constrNames[0])
h.addConstr(x1 + 4*x2 <= 120, name = constrNames[1])

h.setMaximize()

status = h.writeModel('Chip.lp')
print('writeModel(\'Chip.lp\') status =', status)
status = h.writeModel('Chip.mps')
print('writeModel(\'Chip.mps\') status =', status)

h.solve()

for var in vars:
    print('Make', h.varValue(var), h.varName(var), ': Reduced cost', h.varDual(var))
print('Make', h.varValues(vars), 'of', h.varNames(vars))
print('Make', h.allVarValues(), 'of', h.allVarNames())

for name in constrNames:
    print('Constraint', name, 'has value', h.constrValue(name), 'and dual', h.constrDual(name))

print('Constraints have values', h.constrValues(constrNames), 'and duals', h.constrDuals(constrNames))

print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

print('Optimal objective value is', h.getObjectiveValue())


