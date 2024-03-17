# model and solve the problem, first as an LP, then as a MIP
#
# minimize 8 TypeA +  10 TypeB
# s. t.    2 TypeA +   2 TypeB >=  7: Product1
#          3 TypeA +   4 TypeB >= 12: Product2
#          2 TypeA +   1 TypeB >=  6: Product3
#            TypeA >= 0; TypeB >= 0

import highspy

h = highspy.Highs()
h.silent()

varNames = list()
varNames.append('TypeA')
varNames.append('TypeB')

useTypeA = h.addVar(obj =  8, name = varNames[0])
useTypeB = h.addVar(obj = 10, name = varNames[1])

vars = list()
vars.append(useTypeA)
vars.append(useTypeB)

constrNames = list()
constrNames.append('Product1')
constrNames.append('Product2')
constrNames.append('Product3')

h.addConstr(2*useTypeA + 2*useTypeB >=  7, name = constrNames[0])
h.addConstr(3*useTypeA + 4*useTypeB >= 12, name = constrNames[1])
h.addConstr(2*useTypeA +   useTypeB >=  6, name = constrNames[2])

h.setMinimize()

status = h.writeModel('Distillation.lp')
print('writeModel(\'Distillation.lp\') status =', status)
status = h.writeModel('Distillation.mps')
print('writeModel(\'Distillation.mps\') status =', status)

h.solve()

print()
print('Solved as LP')

for var in vars:
    print('Use', h.varValue(var), h.varName(var), ': Reduced cost', h.varDual(var))
print('Use', h.varValues(vars), 'of', h.varNames(vars))
print('Use', h.allVarValues(), 'of', h.allVarNames())

for name in constrNames:
    print('Constraint', name, 'has value', h.constrValue(name), 'and dual', h.constrDual(name))

print('Constraints have values', h.constrValues(constrNames), 'and duals', h.constrDuals(constrNames))

print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

print('Optimal objective value is', h.getObjectiveValue())

for var in vars:
    h.setInteger(var)

h.solve()

print()
print('Solved as MIP')

for var in vars:
    print('Use', h.varValue(var), h.varName(var))

print('Optimal objective value is', h.getObjectiveValue())
