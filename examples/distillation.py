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

print("DEBUG ", varNames)

update_in_addVariable = True
if update_in_addVariable:
    print("DEBUG update_in_addVariable is True")
else:
    print("DEBUG update_in_addVariable is False")
    
useTypeA = h.addVariable(obj =  8, name = varNames[0], update = update_in_addVariable)
useTypeB = h.addVariable(obj = 10, name = varNames[1], update = update_in_addVariable)

# With update_in_addVariable = False (so runs as originally written,
# with self.update() being called only in addConstr) useTypeB.name is
# "TypeB"
#
# With update_in_addVariable = True (so self.update() is called in addVariable) useTypeB.name is
# "TypeA"
#
# Looks as if the internal stack of names isn't being cleared, so
# TypeA is still on top when useTypeB is being defined

print('\nDEBUG Names: useTypeA', useTypeA.name, "; varNames[0]", varNames[0])
print('DEBUG Names: useTypeB', useTypeB.name, "; varNames[1]", varNames[1])

vars = list()
vars.append(useTypeA)
vars.append(useTypeB)

print("\nDEBUG vars: ", vars)

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
    print('Use {0:6f} of {1:6s}: reduced cost {2:6f}'.format(h.varValue(var), h.varName(var), h.varDual(var)))
print('Use', h.varValues(vars), 'of', h.varNames(vars))
print('Use', h.allVarValues(), 'of', h.allVarNames())

for name in constrNames:
    print('Constraint {0:6s} has value {1:6f} and dual {2:6f}'.format(name, h.constrValue(name), h.constrDual(name)))

print('Constraints have values', h.constrValues(constrNames), 'and duals', h.constrDuals(constrNames))

print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

print('Optimal objective value is', h.getObjectiveValue())

for var in vars:
    h.setInteger(var)

h.solve()

print()
print('Solved as MIP')

for var in vars:
    print('Use {0:6f} of {1:6s}'.format(h.varValue(var), h.varName(var)))

print('Optimal objective value is', h.getObjectiveValue())
