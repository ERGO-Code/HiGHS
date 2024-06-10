# model and solve the problem, first as an LP, then as a MIP
#
# minimize 8 TypeA +  10 TypeB
# s. t.    2 TypeA +   2 TypeB >=  7: Product1
#          3 TypeA +   4 TypeB >= 12: Product2
#          2 TypeA +   1 TypeB >=  6: Product3
#            TypeA >= 0; TypeB >= 0

import highspy

# For printing
width = 4
precision = 3

h = highspy.Highs()
h.silent()

variableNames = list()
variableNames.append('TypeA')
variableNames.append('TypeB')

useTypeA = h.addVariable(obj =  8, name = variableNames[0])
useTypeB = h.addVariable(obj = 10, name = variableNames[1])

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

print()
print('Solve as LP')

h.solve()

for var in vars:
    print('Use {0:.1f} of {1:s}: reduced cost {2:.6f}'.format(h.variableValue(var), h.variableName(var), h.variableDual(var)))
print('Use', h.variableValues(vars), 'of', h.variableNames(vars))
print('Use', h.allVariableValues(), 'of', h.allVariableNames())

for name in constrNames:
    print(f"Constraint {name} has value {h.constrValue(name):{width}.{precision}} and dual  {h.constrDual(name):{width}.{precision}}")

print('Constraints have values', h.constrValues(constrNames), 'and duals', h.constrDuals(constrNames))

print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

for var in vars:
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")

print()
print('Solve as MIP')

for var in vars:
    h.setInteger(var)

h.solve()

for var in vars:
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")

print()
print('Solve as LP with Gomory cut')

# Make the variables continuous
for var in vars:
    h.setContinuous(var)

# Add Gomory cut
h.addConstr(useTypeA + useTypeB >= 4, name = "Gomory")

h.solve()

for var in vars:
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")
