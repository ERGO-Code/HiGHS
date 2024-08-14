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

variableNames = ['TypeA', 'TypeB']
x = h.addVariables(variableNames, obj = [8, 10], name = variableNames[0])

constrNames = ['Product1', 'Product2', 'Product3']
cons = h.addConstrs(2*x['TypeA'] + 2*x['TypeB'] >=  7,
                    3*x['TypeA'] + 4*x['TypeB'] >= 12,
                    2*x['TypeA'] +   x['TypeB'] >=  6, name = constrNames)

h.setMinimize()

status = h.writeModel('Distillation.lp')
print('writeModel(\'Distillation.lp\') status =', status)
status = h.writeModel('Distillation.mps')
print('writeModel(\'Distillation.mps\') status =', status)

print()
print('Solve as LP')

h.solve()

for name, var in x.items():
    print('Use {0:.1f} of {1:s}: reduced cost {2:.6f}'.format(h.variableValue(var), h.variableName(var), h.variableDual(var)))
    
print('Use', h.variableValues(x.values()), 'of', h.variableNames(x.values()))
print('Use', h.allVariableValues(), 'of', h.allVariableNames())

for c in cons:
    print(f"Constraint {c.name} has value {h.constrValue(c):{width}.{precision}} and dual  {h.constrDual(c):{width}.{precision}}")

print('Constraints have values', h.constrValues(cons), 'and duals', h.constrDuals(cons))
print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

for var in x.values():
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")

print()
print('Solve as MIP')

for var in x.values():
    h.setInteger(var)

h.solve()

for var in x.values():
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")

print()
print('Solve as LP with Gomory cut')

# Make the variables continuous
for var in x.values():
    h.setContinuous(var)

# Add Gomory cut
h.addConstr(x['TypeA'] + x['TypeB'] >= 4, name = "Gomory")

h.solve()

for var in x.values():
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")
