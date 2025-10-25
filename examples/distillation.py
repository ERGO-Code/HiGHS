# model and solve the problem, first as an LP, then as a MIP
#
# minimize 8 TypeA +  10 TypeB
# s. t.    2 TypeA +   2 TypeB >=  7: Product1
#          3 TypeA +   4 TypeB >= 12: Product2
#          2 TypeA +   1 TypeB >=  6: Product3
#            TypeA >= 0; TypeB >= 0

import highspy
width, precision = (4, 3)   # for printing

h = highspy.Highs()
h.silent()

variableNames = ['TypeA', 'TypeB']
x = h.addVariables(variableNames, obj = [8, 10], name = variableNames)

constrNames = ['Product1', 'Product2', 'Product3']
cons = h.addConstrs(2*x['TypeA'] + 2*x['TypeB'] >=  7,
                    3*x['TypeA'] + 4*x['TypeB'] >= 12,
                    2*x['TypeA'] +   x['TypeB'] >=  6, name = constrNames)

h.setMinimize()

status = h.writeModel('Distillation.lp')
print(f"writeModel('Distillation.lp') status = {status}")

status = h.writeModel('Distillation.mps')
print(f"writeModel('Distillation.mps') status = {status}")


print()
print('# Solve as LP')

h.solve()

for name, var in x.items():
    print(f'Use {h.variableValue(var):.1f} of {name}: reduced cost {h.variableDual(var):.6f}')
    
print()
print('Use', h.vals(x))
print('Use', h.allVariableValues(), 'of', h.allVariableNames())
print()

for c in cons:
    print(f"Constraint {c.name} has value {h.constrValue(c):{width}.{precision}} and dual  {h.constrDual(c):{width}.{precision}}")

print('Constraints have values', h.constrValues(cons), 'and duals', h.constrDuals(cons))
print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")


print()
print()
print('# Solve as MIP')

h.setInteger(x.values())
h.solve()

for var in x.values():
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")

print()
print()
print('# Solve as LP with Gomory cut')

# Make the variables continuous
h.setContinuous(x.values())

# Add Gomory cut
h.addConstr(x['TypeA'] + x['TypeB'] >= 4, name = "Gomory")
h.solve()

for var in x.values():
    print(f"Use {h.variableValue(var):{width}.{precision}} of {h.variableName(var)}")
print(f"Optimal objective value is {h.getObjectiveValue():{width}.{precision}}")
