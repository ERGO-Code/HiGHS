# model and solve the LP
#
# maximize 10 Tables +  25 SetsOfChairs
# s. t.       Tables +   2 SetsOfChairs <=  80: Assembly
#             Tables +   4 SetsOfChairs <= 120: Finishing
#             Tables >= 0; SetsOfChairs >= 0
import highspy

h = highspy.Highs()
h.silent()

items = ['Tables', 'Sets of chairs']
x = h.addVariables(items, obj = [10, 20], name = items)

constrNames = ['Assembly', 'Finishing']
cons = h.addConstrs(x['Tables'] + 2*x['Sets of chairs'] <=  80, 
                    x['Tables'] + 4*x['Sets of chairs'] <= 120, name = constrNames)

h.setMaximize()

status = h.writeModel('Chip.lp')
print('writeModel(\'Chip.lp\') status =', status)
status = h.writeModel('Chip.mps')
print('writeModel(\'Chip.mps\') status =', status)

h.solve()


for n, var in x.items():
    print('Make', h.variableValue(var), h.variableName(var), ': Reduced cost', h.variableDual(var))
    
print('Make', h.variableValues(x.values()), 'of', h.variableNames(x.values()))
print('Make', h.allVariableValues(), 'of', h.allVariableNames())

for c in cons:
    print('Constraint', c.name, 'has value', h.constrValue(c), 'and dual', h.constrDual(c))
    
print('Constraints have values', h.constrValues(cons), 'and duals', h.constrDuals(cons))
print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())

print('Optimal objective value is', h.getObjectiveValue())


