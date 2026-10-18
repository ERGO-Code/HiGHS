# model and solve the LP
#
# maximize 10 Tables +  25 SetsOfChairs
# s. t.       Tables +   2 SetsOfChairs <=  80: Assembly
#             Tables +   4 SetsOfChairs <= 120: Finishing
#             Tables >= 0; SetsOfChairs >= 0
import highspy

h = highspy.Highs()
h.silent()

items = ['Tables', 'SetsOfChairs']
x = h.addVariables(items, obj = [10, 25], name = items)

constrNames = ['Assembly', 'Finishing']
cons = h.addConstrs(x['Tables'] + 2*x['SetsOfChairs'] <=  80, 
                    x['Tables'] + 4*x['SetsOfChairs'] <= 120, name = constrNames)
h.setMaximize()

status = h.writeModel('Chip.lp')
print(f"writeModel('Chip.lp') status = {status}")
status = h.writeModel('Chip.mps')
print(f"writeModel('Chip.mps') status = {status}\n")

h.solve()

for n, var in x.items():
    print('Make', h.variableValue(var), n, ': Reduced cost', h.variableDual(var))

print()
print('Make', h.vals(x))
print('Make', h.allVariableValues(), 'of', h.allVariableNames())
print()

for c in cons:
    print('Constraint', c.name, 'has value', h.constrValue(c), 'and dual', h.constrDual(c))
    
print('Constraints have values', h.constrValues(cons), 'and duals', h.constrDuals(cons))
print('Constraints have values', h.allConstrValues(), 'and duals', h.allConstrDuals())
print()
print('Optimal objective value is', h.getObjectiveValue())


