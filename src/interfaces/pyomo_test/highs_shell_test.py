from pyomo.core.base.piecewise import Bound
from pyomo.environ import *
import highs
import pytest


EPS = 1e-4


def decimal_error_judge(a, b):
    return (abs(a-b)/max(1, b) <= EPS)


def test_ConcreteModel_1():
    # create a model
    model = ConcreteModel()

    # declare decision variables
    model.x = Var(bounds=(-1, 1))  # NonNegativeReals代表非0实数
    model.y = Var(domain=NonNegativeReals)

    # declare objective
    model.profit = Objective(expr=40 * model.x + 30 * model.y, sense=maximize)

    # declare constraints
    model.demand = Constraint(expr=model.x <= 40)
    model.laborA = Constraint(expr=model.x + model.y <= 80)
    model.laborB = Constraint(expr=2 * model.x + model.y <= 100)

    # model.pprint()
    opt = SolverFactory('highs')

    opt.solve(model)
    assert (decimal_error_judge(model.profit(), 2410.0)) == True
    assert (decimal_error_judge(model.x(), 1.0)) == True
    assert (decimal_error_judge(model.y(), 79.0)) == True
    assert (decimal_error_judge(model.demand(), 1)) == True
    assert (decimal_error_judge(model.laborA(), 80)) == True
    assert (decimal_error_judge(model.laborB(), 81)) == True


def test_ConcreteModel_2():
    A = ['hammer', 'wrench', 'screwdriver', 'towel']
    b = {'hammer': 8, 'wrench': 3, 'screwdriver': 6, 'towel': 11}
    w = {'hammer': 5, 'wrench': 7, 'screwdriver': 4, 'towel': 3}
    W_max = 14
    model = ConcreteModel()
    model.x = Var(A, bounds=(0, 1))
    model.value = Objective(
        expr=sum(b[i]*model.x[i] for i in A), sense=maximize)
    model.weight = Constraint(
        expr=sum(w[i]*model.x[i] for i in A) <= W_max)

    opt = SolverFactory('highs')
    result_obj = opt.solve(model, tee=True)

    assert (decimal_error_judge(model.value(), 25.857142857142858)) == True
    assert (decimal_error_judge(model.x['hammer'](), 1)) == True
    assert (decimal_error_judge(
        model.x['wrench'](), 0.2857142857142857)) == True
    assert (decimal_error_judge(model.x['screwdriver'](), 1)) == True
    assert (decimal_error_judge(model.x['towel'](), 1)) == True


def test_AbstractModel_smalldata_1():
    model = AbstractModel()
    model.ITEMS = Set()
    model.v = Param(model.ITEMS, within=PositiveReals)
    model.w = Param(model.ITEMS, within=PositiveReals)
    model.W_max = Param(within=PositiveReals)
    model.x = Var(model.ITEMS, bounds=(0, 1))

    def value_rule(model):
        return sum(model.v[i]*model.x[i] for i in model.ITEMS)
    model.value = Objective(rule=value_rule, sense=maximize)

    def weight_rule(model):
        return sum(model.w[i]*model.x[i] for i in model.ITEMS) <= model.W_max
    model.weight = Constraint(rule=weight_rule)

    instance = model.create_instance("check/instances/knapsack.dat")
    opt = SolverFactory('highs')
    opt.solve(instance)

    assert (decimal_error_judge(
        instance.value(), 25.857142857142858)) == True
    assert (decimal_error_judge(instance.x['hammer'](), 1)) == True
    assert (decimal_error_judge(
        instance.x['wrench'](), 0.2857142857142857)) == True
    assert (decimal_error_judge(instance.x['screwdriver'](), 1)) == True
    assert (decimal_error_judge(instance.x['towel'](), 1)) == True


def test_AbstractModel_smalldata_2():
    # Minimize
    # Obj: 1.7 * x1 + 7.2 * x2 + 9 x3 + 8.3 * x4 + 2.9 * x5 + 6.3 * x6 + 9.8 * x7 + 0.7 * x8 + 4.5 * x9 + 4.8 * x10
    # + 4.2 * x11 + 9.3 * x12
    # Subject To
    # x16:  x1 - x13 <= 0
    # x17:  x2 - x13 <= 0
    # x18:  x3 - x13 <= 0
    # x19:  x4 - x13 <= 0
    # x20:  x5 - x14 <= 0
    # x21:  x6 - x14 <= 0
    # x22:  x7 - x14 <= 0
    # x23:  x8 - x14 <= 0
    # x24:  x9 - x15 <= 0
    # x25:  x10 - x15 <= 0
    # x26:  x11 - x15 <= 0
    # x27:  x12 - x15 <= 0
    # x28:  x13 + x14 + x15 = 2
    # x29:  x1 + x5 + x9 = 1
    # x30:  x2 + x6 + x10 = 1
    # x31:  x3 + x7 + x11 = 1
    # x32:  x4 + x8 + x12 = 1
    # Bounds
    # 0 <= x1 <= 1
    # 0 <= x2 <= 1
    # 0 <= x3 <= 1
    # 0 <= x4 <= 1
    # 0 <= x5 <= 1
    # 0 <= x6 <= 1
    # 0 <= x7 <= 1
    # 0 <= x8 <= 1
    # 0 <= x9 <= 1
    # 0 <= x10 <= 1
    # 0 <= x11 <= 1
    # 0 <= x12 <= 1
    # 0 <= x13 <= 1
    # 0 <= x14 <= 1
    # 0 <= x15 <= 1
    model = AbstractModel()
    model.N = Param(within=PositiveIntegers)
    model.P = Param(within=RangeSet(model.N))
    model.M = Param(within=PositiveIntegers)
    model.Locations = RangeSet(model.N)
    model.Customers = RangeSet(model.M)
    model.d = Param(model.Locations, model.Customers)
    model.x = Var(model.Locations, model.Customers, bounds=(0.0, 1.0))
    model.y = Var(model.Locations, bounds=(0, 1))

    def obj_rule(model):
        return sum(model.d[n, m]*model.x[n, m]
                   for n in model.Locations for m in model.Customers)
    model.obj = Objective(rule=obj_rule)

    def single_x_rule(model, m):
        return sum(model.x[n, m] for n in model.Locations) == 1.0
    model.single_x = Constraint(model.Customers, rule=single_x_rule)

    def bound_y_rule(model, n, m):
        return model.x[n, m] - model.y[n] <= 0.0
    model.bound_y = Constraint(model.Locations, model.Customers,
                               rule=bound_y_rule)

    def num_facilities_rule(model):
        return sum(model.y[n] for n in model.Locations) == model.P
    model.num_facilities = Constraint(rule=num_facilities_rule)

    instance = model.create_instance("check/instances/pmedian.dat")
    opt = SolverFactory('highs')
    opt.solve(instance)

    except_x_1 = [(2, 1), (2, 4), (3, 2), (3, 3)]

    assert (decimal_error_judge(instance.obj(), 12.6))
    for i in range(1, 4):
        for j in range(1, 5):
            if except_x_1.count((i, j)):
                assert (decimal_error_judge(instance.x[i, j](), 1)) == True
            else:
                assert (decimal_error_judge(instance.x[i, j](), 0)) == True
    assert (decimal_error_judge(instance.y[1](), 0)) == True
    assert (decimal_error_judge(instance.y[2](), 1)) == True
    assert (decimal_error_judge(instance.y[3](), 1)) == True


def test_AbstractModel_bigdata_1():
    model = AbstractModel()
    model.N = Param(within=PositiveIntegers)
    model.M = Param(within=PositiveIntegers)
    model.Locations = RangeSet(model.N)
    model.Customers = RangeSet(model.M)

    model.C = Param(model.Locations)
    model.D = Param(model.Customers)
    model.W = Param(model.Customers, model.Locations)
    model.x = Var(model.Locations, within=NonNegativeReals)

    def value_rule(model):
        return sum(model.C[i]*model.x[i] for i in model.Locations)
    model.value = Objective(rule=value_rule, sense=minimize)

    def single_x_rule(model, m):
        return sum(model.W[m, n]*model.x[n] for n in model.Locations) >= model.D[m]
    model.single_x = Constraint(model.Customers, rule=single_x_rule)

    # There is a 20*20 constraint matrix.
    instance = model.create_instance("check/instances/big_data_1.dat")
    opt = SolverFactory('highs')

    results = opt.solve(instance).write()
    assert (decimal_error_judge(instance.value(), 56)) == True


def test_AbstractModel_bigdata_2():
    model = AbstractModel()
    model.N = Param(within=PositiveIntegers)
    model.M = Param(within=PositiveIntegers)
    model.Locations = RangeSet(model.N)
    model.Customers = RangeSet(model.M)

    model.C = Param(model.Locations)
    model.D = Param(model.Customers)
    model.W = Param(model.Customers, model.Locations)
    model.x = Var(model.Locations, within=NonNegativeReals)

    def value_rule(model):
        return sum(model.C[i]*model.x[i] for i in model.Locations)
    model.value = Objective(rule=value_rule, sense=minimize)

    def single_x_rule(model, m):
        return sum(model.W[m, n]*model.x[n] for n in model.Locations) >= model.D[m]
    model.single_x = Constraint(model.Customers, rule=single_x_rule)

    # There is a 90*900 constraint matrix.
    instance = model.create_instance("check/instances/big_data_2.dat")
    opt = SolverFactory('highs')

    results = opt.solve(instance).write()
    assert (decimal_error_judge(instance.value(), 9168492)) == True

def test_AbstractModel_maxflow():
    model = AbstractModel()

    model.nodes = Set()
    model.arcs = Set(within=model.nodes*model.nodes)
    model.sources = Set(within=model.nodes)
    model.sinks = Set(within=model.nodes)
    model.upperBound = Param(model.arcs)
    model.supply = Param(model.sources)
    model.demand = Param(model.sinks)
    model.amount = Var(model.arcs, within=NonNegativeReals)

    def totalRule(model):
        expression = sum(
            model.amount[i, j]
            for (i, j) in model.arcs
            if j in model.sinks
        )
        return expression

    model.maxFlow = Objective(rule=totalRule, sense=maximize)

    def maxRule(model, arcIn, arcOut):
        constraint_equation = (
            model.amount[arcIn, arcOut] <= model.upperBound[arcIn, arcOut])
        return constraint_equation

    model.loadOnArc = Constraint(model.arcs, rule=maxRule)

    def flowRule(model, node):
        if node in model.sources:
            flow_out = sum(
                model.amount[i, j]
                for (i, j) in model.arcs
                if i == node
            )
            constraint_equation = (flow_out <= model.supply[node])

        elif node in model.sinks:
            flow_in = sum(
                model.amount[i, j]
                for (i, j) in model.arcs
                if j == node
            )
            constraint_equation = (flow_in >= model.demand[node])

        else:
            amountIn = sum(
                model.amount[i, j]
                for (i, j) in model.arcs
                if j == node
            )
            amountOut = sum(
                model.amount[i, j]
                for (i, j) in model.arcs
                if i == node
            )
            constraint_equation = (amountIn == amountOut)

        return constraint_equation

    model.flow = Constraint(model.nodes, rule=flowRule)

    instance = model.create_instance("check/instances/maxflow.dat")

    instance.pprint()

    opt = SolverFactory('highs')

    results = opt.solve(instance).write()
    assert (decimal_error_judge(instance.maxFlow(), 16.0)) == True


def test_AbstractModel_networkflow():
    model = AbstractModel()

    model.places = Set()
    model.routes = Set(within=model.places*model.places)
    model.supply = Param(model.places)
    model.demand = Param(model.places)
    model.cost = Param(model.routes)
    model.minimum = Param(model.routes)
    model.maximum = Param(model.routes)
    model.amount = Var(model.routes, within=NonNegativeReals)
    model.excess = Var(model.places, within=NonNegativeReals)

    def costRule(model):
        return sum(model.cost[n]*model.amount[n] for n in model.routes)

    model.costTotal = Objective(rule=costRule)

    def loadRule(model, i, j):
        return (model.minimum[i, j], model.amount[i, j], model.maximum[i, j])

    model.loadOnRoad = Constraint(model.routes, rule=loadRule)

    def supplyDemandRule(model, nn):

        amountIn = sum(model.amount[i, j]
                       for (i, j) in model.routes if j == nn)
        amountOut = sum(model.amount[i, j]
                        for (i, j) in model.routes if i == nn)

        input = amountIn + model.supply[nn]
        output = amountOut + model.demand[nn] + model.excess[nn]

        return input == output

    model.supplyDemand = Constraint(model.places, rule=supplyDemandRule)
    instance = model.create_instance("check/instances/networkflow.dat")

    instance.pprint()

    opt = SolverFactory('highs')

    results = opt.solve(instance).write()
    assert (decimal_error_judge(instance.costTotal(), 11575.0)) == True
