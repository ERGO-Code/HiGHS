import highspy, random, time, math

CG_GAP_REL = 1e-4

random.seed(100)

#Total number of items
N = 60

#Item weight
weights = [round(random.uniform(1, 10),4) for _ in range(N)]

# Bin capacity
capacity = 15

# Maximum number of bins
B = N # Easy upper bound on the number of bins

def solveCompactModel():
    m = highspy.Highs()
    m.setOptionValue('output_flag', True)
    m.setOptionValue('mip_abs_gap', 1-1e-5) # Since all variables in the objective are binary, it should stop when the absolute gap is below 1

    idx = 0
    x = dict()
    for i in range(N):
        for j in range(B):
            m.addVar(0, 1)
            x[i,j] = idx
            idx += 1
    
    y = dict()
    for j in range(B):
        m.addVar(0, 1)
        y[j] = idx
        idx += 1

    m.changeColsIntegrality(
        len(x)+len(y),
        list(range(len(x)+len(y))),
        [highspy.HighsVarType.kInteger]*(len(x)+len(y))
    )

    # Min: sum_{j} y_{j}
    m.changeObjectiveSense(highspy.ObjSense.kMinimize)
    m.changeColsCost(len(y), list(y.values()), [1]*len(y))
    
    for i in range(N):
        # \sum_{j} x_{ij} == 1, \foreach i ---> 1 <= \sum_{j} x_{ij} <= 1
        m.addRow(
            1,                                              # lhs
            1,                                              # rhs
            B,                                              # Number of non-zero variables
            [x[i,j] for j in range(B)],                     # Indexes of variable
            [1] * B                                         # Coefficients of variables
        )
        for j in range(B):
            # x_{ij} <= y_{j}, \foreach i,j ---> -inf <= x_{ij}-y_{j} <= 0
            m.addRow(
                -highspy.kHighsInf,                         # lhs
                0,                                          # rhs
                2,                                          # Number of non-zero variables 
                [x[i,j], y[j]],                             # Indexes of variable
                [1, -1]                                     # Coefficients of variables
            )

    for j in range(B):
        # sum_{i}x_{ij} * w_{i} <= c*y_{j}, \foreach j ---> -inf <= sum_{i}x_{ij} * w_{i} - c*y_{j} <= 0
        m.addRow(
            -highspy.kHighsInf,                             # lhs
            0,                                              # rhs
            N+1,                                            # Number of non-zero variables 
            [x[i,j] for i in range(N)] + [y[j]],            # Indexes of variable
            [weights[i] for i in range(N)] + [-capacity]    # Coefficients of variables  
        )
    
    m.run()

    vals = list(m.getSolution().col_value)

    bins = [
        [i for i in range(N) if vals[x[i,j]] > 0.9]
        for j in range(B)
        if vals[y[j]] > 0.9
    ]
    
    return bins


def createMasterProblem(columns: list, solve=True):
    m = highspy.Highs()
    m.setOptionValue('output_flag', False)
    m.setOptionValue('random_seed', 100)
    
    m.addVars(len(columns), [0]*len(columns), [1]*len(columns))
    
    # Min: \sum_{k} \lambda_{k}
    m.changeObjectiveSense(highspy.ObjSense.kMinimize)
    m.changeColsCost(len(columns), list(range(len(columns))), [1]*len(columns))
    
    for i in range(N):
        # \sum_{k \in K_i} \lambda_{k} == 1, \foreach i ---> 1 <= \sum_{k \in K_i} \lambda_{k} <= 1, \foreach i
        # which K_i is the set of columns where i appears
        columns_with_item_i = [idx for idx, column in enumerate(columns) if i in column]
        m.addRow(
            1,                                  # lhs
            1,                                  # rhs
            len(columns_with_item_i),           # Number of non-zero variables
            columns_with_item_i,                # Indexes of variable
            [1]*len(columns_with_item_i)        # Coefficients of variables
        )

    if solve:
        m.run() 

        solution = m.getSolution()
        vals = list(solution.col_value)
        duals = list(solution.row_dual)

        return m, vals, duals
    else:
        return m, None, None

def solveSubproblemNotExact(duals):

    sorting = sorted(
        [
            (i, -duals[i]/weights[i])
            for i in range(N)
            if duals[i] >= 0
        ],
        key=lambda el: el[1]
    )

    sol = []

    tt_weight = 0
    for i, _ in sorting:
        if weights[i] + tt_weight <= capacity:
            sol += [i]
            tt_weight += weights[i]
    
    obj = 1 - sum(duals[i] for i in sol)

    return obj, sorted(sol)


def solveSubproblemExact(duals):
    m = highspy.Highs()
    m.setOptionValue('output_flag', False)
    m.setOptionValue('random_seed', 100)
    m.setOptionValue('mip_abs_gap', 1-1e-5) # Since all variables are binary, it should stop when the absolute gap is below 1

    m.addVars(N, [0]*N, [1]*N)
    m.changeColsIntegrality(
        N,
        list(range(N)),
        [highspy.HighsVarType.kInteger]*N
    )

    # Min: 1 - \sum_{i} \mu_{i} * x_{i}
    # Where \mu_{i} is the dual variable for the i-th row of the master problem
    m.changeColsCost(N, list(range(N)), [-duals[i] for i in range(N)])
    m.changeObjectiveSense(highspy.ObjSense.kMinimize)

    # sum_{i}x_{i} * w_{i} <= c, \foreach j ---> -inf <= sum_{i}x_{i} * w_{i} <= c, \foreach j
    m.addRow(
        -highspy.kHighsInf,
        capacity,
        N,
        [i for i in range(N)],
        [weights[i] for i in range(N)]
    )
    
    m.run()

    obj = m.getInfoValue("objective_function_value")[1]
    solution = m.getSolution()
    vals = list(solution.col_value)
    
    new_column = sorted([i for i in range(N) if vals[i] > 0.9])

    return 1+obj, new_column

def getUniqueColumns(columns):
    unique_columns = []
    for column in columns:
        if column not in unique_columns:
            unique_columns += [column]

    if len(unique_columns) != len(columns):
        return unique_columns
    
    return columns

def generateColumns(columns: list, m, msg=True, solve_exact = False, start_time = 0):
    vals = list(m.getSolution().col_value)
    duals = list(m.getSolution().row_dual)
    ub = sum(vals)

    best_gap = math.inf

    iter = 0
    while True:
        iter += 1     
        if not solve_exact:
            obj_sub, new_column = solveSubproblemNotExact(duals)
        else:
            obj_sub, new_column = solveSubproblemExact(duals)

        lb = ub+obj_sub
        gap = (ub-lb)/ub
        if gap < best_gap:
            best_gap = gap
            if msg:
                row = [
                    f"{len(columns)}",
                    f"{round(obj_sub, 3) :.4f}",
                    f'{round(lb, 3) :.4f}',
                    f"{round(ub, 3) :.4f}",
                    f"{max(0, gap) :.3%}", 
                    f"{round(time.perf_counter()-start_time, 2) :.2f}"
                ]
                if iter == 1:
                    header = ["Columns", 'Pricing', 'LB', 'UB', 'gap', 'Time']
                    row_format = "".join([
                        "{:>"+str(max(len(row_el), len(header_el))+3)+"}"
                        for row_el, header_el in zip(row, header)
                    ])
                    print(row_format.format(*header))

                print(row_format.format(*row))

        if obj_sub == 0 or gap < CG_GAP_REL:
            break

        # Adds a new column to the master problem
        columns += [new_column]
        m.addCol(
            1,                      # cost
            0,                      # lower bound
            1,                      # upper bound
            len(new_column),        # number of rows
            new_column,             # indexes of rows
            [1]*len(new_column)     # coefficients of rows
        )

        # Solves the master problem
        m.run()

        solution = m.getSolution()
        vals = list(solution.col_value)
        duals = list(solution.row_dual)
        ub = sum(vals)

    return m, columns, vals

def getFractionalColumns(vals: list, columns: list):
    return sorted(
        [
            idx
            for idx, val in enumerate(vals)
            if abs(round(val)-val) > 1e-6
        ],
        key=lambda idx: -len(columns[idx]) # Orders by the most fractional
    )

class Node:
    id: int
    layer: int
    assigned_columns: list[int]

    value: float
    final_columns: list[int]
    final_columns_vals: list[float]
    final_fractional_columns: list[int]

    def __init__(self, id: int, layer: int, assigned_columns: list[int], parent = None):
        self.id = id
        self.layer = layer
        self.parent = parent
        self.assigned_columns = assigned_columns


def solveNode(node: Node, columns: list[list[int]], m):
    lower_bounds = [
        0 if idx not in node.assigned_columns
        else 1
        for idx in range(len(columns))
    ]
    m.changeColsBounds(
        len(columns),                       # Qty columns to change bounds
        list(range(len(columns))),          # Which columns
        lower_bounds,                       # Lower bound (0 if not assigned else 1)
        [1]*len(columns)                    # Upper bound (always 1)
    )

    # Making sure that there is no repeated index
    assert len(node.assigned_columns) == len(set(node.assigned_columns))

    m.run()

    if m.getModelStatus() == highspy.HighsModelStatus.kOptimal:        
        m, columns, vals = generateColumns(columns, m, msg=False, solve_exact=False)
        # Proving optimallity on node
        m, columns, vals = generateColumns(columns, m, msg=False, solve_exact=True)

        node.value = sum(vals)
        node.final_columns = columns
        node.final_columns_vals = vals
        node.final_fractional_columns = getFractionalColumns(vals, columns)

        # Making sure that no column already assigned to 1 in the node
        # is returned back with a fractional value
        assert len(set(node.final_fractional_columns).intersection(set(node.assigned_columns))) == 0

        return True, m, node
    
    return False, m, node


def getKey(columns_idxs):
    return ",".join(sorted([f"{idx}" for idx in columns_idxs]))

def branchAndPrice(m, vals, columns, start_time=0):    
    header = ["NodesExpl", "TreeSize", "CurrCols", "FracCols", "Time"]
    row_format ="{:>12}" * (len(header))
    print(row_format.format(*header))
    
    fractional_columns = getFractionalColumns(vals, columns)

    node_id = 0

    branch_tree: list[Node] = []
    branch_tree_dict = dict()
    for column_idx in fractional_columns:
        node = Node(node_id, 0, [column_idx])
        branch_tree += [node]
        branch_tree_dict[getKey([column_idx])] = node
        node_id += 1

    best_qtd_frac_columns = math.inf

    qty_nodes_visited = 0

    while len(branch_tree) > 0:
        # Searchs for a node in the last layer,
        # with the least quantity of fractional columns of its parent,
        # and the highest value of its original fractional columns
        node = min(branch_tree, key=lambda node: (
            -node.layer,
            len(node.parent.final_fractional_columns) if node.parent is not None else math.inf,
            -node.parent.final_columns_vals[node.assigned_columns[-1]] if node.parent is not None else 0
        ))
        branch_tree.remove(node)

        # Solves the node with column generation
        mp_is_feasible, m, node = solveNode(node, columns, m)
        
        qty_nodes_visited += 1
        if mp_is_feasible:
            
            # It will only add columns if the master problem is feasible
            columns = node.final_columns
            
            qtd_fractional_columns = len(node.final_fractional_columns)

            if qtd_fractional_columns < best_qtd_frac_columns:
                best_qtd_frac_columns = qtd_fractional_columns
                print(row_format.format(*[
                    qty_nodes_visited,
                    len(branch_tree),
                    len(columns),
                    best_qtd_frac_columns,
                    f"{round(time.perf_counter()-start_time, 2) :.2f}"
                ]))

            if qtd_fractional_columns > 0:
                # Adds all fractional columns as new nodes
                for column_idx in node.final_fractional_columns:
                    new_assigned_columns = node.assigned_columns + [column_idx]
                    key = getKey(new_assigned_columns)
                    if key not in branch_tree_dict:
                        node_id += 1
                        new_node = Node(node_id, node.layer+1, new_assigned_columns, node)
                        branch_tree = [new_node] + branch_tree
                        branch_tree_dict[key] = new_node

            else:
                return [column for idx, column in enumerate(node.final_columns) if node.final_columns_vals[idx] > 0.9]

    raise Exception("It should not get to this point")

if __name__ == '__main__':

    start_time = time.perf_counter()
    bins = solveCompactModel()

    print(f"\nSolution by compact model: {len(bins)} bins")
    for bin, items in enumerate(bins):
        tt_weight = round(sum(weights[i] for i in items))
        print(f"Bin {bin+1} ({tt_weight} <= {capacity}): {items}")
        assert tt_weight <= capacity

    print(f"Finished in {round(time.perf_counter()-start_time,2): .2f}s")

    print()

    start_time = time.perf_counter()

    columns = [[i] for i in range(N)]

    m, vals, duals = createMasterProblem(columns)
    print("\nSolving root node:")
    m, columns, vals = generateColumns(columns, m, solve_exact=False, start_time=start_time)
    print("\nProving optimallity on root node:")
    m, columns, vals = generateColumns(columns, m, solve_exact=True, start_time=start_time)

    print("\nBranch-and-price:")
    bins_cg = branchAndPrice(m, vals, columns, start_time)

    print(f"\nSolution by column generation: {len(bins_cg)} bins")
    for bin, items in enumerate(bins_cg):
        tt_weight = round(sum(weights[i] for i in items))
        print(f"Bin {bin+1} ({tt_weight} <= {capacity}): {items}")
        assert tt_weight <= capacity

    print(f"Finished in {round(time.perf_counter()-start_time,2): .2f}s")
