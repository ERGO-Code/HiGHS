# An example branch-and-price algorithm for the bin packing problem
#
from collections import defaultdict
from operator import itemgetter
import highspy, random, time, math

# Relative gap for column generation
CG_GAP_REL = 1e-4

# Random seeds for reproducibility
random.seed(100)
SEED = 100

# Instance parameters
NumberItems = 40
ItemWeights = [round(random.uniform(1, 10), 1) for _ in range(NumberItems)]
BinCapacity = 15


#
# Solve instance with greedy first fit decreasing (FFD) heuristic
#
def solveGreedyModel():
    bins = defaultdict(float)
    solution = defaultdict(list)

    for item, w in sorted(enumerate(ItemWeights), reverse=True, key=itemgetter(1)):
        index = next((i for i, W in bins.items() if W + w <= BinCapacity), len(bins))
        bins[index] += w
        solution[index].append(item)
    
    return list(solution.values())

#
# The compact bin-packing model (for comparison with branch-and-price)
#
# min sum_{j} y_{j}
# s.t. 
#   sum_{j} x_{ij} == 1, \forall i                 // Each item i is assigned to exactly one bin
#   sum_{i} w_{i} * x_{ij} <= c * y_{j}, \forall j // Bin capacity constraint
#
#   x_{ij} = 1 if item i is assigned to bin j, 0 otherwise
#   y_{j}  = 1 if bin j is open, 0 otherwise
#
def solveCompactModel():
    # Estimate maximum number of bins with first fit decreasing (for smaller sized model)
    greedy_bins = solveGreedyModel()
    B = len(greedy_bins)

    m = highspy.Highs()
    m.setOptionValue('output_flag', True)
    m.setOptionValue('mip_abs_gap', 1-1e-5) # Since all variables in the objective are binary, it should stop when the absolute gap is below 1
    m.setOptionValue('random_seed', SEED)

    # add binary variables x_{ij} and y_{j}
    x = {(i,j): i*B + j for i in range(NumberItems) for j in range(B) }
    y = [len(x) + j for j in range(B)]

    m.addVars(len(x), [0]*len(x), [1]*len(x));  # x_{ij} \in {0,1}, \forall i,j
    m.addVars(len(y), [0]*len(y), [1]*len(y));  # y_{j}  \in {0,1}, \forall j
    m.changeColsIntegrality(
        len(x)+len(y),
        list(range(len(x)+len(y))),
        [highspy.HighsVarType.kInteger]*(len(x)+len(y))
    )

    # min sum_{j} y_{j}
    m.changeObjectiveSense(highspy.ObjSense.kMinimize)
    m.changeColsCost(len(y), y, [1]*len(y))
    
    #      \sum_{j} x_{ij} == 1, \foreach i 
    # 1 <= \sum_{j} x_{ij} <= 1
    for i in range(NumberItems):
        m.addRow(
            1,                                              # lhs
            1,                                              # rhs
            B,                                              # Number of non-zero variables
            [x[i,j] for j in range(B)],                     # Indexes of variable
            [1] * B                                         # Coefficients of variables
        )

    #         sum_{i} w_{i} * x_{ij} <= c * y_{j}, \foreach j
    # -inf <= sum_{i} w_{i} * x_{ij}  - c * y_{j} <= 0
    for j in range(B):
        m.addRow(
            -highspy.kHighsInf,                                             # lhs
            0,                                                              # rhs
            NumberItems+1,                                                  # Number of non-zero variables 
            [x[i,j] for i in range(NumberItems)] + [y[j]],                  # Indexes of variable
            [ItemWeights[i] for i in range(NumberItems)] + [-BinCapacity]   # Coefficients of variables  
        )
        
    m.run()
    vals = list(m.getSolution().col_value)

    bins = [
        [i for i in range(NumberItems) if vals[x[i,j]] > 0.9]
        for j in range(B)
        if vals[y[j]] > 0.9
    ]
    
    return bins

#
# Master problem for column generation
# A column represents a set of items packed into a bin
#
# min \sum_{k} \lambda_{k}
# s.t. 
#   \sum_{k \in K_i} \lambda_{k} == 1, \foreach i  (\mu_i duals)  // Each item i is packed exactly once (but may be spread over multiple columns)
#                                                                 //   where K_i is the set of columns where i appears
#   0 <= \lambda_{k} <= 1, \foreach k                             // Can use each column at most once (can be fractional)
#
def createMasterProblem(columns: list):
    m = highspy.Highs()
    m.setOptionValue('output_flag', False)
    m.setOptionValue('random_seed', SEED)
    
    m.addVars(len(columns), [0]*len(columns), [1]*len(columns))
    
    # min \sum_{k} \lambda_{k}
    m.changeObjectiveSense(highspy.ObjSense.kMinimize)
    m.changeColsCost(len(columns), list(range(len(columns))), [1]*len(columns))
    
    #      \sum_{k \in K_i} \lambda_{k} == 1, \foreach i
    # 1 <= \sum_{k \in K_i} \lambda_{k} <= 1
    for i in range(NumberItems):
        K_i = [k for k, column in enumerate(columns) if i in column]
        m.addRow(1, 1, len(K_i), K_i, [1]*len(K_i))

    m.run() 

    solution = m.getSolution()
    vals = list(solution.col_value)
    duals = list(solution.row_dual)

    return m, vals, duals

#
# Solve the knapsack subproblem exactly with HiGHS 
#
# 1 - max \sum_{i} \mu_{i} * z_{i}
# s.t.
#   \sum_{i} w_{i} * z_{i} <= capacity
#
#   z_{i} = 1 if item i is packed, 0 otherwise
#
def solveSubproblemExact(duals):
    m = highspy.Highs()
    m.setOptionValue('output_flag', False)
    m.setOptionValue('random_seed', SEED)

    m.addVars(NumberItems, [0]*NumberItems, [1]*NumberItems)
    m.changeColsIntegrality(NumberItems, list(range(NumberItems)), [highspy.HighsVarType.kInteger]*NumberItems)

    # max \sum_{i} \mu_{i} * z_{i}
    # where \mu_{i} is the dual variable for the i-th row of the master problem
    m.changeColsCost(NumberItems, list(range(NumberItems)), duals)
    m.changeObjectiveSense(highspy.ObjSense.kMaximize)

    #         sum_{i} w_{i} * z_{i} <= capacity, \foreach j
    # -inf <= sum_{i} w_{i} * z_{i} <= capacity
    m.addRow(
        -highspy.kHighsInf,
        BinCapacity,
        NumberItems,
        list(range(NumberItems)),
        ItemWeights
    )
    
    m.run()

    vals = list(m.getSolution().col_value)
    new_column = sorted([i for i in range(NumberItems) if vals[i] > 0.9])

    return 1 - m.getObjectiveValue(), new_column


# Solve the knapsack subproblem with greedy heuristic
def solveSubproblemNotExact(duals):
    total_weight = 0
    new_column = []
    
    for i in sorted(range(NumberItems), key=lambda i: -duals[i]/ItemWeights[i]):
        if duals[i] >= 0 and ItemWeights[i] + total_weight <= BinCapacity:
            total_weight += ItemWeights[i]
            new_column += [i]
    
    return 1 - sum(duals[i] for i in new_column), sorted(new_column)


#
# Generate columns for the master problem
#
def generateColumns(columns: list, m, msg=True, solve_exact = False, start_time = 0):
    best_gap = math.inf

    iter = 0
    while True:
        solution = m.getSolution()
        duals = list(solution.row_dual)
        ub = m.getObjectiveValue()

        # solve sub problem to generate new column
        if not solve_exact:
            obj_sub, new_column = solveSubproblemNotExact(duals)
        else:
            obj_sub, new_column = solveSubproblemExact(duals)

        iter += 1     
        obj = ub + min(0, obj_sub)  # we terminate if obj_sub >= 0
        gap = (ub-obj)/ub
        
        if gap < best_gap:
            best_gap = gap
            if msg:
                row = [
                    f"{len(columns)}",
                    f"{round(obj_sub, 3) :.4f}",
                    f'{round(obj, 3) :.4f}',
                    f"{round(ub, 3) :.4f}",
                    f"{max(0, gap) :.3%}", 
                    f"{round(time.perf_counter()-start_time, 2) :.2f}"
                ]
                if iter == 1:
                    header = ["Columns", 'Pricing', 'Obj', 'UB', 'gap', 'Time']
                    row_format = "".join([
                        "{:>"+str(max(len(row_el), len(header_el))+3)+"}"
                        for row_el, header_el in zip(row, header)
                    ])
                    print(row_format.format(*header))

                print(row_format.format(*row))

        # terminate if no column with good reduced cost is found
        if obj_sub >= 0 or gap < CG_GAP_REL:
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

    return m, columns, list(solution.col_value)


#
# Branching helper structures/functions
#
class Node:
    layer: int
    assigned_columns: list[int]

    value: float
    final_columns: list[int]
    final_columns_vals: list[float]
    final_fractional_columns: list[int]

    def __init__(self, layer: int, assigned_columns: list[int], parent = None):
        self.layer = layer
        self.parent = parent
        self.assigned_columns = assigned_columns


def getFractionalColumns(vals: list, columns: list):
    return sorted(
        [k for k, val in enumerate(vals) if abs(round(val) - val) > 1e-6],
        key=lambda k: -len(columns[k]) # order by most fractional (largest number of items)
    )


def solveNode(node: Node, columns: list[list[int]], m):
    # "enforce" the branch constraints, i.e., selected columns in the current node
    lower_bounds = [int(k in node.assigned_columns) for k in range(len(columns))]

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
        m, columns, vals = generateColumns(columns, m, msg=False, solve_exact=False) # Generate columns quickly
        m, columns, vals = generateColumns(columns, m, msg=False, solve_exact=True)  # Prove optimality on node

        node.value = sum(vals)
        node.final_columns = columns
        node.final_columns_vals = vals
        node.final_fractional_columns = getFractionalColumns(vals, columns)

        # ensure no column already assigned to 1 in the node is returned back with a fractional value
        assert len(set(node.final_fractional_columns).intersection(set(node.assigned_columns))) == 0
        return True, m, node
    
    return False, m, node


#
# Branch-and-price algorithm
#
# Note: Branch selection has been chosen for simple implementation and it's ability to avoid symmetry. 
# This comes at the cost of an unbalanced search tree (other approaches may solve the problem with fewer branches), 
# and it is less likely to be generalizable to other problems.
#
# Specifically, the code "up-branches" columns with fractional value in the RMP solution, i.e., forces specific columns 
# to be selected in RMP.  The subproblems don't need to explicitly enforce this constraint (unlike other branching strategies).
#
def branchAndPrice(m, vals, columns, start_time=0):    
    header = ["NodesExpl", "TreeSize", "CurrCols", "FracCols", "UB", "Time"]
    row_format ="{:>12}" * (len(header))
    print(row_format.format(*header))

    root_node = Node(-1, [])
    root_node.value = sum(vals)
    root_node.final_columns = columns
    root_node.final_columns_vals = vals
    root_node.final_fractional_columns = getFractionalColumns(vals, columns)

    # create initial branches for each column with fractional value in RMP solution
    branch_tree: list[Node] = [Node(0, [column_idx]) for column_idx in root_node.final_fractional_columns]
    branch_tree_dict = {frozenset(node.assigned_columns): node for node in branch_tree}

    # RMP lower bound can be rounded up to nearest integer (avoiding floating point precision issues)
    # so, take integer floor, add one if fractional part greater than tolerance
    rmp_LB = int(root_node.value) + int(root_node.value % 1 > 1e-6)
    best_obj = math.inf
    best_node = None
    best_count_frac_columns = math.inf
    count_nodes_visited = 0

    # if no fractional columns, then root node is an optimal integer solution
    if len(root_node.final_fractional_columns) == 0:
        best_node = root_node
        print(row_format.format(*[0, 0, len(columns), 0, best_node.value, f"{round(time.perf_counter()-start_time, 2) :.2f}" ]))

    # explore branches, looking for optimal integer solution
    while len(branch_tree) > 0:
        # Choose next node in branch to evaluate.  Prioritize nodes that are likely to be 
        # integer (DFS with fewer fractional columns), and likely to be optimal (larger value).
        #
        # i.e., maximize branch depth (layer), minimize number of fractional columns of its parent, 
        # and maximize value of its original fractional columns.
        node = min(branch_tree, key=lambda node: (
            -node.layer,
            len(node.parent.final_fractional_columns) if node.parent is not None else math.inf,
            -node.parent.final_columns_vals[node.assigned_columns[-1]] if node.parent is not None else 0
        ))
        branch_tree.remove(node)

        # Solves the node with column generation
        mp_is_feasible, m, node = solveNode(node, columns, m)
        count_nodes_visited += 1

        # Only add columns if the master problem is feasible
        # Prune nodes that cannot lead to optimal solution        
        if mp_is_feasible and node.value < best_obj:
            columns = node.final_columns            
            count_fractional_columns = len(node.final_fractional_columns)

            # if no fractional columns, then this is an integer solution
            # update the best solution here for console output
            if count_fractional_columns == 0 and node.value < best_obj:
                best_obj = int(round(node.value, 0))  # avoid float precision issues
                best_node = node

            # found a less fractional solution or enough time has elapsed
            if count_fractional_columns < best_count_frac_columns or time.perf_counter()-start_time > 5:
                best_count_frac_columns = count_fractional_columns
                print(row_format.format(*[
                    count_nodes_visited,
                    len(branch_tree),
                    len(columns),
                    best_count_frac_columns,
                    best_obj if best_obj < math.inf else "-",
                    f"{round(time.perf_counter()-start_time, 2) :.2f}"
                ]))

            # add all fractional columns as new branches
            if count_fractional_columns > 0:
                for column_idx in node.final_fractional_columns:
                    new_assigned_columns = node.assigned_columns + [column_idx]
                    key = frozenset(new_assigned_columns)
                    
                    if key not in branch_tree_dict:
                        new_node = Node(node.layer + 1, new_assigned_columns, node)
                        branch_tree = [new_node] + branch_tree
                        branch_tree_dict[key] = new_node

            # if no fractional columns, then this is an integer solution
            # stop if we've found a provably optimal solution (using lower bound from RMP root node)
            elif rmp_LB == best_obj:
                break

    # explored the entire tree, so best found solution is optimal
    if best_node is not None:
        return [column for idx, column in enumerate(best_node.final_columns) if best_node.final_columns_vals[idx] > 0.9]
    
    raise Exception("It should not get to this point")


if __name__ == '__main__':
    if max(ItemWeights) > BinCapacity:
        print(f"Instance is infeasible: item {max(enumerate(ItemWeights), key=itemgetter(1))[0]} has weight {max(ItemWeights)} > {BinCapacity} (bin capacity).")
        exit()

    start_time = time.perf_counter()
    greedy_bins = solveGreedyModel()
    greedy_time = time.perf_counter()-start_time
    print(f"Greedy estimate: {len(greedy_bins)} bins")
    print(f"Finished in {greedy_time: .2f}s\n")

    start_time = time.perf_counter()
    compact_bins = solveCompactModel()
    compact_time = time.perf_counter()-start_time

    print(f"\nSolution by compact model: {len(compact_bins)} bins")
    for bin, items in enumerate(compact_bins):
        tt_weight = round(sum(ItemWeights[i] for i in items))
        print(f"Bin {bin+1} ({tt_weight} <= {BinCapacity}): {items}")
        assert tt_weight <= BinCapacity

    print(f"Finished in {compact_time: .2f}s\n")

    # start with initial set of columns for feasible master problem
    start_time = time.perf_counter()
    columns = [[i] for i in range(NumberItems)] + greedy_bins
    m, vals, duals = createMasterProblem(columns)
    
    print("\nSolving root node:")
    m, columns, vals = generateColumns(columns, m, solve_exact=False, start_time=start_time)
    
    print("\nProving optimality on root node:")
    m, columns, vals = generateColumns(columns, m, solve_exact=True, start_time=start_time)

    print("\nBranch-and-price:")
    cg_bins = branchAndPrice(m, vals, columns, start_time)
    cg_time = time.perf_counter()-start_time

    print(f"\nSolution by column generation: {len(cg_bins)} bins")
    for bin, items in enumerate(cg_bins):
        tt_weight = round(sum(ItemWeights[i] for i in items))
        print(f"Bin {bin+1} ({tt_weight} <= {BinCapacity}): {items}")
        assert tt_weight <= BinCapacity

    print(f"Finished in {cg_time: .2f}s\n")

    print(f"Greedy : {len(greedy_bins)} bins, {greedy_time:6.2f}s")
    print(f"Compact: {len(compact_bins)} bins, {compact_time:6.2f}s")
    print(f"ColGen : {len(cg_bins)} bins, {cg_time:6.2f}s\n")

