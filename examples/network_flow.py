# Example of a shortest path network flow in a graph
# Shows integration of highspy with networkx

import highspy
import networkx as nx

orig, dest = ('A', 'D')

# create directed graph with edge weights (distances)
G = nx.DiGraph()
G.add_weighted_edges_from([('A', 'B', 2.0), ('B', 'C', 3.0), ('A', 'C', 1.5), ('B', 'D', 2.5), ('C', 'D', 1.0)])

h = highspy.Highs()
h.silent()

x = h.addBinaries(G.edges, obj=nx.get_edge_attributes(G, 'weight'))

# add flow conservation constraints
#                        {  1  if n = orig
#   sum(out) - sum(in) = { -1  if n = dest
#                        {  0     otherwise
rhs  = lambda n: 1 if n == orig else -1 if n == dest else 0
flow = lambda E: sum((x[e] for e in E))

h.addConstrs(flow(G.out_edges(n)) - flow(G.in_edges(n)) == rhs(n) for n in G.nodes)
h.minimize()

# Print the solution
print('Shortest path from', orig, 'to', dest, 'is: ', end = '')
sol = h.vals(x)

n = orig
while n != dest:
    print(n, end=' ')
    n = next(e[1] for e in G.out_edges(n) if sol[e] > 0.5)

print(dest)
