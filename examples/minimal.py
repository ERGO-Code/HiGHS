import highspy
import time
import networkx as nx
from random import randint

h = highspy.Highs()

(x1, x2) = h.addVariables(2, lb = -h.inf)

h.addConstrs(x2 - x1 >= 2, 
             x1 + x2 >= 0)

h.minimize(x2)



h = highspy.Highs()


G = nx.circular_ladder_graph(5).to_directed()
nx.set_edge_attributes(G, {e: {'weight': randint(1, 9)} for e in G.edges})

d = h.addBinaries(G.edges, obj=nx.get_edge_attributes(G, 'weight'))

h.addConstrs(sum(d[e] for e in G.in_edges(i)) - sum(d[e] for e in G.out_edges(i)) == 0 for i in G.nodes)

h = highspy.Highs()

ts = time.time()
perf1 = [h.addBinary() for _ in range(1000000)]
t1 = time.time() - ts
print(t1)

h = highspy.Highs()

ts = time.time()
perf2 = h.addVariables(1000000)
t2 = time.time() - ts
print(t2)

