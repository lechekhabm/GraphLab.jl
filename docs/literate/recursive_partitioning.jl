# # Recursive partitioning
#
# This example shows how to partition a graph into more than two parts using
# recursive bisection.
#
# Recursive bisection repeatedly applies a bisection method until the requested
# number of parts is reached.

using GraphLab

#-
# ## Build a graph

A, coords = GraphLab.grid_graph(10, 50, π / 3);

#-
# ## Partition into four parts
#
# Here we use spectral bisection recursively to obtain `k = 4` parts.

k = 4
p = GraphLab.recursive_bisection(GraphLab.part_spectral, k, A);

#-
# ## Evaluate the partition

edge_cut = GraphLab.count_edge_cut(A, p)
balance = GraphLab.compute_partition_balance(p)

(edge_cut=edge_cut, balance=balance)

#-
# ## Visualize the recursive partition

GraphLab.draw_graph(A, coords, p)

