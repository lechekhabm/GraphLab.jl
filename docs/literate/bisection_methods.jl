# # Comparing bisection methods
#
# This example compares several graph bisection methods on the same small graph.
#
# We use a small generated graph instead of loading many mesh files. This keeps
# the example focused on the main ideas.

using GraphLab

#-
# ## Build a test graph

A, coords = GraphLab.grid_graph(10, 50, π / 3);

#-
# ## Run different bisection methods
#
# We compute four different partitions:
#
# - coordinate bisection,
# - inertial bisection,
# - spectral bisection,
# - METIS recursive bisection.

p_coord = GraphLab.part_coordinate(A, coords);
p_inertial = GraphLab.part_inertial(A, coords);
p_spectral = GraphLab.part_spectral(A);
p_metis = GraphLab.part_metis(A, 2, :RECURSIVE);

#-
# ## Compare edge cuts
#
# A smaller edge cut means that fewer edges cross between the two parts.

edge_cuts = (
    coordinate=GraphLab.count_edge_cut(A, p_coord),
    inertial=GraphLab.count_edge_cut(A, p_inertial),
    spectral=GraphLab.count_edge_cut(A, p_spectral),
    metis=GraphLab.count_edge_cut(A, p_metis),
)

#-
# ## Compare balance
#
# The balance measures how evenly the vertices are distributed between the two
# parts.

balances = (
    coordinate=GraphLab.compute_partition_balance(p_coord),
    inertial=GraphLab.compute_partition_balance(p_inertial),
    spectral=GraphLab.compute_partition_balance(p_spectral),
    metis=GraphLab.compute_partition_balance(p_metis),
)

#-
# ## Summary

(edge_cuts=edge_cuts, balances=balances)