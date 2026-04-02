# Example 2
#
# This example compares recursive partitioning methods on the Swiss graph for `k = 8` and `k = 16`.

using SparseArrays, DelimitedFiles, MAT, Plots, PrettyTables
using GraphLab

function benchmark(A::SparseMatrixCSC,
    coords::Matrix,
    k::Int64,
    prefix::String)


    ## Partition methods
    p_coord = GraphLab.recursive_bisection(GraphLab.part_coordinate, k, A, coords)
    p_inertial = GraphLab.recursive_bisection(GraphLab.part_inertial, k, A, coords)
    p_spectral = GraphLab.recursive_bisection(GraphLab.part_spectral, k, A)

    p_metis_rec = GraphLab.part_metis(A, k, :RECURSIVE)
    p_metis_kway = GraphLab.part_metis(A, k, :KWAY)

    ## Save visualizations
    GraphLab.draw_graph(A, coords, p_coord, file_name=prefix * "_coordinate.png")
    GraphLab.draw_graph(A, coords, p_inertial, file_name=prefix * "_inertial.png")
    GraphLab.draw_graph(A, coords, p_spectral, file_name=prefix * "_spectral.png")
    GraphLab.draw_graph(A, coords, p_metis_rec, file_name=prefix * "_metis_rec.png")
    GraphLab.draw_graph(A, coords, p_metis_kway, file_name=prefix * "_metis_kway.png")

    ## Compute edge cuts
    edge_cuts = [
        GraphLab.count_edge_cut(A, p_coord),
        GraphLab.count_edge_cut(A, p_inertial),
        GraphLab.count_edge_cut(A, p_spectral),
        GraphLab.count_edge_cut(A, p_metis_rec),
        GraphLab.count_edge_cut(A, p_metis_kway)
    ]

    ## Compute balance ratios
    balances = [
        GraphLab.compute_partition_balance(p_coord),
        GraphLab.compute_partition_balance(p_inertial),
        GraphLab.compute_partition_balance(p_spectral),
        GraphLab.compute_partition_balance(p_metis_rec),
        GraphLab.compute_partition_balance(p_metis_kway)
    ]

    return edge_cuts, balances
end

## List of .mat files
files = [
    "Swiss_graph.mat"
]

## Directory containing the files
base_dir = "examples/meshes/"

for k in [8, 16]

    final_cuts = Matrix{Any}(undef, length(files), 6)
    final_balances = Matrix{Any}(undef, length(files), 6)

    ## Process each file
    for (i, file_name) in enumerate(files)
        ## Open the .mat file
        file_path = joinpath(base_dir, file_name)
        file = MAT.matopen(file_path)
        data = MAT.read(file)

        ## Extract adjacency matrix and coordinates
        ## note "Swiss_graph.mat" as sliglty different data structure
        if file_name == "Swiss_graph.mat"
            A = sparse(data["CH_adj"])
            coords = data["CH_coords"]
        else
            A = data["Problem"]["A"]
            coords = data["Problem"]["aux"]["coord"]
        end

        ## Symmetrize A
        A = (A + transpose(A)) / 2

        ## Generate a prefix for output files based on the file name
        prefix = "ex2_" * replace(file_name, ".mat" => "")

        ## Run the benchmark function
        println("Processing file: $file_path with prefix: $prefix")
        edge_cuts, balances = benchmark(A, coords, k, prefix)

        final_cuts[i, 1] = file_name
        final_cuts[i, 2] = edge_cuts[1]
        final_cuts[i, 3] = edge_cuts[2]
        final_cuts[i, 4] = edge_cuts[3]
        final_cuts[i, 5] = edge_cuts[4]
        final_cuts[i, 6] = edge_cuts[5]

        final_balances[i, 1] = file_name
        final_balances[i, 2] = round(balances[1], digits=2)
        final_balances[i, 3] = round(balances[2], digits=2)
        final_balances[i, 4] = round(balances[3], digits=2)
        final_balances[i, 5] = round(balances[4], digits=2)
        final_balances[i, 6] = round(balances[5], digits=2)

        ## Close the .mat file
        MAT.close(file)
    end

    ## Pretty print the final results
    header = [
        "Mesh Name",
        "Recursive($k):part_coordinate",
        "Recursive($k):part_inertial",
        "Recursive($k):part_spectral",
        "part_metis_rec($k)",
        "part_metis_kway($k)"
    ]
    pretty_table(final_cuts; column_labels=header)
    pretty_table(final_balances; column_labels=header)
end