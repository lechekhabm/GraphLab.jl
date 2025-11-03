using Pkg
Pkg.activate(@__DIR__)

using SparseArrays,DelimitedFiles, MAT, Plots, PrettyTables
using GraphLab

"""
    benchmark(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, coords::Matrix{Int64}, prefix::String)

Run various graph partitioning methods on `A` and save visualizations.

# Arguments
- `A::SparseMatrixCSC{Float64, Int64}`: Adjacency matrix of the graph.
- `coords::Matrix{Int64}`: Node coordinates for visualization.
- `prefix::String`: Filename prefix for saved images.

# Methods Used
- Coordinate-based (`part_coordinate`)
- Inertial (`part_inertial`)
- Spectral (`part_spectral`)
- METIS recursive bisection (`part_metis`)

# Output
- Saves partitioned graph visualizations as PNG files.
"""
function benchmark( A::SparseMatrixCSC, 
                    coords::Matrix, 
                    prefix::String)
    # Partition methods
    p_coord    = GraphLab.part_coordinate(A, coords)
    p_inertial = GraphLab.part_inertial(A, coords)
    p_spectral = GraphLab.part_spectral(A)
    p_metis    = GraphLab.part_metis(A, 2, :RECURSIVE)

    # Save visualizations
    GraphLab.draw_graph(A, coords, p_coord, file_name=prefix * "_coordinate.png")
    GraphLab.draw_graph(A, coords, p_inertial, file_name=prefix * "_inertial.png")
    GraphLab.draw_graph(A, coords, p_spectral, file_name=prefix * "_spectral.png")
    GraphLab.draw_graph(A, coords, p_metis, file_name=prefix * "_metis.png")

    # Compute edge cuts
    edge_cuts = [
        GraphLab.count_edge_cut(A, p_coord),
        GraphLab.count_edge_cut(A, p_inertial),
        GraphLab.count_edge_cut(A, p_spectral),
        GraphLab.count_edge_cut(A, p_metis)
    ]

    norm_cuts = [
        GraphLab.normalized_cut(A, p_coord),
        GraphLab.normalized_cut(A, p_inertial),
        GraphLab.normalized_cut(A, p_spectral),
        GraphLab.normalized_cut(A, p_metis)
    ]

    ratio_cuts = [
        GraphLab.ratio_cut(A, p_coord),
        GraphLab.ratio_cut(A, p_inertial),
        GraphLab.ratio_cut(A, p_spectral),
        GraphLab.ratio_cut(A, p_metis)
    ]

    # Compute balance ratios
    balances = [
        GraphLab.compute_partition_balance(p_coord),
        GraphLab.compute_partition_balance(p_inertial),
        GraphLab.compute_partition_balance(p_spectral),
        GraphLab.compute_partition_balance(p_metis)
    ]

    return edge_cuts, norm_cuts, ratio_cuts, balances
end

# List of .mat files
files = [
    "3elt.mat",
    "airfoil1.mat",
    "barth4.mat",
    "crack.mat",
    "mesh1e1.mat",
    "mesh2e1.mat",
    "mesh3e1.mat",
    "netz4504_dual.mat",
    "stufe.mat",
    "ukerbe1.mat"
]

# Directory containing the files
base_dir = "meshes/"

final_results = Matrix{Any}(undef,length(files),17)

# Process each file
for (i, file_name) in enumerate(files)
    # Open the .mat file
    file_path = joinpath(base_dir, file_name)
    file = MAT.matopen(file_path)
    data = MAT.read(file)

    # Extract adjacency matrix and coordinates
    # note "Swiss_graph.mat" as sliglty different data structure
    if file_name == "Swiss_graph.mat" 
        A = sparse(data["CH_adj"])
        coords = data["CH_coords"]
    else
        A = data["Problem"]["A"]
        coords = data["Problem"]["aux"]["coord"]
    end

    # Symmetrize A
    A = (A + transpose(A))/2

    # Generate a prefix for output files based on the file name
    prefix = "ex1_" * replace(file_name, ".mat" => "")

    # Run the benchmark function
    println("Processing file: $file_path with prefix: $prefix")
    edge_cuts, norm_cuts, ratio_cuts, balances = benchmark(A, coords, prefix)

    final_results[i, 1] = file_name
    final_results[i, 2] = edge_cuts[1]
    final_results[i, 3] = norm_cuts[1]
    final_results[i, 4] = ratio_cuts[1]
    final_results[i, 5] = balances[1]  # Balance for part_coordinate
    final_results[i, 6] = edge_cuts[2]
    final_results[i, 7] = norm_cuts[2]
    final_results[i, 8] = ratio_cuts[2]
    final_results[i, 9] = balances[2]  # Balance for part_inertial
    final_results[i, 10] = edge_cuts[3]
    final_results[i, 11] = norm_cuts[3]
    final_results[i, 12] = ratio_cuts[3]
    final_results[i, 13] = balances[3]  # Balance for part_spectral
    final_results[i, 14] = edge_cuts[4]
    final_results[i, 15] = norm_cuts[4]
    final_results[i, 16] = ratio_cuts[4]
    final_results[i, 17] = balances[4]  # Balance for part_metis
   
    # Close the .mat file
    MAT.close(file)
end

# Pretty print the final results
header = [
    "Mesh Name",
    "Edge Cut (Coord)", "Norm. Cut (Coord)", "Ratio Cut (Coord)", "Balance (Coord)",
    "Edge Cut (Inertial)", "Norm. Cut (Inertial)", "Ratio Cut (Inertial)", "Balance (Inertial)",
    "Edge Cut (Spectral)", "Norm. Cut (Spectral)", "Ratio Cut (Spectral)", "Balance (Spectral)",
    "Edge Cut (METIS)", "Norm. Cut (METIS)", "Ratio Cut (METIS)", "Balance (METIS)"
]

pretty_table(final_results, header=header)