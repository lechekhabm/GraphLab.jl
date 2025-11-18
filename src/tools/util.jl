using Graphs, Artifacts, MAT, SparseArrays, LinearAlgebra, CSV, NearestNeighbors


"""
    count_edge_cut(A::AbstractMatrix, p::AbstractVector)

Count the number of edges that cross partitions in the graph `A`.

# Arguments
- `A::AbstractMatrix`: Adjacency matrix of the graph.
- `p::AbstractVector`: Partition vector where `p[v]` represents the partition of vertex `v`.

# Returns
- The number of edges that connect nodes in different partitions.

# Example
```julia-repl
julia> count_edge_cut(A, p)
15
```
"""
function count_edge_cut(A::AbstractMatrix, p::AbstractVector)
    # Convert adjacency matrix to Graph object
    g = Graphs.SimpleGraph(A)

    # Initialize the edge cut count
    edge_cut = 0

    # Iterate over edges
    for e in Graphs.edges(g)
        u, v = Graphs.src(e), Graphs.dst(e)
        if p[u] != p[v]  # Check if edge crosses partitions
            edge_cut += 1
        end
    end

    return edge_cut
end


"""
    ratio_cut(A::AbstractMatrix, p::AbstractVector{<:Integer}) -> Float64

Compute the Ratio Cut:
    Rcut(p) = sum_j cut(p_j, ¯p_j) / |p_j|

- `A`: (0/1) adjacency (dense or sparse). Nonzeros are treated as edges.
- `p`: partition labels (length == number of vertices). Labels can be any integers.

Counts each inter-cluster edge once for each endpoint's cluster (as per definition).
"""
function ratio_cut(A::AbstractMatrix, p::AbstractVector{<:Integer})
    g = SimpleGraph(A)
    nv(g) == length(p) || throw(ArgumentError("length(p) must equal number of vertices"))

    # Map labels -> contiguous indices 1..k
    labels = collect(unique(p))
    k = length(labels)
    idx_of = Dict(l => i for (i, l) in enumerate(labels))

    # Sizes and boundary counts per cluster
    sizes = zeros(Int, k)
    for l in p
        sizes[idx_of[l]] += 1
    end
    any(==(0), sizes) && throw(ArgumentError("Empty clusters are not allowed"))

    boundary = zeros(Int, k)

    # Accumulate boundary counts (each crossing edge adds 1 to both endpoint clusters)
    for e in edges(g)
        u, v = src(e), dst(e)
        lu, lv = p[u], p[v]
        if lu != lv
            iu, iv = idx_of[lu], idx_of[lv]
            boundary[iu] += 1
            boundary[iv] += 1
        end
    end

    # Sum_j cut(p_j, ¯p_j)/|p_j|
    return sum(boundary[i] / sizes[i] for i in 1:k)
end


using Graphs

"""
    normalized_cut(A::AbstractMatrix, p::AbstractVector{<:Integer}) -> Float64

Compute Ncut(p) = ∑_j cut(p_j, ¯p_j) / vol(p_j), where vol(p_j) = ∑_{v∈p_j} degree(v).

Assumes an undirected, unweighted adjacency (nonzeros -> edges).
"""
function normalized_cut(A::AbstractMatrix, p::AbstractVector{<:Integer})
    g = SimpleGraph(A)
    nv(g) == length(p) || throw(ArgumentError("length(p) must equal number of vertices"))

    # relabel cluster ids to 1..k
    labs = collect(unique(p))
    k = length(labs)
    idx_of = Dict(l => i for (i,l) in enumerate(labs))

    # volumes: sum of degrees per cluster
    vol = zeros(Int, k)
    @inbounds for v in 1:nv(g)
        vol[idx_of[p[v]]] += degree(g, v)
    end
    any(vol .== 0) && throw(ArgumentError("Cluster with zero volume (all isolated nodes) not allowed"))

    # boundary cuts per cluster (each crossing edge adds 1 to both endpoint clusters)
    bnd = zeros(Int, k)
    @inbounds for e in edges(g)
        u, v = src(e), dst(e)
        lu, lv = p[u], p[v]
        if lu != lv
            iu, iv = idx_of[lu], idx_of[lv]
            bnd[iu] += 1
            bnd[iv] += 1
        end
    end

    return sum(bnd[i] / vol[i] for i in 1:k)
end




"""
    build_adjacency(edges::Matrix{Int}, num_nodes::Int)

Construct the adjacency matrix of a graph from an edge list.

# Arguments
- `edges::Matrix`: Matrix where each row represents an edge `[u, v]`.
- `num_nodes::Int`: Total number of nodes in the graph.

# Returns
- A symmetric sparse adjacency matrix (`SparseMatrixCSC{Int, Int}`).

# Example
```julia-repl
julia> edges = [1 2; 2 3; 3 1]
julia> A = build_adjacency(edges, 3)
```
"""
function build_adjacency(edges::Matrix{Int}, num_nodes::Int)
    adjacency_matrix = spzeros(Int, num_nodes, num_nodes)
    for edge in eachrow(edges)
        adjacency_matrix[edge[1], edge[2]] = 1
        adjacency_matrix[edge[2], edge[1]] = 1  # Ensure symmetry
    end
    return adjacency_matrix
end


"""
    build_adjacency(type::String)

Generate a predefined adjacency matrix and corresponding node coordinates.

# Arguments
- `type::String`: Type of graph to generate.  
  - `"network"`: A predefined network structure.  
  - `"triangles"`: A small triangular mesh structure.

# Returns
- `A::SparseMatrixCSC`: The sparse adjacency matrix of the graph.
- `coords::Matrix`: Node coordinates for visualization.

# Example
```julia-repl
julia> A, coords = build_adjacency("network")
```
"""
function build_adjacency(type::String)
    if type=="network"
        A = [[0. 1. 1. 0. 0. 1. 0. 0. 1. 1.]
            [1. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
            [1. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
            [0. 0. 0. 0. 1. 1. 0. 0. 0. 0.]
            [0. 0. 0. 1. 0. 1. 0. 0. 0. 0.]
            [1. 0. 0. 1. 1. 0. 1. 1. 0. 0.]
            [0. 0. 0. 0. 0. 1. 0. 1. 0. 0.]
            [0. 0. 0. 0. 0. 1. 1. 0. 0. 0.]
            [1. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
            [1. 0. 0. 0. 0. 0. 0. 0. 1. 0.]]
        A = sparse(A)
        x = [3.0 1.0 2.0 4.0 5.0 3.0 2.0 1.0 5.0 4.0]
        y = [2.0 1.0 1.0 5.0 5.0 4.0 5.0 5.0 1.0 1.0]
    elseif type=="triangles"
        A = [   [0. 1. 1. 1. 0. 0. ]
                [1. 0. 1. 0. 1. 0. ]
                [1. 1. 0. 0. 0. 1. ]
                [1. 0. 0. 0. 1. 1. ]
                [0. 1. 0. 1. 0. 1. ]
                [0. 0. 1. 1. 1. 0. ]]
        A = sparse(A)
        x = [2.0 4.0 4.0 0.0 6.0 6.0]
        y = [3.0 2.0 4.0 3.0 0.0 6.0]
    end
    
    coords = hcat(vec(x), vec(y))

    return A, coords
end


"""
    compute_partition_balance(p::AbstractVector) -> Float64

Computes the balance metric of a given graph partitioning.

# Parameters
- `p::AbstractVector`: A vector where `p[i]` represents the partition index assigned to vertex `i`.

# Returns
- `Float64`: The balance metric, defined as the ratio of the largest partition size to the ideal partition size.
  A value close to `1.0` indicates a well-balanced partitioning, while higher values suggest imbalance.

# Example
```julia-repl
julia> p = [1, 1, 2, 2, 2, 3, 3, 3, 3]  # Example partitioning
julia> balance = compute_partition_balance(p)
julia> println(balance)  # Output close to 1 for balanced partitions
```
""" 
function compute_partition_balance(p::AbstractVector)
    # Count the number of vertices in each partition
    partition_sizes = Dict{eltype(p), Int}()
    
    for part in p
        partition_sizes[part] = get(partition_sizes, part, 0) + 1
    end

    # Compute balance metric
    n = length(p)
    k = length(partition_sizes)
    ideal_size = n / k
    max_size = maximum(values(partition_sizes))

    return max_size / ideal_size  # Balance metric (1 is ideal)
end


"""
    grid_graph(n::Int, m::Int, α::Float64)

Returns the adjacency matrix `A::SparseMatrixCSC` and the coordinates `coords::Matrix{Float64}`
of an `n × m` grid graph rotated by angle `α` (in radians).

Vertices are ordered row-wise: vertex `i,j` has index `(i-1)*m + j`.
"""
function grid_graph(n::Int, m::Int, α::Float64)
    N = n * m
    A = spzeros(N, N)
    coords = zeros(2, N)

    # Generate original (unrotated) coordinates and edges
    for i in 1:n
        for j in 1:m
            idx = (i-1)*m + j
            x, y = j, n - i + 1
            coords[:, idx] = [x, y]
            # Right neighbor
            if j < m
                neighbor_idx = idx + 1
                A[idx, neighbor_idx] = 1
                A[neighbor_idx, idx] = 1
            end
            # Bottom neighbor
            if i < n
                neighbor_idx = idx + m
                A[idx, neighbor_idx] = 1
                A[neighbor_idx, idx] = 1
            end
        end
    end

    # Apply rotation by α
    R = [cos(α) -sin(α); sin(α) cos(α)]
    coords .= R * coords

    return A, Matrix(coords')
end


"""
    _partition(coords::Matrix, v::Vector)

Compute a partition based on `coords` using a direction vector `v`.

# Arguments
- `coords::Matrix`: Node coordinates in a 2D space.
- `v::Vector`: Direction vector defining the partitioning line.

# Returns
- A tuple of two vectors: indices of nodes in each partition.

# Example
```julia-repl
julia> _partition(coords, [0, 1])
([1, 2, 3, 9, 10], [4, 5, 6, 7, 8])
```
"""
function _partition(coords::Matrix, v::Vector)
    n, d = size(coords)
    
    v = v[:]
    dotprod = coords * v
    split = median(dotprod)
    a = findall(x -> x < split, dotprod)
    b = findall(x -> x >= split, dotprod)
    c = findall(x -> x == split, dotprod)
    nc = length(c)
    # nc = 0
    
    if nc != 0
        na = length(a)
        nca = Int64(max(ceil(n/2)-na, 0))
        nca = Int64(min(nca, nc))

        if nca > 0
            a = [a; c[1:nca]]
        end
        if nca < nc
            b = [b; c[nca + 1: nc]]
        end
    end

    return a, b
end


"""
    airfoil() -> (A::SparseMatrixCSC, coords::Matrix)

Load the SuiteSparse airfoil example from the `airfoil1` artifact.
"""
function airfoil()
    dir  = artifact"airfoil1"
    file = joinpath(dir, "airfoil", "airfoil1.mat")
    d      = MAT.matread(file)
    A      = sparse(d["Problem"]["A"])
    coords = Matrix(d["Problem"]["aux"]["coord"])
    return A, coords
end

"""
    swiss() -> (A::SparseMatrixCSC, coords::Matrix)

Load the Swiss road graph from the `swiss_graph` artifact.
"""
function swiss()
    dir  = artifact"swiss_graph"                    # ← name in Artifacts.toml
    file = joinpath(dir, "swiss", "Swiss_graph.mat")# ← path inside tarball
    d = MAT.matread(file)
    A = sparse(d["CH_adj"])
    A = (A + transpose(A)) / 2
    coords = Matrix(d["CH_coords"])
    return A, coords
end


"""
    france() -> (A::SparseMatrixCSC, coords::Matrix)

Load the France graph from the `france_graph` artifact.
"""
function france()
    dir  = artifact"france_graph"
    file = joinpath(dir, "france", "france_graph.mat")
    d = MAT.matread(file)
    A = sparse(d["A"])
    A = (A + transpose(A)) / 2
    coords = Matrix(d["coords"])
    return A, coords
end


load(name::Symbol) = name === :airfoil ? airfoil() :
                     name === :swiss   ? swiss()   :
                     name === :france   ? france()   :
                     error("Unknown dataset :$name. Available: :airfoil, :swiss")


load(name::AbstractString) = load(Symbol(name))

function build_knn_adjacency(csvfile::AbstractString; k::Integer=5, header::Bool=false)    
    # 1. Load coordinates from CSV
    # assume columns: x, y (and maybe more dims)
    tbl = CSV.File(csvfile; header=header)

    # turn into coords::Matrix{Float64} of size (n, d)
    coords = (reduce(hcat, (collect(col) for col in tbl)) |> Matrix{Float64})

    d, n = size(coords)  # n points, d dims

    # 2. Build k-d tree (NearestNeighbors expects d×n)
    tree = KDTree(coords)  # d×n

    # 3. Query k+1 neighbors (includes the point itself)
    idxs, dists = knn(tree, coords, k)

    # 4. Build sparse adjacency matrix A (undirected, unweighted = 1)
    I = Int[]
    J = Int[]
    V = Int[]

    for i in 1:n
        nbrs = idxs[i]
        for j in nbrs          
            if j == i
                continue
            end         # neighbor index
            push!(I, i); push!(J, j); push!(V, 1)
            push!(I, j); push!(J, i); push!(V, 1)  # make symmetric
        end
    end

    A = sparse(I, J, V, n, n)

    # (Optional) remove multi-edges >1 to 1.0
    # A .= (A .> 0.0)

    coords = Matrix(coords')

    return A, coords
end
