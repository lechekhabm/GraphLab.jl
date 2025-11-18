# ------------------------------------------------------------------------
# Based on work by:
# - John R. Gilbert (UC Santa Barbara): http://www.cs.ucsb.edu/~gilbert/
# - Shang-Hua Teng (University of Southern California): http://www-bcf.usc.edu/~shanghua/
# - Yingzhou Li (Fudan University): https://yingzhouli.com
#
# Note:
# The following acknowledgement is from Yingzhou Li's original codebase:
# "Thanks to Tim Davis for updating the toolbox to MATLAB 5"
# Source: https://github.com/YingzhouLi/meshpart/

using Graphs, AMD, GraphsMatching, JuMP, Cbc


"""
    wait_for_key(prompt)

Prints a prompt to the standard output and waits for a single keypress from the user.

# Arguments
- `prompt::String`: The message to display before waiting for input.

# Returns
- `nothing`.
"""
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)


"""
    _get_maximum_matching(g::Graph)

Computes a maximum weight matching of the graph `g` using a silent CBC solver.

# Arguments
- `g::Graph`: The input graph.

# Returns
- `match`: The matching object containing mate assignments and total weight.

# Notes
- Assumes unit weights unless otherwise specified.
"""
function _get_maximum_matching(g::Graph)
  match = maximum_weight_matching(g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0, MOI.Silent() => true));
  # match.weight ≈ 1

  return match
end

"""
    _vtxsep(gborder, p1border, p2border, match)

Computes a vertex separator based on alternating paths in a bipartite border graph.

# Arguments
- `gborder::Graph`: Bipartite graph between border nodes.
- `p1border::Vector{Int}`: Border nodes from partition 1.
- `p2border::Vector{Int}`: Border nodes from partition 2.
- `match`: Maximum matching result on `gborder`.

# Returns
- `(sep, new_p1border, new_p2border)`: 
  - `sep`: Separator set of vertices.
  - `new_p1border`: Updated partition 1 border after separator removal.
  - `new_p2border`: Updated partition 2 border after separator removal.

# Notes
- Traverses alternating paths between matching and non-matching edges to construct the separator.
"""
function _vtxsep(gborder, p1border, p2border, match)
  offset = length(p1border)
  total = nv(gborder)
  global_border_nodes = vcat(p1border, p2border)

  visited = falses(total)
  queue = Vector{Tuple{Int, Bool}}()  # (node, phase): false = non-matching edge, true = matching edge

  # Start from unmatched p1 nodes (partition A), phase=false (expect non-matching)
  for i in 1:offset
      if match.mate[i] == -1
          push!(queue, (i, false))
          visited[i] = true
      end
  end

  while !isempty(queue)
      (u, phase) = popfirst!(queue)
      for v in neighbors(gborder, u)
          if !visited[v]
              if phase == false
                  if match.mate[u] != v  # non-matching edge
                      visited[v] = true
                      push!(queue, (v, true))
                  end
              else
                  if match.mate[u] == v  # matching edge
                      visited[v] = true
                      push!(queue, (v, false))
                  end
              end
          end
      end
  end

  # Build separator using the constructive proof
  sepp1 = Set{Int}()
  sepp2 = Set{Int}()

  for i in 1:offset
      if match.mate[i] != -1 && !visited[i]
          push!(sepp1, p1border[i])
      end
  end

  for j in 1:(total - offset)
      v = offset + j
      if visited[v]
          push!(sepp2, p2border[j])
      end
  end

  sep = collect(union(sepp1, sepp2))
  return sep, setdiff(p1border, sepp1), setdiff(p2border, sepp2)
  
end


"""
    nested_dissection(A::SparseMatrixCSC, method::Function; coords=nothing, minsep=5, verbose=false)

Computes a nested dissection ordering of the sparse matrix `A` using a given separator method.

# Arguments
- `A::SparseMatrixCSC`: The sparse matrix.
- `method::Function`: Function to compute graph separators.
- `coords::Union{Matrix, Nothing}`: Optional node coordinates for geometric methods (default is `nothing`).
- `minsep::Int`: Minimum component size to stop recursion (default is 10).
- `verbose::Bool`: Whether to display visualizations and pause for input during recursion (default is `false`).

# Returns
- `pA::Vector{Int}`: The computed permutation vector.
"""
function nested_dissection(A::SparseMatrixCSC, method::Function; coords::Union{Matrix, Nothing}=nothing, minsep::Int=10, verbose::Bool=false)
    n = size(A, 1)
    p = Vector{Int}(undef, n)

    pA = _nested_dissection(A, method; coords = coords, minsep = minsep, verbose = verbose)

    return pA
end


"""
    _nested_dissection(A::SparseMatrixCSC, method::Function; coords=nothing, minsep, verbose)

Recursively computes a nested dissection ordering of a sparse matrix, handling disconnected components separately.

# Arguments
- `A::SparseMatrixCSC`: The adjacency matrix of the graph to partition.
- `method::Function`: Function used to compute separators.
- `coords::Union{Matrix, Nothing}`: Optional node coordinates if the separator method needs them.
- `minsep::Int`: Minimum number of nodes to stop recursion.
- `verbose::Bool`: Whether to visualize partitions and wait for user input during recursion.

# Returns
- `perm::Vector{Int}`: The final permutation vector.

# Notes
- Uses approximate minimum degree ordering for small components.
"""
function _nested_dissection(A::SparseMatrixCSC, method::Function; coords::Union{Matrix, Nothing}=nothing, minsep::Int, verbose::Bool)
    n = size(A, 1)
    perm = zeros(Int, n)

    # Find connected components
    g = SimpleGraph(A)
    components = connected_components(g)
    comp_count = size(components)[1]
    offset = 0

    numbered = 0
    for (i, component) in enumerate(components)

      nC = length(component)
      A_sub = A[component, component]
      if nC <= minsep
        pA = symamd(SparseMatrixCSC{Float64, Int64}(A_sub)) # Approximate Minimum Degree ordering

      else
        sub_coords = coords === nothing ? nothing : coords[component, :]
        if sub_coords === nothing
            p_sub = method(A_sub)
        else
            p_sub = method(A_sub, sub_coords)
        end
        p_sub1 = findall(==(1), p_sub)
        p_sub2 = findall(==(2), p_sub)

        # Find border nodes between p1 and p2
        Abord = A_sub[p_sub1, p_sub2]
        rows, cols, _ = findnz(Abord)
        p1border = unique(p_sub1[rows])
        p2border = unique(p_sub2[cols])

        n1 = length(p1border)
        n2 = length(p2border)

        # Build the graph of border vertices
        # gborder = Graph(offset + length(p2border))
        gborder = Graph(n1 + n2)

        # for i in 1:offset
        #     for j in 1:length(p2border)
        #         if A_sub[p1border[i], p2border[j]] != 0
        #             add_edge!(gborder, i, offset + j)
        #         end
        #     end
        # end

        for (i, v1) in enumerate(p1border)       # i in 1:n1
            for (j, v2) in enumerate(p2border)   # j in 1:n2
                if A_sub[v1, v2] != 0
                    add_edge!(gborder, i, n1 + j)
                end
            end
        end

        # Compute maximum matching
        match =  _get_maximum_matching(gborder)

        sep,_,_ = _vtxsep(gborder, p1border, p2border, match);

        vtx1 = setdiff(p_sub1, sep)
        vtx2 = setdiff(p_sub2, sep)
        
        if verbose
          display(draw_graph(A_sub, sub_coords, p_sub))
          wait_for_key("Press Enter to continue...")
        end
        
        q1 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx1, vtx1], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx1, vtx1], method; coords=sub_coords[vtx1, :], minsep=minsep, verbose=verbose)

        q2 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx2, vtx2], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx2, vtx2], method; coords=sub_coords[vtx2, :], minsep=minsep, verbose=verbose)
        
        if q1 === nothing
          q1 = []
        else
          q1 = vec(Int.(q1))
        end
    
        if q2 === nothing
          q2 = []
        else
          q2 = vec(Int.(q2))
        end

        pA = vcat(vtx1[q1], vtx2[q2], sep)

        # Fallback: if the computed permutation has the wrong size, warn and return the identity permutation
        if length(pA) != size(A_sub, 1)
            @warn "Wrong permutation length" length_pA = length(pA) expected = size(A_sub, 1)
            wait_for_key("Press Enter to continue...")
            return collect(1:size(A_sub, 1))  # Identity permutation
        end
      end
    perm[numbered + 1:numbered + nC] = component[pA]
    numbered = numbered + nC
end
  return perm
end