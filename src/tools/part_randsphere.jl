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

"""
Stereographically lift points to the unit sphere in (d+1)-dimensional space.

coords: n × d matrix of points
Returns: n × (d+1) matrix on the unit sphere
"""
function _stereoup(coords::Matrix{Float64})
    n, d = size(coords)
    lifted = zeros(n, d+1)

    for i in 1:n
        x = coords[i, :]
        r2 = dot(x, x)
        lifted[i, 1:d] .= 2 * x
        lifted[i, d+1] = r2 - 1
        lifted[i, :] ./= r2 + 1
    end

    return lifted
end


"""
Approximate _centerpoint by coordinate-wise median of a random sample.

Arguments:
- X::Matrix{Float64}: n × d data matrix
- sample_size::Int: number of random rows to sample

Returns:
- center::Vector{Float64}: estimated _centerpoint
"""
function _centerpoint(X::Matrix{Float64}, sample_size::Int)
    n, d = size(X)
    idx = rand(1:n, min(sample_size, n))
    sample = X[idx, :]
    c = mapslices(median, sample; dims=1)[:]
    return c
end


# """
# Conformal map that moves `center` to the origin on the unit sphere.

# Arguments:
# - center::Vector{Float64}: _centerpoint on the sphere (length d+1)
# - X::Matrix{Float64}: n × (d+1) matrix of lifted points

# Returns:
# - Y::Matrix{Float64}: transformed points (same size as X)
# """
# function _conmap(c::Vector{Float64}, xyz::Matrix{Float64})
#     d = size(xyz, 2)
  
#     # Compute reflection and stretch
#     Q, r = _reflector(c)
#     alpha = sqrt((1+r)/(1-r))
  
#     # Reflect
#     xyzref = xyz * Q;
  
#     # Handle north pole in the stereographic projection
#     norths = findall(x -> abs(x - 1) < sqrt(eps()), xyzref[:, d]) 
#     if !isempty(norths)
#       xyzn = xyzref[norths, :]
#       # xyzref[norths, :] = zeros(length(norths), d)
#       xyzref[norths, :] .= 0
#     end
  
#     xyref = _stereodown(xyzref)
  
#     # Stretch
#     xymap = xyref ./ alpha
#     xyzmap = _stereoup(xymap)
  
#     if !isempty(norths)
#       xyzmap[norths,:] = xyzn
#     end

#     println(xyzmap)
  
#     return xyzmap, xymap
# end

function _conmap(c::Vector{Float64}, xyz::Matrix{Float64})
    d = size(xyz, 2)

    # 1. reflection / rotation
    Q, r = _reflector(c)
    alpha = sqrt((1 + r) / (1 - r))   # stretch factor (as in the original code)

    # 2. rotate all points
    xyzref = xyz * Q

    # 3. handle north pole (z == 1) specially to avoid stereographic blowup
    norths = findall(x -> abs(x - 1) < sqrt(eps(Float64)), xyzref[:, d])
    xyzn = xyzref[norths, :]
    if !isempty(norths)
        xyzref[norths, :] .= 0
    end

    # 4. stereographic projection down to R^d
    xyref = _stereodown(xyzref)

    # 5. stretch in the plane
    xymap = xyref ./ alpha

    # 6. lift back to the sphere
    xyzmap = _stereoup(xymap)

    # 7. restore exact north-pole points
    if !isempty(norths)
        xyzmap[norths, :] .= xyzn
    end

    return xyzmap, xymap
end


# """
#     _reflector(c::AbstractVector)

# Constructs a Householder reflection matrix `Q` that maps the input vector `c` (reversed)
# to a multiple of the first basis vector. This is useful for aligning a given direction with
# a coordinate axis in geometric transformations.

# # Arguments
# - `c::AbstractVector`: A real vector of length `d`.

# # Returns
# - `Q::Matrix{Float64}`: A `d × d` orthogonal matrix representing the Householder reflection.
# - `r::Float64`: The leading entry of the reflected vector (i.e., norm or signed component).

# # Notes
# - Internally, this performs a QR decomposition of the reversed vector `c[end:-1:1]`
#   and adjusts the result to restore the original ordering.
# """
# function _reflector(c)
#     d = length(c)
#     p = collect(d:-1:1)
#     Q, r = qr(c[collect(d:-1:1)])
#     Q = Q[p, p]
#     r = r[1]
#     return Q, r
# end

"""
    _reflector(c)

Householder reflection matrix Q such that

    Q * c = t,    where t = (0, ..., 0, r),  r = ‖c‖₂

and Q is orthogonal (Q'Q = I). This matches the "rotate so that the
centerpoint becomes a point (0,...,0,T) on the (d+1)-st axis)" step.

Returns:
- Q :: Matrix{Float64} (d × d)
- r :: Float64          (‖c‖₂)
"""
function _reflector(c::AbstractVector)
    d = length(c)
    r = norm(c)

    if r == 0
        return Matrix{Float64}(I, d, d), 0.0
    end

    t = zeros(Float64, d)
    t[end] = r                       # target vector (0,...,0,r)

    v = c .- t
    nv2 = dot(v, v)

    if nv2 == 0
        # c already aligned with last axis
        return Matrix{Float64}(I, d, d), r
    end

    Q = Matrix{Float64}(I, d, d) .- 2.0 * (v * v') / nv2
    return Q, r
end



"""
    _stereodown(xyz::AbstractMatrix)

Performs the inverse of stereographic projection, mapping points from the unit sphere
in ℝ^{d} (embedded in ℝ^{d+1}) back to Euclidean space ℝ^{d}.

# Arguments
- `xyz::AbstractMatrix`: An `n × (d + 1)` matrix where each row is a point on the unit sphere in ℝ^d.

# Returns
- `xy::Matrix{Float64}`: An `n × d` matrix of projected points in Euclidean space.

# Notes
- Assumes the last coordinate in each row of `xyz` corresponds to the vertical axis (pole).
- This operation is the inverse of the stereographic projection performed in `_stereoup`.
"""
function _stereodown(xyz)
    n, dim = size(xyz)
    xy = xyz[:, 1:dim-1] ./ (1 .- xyz[:, dim])  # broadcasting does the repetition
    
  
    return xy
end


"""
  _sepcircle(A::SparseMatrixCSC, X::Matrix{Float64}, ntrials::Int)
  
Try ntrials random great circle cuts on unit-sphere data.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- X::Matrix{Float64}: n × (d+1) unit-sphere coordinates
- ntrials::Int: number of random directions to try

Returns:
- bestdir::Vector{Float64}: direction vector of best cut
- mincut::Float64: minimum edge cut found
"""
function _sepcircle(A::SparseMatrixCSC, X::Matrix{Float64}, ntrials::Int)
    _, d = size(X)
    M = (X' * X) ^ 2
  
    vv = randn(ntrials, d) * M
    
    bestcut = Inf;
    bestdir = vv[1, :]
  
    for i in 1:size(vv, 1)
      v = vv[i, :]
  
      if norm(v) != 0
        currentcut = _sepquality(v, A, X)
      end
  
      if bestcut > currentcut
        bestcut = currentcut
        bestdir = v
      end
    end
  
    return bestdir, bestcut 
  
end


"""
    _sepquality(v, A, xyz)

Evaluates the quality of a geometric separator defined by a direction vector `v`,
by computing the number of graph edges cut by the resulting _partition.

# Arguments
- `v::Vector{Float64}`: A direction vector defining the separating hyperplane (e.g., great circle).
- `A::SparseMatrixCSC`: The adjacency matrix of the (undirected) graph.
- `xyz::Matrix{Float64}`: An `n × (d+1)` matrix of vertex coordinates on the unit sphere (in ℝ^{d+1}).

# Returns
- `cutsize::Int`: The number of edges crossing between the two sides of the _partition.

# Notes
- Vertices are _partitioned based on the sign of their projection onto `v`.
- The number of crossing edges is computed using nonzero entries in `A[a, b]` and `A[b, a]'`.
"""
function _sepquality(v::Vector{Float64}, A::SparseMatrixCSC, xyz::Matrix{Float64})
    a, b = _partition(xyz, v)
    # cutsize = nnz(Bool.(A[a, b]) .| Bool.(A[b, a]'))
    cutsize = nnz((A[a, b] .!= 0) .| (A[b, a]' .!= 0))
  
    return cutsize
end


"""
  _sepline(A, xy, ntrials)

Try ntrials random straight hyperplanes in Euclidean space.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- coords::Matrix{Float64}: n × d original coordinates
- ntrials::Int: number of directions to try

Returns:
- bestdir::Vector{Float64}: best direction
- mincut::Float64: corresponding edge cut
"""
function _sepline(A::SparseMatrixCSC, xy::Matrix{Float64}, ntrials::Int)
    d = size(xy, 2)
  
    _, S, V = svd(xy);
  
    exponent = 2*(d+1)/(ntrials-1)
    s = diagm(S).^exponent
    W = V * s * V'
    vv = randn(ntrials, d) * W
    rownorms = sqrt.(sum(vv .* vv, dims=2))
    vv = Diagonal(vec(1 ./ rownorms)) * vv
  
    bestcut = Inf
    bestdir = vv[1, :]
    for i in 1:ntrials
      v = vv[i, :]
      cut = _sepquality(v, A, xy)
      if cut < bestcut
        bestcut = cut
        bestdir = v
      end
    end
  
    return bestdir, bestcut 
  end


"""
    part_randsphere(A, coords; ntrials)

Geometric _partitioning using random spheres and lines.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- coords::Matrix{Float64}: n × d node positions
- ntries::Int (optional): number of directions to try (default: 30)

Returns:
- part1, part2: vectors of node indices
"""
function part_randsphere(A::SparseMatrixCSC, coords::Matrix{Float64}; ntrials::Int=30)
    n, d = size(coords)

  # How to split the tries
  nlines = floor(Int, (ntrials / 2) ^ (d / (d + 1))) # number of lines to try (fallback)
  nouter = ceil(Int, log(ntrials - nlines + 1) / log(20)) # number of _centerpoints to try
  ninner = floor(Int, (ntrials - nlines) / nouter) # number of random directions (great circles) to test for each _centerpoint.
  nlines = ntrials - nouter * ninner # Rounding

  # Normalize coordinates: center and scale
  xy_shift = mean(coords, dims=1)
  xy = coords .- xy_shift
  scale = maximum(abs.(xy))
  xy ./= scale

  # Project to sphere
  xyz = _stereoup(xy)

  # Spherical sperator search
  csample = min(n, (d + 3)^4)
  bestcut = Inf
  bestdir = bestcpt = fill(NaN, d + 1)

  for _ in 1:nouter
    cpt = _centerpoint(xyz, csample)
    if 1 - norm(cpt) < sqrt(eps())
      dc = randn(size(cpt))
      cpt = 0.9 * cpt + dc / (20 * norm(dc))
    end
    xyzmap, _ = _conmap(cpt, xyz)
    circledir, circlecut = _sepcircle(A, xyzmap, ninner)
    if circlecut < bestcut
      bestcut = circlecut
      bestdir = circledir
      bestcpt = cpt
    end
  end

  linedir, linecut = _sepline(A, xy, nlines)

  if linecut < bestcut
    bestcut = linecut
    bestdir = linedir
    p1, p2 = _partition(xy, linedir)
  else
    bestmap, _  = _conmap(bestcpt, xyz)
    p1, p2 = _partition(bestmap, bestdir)
  end

  p = ones(Int, n)
  p[p1] .= 1
  p[p2] .= 2

  return p
end


"""
    _partition(X, direction)

Splits a set of points into two parts by projecting them onto a given direction vector
and _partitioning at the median projection value.

# Arguments
- `X::Matrix{Float64}`: An `n × d` matrix of points (each row is a point in ℝ^d or ℝ^{d+1}).
- `direction::Vector{Float64}`: A direction vector along which to project the points.

# Returns
- `part1::Vector{Int}`: Indices of points with projection ≤ median (one side of the cut).
- `part2::Vector{Int}`: Indices of points with projection > median (the other side).

# Notes
- This ensures a balanced _partition by cutting at the median of the projected values.
- Commonly used to implement great-circle or hyperplane separators.
"""
function _partition(X::Matrix{Float64}, direction::Vector{Float64})
    proj = X * direction
    threshold = median(proj)
    part1 = findall(proj .<= threshold)
    part2 = findall(proj .> threshold)
    return part1, part2
end
