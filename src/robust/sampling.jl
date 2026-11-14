################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): sampling the parameter box.         #
#                                                                               #
# The robust oracle maximises the constraint violation over the box. Its      #
# starting points come from a dictionary built once, independent of the       #
# hosting decision:                                                            #
#   1. a modal map of the nominal network (resonance mode analysis, Xu et al.  #
#      2005): the parallel resonances, i.e. the local minima of the            #
#      eigenvalue magnitudes |η_m(h)| of the reduced admittance;               #
#   2. their sensitivities to the parameters (Huang, Cui & Xu 2007), here as   #
#      the shift of the resonance order per latent coordinate;                 #
#   3. seeds per (mode, harmonic): where the resonance crosses the harmonic    #
#      inside the box, the point that minimises |η_m(h)|; elsewhere the vertex #
#      that pushes the resonance towards the harmonic;                         #
#   4. a Sobol background, which also covers series resonance;                 #
# with the constraints of every point stored in compact form.                  #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

import Sobol
import LinearAlgebra: eigen, dot, norm

# ──────────────────────────────────────────────────────────────────────────────
# Reduced admittance and its eigen-decomposition
# ──────────────────────────────────────────────────────────────────────────────

"""
    reduced_admittance(hdata, model, p, h; sequence = :positive) -> (Y, buses)

The harmonic admittance at the parameter point `p` and the (real) harmonic order `h`, with
the clean buses (ideal harmonic grounds) eliminated, and the bus ids of its rows.
"""
function reduced_admittance(hdata::Dict, model::AbstractNetworkModel, p, h::Real;
                            sequence::Symbol = :positive)
    Y, idx = admittance_at(hdata, model, p, h; sequence)
    fund   = hdata["nw"]["1"]["bus"]
    buses  = sort([b["bus_i"] for b in values(fund) if !get(b, "is_clean", false)])
    keep   = [idx[b] for b in buses]
    Yr     = Y[keep, keep]
    sequence === :zero && (Yr = first(_earth_floating(Yr)))
    return Yr, buses
end

# Eigenpairs of the complex symmetric Y, with eigenvectors normalised as tᵀt = 1 (Huang,
# Cui & Xu 2007), so that Y = T·diag(η)·Tᵀ and ∂η/∂p = tᵀ(∂Y/∂p)t.
function _modes(Y::AbstractMatrix)
    F = eigen(Matrix(Y))
    T = similar(F.vectors)
    for m in axes(T, 2)
        t = F.vectors[:, m]
        T[:, m] = t / sqrt(transpose(t) * t)
    end
    return F.values, T
end

# similarity of two eigenvectors, insensitive to their scaling: |tᴴs| / (‖t‖‖s‖)
_overlap(t, s) = abs(dot(t, s)) / (norm(t) * norm(s))

# Eigenvalue and eigenvector of Y that continue the mode with eigenvector `tref`.
function _track(Y, tref)
    η, T = _modes(Y)
    m    = argmax([_overlap(T[:, j], tref) for j in axes(T, 2)])
    return η[m], T[:, m]
end

# Minimum of a parabola through three equally spaced samples (x-step `dx`, centre `x0`).
function _parabola_min(x0, dx, y₋, y₀, y₊)
    den = y₋ - 2y₀ + y₊
    den > 0 || return x0
    return x0 + dx * (y₋ - y₊) / (2den)
end

"""
    Resonance

A parallel resonance of the reduced admittance: mode `mode` of the modal scan, harmonic
order `h` (real) at the local minimum of |η|, the eigenvalue `η` there and the eigenvector
`t` (normalised as tᵀt = 1) over the buses `buses`.
"""
struct Resonance
    mode::Int
    h::Float64
    η::ComplexF64
    t::Vector{ComplexF64}
    buses::Vector{Int}
end

"""
    modal_scan(hdata, model, p; hmin = 1.0, hmax, step = 0.02) -> Vector{Resonance}

Resonance mode analysis (Xu et al. 2005) at the parameter point `p`: the eigenvalues of the
reduced admittance over the continuous band `hmin:step:hmax` (between the harmonics too),
each mode tracked by eigenvector continuity, and the interior local minima of |η_m(h)|,
refined by a parabola. `hmax` defaults to 1.25 times the highest harmonic of `hdata`: a
resonance just above the band can still move onto a harmonic inside the parameter box.
"""
function modal_scan(hdata::Dict, model::AbstractNetworkModel, p;
                    hmin::Real = 1.0, hmax::Real = 1.25 * _max_harmonic(hdata),
                    step::Real = 0.02, nc::NetworkCache = NetworkCache(hdata, model))
    hs = collect(hmin:step:hmax)
    buses = nc.buses
    dec = _tmap(h -> _modes(first(reduced_admittance(nc, p, h))), hs)  # in parallel
    η, T = dec[1]
    nm = length(η)
    E  = zeros(ComplexF64, nm, length(hs))      # η of each tracked mode over hs
    V  = Array{ComplexF64}(undef, size(T, 1), nm, length(hs))
    E[:, 1] = η; V[:, :, 1] = T
    for k in 2:length(hs)                        # continuation by eigenvector overlap: serial
        ηk, Tk = dec[k]
        free = collect(1:nm)
        for m in 1:nm                            # greedy continuation by eigenvector overlap
            j = free[argmax([_overlap(Tk[:, j], V[:, m, k-1]) for j in free])]
            E[m, k] = ηk[j]; V[:, m, k] = Tk[:, j]
            deleteat!(free, findfirst(==(j), free))
        end
    end
    res = Resonance[]
    for m in 1:nm, k in 2:length(hs)-1
        a, b, c = abs(E[m, k-1]), abs(E[m, k]), abs(E[m, k+1])
        (b < a && b <= c) || continue
        push!(res, Resonance(m, _parabola_min(hs[k], step, a, b, c), E[m, k], V[:, m, k], buses))
    end
    return sort!(res; by = r -> r.h)
end

_max_harmonic(hdata::Dict) = maximum(parse(Int, nw) for nw in keys(hdata["nw"]))

# Vertex of the parabola through three points (h1 < h2 < h3), or h2 when they are collinear
function _parabola3(h1, f1, h2, f2, h3, f3)
    num = (h2 - h1)^2 * (f2 - f3) - (h2 - h3)^2 * (f2 - f1)
    den = (h2 - h1) * (f2 - f3) - (h2 - h3) * (f2 - f1)
    return den == 0 ? h2 : h2 - num / (2den)
end

"""
    locate_resonance(hdata, model, p, r; window = 0.5, step = 0.01) -> Resonance

Follow the resonance `r` to the parameter point `p`: the minimum of |η| of the mode whose
eigenvector continues `r.t`, within `r.h ± window`. Starting at `r.h`, the search steps
towards decreasing |η| with doubling steps until the minimum is bracketed, then narrows the
bracket by parabolic interpolation to the width `2step`: typically 6–12 frequencies instead
of a scan of the whole window. Each new frequency continues the mode from the nearest point
already evaluated.
"""
function locate_resonance(hdata::Dict, model::AbstractNetworkModel, p, r::Resonance;
                          window::Real = 0.5, step::Real = 0.01,
                          nc::NetworkCache = NetworkCache(hdata, model))
    ev(h, tref) = _track(first(reduced_admittance(nc, p, h)), tref)
    best(pts) = pts[argmin([abs(q[2]) for q in pts])]
    h0 = r.h
    η0, t0 = ev(h0, r.t)
    ηp, tp = ev(h0 + step, t0)
    ηm, tm = ev(h0 - step, t0)
    pts = [(h0 - step, ηm, tm), (h0, η0, t0), (h0 + step, ηp, tp)]
    if !(abs(η0) <= abs(ηp) && abs(η0) <= abs(ηm))
        dir  = abs(ηp) < abs(ηm) ? 1 : -1
        a, b = pts[2], dir > 0 ? pts[3] : pts[1]          # b improves on a
        s = step
        while true
            s *= 2
            hn = b[1] + dir * s
            if abs(hn - h0) > window || hn <= step         # no minimum in the window: its edge
                return Resonance(r.mode, b[1], b[2], b[3], r.buses)
            end
            ηn, tn = ev(hn, b[3])
            if abs(ηn) >= abs(b[2])
                pts = sort([a, b, (hn, ηn, tn)]; by = first)
                break
            end
            a, b = b, (hn, ηn, tn)
        end
    end
    # narrow the bracket (h1, h2, h3), |η(h2)| smallest, by parabolic interpolation
    for _ in 1:12
        (h1, e1, _), (h2, e2, t2), (h3, e3, _) = pts
        h3 - h1 <= 2step && break
        hv = _parabola3(h1, abs(e1), h2, abs(e2), h3, abs(e3))
        # keep the trial strictly inside, away from h2: else bisect the larger side
        if !(h1 < hv < h3) || abs(hv - h2) < step / 2
            hv = (h2 - h1 > h3 - h2) ? (h1 + h2) / 2 : (h2 + h3) / 2
        end
        ηv, tv = ev(hv, t2)
        v = (hv, ηv, tv)
        if abs(ηv) < abs(e2)
            pts = hv < h2 ? [pts[1], v, pts[2]] : [pts[2], v, pts[3]]
        else
            pts = hv < h2 ? [v, pts[2], pts[3]] : [pts[1], pts[2], v]
        end
    end
    (h1, e1, _), (h2, e2, t2), (h3, e3, _) = pts
    hk = clamp(_parabola3(h1, abs(e1), h2, abs(e2), h3, abs(e3)), h1, h3)
    return Resonance(r.mode, hk, e2, t2, r.buses)
end

# ──────────────────────────────────────────────────────────────────────────────
# Sensitivities
# ──────────────────────────────────────────────────────────────────────────────

"""
    resonance_sensitivity(hdata, model, box, r; δ = 0.25) -> Vector{Float64}

The shift of the resonance order of `r` per latent coordinate of `box`, ∂h_res/∂ξ_i, by
central differences of [`locate_resonance`](@ref) at ξ = ±δ·e_i around the nominal point.
Over the whole box (ξ_i ∈ [-1, 1]) the first-order range of the resonance is
h_res ± Σ_i |∂h_res/∂ξ_i|. (Admittance-based analysis has no resonance-frequency
sensitivity of its own, Huang, Cui & Xu 2007; this is its finite-difference counterpart.)
"""
function resonance_sensitivity(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox,
                               r::Resonance; δ::Real = 0.25,
                               nc::NetworkCache = NetworkCache(hdata, model))
    d = length(box)
    return _tmap(1:d) do i                       # the coordinates are independent
        e = zeros(d); e[i] = δ
        h₊ = locate_resonance(hdata, model, latent(box, e), r; nc).h
        h₋ = locate_resonance(hdata, model, latent(box, -e), r; nc).h
        (h₊ - h₋) / (2δ)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Seeds
# ──────────────────────────────────────────────────────────────────────────────

"""
    Seed

A starting point of the oracle from the modal map: latent point `ξ`, the resonance `r` and
harmonic `h` it targets, whether the resonance `crosses` the harmonic inside the box (first
order), and its `score`: the dominant-mode estimate max_n|t_n| · max_u|t_u| / |η_m(h)| at
`ξ` (excitable from the units u, observable at the buses n).
"""
struct Seed
    ξ::Vector{Float64}
    r::Resonance
    h::Int
    crosses::Bool
    score::Float64
end

# |η| of the mode continuing `t` at (ξ, h), and its eigenvector
function _mode_at(hdata, model, box, ξ, h, t; nc::NetworkCache = NetworkCache(hdata, model))
    return _track(first(reduced_admittance(nc, latent(box, ξ), h;
                                           sequence = harmonic_sequence(h))), t)
end

# Projected gradient descent of |η_m(h; ξ)| over the coordinates `active`, central
# finite differences in parallel, serial backtracking; ξ stays in [-1, 1]. The result
# does not depend on the number of threads.
function _min_eta(hdata, model, box, ξ0, h, t; active, iters = 12, fd = 1e-3,
                  nc::NetworkCache = NetworkCache(hdata, model))
    ξ = copy(ξ0)
    f(x) = abs(first(_mode_at(hdata, model, box, x, h, t; nc)))
    fξ = f(ξ)
    for _ in 1:iters
        g = zeros(length(ξ))
        _tforeach(active) do i
            e = zeros(length(ξ)); e[i] = fd
            g[i] = (f(clamp.(ξ .+ e, -1, 1)) - f(clamp.(ξ .- e, -1, 1))) / (2fd)
        end
        norm(g) > 0 || break
        α, improved = 1.0, false
        while α > 1e-4
            x = clamp.(ξ .- α .* g ./ norm(g), -1, 1)
            fx = f(x)
            if fx < fξ
                ξ, fξ, improved = x, fx, true
                break
            end
            α /= 2
        end
        improved || break
    end
    return ξ
end

"""
    resonance_seeds(hdata, model, box; harmonics, top = 3, units) -> Vector{Seed}

Seeds of the oracle from the modal map of the nominal network, for every pair of parallel
resonance and harmonic in `harmonics` (default: the harmonics of `hdata`):

- **crossing** (the harmonic lies in the first-order range of the resonance over the box):
  the minimum-norm latent step that moves the resonance onto the harmonic, refined by
  minimising |η_m(h)| in the `top` most sensitive coordinates;
- **non-crossing**: the vertex that moves the resonance towards the harmonic, and the face
  point that moves only the `top` most sensitive coordinates.

Each seed is scored by the dominant-mode estimate (see [`Seed`](@ref)) over the buses of
the harmonic units `units` (default: the load buses). Sensitivities are computed first,
with their coordinates in parallel, then the (resonance, harmonic) pairs and the central
differences within each descent step run in parallel. Seeds are collected in resonance
then harmonic order before sorting by score; the result does not depend on the number
of threads.
"""
function resonance_seeds(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox;
                         harmonics = sort([parse(Int, nw) for nw in keys(hdata["nw"]) if nw != "1"]),
                         top::Int = 3, units = nothing)
    d    = length(box)
    nc   = NetworkCache(hdata, model)               # one admittance cache for the whole seeding
    res  = modal_scan(hdata, model, nominal_point(box); nc)
    ubus = units === nothing ?
           unique(ld["load_bus"] for ld in values(hdata["nw"]["1"]["load"])) :
           unique(hdata["nw"]["1"]["load"]["$u"]["load_bus"] for u in units)
    sens = [resonance_sensitivity(hdata, model, box, r; nc) for r in res]
    groups = _tmap([(i, h) for i in eachindex(res) for h in harmonics]) do (i, h)
        r, s  = res[i], sens[i]
        reach = sum(abs, s)
        order = sortperm(abs.(s); rev = true)
        act   = order[1:min(top, d)]
        upos  = [findfirst(==(b), r.buses) for b in ubus if b in r.buses]
        group = Seed[]
        crosses = abs(h - r.h) <= reach
        cands = Vector{Float64}[]
        if crosses
            ξ = norm(s) > 0 ? clamp.(s .* (h - r.h) ./ sum(abs2, s), -1, 1) : zeros(d)
            push!(cands, _min_eta(hdata, model, box, ξ, h, r.t; active = act, nc))
        else
            # the frequency coordinates push the resonance towards h; the coordinates
            # that hardly move it (e.g., resistance: damping, not tuning) are then set by
            # minimising |η_m(h)| with the frequency coordinates held
            dir  = sign(h - r.h)
            tune = [i for i in 1:d if abs(s[i]) > 1e-3 * maximum(abs, s)]
            rest = setdiff(1:d, tune)
            v = zeros(d); v[tune] = dir .* sign.(s[tune])
            isempty(rest) || (v = _min_eta(hdata, model, box, v, h, r.t; active = rest, nc))
            push!(cands, v)
            f = zeros(d); f[act] = dir .* sign.(s[act]); push!(cands, f)
        end
        for ξ in unique(cands)
            η, t = _mode_at(hdata, model, box, ξ, h, r.t; nc)
            sc = isempty(upos) ? 0.0 : maximum(abs, t) * maximum(abs.(t[upos])) / abs(η)
            push!(group, Seed(ξ, r, h, crosses, sc))
        end
        return group
    end
    seeds = reduce(vcat, groups; init = Seed[])
    return sort!(seeds; by = s -> -s.score)
end

# ──────────────────────────────────────────────────────────────────────────────
# Background
# ──────────────────────────────────────────────────────────────────────────────

"""
    background_points(box, n; seed = 1) -> Vector{Vector{Float64}}

`n` latent points of a Sobol sequence over [-1, 1]^d, randomised by one uniform shift modulo
1 (Cranley–Patterson rotation) drawn from `seed`, so different seeds give independent sets
for the a-posteriori check.
"""
function background_points(box::ParameterBox, n::Int; seed::Integer = 1)
    d   = length(box)
    seq = Sobol.SobolSeq(d)
    sh  = rand(Random.MersenneTwister(seed), d)
    return [2 .* mod.(Sobol.next!(seq) .+ sh, 1.0) .- 1 for _ in 1:n]
end

# ──────────────────────────────────────────────────────────────────────────────
# Dictionary
# ──────────────────────────────────────────────────────────────────────────────

# Second moments grouped by the number of harmonics of each constraint.
struct _SecondMoments
    single::Vector{Float64}
    full::Array{Float64,3}
    subset::Vector{Float64}
end

const _EMPTY_SECOND_MOMENTS = Float64[]
const _EMPTY_FULL_MOMENTS = zeros(0, 0, 0)

# Shared indices of one constraint in the point's means and second moments.
struct _ConstraintBlock
    n::Int                            # number of harmonics
    position::Int                     # position of a single harmonic
    μ::UnitRange{Int}
    M::Int                            # index in its second-moment group
end

"""
    ParameterDictionary

The points of the oracle with their constraints in compact form. For a constraint `c` and
the hosting variables g (a vector over `harmonics`), E[S] = μᵀg and Var[S] = gᵀMg with
M[h,h'] = Σ_{k≥1} ψ_k·a_{h,k}·a_{h',k}, so

    F_c(g) = μᵀg + λ·√(gᵀMg) ≤ limit.

Each point stores |H_c| means and |H_c|² second moments per constraint. The means
are flat blocks in `μ`; `M` groups single-harmonic scalars, full-harmonic matrices
in a column-major array, and the remaining blocks in a flat vector. `positions[j]`
records the positions of its harmonics in `harmonics`; `blocks[j]` records its mean
range and second-moment index, shared by all points. `keys[j]` is `(kind, key)` of
constraint j, in the order of
[`shhc_constraints`](@ref); `labels[s]` is `:nominal`, `:vertex`, `:corner`, `:seed`, `:background`
or `:refined` (added by the oracle). Points are appended with [`add_point!`](@ref).
"""
struct ParameterDictionary
    box::ParameterBox
    harmonics::Vector{Int}
    keys::Vector{Tuple{Symbol,Any}}
    limit::Vector{Float64}
    λ::Vector{Float64}
    points::Vector{Vector{Float64}}
    labels::Vector{Symbol}
    positions::Vector{Vector{Int}}     # per constraint: positions in harmonics
    blocks::Vector{_ConstraintBlock}   # per constraint: indices in μ and M
    sizes::NTuple{3,Int}               # single scalars, full matrices, subset entries
    μ::Vector{Vector{Float64}}         # per point: means of each constraint
    M::Vector{_SecondMoments}          # per point: grouped second-moment blocks
end

# Ranges and group indices are shared by every point in the dictionary.
function _dictionary_layout(H, positions)
    blocks = _ConstraintBlock[]
    start = 1
    single, full, subset = 0, 0, 0
    for pos in positions
        n = length(pos)
        if n == 1
            single += 1; idx = single
        elseif n == length(H)
            full += 1; idx = full
        else
            idx = subset + 1; subset += n^2
        end
        push!(blocks, _ConstraintBlock(n, n == 1 ? only(pos) : 0, start:start+n-1, idx))
        start += n
    end
    return blocks, (single, full, subset)
end

# Compact means in flat constraint blocks and second moments grouped by harmonic count.
# Full-harmonic blocks retain the original three-dimensional array views for dot.
# The harmonic positions, mean ranges and second-moment indices are shared by all points.
function _compact(cons, moments, H, positions, blocks, sizes)
    (; K, psi_norms) = moments
    μ = zeros(isempty(blocks) ? 0 : last(last(blocks).μ))
    M = _SecondMoments(sizes[1] == 0 ? _EMPTY_SECOND_MOMENTS : zeros(sizes[1]),
                       sizes[2] == 0 ? _EMPTY_FULL_MOMENTS : zeros(length(H), length(H), sizes[2]),
                       sizes[3] == 0 ? _EMPTY_SECOND_MOMENTS : zeros(sizes[3]))
    for (j, c) in enumerate(cons)
        pos = positions[j]
        block = blocks[j]
        n = block.n
        for (b, i2) in enumerate(pos)
            B = c.coef[H[i2]]
            μ[first(block.μ) + b - 1] = B[1]
            for (a, i) in enumerate(pos)
                A = c.coef[H[i]]
                value = sum(psi_norms[k+1] * A[k+1] * B[k+1]
                            for k in K if k != 0; init = 0.0)
                if n == 1
                    M.single[block.M] = value
                elseif n == length(H)
                    M.full[a, b, block.M] = value
                else
                    M.subset[block.M + a - 1 + (b - 1) * n] = value
                end
            end
        end
    end
    return μ, M
end

# compact constraints at the latent point ξ, checked against the dictionary's constraint set
function _compact_point(dict::ParameterDictionary, hdata, model, ξ)
    m = moments_at(hdata, model, latent(dict.box, ξ))
    c = shhc_constraints(hdata, m)
    length(c) == length(dict.keys) ||
        error("add_point!: the constraint set changed at ξ = $ξ")
    for (j, x) in enumerate(c)
        x.kind === first(dict.keys[j]) && x.key == last(dict.keys[j]) ||
            error("add_point!: the constraint set changed at ξ = $ξ")
        pos = dict.positions[j]
        length(x.coef) == length(pos) && all(haskey(x.coef, dict.harmonics[i]) for i in pos) ||
            error("add_point!: the constraint harmonics changed at ξ = $ξ")
    end
    return _compact(c, m, dict.harmonics, dict.positions, dict.blocks, dict.sizes)
end

"""
    add_point!(dict, hdata, model, ξ, label) -> Int

Evaluate the constraints at the latent point `ξ` and append it to `dict` with `label`;
returns its index.
"""
function add_point!(dict::ParameterDictionary, hdata::Dict, model::AbstractNetworkModel,
                    ξ::AbstractVector, label::Symbol)
    return only(add_points!(dict, hdata, model, [ξ], label))
end

"""
    add_points!(dict, hdata, model, points, label) -> Vector{Int}

Evaluate the constraints at the latent `points` in parallel and append them to `dict`, in
the order given, with `label`; returns their indices.
"""
function add_points!(dict::ParameterDictionary, hdata::Dict, model::AbstractNetworkModel,
                     points::AbstractVector, label::Symbol)
    data = _tmap(ξ -> _compact_point(dict, hdata, model, ξ), points)
    idx  = Int[]
    for (ξ, (μ, M)) in zip(points, data)
        push!(dict.points, collect(Float64, ξ)); push!(dict.labels, label)
        push!(dict.μ, μ); push!(dict.M, M)
        push!(idx, length(dict.points))
    end
    return idx
end

# Sign-predicted corners of the mean coefficients, independent of the hosting decision.
# Central differences at ξ = 0 use δ = 1e-4; sensitivities at most 1e-6 times the
# largest magnitude for each (constraint, harmonic) are treated as zero.
function _mean_corners(dict::ParameterDictionary, hdata, model)
    d = length(dict.box)
    δ = 1e-4
    sensitivities = _tmap(1:d) do j
        e = zeros(d); e[j] = δ
        μ₊, _ = _compact_point(dict, hdata, model, e)
        μ₋, _ = _compact_point(dict, hdata, model, -e)
        (μ₊ - μ₋) / (2δ)
    end
    points = Vector{Float64}[]
    for c in eachindex(dict.keys), h in eachindex(dict.harmonics)
        dict.limit[c] > 0.0 || continue
        a = findfirst(==(h), dict.positions[c])
        s = a === nothing ? zeros(d) :
            [sensitivities[j][first(dict.blocks[c].μ) + a - 1] for j in 1:d]
        threshold = 1e-6 * maximum(abs, s; init = 0.0)
        push!(points, [abs(x) <= threshold ? 0.0 : sign(x) for x in s])
    end
    return unique(points)
end

"""
    build_dictionary(hdata, model, box; seeds = Seed[], nbackground = 1024,
                     vertices = length(box) <= 10, corners = !vertices, seed = 1) -> ParameterDictionary

Evaluate the constraints of the robust program at the nominal point, the vertices (when
`vertices`), the sign-predicted corners (when `corners`), the `seeds` and `nbackground` Sobol
points, and store them in compact form. By default, use the vertices for d ≤ 10 and the
sign-predicted corners otherwise. Worst cases of these constraints can sit at corners,
and random points in high dimension rarely reach them.

The corners follow the signs of each mean coefficient μ[h, c] for positive-limit
constraints, independent of the hosting decision. Central differences at ξ = 0 use
δ = 1e-4 in latent coordinates; magnitudes at most 1e-6 times the largest sensitivity
of each (constraint, harmonic) are treated as zero and keep that coordinate nominal.
Duplicates are removed before evaluation: 2d extra compact evaluations plus one per
distinct corner. `hdata` must carry the PCE data (Cantelli or Gaussian bound).
"""
function build_dictionary(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox;
                          seeds = Seed[], nbackground::Int = 1024,
                          vertices::Bool = length(box) <= 10, corners::Bool = !vertices,
                          seed::Integer = 1)
    d  = length(box)
    m0 = moments_at(hdata, model, nominal_point(box))
    c0 = shhc_constraints(hdata, m0)
    H   = sort(m0.harm_nw_ints)
    pos = [findall(h -> haskey(c.coef, h), H) for c in c0]
    blocks, sizes = _dictionary_layout(H, pos)
    dict = ParameterDictionary(box, H, [(c.kind, c.key) for c in c0],
                               [c.limit for c in c0], [c.λ for c in c0],
                               Vector{Float64}[], Symbol[], pos, blocks, sizes,
                               Vector{Float64}[], _SecondMoments[])
    add_point!(dict, hdata, model, zeros(d), :nominal)
    if vertices
        add_points!(dict, hdata, model,
                    [collect(v) for v in Iterators.product(ntuple(_ -> (-1.0, 1.0), d)...)][:],
                    :vertex)
    end
    corners && add_points!(dict, hdata, model, _mean_corners(dict, hdata, model), :corner)
    add_points!(dict, hdata, model, [s.ξ for s in seeds], :seed)
    add_points!(dict, hdata, model, background_points(box, nbackground; seed), :background)
    return dict
end

# Dense arrays preserve the reduction order of harmonic subsets.
function _subset_moments(dict::ParameterDictionary, j, μ, M, gv)
    μfull = zeros(length(gv))
    Mfull = zeros(length(gv), length(gv), 1)
    pos = dict.positions[j]
    block = dict.blocks[j]
    n   = block.n
    for (b, i2) in enumerate(pos)
        μfull[i2] = μ[first(block.μ) + b - 1]
        for (a, i) in enumerate(pos)
            Mfull[i, i2, 1] = M.subset[block.M + a - 1 + (b - 1) * n]
        end
    end
    return dot(μfull, gv), dict.λ[j] <= 0.0 ? 0.0 : dot(gv, view(Mfull, :, :, 1), gv)
end

# normalized violation of constraint j from its compact form at one point
@inline function _compact_violation(dict::ParameterDictionary, j::Int, μ, M, gv)
    block = dict.blocks[j]
    n   = block.n
    if n == 1
        g  = gv[block.position]
        μg = μ[first(block.μ)] * g
        F  = dict.λ[j] <= 0.0 ? μg :
             μg + dict.λ[j] * sqrt(max((M.single[block.M] * g) * g, 0.0))
    elseif n == length(gv)
        μg = dot(view(μ, block.μ), gv)
        F  = dict.λ[j] <= 0.0 ? μg :
             μg + dict.λ[j] * sqrt(max(dot(gv, view(M.full, :, :, block.M), gv), 0.0))
    else
        μg, σ² = _subset_moments(dict, j, μ, M, gv)
        F = dict.λ[j] <= 0.0 ? μg : μg + dict.λ[j] * sqrt(max(σ², 0.0))
    end
    L = dict.limit[j]
    return L > 0.0 ? F / L - 1.0 : F - L
end

"""
    dictionary_violations(dict, g) -> Matrix{Float64}

Normalized violation F_c(g)/limit − 1 of every constraint (rows) at every point of `dict`
(columns), for the hosting variables `g` (h ⇒ g_h = I_h²). Constraints with a non-positive
limit return F_c(g) − limit, as [`normalized_violation`](@ref).
"""
function dictionary_violations(dict::ParameterDictionary, g)
    gv = [g[h] for h in dict.harmonics]
    V  = Matrix{Float64}(undef, length(dict.keys), length(dict.points))
    return _dictionary_violations!(V, dict, gv)
end

# Fill the violations with a decision vector shared by the scan and oracle updates.
@inline function _dictionary_violations!(V, dict::ParameterDictionary, gv)
    for s in eachindex(dict.points)
        μ, M = dict.μ[s], dict.M[s]
        @inbounds for j in eachindex(dict.keys)
            V[j, s] = _compact_violation(dict, j, μ, M, gv)
        end
    end
    return V
end
