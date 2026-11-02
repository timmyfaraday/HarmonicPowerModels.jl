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
    return Y[keep, keep], buses
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
                    step::Real = 0.02)
    hs = collect(hmin:step:hmax)
    Y, buses = reduced_admittance(hdata, model, p, hs[1])
    η, T = _modes(Y)
    nm = length(η)
    E  = zeros(ComplexF64, nm, length(hs))      # η of each tracked mode over hs
    V  = Array{ComplexF64}(undef, size(T, 1), nm, length(hs))
    E[:, 1] = η; V[:, :, 1] = T
    for k in 2:length(hs)
        ηk, Tk = _modes(first(reduced_admittance(hdata, model, p, hs[k])))
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

"""
    locate_resonance(hdata, model, p, r; window = 0.5, step = 0.01) -> Resonance

Follow the resonance `r` to the parameter point `p`: the minimum of |η| of the mode whose
eigenvector continues `r.t`, searched within `r.h ± window`.
"""
function locate_resonance(hdata::Dict, model::AbstractNetworkModel, p, r::Resonance;
                          window::Real = 0.5, step::Real = 0.01)
    hs = collect(max(r.h - window, step):step:r.h + window)
    tr = r.t
    ηs = ComplexF64[]; ts = Vector{ComplexF64}[]
    for h in hs
        η, t = _track(first(reduced_admittance(hdata, model, p, h)), tr)
        push!(ηs, η); push!(ts, t); tr = t
    end
    k = argmin(abs.(ηs))
    hk = 1 < k < length(hs) ? _parabola_min(hs[k], step, abs(ηs[k-1]), abs(ηs[k]), abs(ηs[k+1])) : hs[k]
    return Resonance(r.mode, hk, ηs[k], ts[k], r.buses)
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
                               r::Resonance; δ::Real = 0.25)
    d = length(box)
    s = zeros(d)
    for i in 1:d
        e = zeros(d); e[i] = δ
        h₊ = locate_resonance(hdata, model, latent(box, e), r).h
        h₋ = locate_resonance(hdata, model, latent(box, -e), r).h
        s[i] = (h₊ - h₋) / (2δ)
    end
    return s
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
function _mode_at(hdata, model, box, ξ, h, t)
    return _track(first(reduced_admittance(hdata, model, latent(box, ξ), h;
                                           sequence = harmonic_sequence(h))), t)
end

# Projected gradient descent of |η_m(h; ξ)| over the coordinates `active`, central
# finite differences, backtracking; ξ stays in [-1, 1].
function _min_eta(hdata, model, box, ξ0, h, t; active, iters = 12, fd = 1e-3)
    ξ = copy(ξ0)
    f(x) = abs(first(_mode_at(hdata, model, box, x, h, t)))
    fξ = f(ξ)
    for _ in 1:iters
        g = zeros(length(ξ))
        for i in active
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
the harmonic units `units` (default: the load buses).
"""
function resonance_seeds(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox;
                         harmonics = sort([parse(Int, nw) for nw in keys(hdata["nw"]) if nw != "1"]),
                         top::Int = 3, units = nothing)
    d    = length(box)
    res  = modal_scan(hdata, model, nominal_point(box))
    ubus = units === nothing ?
           unique(ld["load_bus"] for ld in values(hdata["nw"]["1"]["load"])) :
           unique(hdata["nw"]["1"]["load"]["$u"]["load_bus"] for u in units)
    seeds = Seed[]
    for r in res
        s     = resonance_sensitivity(hdata, model, box, r)
        reach = sum(abs, s)
        order = sortperm(abs.(s); rev = true)
        act   = order[1:min(top, d)]
        upos  = [findfirst(==(b), r.buses) for b in ubus if b in r.buses]
        for h in harmonics
            crosses = abs(h - r.h) <= reach
            cands = Vector{Float64}[]
            if crosses
                ξ = norm(s) > 0 ? clamp.(s .* (h - r.h) ./ sum(abs2, s), -1, 1) : zeros(d)
                push!(cands, _min_eta(hdata, model, box, ξ, h, r.t; active = act))
            else
                # the frequency coordinates push the resonance towards h; the coordinates
                # that hardly move it (e.g., resistance: damping, not tuning) are then set by
                # minimising |η_m(h)| with the frequency coordinates held
                dir  = sign(h - r.h)
                tune = [i for i in 1:d if abs(s[i]) > 1e-3 * maximum(abs, s)]
                rest = setdiff(1:d, tune)
                v = zeros(d); v[tune] = dir .* sign.(s[tune])
                isempty(rest) || (v = _min_eta(hdata, model, box, v, h, r.t; active = rest))
                push!(cands, v)
                f = zeros(d); f[act] = dir .* sign.(s[act]); push!(cands, f)
            end
            for ξ in unique(cands)
                η, t = _mode_at(hdata, model, box, ξ, h, r.t)
                sc = isempty(upos) ? 0.0 : maximum(abs, t) * maximum(abs.(t[upos])) / abs(η)
                push!(seeds, Seed(ξ, r, h, crosses, sc))
            end
        end
    end
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

"""
    ParameterDictionary

The points of the oracle with their constraints in compact form. For a constraint `c` and
the hosting variables g (a vector over `harmonics`), E[S] = μᵀg and Var[S] = gᵀMg with
M[h,h'] = Σ_{k≥1} ψ_k·a_{h,k}·a_{h',k}, so

    F_c(g) = μᵀg + λ·√(gᵀMg) ≤ limit,

which needs |H| + |H|² numbers per constraint and point instead of the full PCE
coefficients. `keys[j]` is `(kind, key)` of constraint j, in the order of
[`shhc_constraints`](@ref); `labels[s]` is `:nominal`, `:vertex`, `:seed`, `:background`
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
    μ::Vector{Matrix{Float64}}          # per point: (|H|, constraints)
    M::Vector{Array{Float64,3}}         # per point: (|H|, |H|, constraints)
end

# compact (μ, M) of the constraints at one point, in the harmonic order `H`
function _compact(cons, moments, H)
    (; K, psi_norms) = moments
    nc = length(cons)
    μ  = zeros(length(H), nc)
    M  = zeros(length(H), length(H), nc)
    for (j, c) in enumerate(cons)
        for (a, h) in enumerate(H)
            haskey(c.coef, h) || continue
            μ[a, j] = c.coef[h][1]
            for (b, h2) in enumerate(H)
                haskey(c.coef, h2) || continue
                M[a, b, j] = sum(psi_norms[k+1] * c.coef[h][k+1] * c.coef[h2][k+1]
                                 for k in K if k != 0; init = 0.0)
            end
        end
    end
    return μ, M
end

"""
    add_point!(dict, hdata, model, ξ, label) -> Int

Evaluate the constraints at the latent point `ξ` and append it to `dict` with `label`;
returns its index.
"""
function add_point!(dict::ParameterDictionary, hdata::Dict, model::AbstractNetworkModel,
                    ξ::AbstractVector, label::Symbol)
    m = moments_at(hdata, model, latent(dict.box, ξ))
    c = shhc_constraints(hdata, m)
    [(x.kind, x.key) for x in c] == dict.keys ||
        error("add_point!: the constraint set changed at ξ = $ξ")
    μ, M = _compact(c, m, dict.harmonics)
    push!(dict.points, collect(Float64, ξ)); push!(dict.labels, label)
    push!(dict.μ, μ); push!(dict.M, M)
    return length(dict.points)
end

"""
    build_dictionary(hdata, model, box; seeds = Seed[], nbackground = 1024,
                     vertices = length(box) <= 10, seed = 1) -> ParameterDictionary

Evaluate the constraints of the robust program at the nominal point, the vertices (when
`vertices`), the `seeds` and `nbackground` Sobol points, and store them in compact form.
`hdata` must carry the PCE data (Cantelli or Gaussian bound).
"""
function build_dictionary(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox;
                          seeds = Seed[], nbackground::Int = 1024,
                          vertices::Bool = length(box) <= 10, seed::Integer = 1)
    d  = length(box)
    m0 = moments_at(hdata, model, nominal_point(box))
    c0 = shhc_constraints(hdata, m0)
    dict = ParameterDictionary(box, sort(m0.harm_nw_ints), [(c.kind, c.key) for c in c0],
                               [c.limit for c in c0], [c.λ for c in c0],
                               Vector{Float64}[], Symbol[], Matrix{Float64}[], Array{Float64,3}[])
    add_point!(dict, hdata, model, zeros(d), :nominal)
    if vertices
        for v in Iterators.product(ntuple(_ -> (-1.0, 1.0), d)...)
            add_point!(dict, hdata, model, collect(v), :vertex)
        end
    end
    for s in seeds
        add_point!(dict, hdata, model, s.ξ, :seed)
    end
    for ξ in background_points(box, nbackground; seed)
        add_point!(dict, hdata, model, ξ, :background)
    end
    return dict
end

# normalized violation of constraint j from its compact form at one point
function _compact_violation(dict::ParameterDictionary, j::Int, μ, M, gv)
    μg = dot(view(μ, :, j), gv)
    F  = dict.λ[j] <= 0.0 ? μg : μg + dict.λ[j] * sqrt(max(dot(gv, view(M, :, :, j), gv), 0.0))
    L  = dict.limit[j]
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
    nc, np = length(dict.keys), length(dict.points)
    V = zeros(nc, np)
    for s in 1:np, j in 1:nc
        V[j, s] = _compact_violation(dict, j, dict.μ[s], dict.M[s], gv)
    end
    return V
end
