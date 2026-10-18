# HarmonicPowerModels.jl                                                       #
# Unified PCE workspace built on PolyChaos.jl.                                  #
#                                                                               #
# Convention — monic Beta01 basis w.r.t. the PROBABILITY measure               #
#   ω_u ~ Beta(α̂, β̂) on [0,1],  ψ_k = monic orthogonal polys of PolyChaos  #
#   Ψ_κ(ω) = ∏_u ψ_{κ_u}(ω_u),   ⟨Ψ_0, Ψ_0⟩ = 1                          #
#                                                                               #
# Because the measure is the Beta pdf itself (PolyChaos Beta01Measure):        #
#   E[X]   = c_0                              (first coefficient)              #
#   Var[X] = Σ_{k≥1} ⟨Ψ_k,Ψ_k⟩ · c_k²      (no total-measure correction)   #
# This removes the ξ = 2ω−1 affine map, the a = β̂−1 / b = α̂−1 parameter   #
# swap, and the norm0 = √⟨Ψ_0,Ψ_0⟩ corrections of the legacy Jacobi code.  #
#                                                                               #
# Galerkin tensor — stored SPARSE (99.9 % of entries vanish at IEEE-30 size): #
#   M[j,l,k] = ⟨Ψ_j Ψ_l, Ψ_k⟩ / ⟨Ψ_k, Ψ_k⟩    (1-indexed)                #
# kept as nonzero triples with j ≤ l; `galerkin_product` restores symmetry.   #
#                                                                               #
# PolyChaos pitfalls handled here:                                              #
# - computeSP2(::MultiOrthoPoly) is broken in PolyChaos 1.1 (accesses a       #
#   nonexistent field); multivariate norms are built from the univariate       #
#   ones and the multi-index set instead.                                       #
# - Tensor(3, ·) needs ⌈(3δ+1)/2⌉ quadrature nodes; the default              #
#   Nrec = δ+1 throws "not enough nodes". The constructor uses Nrec = 3δ+2.  #
# - Multi-index ordering is PolyChaos's degree-graded `mop.ind` (NOT the      #
#   lexicographic ordering of the legacy code); row 1 is the constant term.    #
#                                                                               #
# Trigonometric coefficients — cos(φ)/sin(φ) for an angle affine in one germ  #
# are projected onto the basis EXACTLY, by the Gauss rule already carried by   #
# `uni.quad` (`affine_trig_pce`). The superseded Taylor route (`trig_pce` plus #
# a π/2 support rotation) is kept only for angles that are not affine in a     #
# single germ; it is 1–2 orders of magnitude less accurate in the variance.    #
# An angle affine in two germs (own + shared type germ) is projected exactly by #
# `affine2_trig_pce`, as the product of the two univariate projections.         #
# Authors: Tom Van Acker                                                        #

using PolyChaos

"""
    PCEBasis

Immutable workspace for a multivariate PCE in `N` independent Beta(α̂, β̂)
germs at total degree `deg`. All quantities refer to the monic Beta01 basis
w.r.t. the probability measure (see file header).

# Fields
- `N`, `deg`, `dim` — germ count, total degree δ, basis size P = C(N+δ, δ).
- `uni`             — shared univariate `Beta01OrthoPoly` (with quadrature).
- `mop`             — `MultiOrthoPoly` product basis.
- `ind`             — P×N multi-index matrix (`mop.ind`, degree-graded).
- `sp2`             — ⟨Ψ_k,Ψ_k⟩, length P; `sp2[1] == 1`.
- `M_nz`            — nonzero Galerkin-tensor triples `(j, l, k, M[j,l,k])`
                      with j ≤ l, 1-indexed.
"""
struct PCEBasis
    N    ::Int
    deg  ::Int
    dim  ::Int
    uni  ::PolyChaos.Beta01OrthoPoly
    mop  ::PolyChaos.MultiOrthoPoly
    ind  ::Matrix{Int}
    sp2  ::Vector{Float64}
    M_nz ::Vector{Tuple{Int,Int,Int,Float64}}
end

"""
    PCEBasis(N, α̂, β̂, δ; nrec = 3δ + 2) -> PCEBasis

Build the PCE workspace for `N` independent Beta(α̂, β̂) germs at total
degree `δ`. `nrec` controls the number of recurrence coefficients and thereby
the quadrature size (`nrec − 1` nodes); the default integrates the degree-3δ
triple products of the Galerkin tensor exactly.
"""
function PCEBasis(N::Int, α̂::Real, β̂::Real, δ::Int; nrec::Int = 3δ + 2)
    N ≥ 1 || error("PCEBasis: N must be ≥ 1")
    uni = Beta01OrthoPoly(δ, α̂, β̂; Nrec = nrec)
    mop = MultiOrthoPoly([uni for _ in 1:N], δ)
    ind = mop.ind
    P   = mop.dim

    # multivariate norms ⟨Ψ_k,Ψ_k⟩ = ∏_u ⟨ψ_{κ_u},ψ_{κ_u}⟩
    # (computeSP2(::MultiOrthoPoly) is broken in PolyChaos 1.1 — see header)
    sp2_uni = computeSP2(uni)
    sp2     = [prod(sp2_uni[ind[k, u] + 1] for u in 1:N) for k in 1:P]

    # sparse Galerkin tensor: M[j,l,k] = ⟨Ψ_jΨ_l,Ψ_k⟩/⟨Ψ_k,Ψ_k⟩, j ≤ l
    t3   = Tensor(3, mop)
    M_nz = Tuple{Int,Int,Int,Float64}[]
    idx  = Vector{Int}(undef, 3)
    for j in 1:P, l in j:P, k in 1:P
        idx[1] = j - 1; idx[2] = l - 1; idx[3] = k - 1
        v = t3.get(idx)
        abs(v) > 1e-14 && push!(M_nz, (j, l, k, v / sp2[k]))
    end

    return PCEBasis(N, δ, P, uni, mop, ind, sp2, M_nz)
end

"""
    dense_tensor(pb::PCEBasis) -> Array{Float64,3}

Materialise the Galerkin tensor M[j,l,k] (1-indexed, dim³). Intended for the
small univariate basis stored in `hdata["pce"]`; prefer `galerkin_product`
(sparse) for computations.
"""
function dense_tensor(pb::PCEBasis)
    M = zeros(pb.dim, pb.dim, pb.dim)
    for (j, l, k, v) in pb.M_nz
        M[j, l, k] = v
        M[l, j, k] = v
    end
    return M
end

"""
    galerkin_product(pb::PCEBasis, x, y) -> Vector{Float64}

PCE coefficients of the product X·Y from the coefficient vectors `x` and `y`,
via the sparse Galerkin tensor:  z_k = Σ_{j,l} M[j,l,k] x_j y_l.
"""
function galerkin_product(pb::PCEBasis, x::AbstractVector{<:Real},
                          y::AbstractVector{<:Real})
    z = zeros(pb.dim)
    for (j, l, k, v) in pb.M_nz
        if j == l
            z[k] += v * x[j] * y[j]
        else
            z[k] += v * (x[j] * y[l] + x[l] * y[j])
        end
    end
    return z
end

"""
    affine_pce(pb::PCEBasis, u, v_min, v_max) -> Vector{Float64}

PCE coefficients (length `pb.dim`) of the affine function

    φ(ω) = v_min + (v_max − v_min) · ω_u,   ω_u ~ Beta(α̂, β̂) on [0,1],

where `u` selects the germ dimension. Exact (no quadrature): uses
`convert2affinePCE` and scatters into the multivariate index set.

A degenerate range (`v_min == v_max`) is a deterministic (point-mass) angle:
the function is the constant `v_min`, so only the k=0 coefficient is non-zero.
This is the mechanism for reducing the stochastic models to the deterministic
one (variance zero ⇒ the Cantelli bound collapses to the deterministic limit).
"""
function affine_pce(pb::PCEBasis, u::Int, v_min::Real, v_max::Real)
    if v_min == v_max
        c = zeros(pb.dim)
        c[1] = Float64(v_min)          # constant term at the k=0 (constant) index
        return c
    end
    c_uni = convert2affinePCE(Float64(v_min), Float64(v_max), pb.uni)  # [c₀, c₁]
    return Vector{Float64}(assign2multi(c_uni, u, pb.ind))
end

"""
    trig_pce(pb::PCEBasis, phi; n_terms = 4) -> (cos_coefs, sin_coefs)

PCE coefficients of cos(φ) and sin(φ) via truncated Taylor series and
Galerkin products:

    cos(φ) ≈ Σ_{m=0}^{n_terms}   (−1)^m φ^{2m}   / (2m)!
    sin(φ) ≈ Σ_{m=0}^{n_terms−1} (−1)^m φ^{2m+1} / (2m+1)!

General-purpose fallback for an angle given by an arbitrary PCE coefficient
vector. Accuracy degrades for |φ| ≫ π/2, and each Galerkin product truncates
back to total degree `pb.deg`, so the error accumulates with `n_terms`.

!!! note "Prefer `affine_trig_pce` for an affine angle"
    When φ is affine in a single germ — which is the case for every harmonic
    unit in the SHHC models — [`affine_trig_pce`](@ref) projects cos/sin onto
    the basis exactly by quadrature and is one to two orders of magnitude more
    accurate at no extra cost. This function is retained only for angles that
    are not affine in one germ.
"""
function trig_pce(pb::PCEBasis, phi::AbstractVector{<:Real}; n_terms::Int = 4)
    maxpow = 2n_terms
    powers = Vector{Vector{Float64}}(undef, maxpow)
    powers[1] = Vector{Float64}(phi)
    for n in 2:maxpow
        powers[n] = galerkin_product(pb, powers[n-1], phi)
    end

    cos_c    = zeros(pb.dim)
    cos_c[1] = 1.0
    for m in 1:n_terms
        cos_c .+= (-1)^m / factorial(2m) .* powers[2m]
    end

    sin_c = zeros(pb.dim)
    for m in 0:n_terms-1
        sin_c .+= (-1)^m / factorial(2m + 1) .* powers[2m + 1]
    end

    return cos_c, sin_c
end

"""
    affine_trig_pce(pb::PCEBasis, u, v_min, v_max) -> (cos_coefs, sin_coefs)

PCE coefficients of cos(φ) and sin(φ) for the affine angle

    φ(ω) = v_min + (v_max − v_min)·ω_u ,   ω_u ~ Beta(α̂, β̂) on [0,1],

with `v_min`, `v_max` in **radians** and `u` the germ dimension.

Because φ depends on the single germ ω_u, cos(φ) and sin(φ) do too, so their
PCE is supported on that germ's univariate indices. The coefficients are the
exact Galerkin projections, evaluated with the Gauss rule that already lives
in `pb.uni`:

    c_j = Σ_q w_q · f(φ(ξ_q)) · ψ_j(ξ_q) / ⟨ψ_j, ψ_j⟩ ,   f ∈ {cos, sin},

then scattered onto the multivariate index set. The only error left is the
degree-`pb.deg` truncation itself; there is no Taylor truncation and no
expansion point, so accuracy does not depend on where the support sits or how
wide it is. Measured against the exact moments, the relative error in the
standard deviation is ~1e-5 over the angle ranges used for LCC converters,
against ~1e-3 for the Taylor route it replaces.

A degenerate range (`v_min == v_max`) is a deterministic angle: only the k=0
coefficient is non-zero, which is the mechanism that reduces the stochastic
models to the deterministic one.
"""
function affine_trig_pce(pb::PCEBasis, u::Int, v_min::Real, v_max::Real)
    if v_min == v_max                       # point mass ⇒ constant term only
        cos_c = zeros(pb.dim); cos_c[1] = cos(Float64(v_min))
        sin_c = zeros(pb.dim); sin_c[1] = sin(Float64(v_min))
        return cos_c, sin_c
    end

    nodes   = pb.uni.quad.nodes
    weights = pb.uni.quad.weights          # Beta01 PROBABILITY measure: Σ w = 1
    sp2_uni = computeSP2(pb.uni)

    φ  = v_min .+ (v_max - v_min) .* nodes
    cφ = cos.(φ)
    sφ = sin.(φ)

    cos_uni = zeros(pb.deg + 1)
    sin_uni = zeros(pb.deg + 1)
    for j in 0:pb.deg
        ψ = evaluate(j, nodes, pb.uni)
        cos_uni[j+1] = sum(weights .* cφ .* ψ) / sp2_uni[j+1]
        sin_uni[j+1] = sum(weights .* sφ .* ψ) / sp2_uni[j+1]
    end

    return Vector{Float64}(assign2multi(cos_uni, u, pb.ind)),
           Vector{Float64}(assign2multi(sin_uni, u, pb.ind))
end

"Univariate Galerkin coefficients (length deg+1) of cos/sin of φ = lo + (hi − lo)·ω."
function _uni_trig(pb::PCEBasis, lo::Real, hi::Real)
    nodes, weights = pb.uni.quad.nodes, pb.uni.quad.weights
    sp2_uni = computeSP2(pb.uni)
    φ = lo .+ (hi - lo) .* nodes
    cφ, sφ = cos.(φ), sin.(φ)
    c = zeros(pb.deg + 1); s = zeros(pb.deg + 1)
    for j in 0:pb.deg
        ψ = evaluate(j, nodes, pb.uni)
        c[j+1] = sum(weights .* cφ .* ψ) / sp2_uni[j+1]
        s[j+1] = sum(weights .* sφ .* ψ) / sp2_uni[j+1]
    end
    return c, s
end

"""
    affine2_trig_pce(pb::PCEBasis, u, v_min, v_max, t, s_min, s_max)
        -> (cos_coefs, sin_coefs)

PCE coefficients of cos(φ) and sin(φ) for an angle affine in two independent
germs — the unit's own germ `u` and a germ `t` it shares with other units, e.g.,
a converter type's firing delay:

    φ = φ₁(ω_u) + φ₂(ω_t),   φ₁ = v_min + (v_max − v_min)·ω_u,
                              φ₂ = s_min + (s_max − s_min)·ω_t

(radians; `s_max < s_min` is allowed). Since cos φ = cos φ₁ cos φ₂ − sin φ₁ sin φ₂
and sin φ = sin φ₁ cos φ₂ + cos φ₁ sin φ₂ are sums of products of functions of
independent germs, and the basis is the tensor product of the univariate ones,
the projection is the product of the univariate projections (each exact by the
Gauss rule in `pb.uni`), truncated to total degree `pb.deg`.
"""
function affine2_trig_pce(pb::PCEBasis, u::Int, v_min::Real, v_max::Real,
                          t::Int, s_min::Real, s_max::Real)
    u == t && error("affine2_trig_pce: the two germs must differ")
    c1, s1 = _uni_trig(pb, v_min, v_max)
    c2, s2 = _uni_trig(pb, s_min, s_max)
    cos_c = zeros(pb.dim); sin_c = zeros(pb.dim)
    for k in 1:pb.dim
        row = view(pb.ind, k, :)
        any(row[g] != 0 for g in eachindex(row) if g != u && g != t) && continue
        j, l = row[u] + 1, row[t] + 1
        cos_c[k] = c1[j] * c2[l] - s1[j] * s2[l]
        sin_c[k] = s1[j] * c2[l] + c1[j] * s2[l]
    end
    return cos_c, sin_c
end

"""
    rotated_trig_pce(pb::PCEBasis, u, v_min, v_max; n_terms = 4)
        -> (cos_coefs, sin_coefs)

Deprecated alias for [`affine_trig_pce`](@ref), kept so that existing callers
keep working. `n_terms` is accepted and ignored: the quadrature projection has
no Taylor series to truncate, and no support rotation is needed.
"""
rotated_trig_pce(pb::PCEBasis, u::Int, v_min::Real, v_max::Real;
                 n_terms::Int = 4) = affine_trig_pce(pb, u, v_min, v_max)

"""
    recentre_range(v_min, v_max, centre, Eω) -> (lo, hi)

Shift an angle support so that its expectation lands on `centre`, preserving its
width and shape.

For φ = v_min + (v_max − v_min)·ω with ω ~ Beta(α̂, β̂) on [0,1] and mean fraction
`Eω = α̂/(α̂+β̂)`, the support's own mean is `Eᵛ = v_min + (v_max − v_min)·Eω`, so

    (lo, hi) = (centre + v_min − Eᵛ,  centre + v_max − Eᵛ)

gives `E[φ] = centre` with the same dispersion.

A degenerate support (`v_min == v_max`) collapses onto `centre` exactly, which is
what makes a point-mass angle reduce the stochastic models to the deterministic
one *per unit* rather than only on average.
"""
function recentre_range(v_min::Real, v_max::Real, centre::Real, Eω::Real)
    Eᵛ = v_min + (v_max - v_min) * Eω
    return (Float64(centre + v_min - Eᵛ), Float64(centre + v_max - Eᵛ))
end

"""
    angle_pce_coefficients(pb::PCEBasis, units, harmonics; n_terms = 4)
        -> (cos_coefs, sin_coefs)

For each unit and harmonic h, the length-`pb.dim` PCE coefficient vectors of
cos(φ_{u,h}) and sin(φ_{u,h}), where φ_{u,h} is affine in germ dimension u.

`units` is an `AbstractVector` of `NamedTuple`s with fields `id` and
`h_angles::Dict{Int,Tuple}` mapping harmonic order → (v_min, v_max) in
**degrees**; the position of each unit in the vector is its germ dimension.
Returns two `Dict`s keyed by `(unit_id, h)`.

Each pair is the exact Galerkin projection of cos/sin onto the germ's
univariate indices (see [`affine_trig_pce`](@ref)). `n_terms` is accepted and
ignored; it belonged to the superseded Taylor expansion.
"""
function angle_pce_coefficients(pb::PCEBasis,
                                units::AbstractVector,
                                harmonics::AbstractVector{Int};
                                n_terms::Int = 4)
    cos_coefs = Dict{Tuple{Any,Int}, Vector{Float64}}()
    sin_coefs = Dict{Tuple{Any,Int}, Vector{Float64}}()

    for (u_idx, unit) in enumerate(units)
        for h in harmonics
            haskey(unit.h_angles, h) || continue
            v_min_deg, v_max_deg = unit.h_angles[h]
            cos_c, sin_c = affine_trig_pce(pb, u_idx, deg2rad(v_min_deg),
                                           deg2rad(v_max_deg))
            cos_coefs[(unit.id, h)] = cos_c
            sin_coefs[(unit.id, h)] = sin_c
        end
    end

    return cos_coefs, sin_coefs
end

"""
    cantelli_lambda(ε) -> Float64

One-sided Cantelli (Chebyshev) safety factor for risk level ε: √((1 − ε)/ε).
Distribution-free: P(X > μ + λσ) ≤ ε for any distribution.
"""
cantelli_lambda(ε::Real) = sqrt((1 - ε) / ε)

"""
    gaussian_lambda(ε) -> Float64

Gaussian-assumption safety factor for risk level ε: λ_Φ(ε) = Φ⁻¹(1 − ε), the
standard-normal quantile. Exact (P(X > μ + λσ) = ε) when X ~ N(μ, σ²); much less
conservative than Cantelli (ε=0.05: 1.645 vs 4.359).
"""
gaussian_lambda(ε::Real) = Distributions.quantile(Distributions.Normal(), 1 - ε)

"""
    risk_factor(ε, bound::Symbol) -> Float64

Chance-constraint safety factor λ(ε) for the requested `bound`:
`:cantelli` (default, distribution-free) or `:gaussian` (Gaussian assumption).
"""
risk_factor(ε::Real, bound::Symbol) =
    bound === :gaussian  ? gaussian_lambda(ε) :
    bound === :cantelli  ? cantelli_lambda(ε) :
    bound === :empirical ? 1.0 :  # scalar fallback; per-constraint λ* set separately
    error("risk_factor: unknown bound $bound (use :cantelli, :gaussian or :empirical)")

"""
    pce_mean(coefs) -> Float64

E[X] = c₀ (probability-measure basis, `sp2[1] == 1`).
"""
pce_mean(coefs::AbstractVector{<:Real}) = coefs[1]

"""
    pce_variance(coefs, pb::PCEBasis) -> Float64

Var[X] = Σ_{k≥2} ⟨Ψ_k,Ψ_k⟩ · c_k²  (probability-measure basis; no
total-measure correction needed).
"""
pce_variance(coefs::AbstractVector{<:Real}, pb::PCEBasis) =
    sum(pb.sp2[k] * coefs[k]^2 for k in 2:pb.dim; init = 0.0)
