################################################################################
# HarmonicPowerModels.jl                                                       #
# Exact mean and covariance of squared harmonic magnitudes                     #
#                                                                               #
# A bus voltage or branch current at harmonic h, per unit of hosting current,   #
#   Y_h = Σ_u r_{h,u} · exp(j φ_{h,u}(ω_u)),                                    #
# is a sum of independent unit contributions: unit u has its own germ ω_u, and  #
# that germ drives all harmonics of u, so φ_{h,u} and φ_{h',u} are dependent.   #
# The chance constraints of sHHC_SOC need E|Y_h|² and Cov(|Y_h|², |Y_h'|²).     #
# A polynomial chaos expansion of total degree δ truncates the angle functions, #
# which loses variance for wide angle supports (e.g., 39 % kept at δ = 2 for   #
# the 25th harmonic of a phase-controlled bridge). These moments are instead    #
# computed exactly from the Hoeffding decomposition of the quadratic form:      #
#                                                                               #
#   a_u = r_u e^{jφ_u} = μ_u + ã_u,   M = Σ_u μ_u,                              #
#   |Y|² = c + Σ_u f_u(ω_u) + Σ_{u<v} 2 Re(ã_u conj(ã_v)),                      #
#   c    = |M|² + Σ_u E|ã_u|²,                                                  #
#   f_u  = 2 Re(conj(M) ã_u) + |ã_u|² − E|ã_u|²,                                #
#                                                                               #
# whose terms are mutually orthogonal. Hence                                    #
#   E|Y_h|²               = |M_h|² + Σ_u |r_{h,u}|² (1 − |E e^{jφ_{h,u}}|²),    #
#   Cov(|Y_h|², |Y_h'|²)  = Σ_u E[f_{h,u} f_{h',u}]                              #
#                           + 2 (‖Σ_u S_u‖²_F − Σ_u ‖S_u‖²_F),                  #
# with S_u the 2×2 matrix E[(Re, Im) ã_{h,u} (Re, Im) ã_{h',u}ᵀ]. Only moments  #
# of a single germ up to the fourth order enter; they are integrated by a      #
# Gauss rule for the unit's own Beta measure, whose error is at machine         #
# precision for the angle supports of interest (it does not truncate).          #
# Authors: Tom Van Acker                                                        #
################################################################################

"""
    exact_quadratic_moments(R, ranges, shapes; nquad = 60)
        -> (mean, cov)

Exact mean and covariance over the harmonics of the squared magnitude |Y_h|² of

    Y_h = Σ_u R[h][u] · exp(j φ_{h,u}),   φ_{h,u} = lo_{h,u} + (hi_{h,u} − lo_{h,u}) ω_u,
    ω_u ~ Beta(shapes[u]...) independent across units u,

for one quantity (a bus voltage or a branch current).

# Arguments
- `R`      — `Vector` over the harmonics of `Vector{ComplexF64}` over the units.
- `ranges` — `Vector` over the harmonics of `Vector{Tuple}` (lo, hi) [rad] per unit.
- `shapes` — `Vector{Tuple}` of Beta shapes (α̂, β̂) per unit.
- `nquad`  — Gauss nodes per germ.

Returns `mean::Vector{Float64}` (length H) and `cov::Matrix{Float64}` (H × H).
See the file header for the derivation.
"""
function exact_quadratic_moments(R::AbstractVector, ranges::AbstractVector,
                                 shapes::AbstractVector; nquad::Int = 60)
    rules = _beta_gauss_rules(shapes, nquad)
    return _exact_quadratic_moments([R], ranges, shapes, rules)[1]
end

# Gauss rule for each distinct Beta shape (nodes on [0, 1], weights sum to one)
function _beta_gauss_rules(shapes, nquad)
    rules = Dict{Tuple{Float64,Float64}, Tuple{Vector{Float64},Vector{Float64}}}()
    for s in unique(shapes)
        op = PolyChaos.Beta01OrthoPoly(1, s...; Nrec = nquad + 1)
        rules[s] = (Vector{Float64}(op.quad.nodes), Vector{Float64}(op.quad.weights))
    end
    return rules
end

# Several quantities sharing the units' angle phasors: Rs[q][h][u]. The phasors,
# their means and centred values are computed once per (harmonic, unit).
function _exact_quadratic_moments(Rs::AbstractVector, ranges::AbstractVector,
                                  shapes::AbstractVector, rules::AbstractDict)
    H, U = length(ranges), length(shapes)
    # centred phasors P̃[h][u] (vector over nodes) and means μ[h][u]
    w  = [rules[shapes[u]][2] for u in 1:U]
    Pc = [Vector{Vector{ComplexF64}}(undef, U) for _ in 1:H]
    μ  = [zeros(ComplexF64, U) for _ in 1:H]
    for h in 1:H, u in 1:U
        x, wu = rules[shapes[u]]
        lo, hi = ranges[h][u]
        P = cis.(lo .+ (hi - lo) .* x)
        μ[h][u] = sum(wu .* P)
        Pc[h][u] = P .- μ[h][u]
    end
    out = Vector{Tuple{Vector{Float64},Matrix{Float64}}}(undef, length(Rs))
    for (q, R) in enumerate(Rs)
        M  = [sum(R[h][u] * μ[h][u] for u in 1:U) for h in 1:H]
        va = [[abs2(R[h][u]) * sum(w[u] .* abs2.(Pc[h][u])) for u in 1:U] for h in 1:H]  # E|ã|²
        mean = [abs2(M[h]) + sum(va[h]) for h in 1:H]
        # first-order Hoeffding terms f_{h,u} on the nodes, and centred ã_{h,u}
        A = [[R[h][u] .* Pc[h][u] for u in 1:U] for h in 1:H]
        F = [[2 .* real.(conj(M[h]) .* A[h][u]) .+ abs2.(A[h][u]) .- va[h][u] for u in 1:U] for h in 1:H]
        C = zeros(H, H)
        for h in 1:H, k in h:H
            c1 = 0.0
            Ssum = zeros(2, 2); Snorm = 0.0
            for u in 1:U
                c1 += sum(w[u] .* F[h][u] .* F[k][u])
                a, b = A[h][u], A[k][u]
                S = [sum(w[u] .* real.(a) .* real.(b)) sum(w[u] .* real.(a) .* imag.(b));
                     sum(w[u] .* imag.(a) .* real.(b)) sum(w[u] .* imag.(a) .* imag.(b))]
                Ssum .+= S; Snorm += sum(abs2, S)
            end
            C[h, k] = C[k, h] = c1 + 2 * (sum(abs2, Ssum) - Snorm)
        end
        out[q] = (mean, C)
    end
    return out
end

"""
    _moment_coefficients(mean, cov) -> Vector{Vector{Float64}}

Coefficients over the harmonics, in the layout sHHC_SOC uses for its expansions
with unit basis norms: entry 1 is the mean of |Y_h|², entries 2:H+1 a factor L
with L Lᵀ = cov (symmetric square root, negative eigenvalues of round-off
clipped). Then Var(Σ_h g_h |Y_h|²) = Σ_k (Σ_h L[h,k] g_h)², as for an expansion.
"""
function _moment_coefficients(mean::AbstractVector, cov::AbstractMatrix)
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(cov))
    L = F.vectors * LinearAlgebra.Diagonal(sqrt.(max.(F.values, 0.0)))
    return [vcat(mean[h], L[h, :]) for h in eachindex(mean)]
end
