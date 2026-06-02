################################################################################
# HarmonicPowerModels.jl                                                       #
# PCE utilities for stochastic harmonic hosting capacity (SHHC).               #
# Pure-Julia; no PowerModels dependency.                                        #
################################################################################

using LinearAlgebra
using PolyChaos

# ──────────────────────────────────────────────────────────────────────────────
# 1. build_jacobi_basis
# ──────────────────────────────────────────────────────────────────────────────

"""
    build_jacobi_basis(α, β, δ)

Build a univariate Jacobi PCE basis of degree `δ` for Beta(α,β) on [0,1].

Beta(α,β) on [0,1] is mapped to Jacobi(β-1, α-1) on [-1,1] via t = 2u-1:
    u^(α-1)(1-u)^(β-1) du  →  (1+t)^(α-1)(1-t)^(β-1) dt/2

PolyChaos stores the Jacobi weight as (1-t)^a * (1+t)^b, so
    ashapeParameter = β-1,  bshapeParameter = α-1.
"""
function build_jacobi_basis(α::Real, β::Real, δ::Integer)
    return JacobiOrthoPoly(δ, β - 1, α - 1)
end

# ──────────────────────────────────────────────────────────────────────────────
# 2. compute_multiplication_tensor
# ──────────────────────────────────────────────────────────────────────────────

"""
    compute_multiplication_tensor(basis)

Return the (δ+1)×(δ+1)×(δ+1) tensor M[j,l,k] = ⟨ψ_j ψ_l, ψ_k⟩ / ⟨ψ_k, ψ_k⟩.

Index convention: Julia index (j+1, l+1, k+1) for 0-based degrees j, l, k.

A fresh Gauss-Jacobi rule with ⌈(3δ+2)/2⌉ nodes is built to integrate the
degree-3δ triple products exactly.
"""
function compute_multiplication_tensor(basis::JacobiOrthoPoly)
    δ  = basis.deg
    a  = basis.measure.ashapeParameter
    b  = basis.measure.bshapeParameter
    nq = max(3δ + 2, δ + 2)           # exact for polynomials of degree ≤ 2nq-1
    quad_op = JacobiOrthoPoly(nq, a, b)
    nodes, weights = quad_op.quad.nodes, quad_op.quad.weights
    sp2 = computeSP2(basis)

    M = zeros(δ+1, δ+1, δ+1)
    for j in 0:δ, l in 0:δ, k in 0:δ
        vj = evaluate.(j, nodes, Ref(basis))
        vl = evaluate.(l, nodes, Ref(basis))
        vk = evaluate.(k, nodes, Ref(basis))
        M[j+1, l+1, k+1] = dot(vj .* vl .* vk, weights) / sp2[k+1]
    end
    return M
end

# ──────────────────────────────────────────────────────────────────────────────
# 3. cantelli_lambda
# ──────────────────────────────────────────────────────────────────────────────

"""
    cantelli_lambda(ε)

One-sided Cantelli safety factor for risk level ε: √((1-ε)/ε).
"""
cantelli_lambda(ε::Real) = sqrt((1 - ε) / ε)

# ──────────────────────────────────────────────────────────────────────────────
# 4. pce_expectation
# ──────────────────────────────────────────────────────────────────────────────

"""
    pce_expectation(coeffs)

Return the k=0 PCE coefficient, which equals E[X̂].
"""
pce_expectation(coeffs) = coeffs[1]

# ──────────────────────────────────────────────────────────────────────────────
# 5. pce_variance
# ──────────────────────────────────────────────────────────────────────────────

"""
    pce_variance(coeffs, basis)

Return Var[X̂] = Σ_{k=1}^{δ} ⟨ψ_k, ψ_k⟩ · coeffs[k+1]².
"""
function pce_variance(coeffs, basis::JacobiOrthoPoly)
    sp2 = computeSP2(basis)
    return sum(sp2[k+1] * coeffs[k+1]^2 for k in 1:length(coeffs)-1)
end

# ──────────────────────────────────────────────────────────────────────────────
# 6. affine_angle_pce
# ──────────────────────────────────────────────────────────────────────────────

"""
    affine_angle_pce(v_min, v_max, basis)

PCE coefficients of the affine angle φ = v_min + (v_max - v_min)·ω, where
ω ~ Beta on [0,1] is represented in the Jacobi basis.

Uses Gauss quadrature to project ω onto each basis function:
    ω_k = ⟨ω, ψ_k⟩ / ⟨ψ_k, ψ_k⟩

Then  φ_0 = v_min + (v_max - v_min)·ω_0  (= v_min + (v_max-v_min)·E[ω])
      φ_k = (v_max - v_min)·ω_k  for k ≥ 1  (≈ 0 for k ≥ 2 since ω is linear)
"""
function affine_angle_pce(v_min::Real, v_max::Real, basis::JacobiOrthoPoly)
    δ  = basis.deg
    a  = basis.measure.ashapeParameter
    b  = basis.measure.bshapeParameter
    nq = max(2δ + 2, δ + 2)
    quad_op = JacobiOrthoPoly(nq, a, b)
    nodes, weights = quad_op.quad.nodes, quad_op.quad.weights
    sp2 = computeSP2(basis)

    ω_nodes = @. (nodes + 1) / 2   # affine map [-1,1] → [0,1]

    phi_coeffs = zeros(δ + 1)
    for k in 0:δ
        vk    = evaluate.(k, nodes, Ref(basis))
        ω_k   = dot(ω_nodes .* vk, weights) / sp2[k+1]
        phi_coeffs[k+1] = (k == 0 ? v_min : 0.0) + (v_max - v_min) * ω_k
    end
    return phi_coeffs
end

# ──────────────────────────────────────────────────────────────────────────────
# 7. trig_pce_coefficients
# ──────────────────────────────────────────────────────────────────────────────

"""
    trig_pce_coefficients(phi_coeffs, M; n_terms=4)

PCE coefficients of cos(φ) and sin(φ) via truncated Taylor series and Galerkin
projection using the multiplication tensor M.

Powers of φ are built iteratively:
    (φ^1)_k = phi_coeffs[k+1]
    (φ^n)_k = Σ_{j,l} M[j,l,k] · (φ^{n-1})_j · phi_coeffs[l+1]

Taylor truncation:
    cos(φ) ≈ Σ_{m=0}^{n_terms}  (-1)^m φ^{2m} / (2m)!
    sin(φ) ≈ Σ_{m=0}^{n_terms-1} (-1)^m φ^{2m+1} / (2m+1)!

Returns `(cos_coeffs, sin_coeffs)`, each a vector of length δ+1.
"""
function trig_pce_coefficients(phi_coeffs::AbstractVector, M::Array{<:Real,3};
                                n_terms::Int=4)
    δp1     = length(phi_coeffs)
    maxpow  = 2 * n_terms          # highest power needed (for cos)

    # Build all powers φ^1 … φ^maxpow
    powers = Vector{Vector{Float64}}(undef, maxpow)
    powers[1] = Float64.(phi_coeffs)
    for n in 2:maxpow
        prev = powers[n-1]
        nxt  = zeros(δp1)
        for k in 1:δp1, j in 1:δp1, l in 1:δp1
            nxt[k] += M[j, l, k] * prev[j] * phi_coeffs[l]
        end
        powers[n] = nxt
    end

    # cos(φ) = 1 - φ²/2! + φ⁴/4! - …
    cos_coeffs    = zeros(δp1)
    cos_coeffs[1] = 1.0               # φ^0 contributes only to k=0
    for m in 1:n_terms
        cos_coeffs .+= ((-1)^m / factorial(2m)) .* powers[2m]
    end

    # sin(φ) = φ - φ³/3! + φ⁵/5! - …
    sin_coeffs = zeros(δp1)
    for m in 0:n_terms-1
        sin_coeffs .+= ((-1)^m / factorial(2m + 1)) .* powers[2m + 1]
    end

    return cos_coeffs, sin_coeffs
end
