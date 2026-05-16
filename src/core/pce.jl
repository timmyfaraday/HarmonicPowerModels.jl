################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.1 - added PCE pre-computation layer for PolyChaos v1.1.0                #
################################################################################

"""
Pre-computed PCE data for a Beta(α,β)-distributed uncertainty source.

Fields:
- `deg`    : polynomial degree
- `P`      : index vector 0:deg
- `P_size` : deg + 1
- `opq`    : Beta01OrthoPoly object (includes quadrature)
- `norms`  : squared norms ⟨ψₖ, ψₖ⟩ for k = 0:deg
- `T`      : 3-tensor T[k1,k2,k] = ⟨ψ_{k1-1} ψ_{k2-1} ψ_{k-1}⟩ / ⟨ψ_{k-1}²⟩
             (1-based indexing, so T[k1+1, k2+1, k+1] for 0-based k1,k2,k)
"""
struct PCEData
    deg    :: Int
    P      :: Vector{Int}
    P_size :: Int
    opq    :: _PCE.Beta01OrthoPoly
    norms  :: Vector{Float64}
    T      :: Array{Float64, 3}
end

"""
    build_pce_data(alpha, beta; deg=2) -> PCEData

Construct PCEData for a Beta(alpha, beta) distribution on [0,1].
Uses PolyChaos v1.1.0 API.
"""
function build_pce_data(alpha::Float64, beta::Float64; deg::Int = 2)::PCEData
    @assert alpha > 0.0 "alpha must be positive, got $alpha"
    @assert beta  > 0.0 "beta must be positive, got $beta"
    @assert deg   >= 1  "deg must be at least 1, got $deg"

    opq = _PCE.Beta01OrthoPoly(deg, alpha, beta;
                               Nrec = 2 * deg + 2,
                               addQuadrature = true)

    norms = _PCE.computeSP2(opq)

    t2 = _PCE.Tensor(2, opq)
    t3 = _PCE.Tensor(3, opq)

    T = Array{Float64, 3}(undef, deg + 1, deg + 1, deg + 1)
    for k in 0:deg, k1 in 0:deg, k2 in 0:deg
        T[k1+1, k2+1, k+1] = t3.get([k1, k2, k]) / t2.get([k, k])
    end

    # sanity: T[1,1,1] should equal 1 (ψ₀ = 1, ⟨1·1·1⟩/⟨1⟩ = 1)
    @assert abs(T[1, 1, 1] - 1.0) < 1e-10 "T[1,1,1] != 1: got $(T[1,1,1])"

    return PCEData(deg, collect(0:deg), deg + 1, opq, norms, T)
end

"""
    pce_angle_coefficients(pce, v_min_deg, v_max_deg) -> Vector{Float64}

Compute PCE coefficients for the phase angle φ ~ Uniform(v_min_deg, v_max_deg) [degrees]
that is modelled as a linear function of a Beta(α,β)-distributed variable ξ ∈ [0,1]:

    φ(ξ) = v_min_deg + (v_max_deg - v_min_deg) * ξ

Returns the PCE coefficient vector [φ₀, φ₁, …, φ_deg] in radians.
"""
function pce_angle_coefficients(pce::PCEData,
                                v_min_deg::Float64,
                                v_max_deg::Float64)::Vector{Float64}
    v_min_rad = v_min_deg * π / 180.0
    v_max_rad = v_max_deg * π / 180.0
    Δ = v_max_rad - v_min_rad

    nodes   = pce.opq.quad.nodes
    weights = pce.opq.quad.weights
    nq      = length(nodes)

    # α and β from the Beta01OrthoPoly measure (PolyChaos v1.1.0 field names)
    alpha = pce.opq.measure.ashapeParameter
    beta  = pce.opq.measure.bshapeParameter

    # weight function w(ξ) = ξ^(α-1) * (1-ξ)^(β-1) / B(α,β)  (already in weights)
    # φ(ξ) = v_min_rad + Δ * ξ
    coeffs = zeros(Float64, pce.P_size)
    for k in pce.P
        integral = 0.0
        for qi in 1:nq
            xi  = nodes[qi]
            w   = weights[qi]
            phi = v_min_rad + Δ * xi
            psi = _PCE.evaluate(k, xi, pce.opq)
            integral += w * phi * psi
        end
        coeffs[k+1] = integral / pce.norms[k+1]
    end

    return coeffs
end

"""
    galerkin_cos(pce, phi_coeffs) -> Vector{Float64}

Compute the PCE coefficients of cos(φ) given the PCE coefficients of φ.

Uses the Taylor expansion:
    cos(φ) ≈ cos(φ₀) - sin(φ₀)(φ - φ₀) - cos(φ₀)/2 (φ - φ₀)² + …

truncated to degree `pce.deg`, where products are computed via the Galerkin
multiplication tensor T stored in `pce`.
"""
function galerkin_cos(pce::PCEData, phi_coeffs::Vector{Float64})::Vector{Float64}
    phi0 = phi_coeffs[1]  # mean (coefficient of ψ₀)

    # fluctuation coefficients: δφ = φ - φ₀ψ₀  → same coeffs but zero-th set to 0
    delta = copy(phi_coeffs)
    delta[1] = 0.0

    # Taylor: cos(φ₀ + δφ) = cos(φ₀)·T_even - sin(φ₀)·T_odd
    # where T_even = 1 - δφ²/2 + …  and T_odd = δφ - δφ³/6 + …
    # We keep terms up to δφ^deg using polynomial multiplication via T.

    # Helper: multiply two PCE coefficient vectors using T
    function pce_mul(a::Vector{Float64}, b::Vector{Float64})::Vector{Float64}
        c = zeros(Float64, pce.P_size)
        for k in pce.P, k1 in pce.P, k2 in pce.P
            c[k+1] += pce.T[k1+1, k2+1, k+1] * a[k1+1] * b[k2+1]
        end
        return c
    end

    # Build powers of delta up to pce.deg
    powers = Vector{Vector{Float64}}(undef, pce.deg + 1)
    powers[1] = [1.0; zeros(Float64, pce.P_size - 1)]   # δφ^0 = 1 (ψ₀ basis)
    if pce.deg >= 1
        powers[2] = copy(delta)                           # δφ^1
    end
    for n in 2:pce.deg
        powers[n+1] = pce_mul(powers[n], delta)
    end

    # Accumulate Taylor series
    cos_coeffs = zeros(Float64, pce.P_size)

    cosphi0 = cos(phi0)
    sinphi0 = sin(phi0)

    factorial_n = 1.0
    for n in 0:pce.deg
        factorial_n = (n == 0) ? 1.0 : factorial_n * n
        term = powers[n+1] ./ factorial_n
        if n % 4 == 0
            cos_coeffs .+= cosphi0 .* term
        elseif n % 4 == 1
            cos_coeffs .-= sinphi0 .* term
        elseif n % 4 == 2
            cos_coeffs .-= cosphi0 .* term
        else
            cos_coeffs .+= sinphi0 .* term
        end
    end

    return cos_coeffs
end

"""
    galerkin_sin(pce, phi_coeffs) -> Vector{Float64}

Compute the PCE coefficients of sin(φ) given the PCE coefficients of φ.
Mirrors `galerkin_cos` using the sine Taylor expansion.
"""
function galerkin_sin(pce::PCEData, phi_coeffs::Vector{Float64})::Vector{Float64}
    phi0 = phi_coeffs[1]

    delta = copy(phi_coeffs)
    delta[1] = 0.0

    function pce_mul(a::Vector{Float64}, b::Vector{Float64})::Vector{Float64}
        c = zeros(Float64, pce.P_size)
        for k in pce.P, k1 in pce.P, k2 in pce.P
            c[k+1] += pce.T[k1+1, k2+1, k+1] * a[k1+1] * b[k2+1]
        end
        return c
    end

    powers = Vector{Vector{Float64}}(undef, pce.deg + 1)
    powers[1] = [1.0; zeros(Float64, pce.P_size - 1)]
    if pce.deg >= 1
        powers[2] = copy(delta)
    end
    for n in 2:pce.deg
        powers[n+1] = pce_mul(powers[n], delta)
    end

    sin_coeffs = zeros(Float64, pce.P_size)

    cosphi0 = cos(phi0)
    sinphi0 = sin(phi0)

    factorial_n = 1.0
    for n in 0:pce.deg
        factorial_n = (n == 0) ? 1.0 : factorial_n * n
        term = powers[n+1] ./ factorial_n
        if n % 4 == 0
            sin_coeffs .+= sinphi0 .* term
        elseif n % 4 == 1
            sin_coeffs .+= cosphi0 .* term
        elseif n % 4 == 2
            sin_coeffs .-= sinphi0 .* term
        else
            sin_coeffs .-= cosphi0 .* term
        end
    end

    return sin_coeffs
end

"""
    pce_mean(coeffs) -> Float64

Return the mean of a PCE expansion: E[X] = coeffs[1] (the ψ₀ coefficient).
"""
pce_mean(coeffs::Vector{Float64}) = coeffs[1]

"""
    pce_variance(coeffs, pce) -> Float64

Return the variance of a PCE expansion: Var[X] = Σ_{k≥1} coeffs[k+1]² · ⟨ψₖ²⟩.
"""
function pce_variance(coeffs::Vector{Float64}, pce::PCEData)::Float64
    return sum(coeffs[k+1]^2 * pce.norms[k+1] for k in 1:pce.deg)
end
