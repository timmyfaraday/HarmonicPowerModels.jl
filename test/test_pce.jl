################################################################################
# test/test_pce.jl                                                              #
# Comprehensive tests for src/util/pce.jl (SHHC polynomial chaos utilities).  #
# Runnable standalone: julia --project=. test/test_pce.jl                      #
################################################################################

using Test
using Statistics
using Random
using PolyChaos
using Distributions

include(joinpath(@__DIR__, "..", "src", "util", "pce.jl"))

# ── Parameters (§II, Table III of the SHHC paper) ─────────────────────────────
const _α̂ = 4.0414          # Beta shape parameter α̂
const _β̂ = 4.7762          # Beta shape parameter β̂
const _δ  = 2               # PCE degree

const _d2r = π / 180        # degrees → radians

# Harmonic angle bounds in radians (Table III, converted from degrees)
const _H_BOUNDS = Dict(
     5 => (  4.75 * _d2r,  40.05 * _d2r),
     7 => (-22.77 * _d2r,  27.51 * _d2r),
    11 => ( 79.69 * _d2r, 164.01 * _d2r),
    13 => ( 49.38 * _d2r, 147.67 * _d2r),
)

# ── Shared objects (built once) ────────────────────────────────────────────────
const _basis = build_jacobi_basis(_α̂, _β̂, _δ)
const _M     = compute_multiplication_tensor(_basis)

# Use isapprox directly throughout to avoid conflicts with any ≈ overrides in
# the parent test suite (runtests.jl redefines ≈ without keyword support).

# ══════════════════════════════════════════════════════════════════════════════
# A. Multiplication tensor
# ══════════════════════════════════════════════════════════════════════════════
@testset "A. Multiplication tensor" begin

    @test size(_M) == (_δ+1, _δ+1, _δ+1)

    @test all(isfinite, _M)

    # M[1,1,1] = ⟨ψ_0·ψ_0, ψ_0⟩ / ⟨ψ_0,ψ_0⟩ = 1  (ψ_0 = 1)
    @test isapprox(_M[1, 1, 1], 1.0; atol=1e-10)

    # symmetry: M[j,l,k] = M[l,j,k]
    for j in 1:_δ+1, l in 1:_δ+1, k in 1:_δ+1
        @test isapprox(_M[j, l, k], _M[l, j, k]; atol=1e-10)
    end

    # k=0 slice is diagonal: ⟨ψ_j, ψ_l⟩ = 0 for j≠l (orthogonality)
    for j in 1:_δ+1, l in 1:_δ+1
        j == l && continue
        @test abs(_M[j, l, 1]) < 1e-10
    end

end

# ══════════════════════════════════════════════════════════════════════════════
# B. Angle PCE coefficients
# ══════════════════════════════════════════════════════════════════════════════
@testset "B. Angle PCE coefficients" begin

    v_min, v_max = _H_BOUNDS[5]
    phi = affine_angle_pce(v_min, v_max, _basis)

    # k=0 coefficient is the mean of the scaled Beta, lying in (v_min, v_max)
    @test v_min < phi[1] < v_max

    # mean matches the analytical formula: v_min + (v_max-v_min)·α̂/(α̂+β̂)
    E_ω = _α̂ / (_α̂ + _β̂)
    @test isapprox(phi[1], v_min + (v_max - v_min) * E_ω; atol=1e-10)

    # k=2 coefficient is zero: φ is linear in ω, so ⟨φ, ψ_{k≥2}⟩ = 0
    @test abs(phi[3]) < 1e-10

end

# ══════════════════════════════════════════════════════════════════════════════
# C. Trigonometric PCE coefficients
# ══════════════════════════════════════════════════════════════════════════════
@testset "C. Trigonometric PCE coefficients" begin

    # ------------------------------------------------------------------
    @testset "deterministic angle  (φ_det = 20°)" begin
        φ_det = 20.0 * _d2r
        phi_coeffs = [φ_det; zeros(_δ)]   # [φ_det, 0, 0]

        cos_c, sin_c = trig_pce_coefficients(phi_coeffs, _M; n_terms=4)

        # k=0 coefficient must recover the scalar trig value
        @test isapprox(cos_c[1], cos(φ_det); atol=1e-4)
        @test isapprox(sin_c[1], sin(φ_det); atol=1e-4)

        # no stochastic variation for deterministic input
        for k in 2:_δ+1
            @test abs(cos_c[k]) < 1e-10
            @test abs(sin_c[k]) < 1e-10
        end
    end

    # ------------------------------------------------------------------
    @testset "stochastic h=5 — positive variance" begin
        v_min, v_max = _H_BOUNDS[5]
        phi_c         = affine_angle_pce(v_min, v_max, _basis)
        cos_c, sin_c  = trig_pce_coefficients(phi_c, _M; n_terms=4)

        @test pce_variance(cos_c, _basis) > 0
        @test pce_variance(sin_c, _basis) > 0
    end

    # ------------------------------------------------------------------
    @testset "stochastic h=5 — Monte Carlo moment comparison" begin
        Random.seed!(42)
        N = 100_000

        v_min, v_max  = _H_BOUNDS[5]
        ω_samples     = rand(Beta(_α̂, _β̂), N)
        phi_samples   = @. v_min + (v_max - v_min) * ω_samples
        cos_samples   = cos.(phi_samples)

        phi_c         = affine_angle_pce(v_min, v_max, _basis)
        cos_c, _      = trig_pce_coefficients(phi_c, _M; n_terms=4)

        mc_mean  = Statistics.mean(cos_samples)
        mc_var   = Statistics.var(cos_samples; corrected=false)
        pce_mean = pce_expectation(cos_c)
        pce_var  = pce_variance(cos_c, _basis)

        # expectation: degree-2 PCE is highly accurate for the mean (< 5 %)
        @test abs(pce_mean - mc_mean) / abs(mc_mean) < 0.05

        # variance: degree-2 PCE has an inherent ~10 % truncation bias relative
        # to the true variance; 15 % covers both truncation and MC noise
        @test abs(pce_var - mc_var) / abs(mc_var) < 0.15

    end

end

# ══════════════════════════════════════════════════════════════════════════════
# D. Cantelli bound
# ══════════════════════════════════════════════════════════════════════════════
@testset "D. Cantelli bound" begin

    # √(0.95/0.05) = √19
    @test isapprox(cantelli_lambda(0.05), sqrt(19.0); atol=1e-12)
    @test isapprox(cantelli_lambda(0.05), 4.358898944; atol=1e-6)

    # √(0.90/0.10) = √9 = 3
    @test isapprox(cantelli_lambda(0.10), 3.0; atol=1e-12)

    # √(0.50/0.50) = 1
    @test isapprox(cantelli_lambda(0.50), 1.0; atol=1e-12)

end
