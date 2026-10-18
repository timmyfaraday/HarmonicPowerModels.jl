################################################################################
# test/test_pce.jl                                                              #
# Tests for src/core/pce_basis.jl (SHHC polynomial chaos workspace).           #
# Runnable standalone: julia --project=. test/test_pce.jl                      #
################################################################################

using Test
using Statistics
using Random
using Distributions
using HarmonicPowerModels

const _HPMt = HarmonicPowerModels

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
const _pb = _HPMt.PCEBasis(1, _α̂, _β̂, _δ)
const _Mt = _HPMt.dense_tensor(_pb)

# ══════════════════════════════════════════════════════════════════════════════
# A. Multiplication tensor
# ══════════════════════════════════════════════════════════════════════════════
@testset "A. Multiplication tensor" begin

    @test size(_Mt) == (_δ+1, _δ+1, _δ+1)

    @test all(isfinite, _Mt)

    # M[1,1,1] = ⟨ψ_0·ψ_0, ψ_0⟩ / ⟨ψ_0,ψ_0⟩ = 1  (ψ_0 = 1)
    @test isapprox(_Mt[1, 1, 1], 1.0; atol=1e-10)

    # probability measure: ⟨ψ_0,ψ_0⟩ = 1
    @test isapprox(_pb.sp2[1], 1.0; atol=1e-13)

    # symmetry: M[j,l,k] = M[l,j,k]
    for j in 1:_δ+1, l in 1:_δ+1, k in 1:_δ+1
        @test isapprox(_Mt[j, l, k], _Mt[l, j, k]; atol=1e-10)
    end

    # k=0 slice is diagonal: ⟨ψ_j, ψ_l⟩ = 0 for j≠l (orthogonality)
    for j in 1:_δ+1, l in 1:_δ+1
        j == l && continue
        @test abs(_Mt[j, l, 1]) < 1e-10
    end

end

# ══════════════════════════════════════════════════════════════════════════════
# B. Angle PCE coefficients
# ══════════════════════════════════════════════════════════════════════════════
@testset "B. Angle PCE coefficients" begin

    v_min, v_max = _H_BOUNDS[5]
    phi = _HPMt.affine_pce(_pb, 1, v_min, v_max)

    # k=0 coefficient is the mean of the scaled Beta, lying in (v_min, v_max)
    @test v_min < phi[1] < v_max

    # mean matches the analytical formula: v_min + (v_max-v_min)·α̂/(α̂+β̂)
    E_ω = _α̂ / (_α̂ + _β̂)
    @test isapprox(phi[1], v_min + (v_max - v_min) * E_ω; atol=1e-10)

    # k=2 coefficient is zero: φ is linear in ω, so ⟨φ, ψ_{k≥2}⟩ = 0
    @test abs(phi[3]) < 1e-10

    # variance matches the analytical Beta variance
    σ²_ω = _α̂ * _β̂ / ((_α̂ + _β̂)^2 * (_α̂ + _β̂ + 1))
    @test isapprox(_HPMt.pce_variance(phi, _pb), (v_max - v_min)^2 * σ²_ω;
                   rtol=1e-10)

end

# ══════════════════════════════════════════════════════════════════════════════
# C. Trigonometric PCE coefficients
# ══════════════════════════════════════════════════════════════════════════════
@testset "C. Trigonometric PCE coefficients" begin

    # ------------------------------------------------------------------
    @testset "deterministic angle  (φ_det = 20°)" begin
        φ_det = 20.0 * _d2r
        phi_coeffs = [φ_det; zeros(_δ)]   # [φ_det, 0, 0]

        cos_c, sin_c = _HPMt.trig_pce(_pb, phi_coeffs; n_terms=4)

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
        cos_c, sin_c = _HPMt.affine_trig_pce(_pb, 1, v_min, v_max)

        @test _HPMt.pce_variance(cos_c, _pb) > 0
        @test _HPMt.pce_variance(sin_c, _pb) > 0

        # rotated_trig_pce is a deprecated alias and must stay bit-identical
        cos_r, sin_r = _HPMt.rotated_trig_pce(_pb, 1, v_min, v_max; n_terms=4)
        @test cos_r == cos_c
        @test sin_r == sin_c
    end

    # ------------------------------------------------------------------
    @testset "affine_trig_pce — point mass and germ support" begin
        # degenerate range ⇒ deterministic angle: constant term only
        φ = 0.3
        c, s = _HPMt.affine_trig_pce(_pb, 1, φ, φ)
        @test c[1] ≈ cos(φ)
        @test s[1] ≈ sin(φ)
        @test _HPMt.pce_variance(c, _pb) == 0.0
        @test _HPMt.pce_variance(s, _pb) == 0.0

        # a two-germ basis: coefficients must be supported on the named germ only
        pb2 = _HPMt.PCEBasis(2, _α̂, _β̂, _δ)
        c2, s2 = _HPMt.affine_trig_pce(pb2, 2, 0.1, 0.5)
        for k in findall(!iszero, c2) ∪ findall(!iszero, s2)
            @test pb2.ind[k, 1] == 0        # nothing on germ 1
        end
    end

    # ------------------------------------------------------------------
    @testset "stochastic — Monte Carlo moment comparison" begin
        Random.seed!(42)
        N = 1_000_000

        for h in (5, 7)
            v_min, v_max = _H_BOUNDS[h]
            ω_samples   = rand(Beta(_α̂, _β̂), N)
            phi_samples = @. v_min + (v_max - v_min) * ω_samples
            cos_samples = cos.(phi_samples)

            cos_c, _ = _HPMt.affine_trig_pce(_pb, 1, v_min, v_max)

            mc_mean  = Statistics.mean(cos_samples)
            mc_var   = Statistics.var(cos_samples; corrected=false)
            pce_mean = _HPMt.pce_mean(cos_c)
            pce_var  = _HPMt.pce_variance(cos_c, _pb)

            # affine_trig_pce projects cos/sin onto the basis exactly by
            # quadrature, so the only error left is the degree-2 truncation: the
            # mean is essentially exact (~1e-5 relative, MC-noise limited) and the
            # variance lands within ~2e-3. These tolerances lock in the accuracy
            # gained over the superseded Taylor route, which needed 1 % and 2 %.
            @test abs(pce_mean - mc_mean) / abs(mc_mean) < 1e-4
            @test abs(pce_var - mc_var) / abs(mc_var)    < 5e-3
        end

    end

end

# ══════════════════════════════════════════════════════════════════════════════
# D. Cantelli bound
# ══════════════════════════════════════════════════════════════════════════════
@testset "D. Cantelli bound" begin

    # √(0.95/0.05) = √19
    @test isapprox(_HPMt.cantelli_lambda(0.05), sqrt(19.0); atol=1e-12)
    @test isapprox(_HPMt.cantelli_lambda(0.05), 4.358898944; atol=1e-6)

    # √(0.90/0.10) = √9 = 3
    @test isapprox(_HPMt.cantelli_lambda(0.10), 3.0; atol=1e-12)

    # √(0.50/0.50) = 1
    @test isapprox(_HPMt.cantelli_lambda(0.50), 1.0; atol=1e-12)

end

# ══════════════════════════════════════════════════════════════════════════════
# E. Multivariate basis (independent germs)
# ══════════════════════════════════════════════════════════════════════════════
@testset "E. Multivariate basis" begin

    @testset "invariants and affine exactness (N=4)" begin
        pb = _HPMt.PCEBasis(4, _α̂, _β̂, _δ)
        @test pb.dim == binomial(4 + _δ, _δ)
        @test isapprox(pb.sp2[1], 1.0; atol=1e-13)   # probability measure
        @test all(pb.ind[1, :] .== 0)                 # row 1 = constant term

        M = _HPMt.dense_tensor(pb)
        @test isapprox(M[1, 1, 1], 1.0; atol=1e-12)
        @test M ≈ permutedims(M, (2, 1, 3))           # j↔l symmetry

        # affine PCE in germ dimension 2 reproduces analytic Beta moments
        v_min, v_max = -0.3, 0.7
        c    = _HPMt.affine_pce(pb, 2, v_min, v_max)
        μ_ω  = _α̂ / (_α̂ + _β̂)
        σ²_ω = _α̂ * _β̂ / ((_α̂ + _β̂)^2 * (_α̂ + _β̂ + 1))
        @test isapprox(_HPMt.pce_mean(c), v_min + (v_max - v_min) * μ_ω;
                       rtol=1e-12)
        @test isapprox(_HPMt.pce_variance(c, pb), (v_max - v_min)^2 * σ²_ω;
                       rtol=1e-12)
    end

    @testset "sparse Galerkin product vs dense tensor (N=21)" begin
        pb = _HPMt.PCEBasis(21, _α̂, _β̂, _δ)
        @test pb.dim == 253

        Random.seed!(7)
        M    = _HPMt.dense_tensor(pb)
        x, y = randn(pb.dim), randn(pb.dim)
        z_sp = _HPMt.galerkin_product(pb, x, y)
        z_dn = [sum(M[j, l, k] * x[j] * y[l]
                    for j in 1:pb.dim, l in 1:pb.dim) for k in 1:pb.dim]
        @test isapprox(z_sp, z_dn; rtol=1e-10)
    end

end
