################################################################################
# Exact mean and covariance of squared harmonic magnitudes                     #
# (`exact_quadratic_moments`, src/core/exact_moments.jl)                        #
################################################################################

using Random, Distributions, Statistics

@testset "Exact moments of squared magnitudes" begin

    @testset "point masses reduce to the deterministic magnitude" begin
        R      = [ComplexF64[1.0 + 0.5im, -0.3 + 2.0im], ComplexF64[0.7, 1.0im]]
        φ      = [[0.3, 1.1], [2.0, -0.4]]
        ranges = [[(φ[h][u], φ[h][u]) for u in 1:2] for h in 1:2]
        m, C   = HarmonicPowerModels.exact_quadratic_moments(R, ranges, [(3.7, 4.4), (4.2, 3.6)])
        @test isapprox(m, [abs2(sum(R[h] .* cis.(φ[h]))) for h in 1:2]; atol = 1e-12)
        @test all(abs.(C) .< 1e-12)
    end

    @testset "a single unit has a deterministic magnitude" begin
        R      = [ComplexF64[2.0 - 1.0im]]
        m, C   = HarmonicPowerModels.exact_quadratic_moments(R, [[(0.0, 11.5)]], [(4.16, 3.55)])
        @test isapprox(m[1], 5.0; atol = 1e-12)
        @test abs(C[1, 1]) < 1e-12
    end

    @testset "agrees with Monte Carlo, supports beyond 360°" begin
        Random.seed!(7)
        U, H   = 4, 3
        shapes = [isodd(u) ? (3.731, 4.439) : (4.1612, 3.5541) for u in 1:U]
        R      = [randn(ComplexF64, U) for _ in 1:H]
        w      = deg2rad.([36.0, 294.7, 660.2])
        ranges = [[(2π * rand(), 0.0) for u in 1:U] for h in 1:H]
        ranges = [[(ranges[h][u][1], ranges[h][u][1] + w[h]) for u in 1:U] for h in 1:H]
        m, C   = HarmonicPowerModels.exact_quadratic_moments(R, ranges, shapes)
        N = 400_000
        ω = hcat([rand(Beta(shapes[u]...), N) for u in 1:U]...)
        Q = hcat([abs2.(sum(R[h][u] .* cis.(ranges[h][u][1] .+ w[h] .* ω[:, u]) for u in 1:U)) for h in 1:H]...)
        @test isapprox(m, vec(mean(Q, dims = 1)); rtol = 1e-2)
        @test isapprox(C, cov(Q); atol = 2e-2 * maximum(abs.(cov(Q))))
    end

    @testset "agrees with a converged expansion" begin
        # two units, one harmonic: at a high degree the expansion of |Y|² is exact
        shapes = [(3.731, 4.439), (4.1612, 3.5541)]
        R      = [ComplexF64[1.0 + 0.2im, 0.6 - 0.8im]]
        ranges = [[(0.2, 0.2 + deg2rad(99.4)), (-1.0, -1.0 + deg2rad(133.8))]]
        m, C   = HarmonicPowerModels.exact_quadratic_moments(R, ranges, shapes)
        pb = HarmonicPowerModels.PCEBasis(shapes, 5)
        cs = [HarmonicPowerModels.affine_trig_pce(pb, u, ranges[1][u]...) for u in 1:2]
        yr = sum(real(R[1][u]) .* cs[u][1] .- imag(R[1][u]) .* cs[u][2] for u in 1:2)
        yi = sum(imag(R[1][u]) .* cs[u][1] .+ real(R[1][u]) .* cs[u][2] for u in 1:2)
        q  = HarmonicPowerModels.galerkin_product(pb, yr, yr) .+ HarmonicPowerModels.galerkin_product(pb, yi, yi)
        @test isapprox(m[1], q[1]; rtol = 1e-6)
        @test isapprox(C[1, 1], sum(pb.sp2[2:end] .* q[2:end] .^ 2); rtol = 1e-4)
    end

    @testset "moment coefficients reproduce the covariance" begin
        C = [2.0 0.5 0.1; 0.5 1.0 0.2; 0.1 0.2 0.5]
        c = HarmonicPowerModels._moment_coefficients([1.0, 2.0, 3.0], C)
        L = reduce(vcat, [ci[2:end]' for ci in c])
        @test [ci[1] for ci in c] == [1.0, 2.0, 3.0]
        @test L * L' ≈ C
    end
end
