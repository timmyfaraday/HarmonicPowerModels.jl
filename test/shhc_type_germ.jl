################################################################################
# test/shhc_type_germ.jl                                                       #
# Shared type germs: a group of units (a converter type) shares one extra germ  #
# that adds a common-mode angle shift, e.g., a common firing delay.             #
#                                                                               #
# Covers the bivariate trig projection (affine2_trig_pce), the data layer       #
# (basis size, physical index set, unit→type map) and the sHHC_SOC / sHHC_QCQP  #
# models on the N=4 industrial network.                                         #
################################################################################

@testset "SHHC — Shared type germs" begin

    α̂, β̂ = 4.0414, 4.7762
    d2r   = π / 180

    # ── 1. Bivariate projection ──────────────────────────────────────────────
    @testset "affine2_trig_pce: moments vs Monte Carlo" begin
        pb = HPM.PCEBasis(3, α̂, β̂, 2)
        v_min, v_max = 4.75d2r, 40.05d2r          # unit's own germ (dimension 1)
        s_min, s_max = -5 * 16.0d2r, -5 * 46.0d2r # type germ (dimension 3), reversed
        cc, sc = HPM.affine2_trig_pce(pb, 1, v_min, v_max, 3, s_min, s_max)
        # only indices supported on germs 1 and 3 are non-zero
        for k in 1:pb.dim
            pb.ind[k, 2] != 0 && @test cc[k] == 0.0 && sc[k] == 0.0
        end
        Random.seed!(7); B = Beta(α̂, β̂); n = 400_000
        φ = [v_min + (v_max - v_min) * rand(B) + s_min + (s_max - s_min) * rand(B) for _ in 1:n]
        for (c, f) in ((cc, cos), (sc, sin))
            x = f.(φ)
            @test isapprox(c[1], Statistics.mean(x); atol = 3e-3)
            var_pce = sum(pb.sp2[k] * c[k]^2 for k in 2:pb.dim)
            @test isapprox(var_pce, Statistics.var(x); rtol = 5e-2)    # degree-2 truncation
        end
    end

    @testset "affine2_trig_pce: degenerate type shift = shifted univariate" begin
        pb = HPM.PCEBasis(2, α̂, β̂, 2)
        v_min, v_max, s = -22.77d2r, 27.51d2r, -7 * 30.0d2r
        c2, s2 = HPM.affine2_trig_pce(pb, 1, v_min, v_max, 2, s, s)
        c1, s1 = HPM.affine_trig_pce(pb, 1, v_min + s, v_max + s)
        @test maximum(abs.(c2 .- c1)) < 1e-12
        @test maximum(abs.(s2 .- s1)) < 1e-12
    end

    # ── 2. Data layer and models on the industrial network ───────────────────
    path  = joinpath(HPM.BASE_DIR, "test", "data", "matpower", "industrial_network_hhc.m")
    _data = HPM.parse_file(path)
    for (_, bus) in _data["bus"]
        bus["ref_angle"] = 0.0
    end
    _data["principle"] = "absolute equality"
    H_TEST = [1, 5, 7]
    fresh() = HPM.replicate(deepcopy(_data); H = H_TEST)
    load_ids = sort(parse.(Int, collect(keys(_data["load"]))))

    harm = Dict{Int,Any}(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r))
    params(; types = nothing, h = harm) = begin
        p = Dict{String,Any}("alpha_hat" => α̂, "beta_hat" => β̂, "delta" => 2,
                             "epsilon" => 0.10, "bound" => "gaussian",
                             "harmonics" => deepcopy(h))
        types === nothing || (p["type_germs"] = types)
        p
    end
    # firing delay 30° mean, 30° wide: harmonic h shifts by −h·α_f
    shift(αc, w) = Dict(h => (-h * (αc - w * α̂ / (α̂ + β̂)) * d2r,
                              -h * (αc + w * (1 - α̂ / (α̂ + β̂))) * d2r) for h in (5, 7))
    typeB(αc, w) = Dict("B" => Dict("units" => load_ids[[2, 4]], "shift" => shift(αc, w)))

    @testset "data layer: basis, index set, type map" begin
        hd = HPM.update_hdata_with_pce_data!(fresh(), params(types = typeB(30.0, 30.0)))
        pce = hd["pce"]
        @test pce["basis"].N == 5                      # 4 own germs + 1 type germ
        @test pce["basis"].dim == binomial(5 + 2, 2)
        @test pce["type_dim"]["B"] == 5
        @test pce["unit_type"] == Dict(load_ids[2] => "B", load_ids[4] => "B")
        ind = pce["basis"].ind
        mixed = [k - 1 for k in 1:pce["basis"].dim
                 if count(!iszero, view(ind, k, :)) == 2 && ind[k, 5] != 0]
        own2  = pce["germ_dim"][load_ids[2]]
        @test any(ind[k+1, own2] != 0 for k in mixed if k in pce["K_phys"])   # (own, type) kept
        own1  = pce["germ_dim"][load_ids[1]]
        @test all(!(k in pce["K_phys"]) for k in mixed if ind[k+1, own1] != 0) # untyped unit: dropped
    end

    @testset "sHHC_SOC: degenerate type shift = per-unit shifted ranges" begin
        # a point-mass type shift must reproduce the same ranges given per unit
        types = Dict("B" => Dict("units" => load_ids[[2, 4]],
                                 "shift" => Dict(5 => (-5 * 30d2r, -5 * 30d2r), 7 => (-7 * 30d2r, -7 * 30d2r))))
        sol_t = HPM.solve_shhc(fresh(), sHHC_SOC, solver_soc; pce_params = params(types = types),
                               hpf_optimizer = solver_nlp)
        hp = Dict{Int,Any}(h => Dict(d => (d in load_ids[[2, 4]] ? (harm[h][1] - h * 30d2r, harm[h][2] - h * 30d2r)
                                                                : harm[h]) for d in load_ids) for h in (5, 7))
        sol_u = HPM.solve_shhc(fresh(), sHHC_SOC, solver_soc; pce_params = params(h = hp),
                               hpf_optimizer = solver_nlp)
        # the larger basis (one more germ) leaves Clarabel at ALMOST_OPTIMAL here
        @test sol_t["termination_status"] in (OPTIMAL, ALMOST_OPTIMAL)
        for h in (5, 7)
            @test isapprox(sol_t["solution"][h], sol_u["solution"][h]; rtol = 1e-4)
        end
    end

    @testset "sHHC_SOC ≈ sHHC_QCQP with a type germ" begin
        p = params(types = typeB(30.0, 30.0))
        sol_s = HPM.solve_shhc(fresh(), sHHC_SOC, solver_soc; pce_params = deepcopy(p),
                               hpf_optimizer = solver_nlp)
        sol_q = HPM.solve_shhc(fresh(), sHHC_QCQP, solver_nlp; pce_params = deepcopy(p))
        @test sol_s["termination_status"] in (OPTIMAL, ALMOST_OPTIMAL)
        tot_s = sum(values(sol_s["solution"]))
        tot_q = sum(sol_q["solution"]["nw"]["$h"]["load"]["$(load_ids[1])"]["cmd"] for h in (5, 7))
        @test isapprox(tot_s, tot_q; rtol = 1e-2)
    end
end
