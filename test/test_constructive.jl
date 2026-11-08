################################################################################
# test/test_constructive.jl                                                    #
# Unit tests for the constructive worst-case angle utilities:                  #
#   reduced_harmonic_impedance, constructive_worst_case_angles,                #
#   recentre_angle_ranges                                                       #
################################################################################

@testset "Constructive worst-case angle utilities" begin

    path = joinpath(HPM.BASE_DIR, "examples", "shhc", "data",
                    "industrial_network_shhc.m")
    data = HPM.parse_file(path)
    H    = [1, 5, 7, 11, 13]
    harm = [5, 7, 11, 13]
    hd   = HPM.replicate(deepcopy(data), H = H)

    Z, bus_pos = reduced_harmonic_impedance(hd)

    # ── A. Z built for every harmonic, square over non-reference buses ──────────
    @test sort(collect(keys(Z))) == harm
    nref = length(bus_pos)
    for h in harm
        @test size(Z[h]) == (nref, nref)
    end
    # slack bus (id 0, bus_type 3) is eliminated → not in bus_pos
    @test !haskey(bus_pos, 0)

    # ── B. Constructive angles align every unit IN PHASE at the binding bus ─────
    nstar = 8
    ang   = constructive_worst_case_angles(Z, bus_pos, hd, nstar)
    @test sort(collect(keys(ang))) == harm
    np = bus_pos[nstar]
    for h in harm, (d_str, load) in hd["nw"]["$h"]["load"]
        d  = parse(Int, d_str)
        bp = get(bus_pos, load["load_bus"], 0)
        bp == 0 && continue
        # Z[n*,d]·exp(i·φ_d) must be real & positive (all contributions parallel)
        contrib = Z[h][np, bp] * exp(im * ang[h][d])
        @test isapprox(imag(contrib), 0.0; atol = 1e-9)
        @test real(contrib) > 0.0
    end

    # ── C. recentre_angle_ranges puts the Beta mean exactly on the centre ───────
    α̂, β̂ = 4.0414, 4.7762
    Eω    = α̂ / (α̂ + β̂)
    base  = Dict(5 => (deg2rad(4.75), deg2rad(40.05)),
                 7 => (deg2rad(-22.77), deg2rad(27.51)),
                 11 => (deg2rad(79.69), deg2rad(164.01)),
                 13 => (deg2rad(49.38), deg2rad(147.67)))
    rr = recentre_angle_ranges(base, ang, α̂, β̂)
    for h in harm, (d, (lo, hi)) in rr[h]
        @test isapprox(lo + (hi - lo) * Eω, ang[h][d]; atol = 1e-9)
        # width is preserved from the section-II range
        @test isapprox(hi - lo, base[h][2] - base[h][1]; atol = 1e-9)
    end

end

@testset "Empirical λ* chance constraint" begin
    path = joinpath(HPM.BASE_DIR, "examples", "shhc", "data",
                    "industrial_network_shhc.m")
    data = HPM.parse_file(path)
    H = [1, 5, 7, 11, 13]; harm = [5, 7, 11, 13]
    d2r = π/180
    pce_params = Dict{String,Any}(
        "alpha_hat"=>4.0414, "beta_hat"=>4.7762, "delta"=>2, "n_terms"=>4, "epsilon"=>0.05,
        "harmonics"=>Dict(5=>(4.75d2r,40.05d2r), 7=>(-22.77d2r,27.51d2r),
                          11=>(79.69d2r,164.01d2r), 13=>(49.38d2r,147.67d2r)))

    # empirical_lambdas needs pce data on hdata
    hd = HPM.replicate(deepcopy(data), H = H)
    HPM.update_hdata_with_pce_data!(hd, pce_params)
    g  = Dict(h => 1e-3 for h in harm)
    λs = empirical_lambdas(hd, g; nsamp = 3000, seed = 1)
    @test Set(keys(λs)) == Set([:ihd, :agg, :rmsi, :q, :g])
    @test all(Set(keys(λs[:q][k])) == Set(keys(λs[k])) for k in (:ihd, :agg, :rmsi))
    @test all(0.0 .≤ values(λs[:ihd]))             # λ* ≥ 0
    @test all(v -> v < HPM.cantelli_lambda(0.05), values(λs[:ihd]))   # tighter than Cantelli
    @test all(isfinite, values(λs[:agg]))

    # iterative calibrated solve converges and is feasible
    sol = HPM.solve_shhc_calibrated(HPM.replicate(deepcopy(data), H=H), sHHC_SOC, solver_soc;
                                    pce_params = deepcopy(pce_params), hpf_optimizer = solver_nlp,
                                    nsamp = 3000, maxiter = 4)
    @test string(sol["termination_status"]) in ("OPTIMAL", "ALMOST_OPTIMAL")
    @test isfinite(sol["objective_value"]) && sol["objective_value"] > 0
    @test sol["iterations"] ≥ 1
    @test haskey(sol, "empirical_lambdas")

    # the default calibrates on the moments of the expansion, so that every IHD bound
    # equals the sampled (1−ε) quantile of |V|²/I²; the sampled moments stay close
    @test !isempty(sol["lambdas_used"][:ihd])
    sol_p = HPM.solve_shhc_calibrated(HPM.replicate(deepcopy(data), H=H), sHHC_SOC, solver_soc;
                                      pce_params = merge(deepcopy(pce_params), Dict("empirical_moments" => "sample")),
                                      hpf_optimizer = solver_nlp, nsamp = 3000, maxiter = 4)
    @test string(sol_p["termination_status"]) in ("OPTIMAL", "ALMOST_OPTIMAL")
    @test isapprox(sol_p["objective_value"], sol["objective_value"]; rtol = 0.05)
end

@testset "Per-unit shapes in PCE data and empirical sampling" begin
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "shhc", "data",
                                   "industrial_network_shhc.m"))
    units = sort(parse.(Int, collect(keys(data["load"]))))[[1, end]]
    a, b = 4.0414, 4.7762
    params = Dict{String,Any}("alpha_hat" => a, "beta_hat" => b, "delta" => 2,
        "epsilon" => 0.05, "harmonic_units" => units, "angle_centre" => "ref_angle",
        "harmonics" => Dict(5 => Dict(units[1] => (0.1, 0.8), units[2] => (-1.2, 0.2))),
        "unit_shapes" => Dict(units[2] => (4.139, 3.659)),
        "type_germs" => Dict("B" => Dict("units" => [units[2]], "shift" => Dict(5 => (-0.2, 0.3)))))
    hd = HPM.replicate(deepcopy(data); H = [1, 5])
    HPM.update_hdata_with_pce_data!(hd, params)
    pce = hd["pce"]
    @test pce["germ_shapes"] == [(a, b), (4.139, 3.659), (a, b)]
    for d in units
        u = pce["germ_dim"][d]
        α, β = pce["germ_shapes"][u]
        lo, hi = pce["angle_ranges"][5][d]
        load = hd["nw"]["5"]["load"][string(d)]
        centre = deg2rad(get(load, "ref_angle", 0.0) + get(load, "angle_offset", 0.0))
        @test isapprox(lo + (hi - lo) * α / (α + β), centre; atol = 1e-14)
        @test isapprox(hi - lo, params["harmonics"][5][d][2] - params["harmonics"][5][d][1]; atol = 1e-14)
        cc, sc = if d == units[2]
            HPM.affine2_trig_pce(pce["basis"], u, lo, hi, pce["type_dim"]["B"], -0.2, 0.3)
        else
            HPM.affine_trig_pce(pce["basis"], u, lo, hi)
        end
        @test pce["cos_phi"][5][d] == cc
        @test pce["sin_phi"][5][d] == sc
    end
    # Explicit default unit shapes and the legacy absence of germ_shapes must
    # preserve the sampled empirical factors, including the shared type stream.
    params["unit_shapes"] = Dict(d => (a, b) for d in units)
    HPM.update_hdata_with_pce_data!(hd, params)
    λs = empirical_lambdas(hd, Dict(5 => 1e-3); nsamp = 100, seed = 19)
    delete!(hd["pce"], "germ_shapes")
    @test empirical_lambdas(hd, Dict(5 => 1e-3); nsamp = 100, seed = 19) == λs
    delete!(params, "unit_shapes")
    HPM.update_hdata_with_pce_data!(hd, params)
    @test empirical_lambdas(hd, Dict(5 => 1e-3); nsamp = 100, seed = 19) == λs
end

@testset "Empirical sampling kernel: per-unit Beta shapes" begin
    # Compare the typed real arithmetic with direct complex network responses,
    # using the same scalar draw stream (unit germs first, then shared type germs).
    nsamp, npos = 40, 2
    Bd = Beta(4.0414, 4.7762)
    Bu = [Bd, Beta(4.139, 3.659), Bd]
    dp, w, tix = [1, 2, 0], [1.0, 0.7, 0.5], [1, 1, 2]
    lo = [0.1 -0.5 0.0; -0.3 0.8 0.0]
    hi = [0.8 0.2 1.0; 0.4 -0.2 1.0]   # also a reversed support
    slo = [-0.2 0.1; -0.4 0.2]; shi = [0.3 0.5; 0.2 0.6]
    Zs = [ComplexF64[1+2im 0.3-0.2im; 0.3-0.2im 2+im],
          ComplexF64[0.4+im 0.1+0.2im; 0.1+0.2im 1+3im]]
    brs = [[(1, 1, 2, 0.3, -0.7, 0.9, 0.1),
            (2, 0, 1, 0.2, -0.5, 1.0, 0.0),
            (3, 2, 0, 0.1, -0.4, 1.0, 0.0)] for _ in Zs]
    aV = [zeros(nsamp, npos) for _ in Zs]
    aI = [zeros(nsamp, 3) for _ in Zs]
    rng = MersenneTwister(17); ref_rng = MersenneTwister(17)
    HPM._sample_aggregates!(aV, aI, rng, Bu, Bd, dp, lo, hi, w, tix,
                            slo, shi, Zs, brs, nsamp, npos)
    refV = [zeros(nsamp, npos) for _ in Zs]
    refI = [zeros(nsamp, 3) for _ in Zs]
    for m in 1:nsamp
        ωu = [rand(ref_rng, B) for B in Bu]
        ωt = [rand(ref_rng, Bd) for _ in 1:2]
        for (i, Z) in enumerate(Zs)
            V = zeros(ComplexF64, npos)
            for u in eachindex(dp)
                dp[u] == 0 && continue
                t = tix[u]
                φ = lo[i, u] + (hi[i, u] - lo[i, u]) * ωu[u] +
                    slo[i, t] + (shi[i, t] - slo[i, t]) * ωt[t]
                V .+= Z[:, dp[u]] .* (w[u] * cis(φ))
            end
            refV[i][m, :] .= abs2.(V)
            for (j, fp, tp, gy, by, cr, ci) in brs[i]
                vf = fp == 0 ? 0.0im : V[fp]
                vt = tp == 0 ? 0.0im : V[tp]
                refI[i][m, j] = abs2(complex(gy, by) * (complex(cr, ci) * vf - vt))
            end
        end
    end
    @test all(isapprox(a, b; rtol = 1e-13, atol = 1e-14) for (a, b) in zip(aV, refV))
    @test all(isapprox(a, b; rtol = 1e-13, atol = 1e-14) for (a, b) in zip(aI, refI))
    @test rand(rng) == rand(ref_rng)
end
