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
    @test Set(keys(λs)) == Set([:ihd, :agg, :rmsi])
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
end
