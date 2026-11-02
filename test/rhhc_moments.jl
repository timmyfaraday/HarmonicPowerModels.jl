################################################################################
# test/rhhc_moments.jl                                                         #
# Moments and chance constraints at a parameter point: the evaluation used by  #
# the robust oracle agrees with the solved program, and the empirical bound is #
# rejected.                                                                    #
################################################################################

@testset "RHHC — Moments and Constraints" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "shhc", "data",
                                   "industrial_network_shhc.m"))
    data["principle"] = "absolute equality"
    pce(bound) = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762,
        "delta" => 2, "n_terms" => 4, "epsilon" => 0.10, "bound" => bound,
        "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r)))

    for bound in ("cantelli", "gaussian")
        @testset "$bound: evaluation matches the solved program" begin
            hd  = HPM.replicate(deepcopy(data); H = [1, 5, 7])
            sol = HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = pce(bound),
                                 hpf_optimizer = solver_nlp)
            @test sol["termination_status"] == OPTIMAL

            m    = shhc_moments(hd)
            cons = shhc_constraints(hd, m)
            @test !isempty(cons)
            @test Set(c.kind for c in cons) ⊆ Set([:ihd, :thd, :rmsv, :rmsi])

            # at the optimum, g_h = I_h²: every constraint holds and at least one binds
            g = Dict(h => I^2 for (h, I) in sol["solution"])
            v = [normalized_violation(c, g, m) for c in cons if c.limit > 0.0]
            @test maximum(v) ≤ 1e-6
            @test maximum(v) ≥ -1e-4

            # homogeneity: scaling g scales every constraint value
            g2 = Dict(h => 2.0 * x for (h, x) in g)
            @test all(isapprox(constraint_value(c, g2, m), 2.0 * constraint_value(c, g, m);
                               rtol = 1e-12) for c in cons)

            # the nominal point of a network model gives the same moments and constraints
            brs   = sort(parse.(Int, collect(keys(hd["nw"]["1"]["branch"]))))
            model = ScaledBranchModel(hd, Dict(b => :cable for b in brs))
            mp    = moments_at(hd, model, NamedTuple())
            @test mp.alpha_V == m.alpha_V && mp.alpha_I == m.alpha_I
            cp = shhc_constraints(hd, mp)
            @test [c.limit for c in cp] == [c.limit for c in cons]
            @test all(normalized_violation(a, g, mp) == normalized_violation(b, g, m)
                      for (a, b) in zip(cp, cons))

            # another point moves the constraint values at the same g
            mc = moments_at(hd, model, (kC_cable = 1.2, kL_cable = 1.1))
            cc = shhc_constraints(hd, mc)
            @test any(normalized_violation(a, g, mc) != normalized_violation(b, g, m)
                      for (a, b) in zip(cc, cons))
        end
    end

    @testset "empirical bound is rejected" begin
        hd  = HPM.replicate(deepcopy(data); H = [1, 5, 7])
        HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = pce("gaussian"),
                       hpf_optimizer = solver_nlp)
        brs   = sort(parse.(Int, collect(keys(hd["nw"]["1"]["branch"]))))
        model = ScaledBranchModel(hd, Dict(b => :cable for b in brs))
        hd["pce"]["bound"] = :empirical
        @test_throws ArgumentError moments_at(hd, model, NamedTuple())
        hd["pce"]["bound"] = :gaussian
        hd["pce"]["lambdas"] = Dict(:ihd => Dict(), :agg => Dict(), :rmsi => Dict())
        @test_throws ArgumentError moments_at(hd, model, NamedTuple())
    end
end
