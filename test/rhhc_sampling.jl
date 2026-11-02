################################################################################
# test/rhhc_sampling.jl                                                        #
# Sampling the parameter box on the hello-world network (one parallel         #
# resonance near the 7th harmonic): modal scan, sensitivities, seeds,         #
# background and the compact dictionary.                                      #
################################################################################

@testset "RHHC — Sampling" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    hd   = HPM.replicate(data; H = [1, 5, 7, 11, 13])
    pp   = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli",
        "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r),
                            11 => (79.69d2r, 164.01d2r), 13 => (49.38d2r, 147.67d2r)))
    sol  = HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = pp, hpf_optimizer = solver_nlp)
    @test sol["termination_status"] == OPTIMAL

    model = ScaledBranchModel(hd, Dict(2 => :cable, 3 => :cable, 4 => :line))
    box   = ParameterBox(:kC_cable => (1.0, 0.9, 1.1), :kL_cable => (1.0, 0.95, 1.05),
                         :kR_cable => (1.0, 0.8, 1.2))

    @testset "modal scan and sensitivities" begin
        res = modal_scan(hd, model, nominal_point(box))
        @test length(res) == 1
        r = only(res)
        @test 7.0 < r.h < 7.5
        @test transpose(r.t) * r.t ≈ 1.0

        # h ∝ 1/√C: ±10 % cable capacitance shifts the resonance by ∓h/2·10 %
        s = resonance_sensitivity(hd, model, box, r)
        @test isapprox(s[1], -r.h / 2 * 0.1; rtol = 0.1)
        @test abs(s[2]) < abs(s[1])                   # the source inductance dominates
        @test abs(s[3]) < 1e-2                        # resistance damps, it does not tune
        @test locate_resonance(hd, model, (kC_cable = 1.1,), r).h < r.h
    end

    seeds = resonance_seeds(hd, model, box)

    @testset "seeds" begin
        @test first(seeds).h == 7 && first(seeds).crosses
        @test all(!s.crosses for s in seeds if s.h != 7)
        @test all(all(-1 .<= s.ξ .<= 1) for s in seeds)
        # the crossing seed lowers the critical eigenvalue at the 7th harmonic
        s7 = first(seeds)
        η0 = abs(first(HPM._mode_at(hd, model, box, zeros(3), 7, s7.r.t)))
        η7 = abs(first(HPM._mode_at(hd, model, box, s7.ξ, 7, s7.r.t)))
        @test η7 < η0
    end

    @testset "background" begin
        b1 = background_points(box, 64; seed = 1)
        @test length(b1) == 64 && all(all(-1 .<= ξ .<= 1) for ξ in b1)
        @test b1 == background_points(box, 64; seed = 1)
        @test b1 != background_points(box, 64; seed = 2)
    end

    @testset "dictionary" begin
        dict = build_dictionary(hd, model, box; seeds, nbackground = 128)
        @test count(==(:vertex), dict.labels) == 8
        @test count(==(:seed), dict.labels) == length(seeds)
        @test dict.labels[1] == :nominal && dict.points[1] == zeros(3)

        g = Dict(h => I^2 for (h, I) in sol["solution"])
        V = dictionary_violations(dict, g)
        pos = dict.limit .> 0
        # the nominal solution is feasible (and binding) at the nominal point …
        @test maximum(V[pos, 1]) ≤ 1e-6
        @test maximum(V[pos, 1]) ≥ -1e-4
        # … but not over the box
        @test maximum(V[pos, :]) > 0.1

        # compact form = full evaluation
        for s in (2, length(dict.points))
            m = moments_at(hd, model, latent(box, dict.points[s]))
            full = [normalized_violation(c, g, m) for c in shhc_constraints(hd, m)]
            @test maximum(abs.(V[:, s] .- full)) < 1e-12
        end

        # the interior worst case: for the 7th harmonic at the end of the cable feeder the
        # crossing seed beats every vertex
        j  = findfirst(==((:ihd, (7, 3))), dict.keys)
        sv = dict.labels .== :seed; vv = dict.labels .== :vertex
        @test maximum(V[j, sv]) > maximum(V[j, vv])
    end
end
