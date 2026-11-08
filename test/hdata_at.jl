@testset "Network model solver data" begin
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    hd = HPM.replicate(data; H = [1, 5, 7, 11, 13])
    original = deepcopy(hd)
    box = ParameterBox(:kR_cable => (1.0, 0.9, 1.1), :kC_cable => (1.0, 0.9, 1.1))
    scaled = ScaledBranchModel(hd, Dict(2 => :cable, 3 => :cable))
    # Sequence-dependent per-meter data exercise the exact π-model without LCM.
    zy(p, f, seq) = (p.kR_cable * 2e-4 * sqrt(f / 50) + im * 3e-4 * f / 50,
                     im * p.kC_cable * (seq == :zero ? 1e-9 : seq == :negative ? 3e-9 : 2e-9) * f / 50)
    exact = ExactPiModel(Dict(2 => (zy = zy, len = 1000.0, zbase = 10.0)))
    pp = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "gaussian",
        "harmonics" => Dict(5 => (4.75, 40.05) .* (π / 180),
                            7 => (-22.77, 27.51) .* (π / 180),
                            11 => (79.69, 164.01) .* (π / 180),
                            13 => (49.38, 147.67) .* (π / 180)))

    @testset "All networks and untouched data" for model in (exact, scaled)
        for p in (nominal_point(box), (kR_cable = 1.1, kC_cable = 0.9))
            result = hdata_at(hd, model, p)
            expected = deepcopy(hd)
            for (nw, ntw) in expected["nw"]
                h = parse(Int, nw)
                seq = nw == "1" ? :positive : harmonic_sequence(h)
                for b in uncertain_branches(model)
                    for (k, v) in pairs(branch_data(model, b, p, h, seq))
                        ntw["branch"][string(b)][string(k)] = v
                    end
                    ntw["branch"][string(b)]["sequence_data"] = seq
                end
            end
            @test isequal(result, expected)
            # Exercise nested data outside the overridden branch table as well.
            result["nw"]["1"]["bus"]["1"]["vm"] = -1.0
            result["nw"]["5"]["branch"]["4"]["br_r"] = -1.0
            @test isequal(hd, original)
        end
    end

    @testset "Solver writes and zero-width RHHC" for model in (exact, scaled)
        for p in (nominal_point(box), (kR_cable = 1.1, kC_cable = 0.9))
            input = hdata_at(hd, model, p)
            hhc = solve_hhc(input, dHHC_SOC, solver_soc, solver_nlp)
            @test hhc["termination_status"] == OPTIMAL
            @test isequal(hd, original)
            hpf = solve_hpf(hdata_at(hd, model, p), PMs.IVRPowerModel, solver_nlp)
            @test hpf["termination_status"] == LOCALLY_SOLVED
            @test isequal(hd, original)
            stochastic = hdata_at(hd, model, p)
            shhc = solve_shhc(stochastic, sHHC_SOC, solver_soc;
                pce_params = deepcopy(pp), hpf_optimizer = solver_nlp)
            @test shhc["termination_status"] == OPTIMAL
            @test haskey(stochastic, "pce")
            @test isequal(hd, original)
            values = [getproperty(p, name) for name in box.names]
            fixed = ParameterBox(box.names, values, values, values, box.scale)
            # RHHC holds its supplied fundamental fixed at this point.
            robust = solve_rhhc(hdata_at(hd, model, p), model, fixed, solver_soc;
                pce_params = deepcopy(pp), hpf_optimizer = solver_nlp,
                nbackground = 0, use_seeds = false)
            @test robust["signed_off"]
            @test isapprox(shhc["objective_value"], robust["upper_bound"]; rtol = 1e-6)
            @test isapprox(shhc["objective_value"], robust["lower_bound"]; rtol = 1e-6)
        end
    end

    @testset "Unit scaling preserves deterministic capacity exactly" begin
        reference = solve_hhc(deepcopy(hd), dHHC_SOC, solver_soc, solver_nlp)
        result = solve_hhc(hdata_at(hd, scaled, nominal_point(box)), dHHC_SOC,
                           solver_soc, solver_nlp)
        @test result["termination_status"] == reference["termination_status"] == OPTIMAL
        @test result["objective"] == reference["objective"]
        @test isequal(hd, original)
    end

    @testset "Zero sequence" begin
        zero = HPM.replicate(deepcopy(data); H = [1, 3])
        p = nominal_point(box)
        result = hdata_at(zero, exact, p)
        for (k, v) in pairs(branch_data(exact, 2, p, 3, :zero))
            @test result["nw"]["3"]["branch"]["2"][string(k)] == v
        end
        @test result["nw"]["3"]["branch"]["2"]["sequence_data"] == :zero
    end
end
