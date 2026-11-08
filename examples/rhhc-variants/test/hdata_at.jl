@testset "Other solvers on the physical IEEE30 network" begin
    s = study("cases/ieee30-physical.toml")
    p = nominal_point(s.box)
    hd = hdata_at(s.hd, s.model, p)
    original = deepcopy(s.hd)
    stochastic = solve_shhc(deepcopy(hd), sHHC_SOC, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt)
    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    nominal = solve_rhhc(deepcopy(s.hd), s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        nbackground = 0, use_seeds = false)
    @test stochastic["termination_status"] == OPTIMAL
    # run.jl reports six decimals. Compare its unrounded result as well.
    @test round(nominal["upper_bound"]; digits = 6) == 0.055640
    @test stochastic["objective_value"] ≈ nominal["upper_bound"] rtol = 1e-6
    for solver in (:hhc, :hpf)
        input = deepcopy(hd)
        result = solver == :hhc ? solve_hhc(input, dHHC_SOC, s.conic, s.ipopt) :
                 solve_hpf(input, PowerModels.IVRPowerModel, s.ipopt)
        @test result["termination_status"] == (solver == :hhc ? OPTIMAL : LOCALLY_SOLVED)
        for (nw, ntw) in input["nw"], b in uncertain_branches(s.model)

            h = parse(Int, nw)
            seq = nw == "1" ? :positive : harmonic_sequence(h)
            for (key, value) in pairs(branch_data(s.model, b, p, h, seq))
                @test ntw["branch"][string(b)][string(key)] == value
            end
        end
    end
    @test isequal(s.hd, original)
end
