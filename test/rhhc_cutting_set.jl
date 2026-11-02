################################################################################
# test/rhhc_cutting_set.jl                                                     #
# The cutting-set method on the hello-world network: a zero-width box gives   #
# sHHC_SOC, the method converges with valid bounds, and the robust solution   #
# holds on an independent grid of the box.                                     #
################################################################################

@testset "RHHC — Cutting Set" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    H    = [1, 5, 7, 11, 13]
    pp   = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli",
        "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r),
                            11 => (79.69d2r, 164.01d2r), 13 => (49.38d2r, 147.67d2r)))
    nom = HPM.solve_shhc(HPM.replicate(deepcopy(data); H), sHHC_SOC, solver_soc;
                         pce_params = deepcopy(pp), hpf_optimizer = solver_nlp)
    classes = Dict(2 => :cable, 3 => :cable, 4 => :line)

    @testset "zero-width box reproduces sHHC_SOC" begin
        hd  = HPM.replicate(deepcopy(data); H)
        box = ParameterBox(:kC_cable => (1.0, 1.0, 1.0), :kR_cable => (1.0, 1.0, 1.0))
        res = solve_rhhc(hd, ScaledBranchModel(hd, classes), box, solver_soc;
                         pce_params = deepcopy(pp), hpf_optimizer = solver_nlp, nbackground = 8)
        @test res["signed_off"]
        @test res["iterations"] == 1
        @test isapprox(res["upper_bound"], nom["objective_value"]; rtol = 1e-6)
        @test isapprox(res["lower_bound"], res["upper_bound"]; rtol = 1e-6)
    end

    @testset "cable box: convergence, bounds and an independent check" begin
        hd    = HPM.replicate(deepcopy(data); H)
        model = ScaledBranchModel(hd, classes)
        box   = ParameterBox(:kC_cable => (1.0, 0.9, 1.1), :kL_cable => (1.0, 0.95, 1.05),
                             :kR_cable => (1.0, 0.8, 1.2))
        res = solve_rhhc(hd, model, box, solver_soc; pce_params = deepcopy(pp),
                         hpf_optimizer = solver_nlp, nbackground = 64)
        @test res["signed_off"]
        @test res["iterations"] > 1                          # the nominal solution was not robust
        @test res["lower_bound"] ≤ res["upper_bound"] + 1e-12
        @test res["gap"] < 1e-3
        @test res["upper_bound"] < nom["objective_value"]    # robustness costs capacity
        hist = res["history"]
        @test first(hist).max_violation > 0.1
        @test all(hist[k+1].upper_bound ≤ hist[k].upper_bound + 1e-9 for k in 1:length(hist)-1)

        # the robust solution on an independent grid of the box (not used by the method)
        g    = Dict(h => I^2 for (h, I) in res["solution"])
        dict = res["dictionary"]; pos = dict.limit .> 0
        grid = range(-1, 1; length = 7)
        vmax = maximum(maximum(HPM._violations_at(hd, model, box, [a, b, c], g)[pos])
                       for a in grid, b in grid, c in grid)
        @test vmax ≤ 1e-6
    end
end
