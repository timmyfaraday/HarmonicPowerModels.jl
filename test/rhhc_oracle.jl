################################################################################
# test/rhhc_oracle.jl                                                          #
# The pessimizing oracle on the hello-world network, against a dense grid of  #
# the parameter box.                                                           #
################################################################################

@testset "RHHC — Oracle" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    hd   = HPM.replicate(data; H = [1, 5, 7, 11, 13])
    pp   = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli",
        "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r),
                            11 => (79.69d2r, 164.01d2r), 13 => (49.38d2r, 147.67d2r)))
    sol  = HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = pp, hpf_optimizer = solver_nlp)
    g    = Dict(h => I^2 for (h, I) in sol["solution"])

    model = ScaledBranchModel(hd, Dict(2 => :cable, 3 => :cable, 4 => :line))
    box   = ParameterBox(:kC_cable => (1.0, 0.9, 1.1), :kL_cable => (1.0, 0.95, 1.05),
                         :kR_cable => (1.0, 0.8, 1.2))
    seeds = resonance_seeds(hd, model, box)

    # dense grid of the box for the reference maxima
    grid  = range(-1, 1; length = 11)
    G     = [HPM._violations_at(hd, model, box, [a, b, c], g) for a in grid, b in grid, c in grid]
    gridmax(j) = maximum(v[j] for v in G)

    @testset "default effort" begin
        dict = build_dictionary(hd, model, box; seeds, nbackground = 64)
        pre  = vec(maximum(dictionary_violations(dict, g); dims = 2))
        r    = pessimize!(dict, hd, model, g)
        @test all(r.violation .>= pre)                       # never worse than the dictionary
        @test max_violation(r) > 0.1                         # the nominal capacity is not robust
        @test max_violation(r) ≥ gridmax(r.worst) - 1e-9     # the worst constraint is found
        @test all(r.violation[j] ≥ maximum(dictionary_violations(dict, g)[j, :]) - 1e-12
                  for j in eachindex(r.violation))
    end

    @testset "full effort refines interior maxima" begin
        dict = build_dictionary(hd, model, box; seeds, nbackground = 64)
        j    = findfirst(==((:ihd, (7, 3))), dict.keys)
        pre  = maximum(dictionary_violations(dict, g)[j, :])
        r    = pessimize!(dict, hd, model, g; ncons = length(dict.keys), starts = 3, iters = 20)
        @test count(==(:refined), dict.labels) > 0
        @test r.violation[j] > pre                           # the ascent improves on the dictionary
        @test r.violation[j] ≥ gridmax(j) - 1e-9             # and reaches the grid maximum
        ξ = dict.points[r.point[j]]
        @test -1 < ξ[1] < 1                                  # interior in the capacitance
        # the stored point reproduces the reported value
        m = moments_at(hd, model, latent(box, ξ))
        @test normalized_violation(shhc_constraints(hd, m)[j], g, m) ≈ r.violation[j]
    end
end
