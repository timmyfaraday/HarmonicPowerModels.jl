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

    @testset "targeted evaluation equals the full evaluation" begin
        ce = ConstraintEvaluator(hd, model, box)
        for ξ in ([0.0, 0.0, 0.0], [0.7, -0.3, 1.0], [-1.0, 1.0, -0.4])
            full = HPM._violations_at(hd, model, box, ξ, g)
            @test maximum(abs.(full .- [fast_violation(ce, j, ξ, g) for j in eachindex(full)])) < 1e-10
        end
        # the adjoint gradient equals central finite differences, for every constraint kind
        for kind in (:ihd, :thd, :rmsv, :rmsi), ξ in ([0.2, -0.4, 0.6], [-0.7, 0.5, -0.1])
            j = findfirst(c -> c.kind === kind, ce.cons)
            j === nothing && continue
            v, gr = fast_violation_gradient(ce, j, ξ, g)
            @test v ≈ fast_violation(ce, j, ξ, g)
            fdg = [(fast_violation(ce, j, ξ .+ 1e-6 .* (1:3 .== i), g) -
                    fast_violation(ce, j, ξ .- 1e-6 .* (1:3 .== i), g)) / 2e-6 for i in 1:3]
            @test isapprox(gr, fdg; rtol = 1e-4, atol = 1e-8)
        end
        # the cached admittance equals the assembled one, between harmonics too
        nc = NetworkCache(hd, model)
        for h in (5, 6.37)
            seq = isinteger(h) ? harmonic_sequence(Int(h)) : :positive
            p   = latent(box, [0.5, -0.5, 0.2])
            @test fast_admittance(nc, p, h; sequence = seq) ≈ first(admittance_at(hd, model, p, h; sequence = seq))
        end
    end

    @testset "dual branch derivatives" begin
        zy = (p, f, seq) -> (0.0001 * sqrt(f / 50) * p.kR_cable +
                             im * 2π * f * 3e-7 * p.kL_cable,
                             im * 2π * f * 3e-10 * p.kC_cable)
        branches = Dict(b => (zy = zy, len = b == 2 ? p -> 15000p.kL_cable : 10000.0,
                               zbase = 10.0) for b in (2, 3))
        exact = ExactPiModel(branches)
        points = ([0.2, -0.4, 0.6], [-0.7, 0.5, -0.1], [0.0, 0.0, 0.0],
                  [0.7, -0.3, 1.0], [-1.0, 1.0, -0.4])
        for m in (model, exact)
            ce = ConstraintEvaluator(hd, m, box)
            for kind in (:ihd, :thd, :rmsv, :rmsi), ξ in points
                j = findfirst(c -> c.kind === kind, ce.cons)
                @test j !== nothing
                v, gd = fast_violation_gradient(ce, j, ξ, g; branch_derivative = :dual)
                vc, gc = fast_violation_gradient(ce, j, ξ, g)
                @test v === vc
                @test isapprox(gd, gc; rtol = 1e-6)
                fdg = [(fast_violation(ce, j, ξ .+ 1e-6 .* (1:3 .== i), g) -
                        fast_violation(ce, j, ξ .- 1e-6 .* (1:3 .== i), g)) / 2e-6 for i in 1:3]
                @test isapprox(gd, fdg; rtol = 1e-4, atol = 1e-8)
            end
            @test_throws ArgumentError fast_violation_gradient(ce, 1, zeros(3), g;
                                                                branch_derivative = :other)
            # used coordinates in both chunks, with unused padding in between
            entries = [box.names[i] => (box.nominal[i], box.lower[i], box.upper[i]) for i in 1:2]
            append!(entries, [Symbol("unused_", i) => (1.0, 0.5, 1.5) for i in 3:32])
            push!(entries, box.names[3] => (box.nominal[3], box.lower[3], box.upper[3]))
            wide = ParameterBox(entries...)
            cw = ConstraintEvaluator(hd, m, wide)
            ξ = collect(range(-0.7, 0.6; length = 33))
            for kind in (:ihd, :thd, :rmsv, :rmsi)
                j = findfirst(c -> c.kind === kind, cw.cons)
                @test isapprox(last(fast_violation_gradient(cw, j, ξ, g; branch_derivative = :dual)),
                               last(fast_violation_gradient(cw, j, ξ, g)); rtol = 1e-6)
            end
        end

        calls = Ref(0)
        counted = function (p, f, seq)
            any(x -> x isa HPM.ForwardDiff.Dual, values(p)) && (calls[] += 1)
            zy(p, f, seq)
        end
        shared = ExactPiModel(Dict(b => merge(br, (zy = counted,)) for (b, br) in branches))
        zy2 = (p, f, seq) -> zy(p, f, seq)
        separate = ExactPiModel(Dict(b => merge(br, (zy = b == 3 ? zy2 : zy,))
                                     for (b, br) in branches))
        @test shared.branches[2].zy === shared.branches[3].zy
        @test separate.branches[2].zy !== separate.branches[3].zy
        cs, ct = ConstraintEvaluator(hd, shared, box), ConstraintEvaluator(hd, separate, box)
        for kind in (:ihd, :thd, :rmsv, :rmsi)
            j = findfirst(c -> c.kind === kind, cs.cons)
            calls[] = 0
            vs, gs = fast_violation_gradient(cs, j, first(points), g; branch_derivative = :dual)
            vt, gt = fast_violation_gradient(ct, j, first(points), g; branch_derivative = :dual)
            @test calls[] == length(cs.cons[j].coef)
            @test vs === vt
            @test gs == gt
        end

        float_zy = (p::NamedTuple{(:kC_cable, :kL_cable, :kR_cable),Tuple{Float64,Float64,Float64}},
                    f, seq) -> zy(p, f, seq)
        plain = ExactPiModel(Dict(2 => merge(branches[2], (zy = float_zy,))))
        cp = ConstraintEvaluator(hd, plain, box)
        j = findfirst(c -> c.kind === :ihd, cp.cons)
        @test_throws r"^branch 2: the data do not accept dual numbers; use branch_derivative = :central\. Original error: MethodError:.*$"s fast_violation_gradient(
            cp, j, first(points), g; branch_derivative = :dual)
        @test all(isfinite, last(fast_violation_gradient(cp, j, first(points), g)))
    end

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
