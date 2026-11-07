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

        # Recompute the former dense storage independently, for every constraint kind.
        @test Set(first.(dict.keys)) == Set((:ihd, :thd, :rmsv, :rmsi))
        gv = [g[h] for h in dict.harmonics]
        for s in eachindex(dict.points)
            m = moments_at(hd, model, latent(box, dict.points[s]))
            cons = shhc_constraints(hd, m)
            μ = zeros(length(gv), length(cons))
            M = zeros(length(gv), length(gv), length(cons))
            for (j, c) in enumerate(cons), (a, h) in enumerate(dict.harmonics)
                haskey(c.coef, h) || continue
                μ[a, j] = c.coef[h][1]
                for (b, h2) in enumerate(dict.harmonics)
                    haskey(c.coef, h2) || continue
                    M[a, b, j] = sum(m.psi_norms[k+1] * c.coef[h][k+1] * c.coef[h2][k+1]
                                     for k in m.K if k != 0; init = 0.0)
                end
            end
            for (j, c) in enumerate(cons)
                μg = HPM.dot(view(μ, :, j), gv)
                F = c.λ <= 0.0 ? μg :
                    μg + c.λ * sqrt(max(HPM.dot(gv, view(M, :, :, j), gv), 0.0))
                dense = c.limit > 0.0 ? F / c.limit - 1.0 : F - c.limit
                @test V[j, s] === dense
            end
        end
        @test length(first(dict.μ)) == sum(length, dict.positions)
        @test sum(length, (dict.M[1].single, dict.M[1].full, dict.M[1].subset)) ==
              sum(p -> length(p)^2, dict.positions)
        changed = deepcopy(dict)
        pop!(changed.positions[1])
        @test_throws ErrorException add_point!(changed, hd, model, zeros(3), :refined)

        # Reading each constraint and filling a supplied matrix allocate no storage.
        read_allocations(dict, gv, j) =
            @allocated HPM._compact_violation(dict, j, dict.μ[1], dict.M[1], gv)
        scan_allocations(V, dict, gv) = @allocated HPM._dictionary_violations!(V, dict, gv)
        W = similar(V)
        HPM._dictionary_violations!(W, dict, gv)
        read_allocations(dict, gv, 1); scan_allocations(W, dict, gv)
        @test all(read_allocations(dict, gv, j) == 0 for j in eachindex(dict.keys))
        @test all(W[i] === V[i] for i in eachindex(V))
        @test scan_allocations(W, dict, gv) == 0

        # A proper harmonic subset keeps the dense dot path through local arrays.
        m = moments_at(hd, model, nominal_point(box))
        c = first(c for c in shhc_constraints(hd, m) if c.kind === :thd)
        H = dict.harmonics
        pos = [[1, 3]]
        subset = HPM.ShhcConstraint(c.kind, c.key,
            Dict(h => c.coef[h] for h in H[only(pos)]), c.limit, c.λ)
        blocks, sizes = HPM._dictionary_layout(H, pos)
        μ, M = HPM._compact([subset], m, H, pos, blocks, sizes)
        partial = HPM.ParameterDictionary(box, H, [(c.kind, c.key)], [c.limit], [c.λ],
            [zeros(length(box))], [:nominal], pos, blocks, sizes,
            [μ], [M])
        μfull = zeros(length(H), 1); Mfull = zeros(length(H), length(H), 1)
        μfull[only(pos), 1] = μ
        Mfull[only(pos), only(pos), 1] = reshape(M.subset, 2, 2)
        F = HPM.dot(view(μfull, :, 1), gv) +
            c.λ * sqrt(max(HPM.dot(gv, view(Mfull, :, :, 1), gv), 0.0))
        dense = c.limit > 0.0 ? F / c.limit - 1.0 : F - c.limit
        @test only(dictionary_violations(partial, g)) === dense

        # the interior worst case: for the 7th harmonic at the end of the cable feeder the
        # crossing seed beats every vertex
        j  = findfirst(==((:ihd, (7, 3))), dict.keys)
        sv = dict.labels .== :seed; vv = dict.labels .== :vertex
        @test maximum(V[j, sv]) > maximum(V[j, vv])
    end

    @testset "sign-predicted corners" begin
        dict = build_dictionary(hd, model, box; vertices = false, seeds, nbackground = 8)
        corners = dict.points[dict.labels .== :corner]
        @test !isempty(corners)
        @test length(corners) == length(unique(corners))
        @test !(:vertex in dict.labels)
        @test findlast(==(:corner), dict.labels) < findfirst(==(:seed), dict.labels)

        # Independently difference the full constraint records, including every harmonic.
        δ = 1e-4
        slopes = map(1:length(box)) do j
            e = zeros(length(box)); e[j] = δ
            left = shhc_constraints(hd, moments_at(hd, model, latent(box, -e)))
            right = shhc_constraints(hd, moments_at(hd, model, latent(box, e)))
            [(get(b.coef, h, [0.0])[1] - get(a.coef, h, [0.0])[1]) / (2δ)
             for h in dict.harmonics, (a, b) in zip(left, right)]
        end
        expected = unique([begin
            sensitivity = [s[h, c] for s in slopes]
            threshold = 1e-6 * maximum(abs, sensitivity)
            [abs(x) <= threshold ? 0.0 : sign(x) for x in sensitivity]
        end for c in eachindex(dict.keys) for h in eachindex(dict.harmonics)
            if dict.limit[c] > 0.0])
        @test Set(Tuple.(corners)) == Set(Tuple.(expected))
        vertices = [collect(v) for v in Iterators.product(ntuple(_ -> (-1.0, 1.0), 3)...)][:]
        matched = filter(v -> v in expected, vertices)
        @test !isempty(matched)
        @test all(v in corners for v in matched)

        disabled = build_dictionary(hd, model, box; vertices = false, corners = false,
                                    nbackground = 8)
        @test disabled.labels == [:nominal; fill(:background, 8)]

        # Reconstruct the dictionary before corners were introduced, in the original order.
        legacy = HPM.ParameterDictionary(box, dict.harmonics, dict.keys, dict.limit, dict.λ,
            Vector{Float64}[], Symbol[], dict.positions, dict.blocks, dict.sizes,
            Vector{Float64}[], HPM._SecondMoments[])
        add_point!(legacy, hd, model, zeros(3), :nominal)
        add_points!(legacy, hd, model, vertices, :vertex)
        add_points!(legacy, hd, model, [s.ξ for s in seeds], :seed)
        add_points!(legacy, hd, model, background_points(box, 8), :background)
        default = build_dictionary(hd, model, box; seeds, nbackground = 8)
        @test all(isequal(getfield(default, f), getfield(legacy, f))
                  for f in fieldnames(HPM.ParameterDictionary) if f !== :M)
        @test all(isequal(a.single, b.single) && isequal(a.full, b.full) &&
                  isequal(a.subset, b.subset) for (a, b) in zip(default.M, legacy.M))

        # At h5 the box stays below resonance. Two coordinates vary; nine fixed
        # coordinates exercise the d > 10 default without enlarging the physical box.
        small = ParameterBox(:kC_cable => (1.0, 0.9, 1.1), :kR_cable => (1.0, 0.8, 1.2))
        padded = ParameterBox([small.names; [Symbol("fixed_", j) for j in 1:9]],
            [small.nominal; ones(9)], [small.lower; ones(9)], [small.upper; ones(9)],
            [small.scale; fill(:linear, 9)])
        prepared = HPM.replicate(deepcopy(data); H = [1, 5])
        scaled = ScaledBranchModel(prepared, Dict(2 => :cable, 3 => :cable, 4 => :line))
        tol = 1e-3
        res = solve_rhhc(prepared, scaled, padded, solver_soc; pce_params = deepcopy(pp),
                         hpf_optimizer = solver_nlp, nbackground = 8, use_seeds = false, tol)
        @test res["signed_off"]
        @test :corner in res["dictionary"].labels
        g = Dict(h => I^2 for (h, I) in res["solution"])
        pos = res["dictionary"].limit .> 0.0
        worst = maximum(maximum(HPM._violations_at(prepared, scaled, small,
                            collect(v), g)[pos]) for v in Iterators.product((-1.0, 1.0), (-1.0, 1.0)))
        @test worst <= tol
    end

end
