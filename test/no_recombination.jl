################################################################################
# test/no_recombination.jl                                                     #
# Hosting capacity without recombination: the unit contributions align in the #
# worst case at every bus and branch, so the squared quantities are the square #
# of the arithmetic sum of their magnitudes, with zero variance.               #
################################################################################

@testset "No Recombination" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    for br in values(data["branch"])        # zero-sequence data for the 3rd harmonic
        br["br_r0"] = 3 * br["br_r"]; br["br_x0"] = 3 * br["br_x"]
        br["b0_fr"] = 0.6 * br["b_fr"]; br["b0_to"] = 0.6 * br["b_to"]
    end
    # no angle data: without recombination the angles do not enter
    pp = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli", "recombination" => false)

    function prepared(H)
        hd = HPM.replicate(deepcopy(data); H)
        HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, solver_nlp)
        HPM.update_hdata_with_pce_data!(hd, deepcopy(pp))
        return hd
    end

    @testset "moments are the squared arithmetic sum" begin
        hd = prepared([1, 3, 5])
        @test hd["pce"]["recombination"] == false
        @test haskey(hd["pce"]["cos_phi"], 3) && haskey(hd["pce"]["cos_phi"], 5)
        m  = HPM.shhc_moments(hd)
        ids = sort([b["bus_i"] for b in values(hd["nw"]["1"]["bus"])])
        bi  = Dict(id => i for (i, id) in enumerate(ids))
        keep = [bi[id] for id in ids if !get(hd["nw"]["1"]["bus"]["$id"], "is_clean", false)]
        units = [bi[ld["load_bus"]] for ld in values(hd["nw"]["1"]["load"])]
        for h in (3, 5)
            Y = HPM._harmonic_admittance(hd["nw"]["$h"], bi)[keep, keep]
            Z = inv(Y); pos = Dict(k => j for (j, k) in enumerate(keep))
            for n in keep
                ref = sum(abs(Z[pos[n], pos[u]]) for u in units)^2
                @test isapprox(m.alpha_V[h][n, 1], ref; rtol = 1e-12)
                @test all(iszero, m.alpha_V[h][n, 2:end])          # no variance
            end
            # a branch current: the arithmetic sum of |y_s (Z_f − Z_t)|
            br = hd["nw"]["$h"]["branch"]["2"]
            v  = HPM._zero_sequence_view(hd["nw"]["$h"])["branch"]["2"]
            ys = 1 / (v["br_r"] + im * v["br_x"])
            f, t = pos[bi[br["f_bus"]]], pos[bi[br["t_bus"]]]
            ref = sum(abs(ys * (Z[f, pos[u]] - Z[t, pos[u]])) for u in units)^2
            @test isapprox(m.alpha_I[h][2][1], ref; rtol = 1e-12)
        end
    end

    @testset "worst-case alignment bounds any fixed angles" begin
        hd  = HPM.replicate(deepcopy(data); H = [1, 3, 5, 7])
        sol = HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = deepcopy(pp),
                             hpf_optimizer = solver_nlp)
        @test sol["termination_status"] == OPTIMAL
        # fixed (zero-width) angles at their mean: contributions may partly cancel
        fixed = Dict{String,Any}(pp); delete!(fixed, "recombination")
        fixed["harmonics"] = Dict(3 => (20d2r, 20d2r), 5 => (30d2r, 30d2r), 7 => (-10d2r, -10d2r))
        hd2  = HPM.replicate(deepcopy(data); H = [1, 3, 5, 7])
        sol2 = HPM.solve_shhc(hd2, sHHC_SOC, solver_soc; pce_params = fixed,
                              hpf_optimizer = solver_nlp)
        @test sol["objective_value"] <= sol2["objective_value"] + 1e-9
    end

    @testset "robust evaluation without recombination" begin
        hd  = prepared([1, 3, 5])
        cables = (2, 3); len = 1.0e3
        zy(b) = (p, f, seq) -> begin
            br = data["branch"]["$b"]; k = seq === :zero ? 3.0 : 1.0; ky = seq === :zero ? 0.6 : 1.0
            ((k * br["br_r"] * sqrt(f / 50) + im * k * br["br_x"] * f / 50) / len,
             im * ky * p.kC * (br["b_fr"] + br["b_to"]) * f / 50 / len)
        end
        model = ExactPiModel(Dict(b => (zy = zy(b), len = len, zbase = 1.0) for b in cables))
        box   = ParameterBox(:kC => (1.0, 0.7, 1.3))
        ce    = ConstraintEvaluator(hd, model, box)
        @test !ce.recomb
        g     = Dict(3 => 2e-4, 5 => 1e-4)
        for ξ in ([0.0], [0.6], [-1.0])
            m    = moments_at(hd, model, latent(box, ξ))
            full = [normalized_violation(c, g, m) for c in shhc_constraints(hd, m)]
            fast = [fast_violation(ce, j, ξ, g) for j in eachindex(ce.cons)]
            @test maximum(abs.(fast .- full)) < 1e-10
        end
        δ = 1e-5
        for kind in (:thd, :ihd, :rmsi)
            for j in findall(c -> c.kind === kind, ce.cons)
                v, gr = fast_violation_gradient(ce, j, [0.25], g)
                @test isapprox(v, fast_violation(ce, j, [0.25], g); atol = 1e-12)
                fd = (fast_violation(ce, j, [0.25 + δ], g) - fast_violation(ce, j, [0.25 - δ], g)) / (2δ)
                @test isapprox(gr[1], fd; rtol = 1e-4, atol = 1e-8)
            end
        end

        # the cutting-set method end to end: bounds, sign-off and the cost of robustness
        hd3 = HPM.replicate(deepcopy(data); H = [1, 3, 5])
        res = solve_rhhc(hd3, model, box, solver_soc; pce_params = deepcopy(pp),
                         hpf_optimizer = solver_nlp, nbackground = 32, posterior = 64)
        @test res["signed_off"]
        @test res["lower_bound"] <= res["upper_bound"] + 1e-12 && res["gap"] < 1e-3
        @test res["posterior"].max_violation <= 1e-6
        hd0 = HPM.replicate(deepcopy(data); H = [1, 3, 5])
        nom = solve_rhhc(hd0, model, ParameterBox(:kC => (1.0, 1.0, 1.0)), solver_soc;
                         pce_params = deepcopy(pp), hpf_optimizer = solver_nlp, nbackground = 0,
                         use_seeds = false)
        @test res["upper_bound"] <= nom["upper_bound"] + 1e-9
    end
end
