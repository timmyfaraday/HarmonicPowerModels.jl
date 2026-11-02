################################################################################
# test/rhhc_network_model.jl                                                   #
# Parametric network model of the robust HHC: parameter box, scaling between   #
# harmonics, network models, and the nominal point reproducing replicate.      #
################################################################################

@testset "RHHC — Network Model" begin

    path = joinpath(HPM.BASE_DIR, "test", "data", "matpower", "industrial_network_hhc.m")
    H    = [1, 5, 7, 11, 13]
    hd   = HPM.replicate(PMs.parse_file(path); H = H)
    brs  = sort(parse.(Int, collect(keys(hd["nw"]["1"]["branch"]))))

    @testset "ParameterBox" begin
        box = ParameterBox(:kC_cable => (1.0, 0.95, 1.10),
                           :rho_earth => (100.0, 10.0, 1000.0, :log))
        @test length(box) == 2
        @test latent(box, [-1.0, -1.0]).kC_cable == 0.95
        @test latent(box, [-1.0, -1.0]).rho_earth ≈ 10.0
        @test latent(box, [1.0, 1.0]).kC_cable == 1.10
        @test latent(box, [0.0, 0.0]).kC_cable ≈ 1.025
        @test latent(box, [1.0, 1.0]).rho_earth ≈ 1000.0
        @test latent(box, [0.0, 0.0]).rho_earth ≈ 100.0          # geometric midpoint
        @test nominal_point(box) == (kC_cable = 1.0, rho_earth = 100.0)
        @test_throws DomainError ParameterBox(:a => (1.0, 1.1, 2.0))
        @test_throws DomainError ParameterBox(:a => (1.0, 0.0, 2.0, :log))
        @test_throws ArgumentError ParameterBox(:a => (1.0, 0.0, 2.0, :cubic))
        @test_throws DimensionMismatch latent(box, [0.0])
    end

    @testset "harmonic_network equals replicate at integer h" begin
        for h in H[2:end]
            ntw = harmonic_network(hd, h; sequence = harmonic_sequence(h))
            ref = hd["nw"]["$h"]
            for b in keys(ref["branch"]), k in ("br_r", "br_x", "b_fr", "b_to")
                @test ntw["branch"][b][k] == ref["branch"][b][k]
            end
            for g in keys(ref["gen"]), k in ("gsc", "bsc")
                haskey(ref["gen"][g], k) && @test ntw["gen"][g][k] == ref["gen"][g][k]
            end
            for x in keys(ref["xfmr"]), k in ("xsc", "r1", "r2", "tr", "ti")
                @test ntw["xfmr"][x][k] == ref["xfmr"][x][k]
            end
        end
    end

    @testset "ScaledBranchModel" begin
        model = ScaledBranchModel(hd, Dict(b => :cable for b in brs); lengths = [brs[end]])
        @test uncertain_branches(model) == brs

        # nominal: every factor 1 ⇒ replicate's data exactly
        nws = networks_at(hd, model, NamedTuple())
        for h in H[2:end], b in brs, k in ("br_r", "br_x", "g_fr", "b_fr", "g_to", "b_to")
            @test nws["$h"]["branch"]["$b"][k] == get(hd["nw"]["$h"]["branch"]["$b"], k, 0.0)
        end
        # the input data are not modified
        @test nws["5"]["branch"] !== hd["nw"]["5"]["branch"]

        # factors act on their own quantity; the length factor on all three
        b, h = brs[2], 7
        ref  = hd["nw"]["$h"]["branch"]["$b"]
        d    = branch_data(model, b, (kR_cable = 1.2, kL_cable = 0.9, kC_cable = 1.1), h, :positive)
        @test d.br_r ≈ 1.2 * ref["br_r"]
        @test d.br_x ≈ 0.9 * ref["br_x"]
        @test d.b_fr ≈ 1.1 * ref["b_fr"]
        bl  = brs[end]
        dl  = branch_data(model, bl, (Symbol("len_", bl) => 2.0,) |> NamedTuple, h, :positive)
        rl  = hd["nw"]["$h"]["branch"]["$bl"]
        @test dl.br_r ≈ 2.0 * rl["br_r"] && dl.br_x ≈ 2.0 * rl["br_x"] && dl.b_fr ≈ 2.0 * rl["b_fr"]
        @test_throws ArgumentError branch_data(model, b, NamedTuple(), 3, :zero)
    end

    @testset "ExactPiModel tends to the lumped model for a short line" begin
        b, zbase, len = brs[2], 1.0, 1.0e3
        br = hd["nw"]["1"]["branch"]["$b"]
        # per-unit-length data reproducing the fundamental lumped values, with r·√h scaling
        zy(p, f, seq) = ((br["br_r"] * sqrt(f / 50) + im * br["br_x"] * f / 50) / len,
                         im * (br["b_fr"] + br["b_to"]) * f / 50 / len)
        model = ExactPiModel(Dict(b => (zy = zy, len = len, zbase = zbase)))
        d   = branch_data(model, b, NamedTuple(), 5, :positive)
        ref = hd["nw"]["5"]["branch"]["$b"]
        @test isapprox(d.br_x, ref["br_x"]; rtol = 1e-3)
        @test isapprox(d.br_r, ref["br_r"]; rtol = 1e-3)
        @test isapprox(d.b_fr, ref["b_fr"]; rtol = 1e-3)
        # an uncertain length as a function of p
        model2 = ExactPiModel(Dict(b => (zy = zy, len = p -> p.len, zbase = zbase)))
        @test branch_data(model2, b, (len = len,), 5, :positive) == d
    end

    @testset "admittance_at" begin
        model = ScaledBranchModel(hd, Dict(b => :cable for b in brs))
        for h in (5, 7)
            Y, idx = admittance_at(hd, model, NamedTuple(), h; sequence = harmonic_sequence(h))
            @test Y == HPM._harmonic_admittance(hd["nw"]["$h"], idx)
        end
        # continuous in h between harmonics
        Ya, _ = admittance_at(hd, model, NamedTuple(), 6.0)
        Yb, _ = admittance_at(hd, model, NamedTuple(), 6.0 + 1e-6)
        @test maximum(abs.(Ya .- Yb)) < 1e-3 * maximum(abs.(Ya))
    end

    @testset "nominal point reproduces sHHC_SOC" begin
        d2r  = π / 180
        data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "shhc", "data",
                                       "industrial_network_shhc.m"))
        data["principle"] = "absolute equality"
        pce_params = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762,
            "delta" => 2, "n_terms" => 4, "epsilon" => 0.10,
            "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r)))
        hd2 = HPM.replicate(deepcopy(data); H = [1, 5, 7])
        ref = HPM.solve_shhc(hd2, sHHC_SOC, solver_soc; pce_params = pce_params,
                             hpf_optimizer = solver_nlp)

        b2    = sort(parse.(Int, collect(keys(hd2["nw"]["1"]["branch"]))))
        model = ScaledBranchModel(hd2, Dict(b => :cable for b in b2))
        m_nom = shhc_moments(hd2)
        m_p   = shhc_moments(hd2; networks = networks_at(hd2, model, NamedTuple()))
        @test m_p.alpha_V == m_nom.alpha_V
        @test m_p.alpha_I == m_nom.alpha_I

        master = shhc_master(hd2, m_nom, solver_soc)
        add_shhc_constraints!(master, hd2, m_p)
        JuMP.optimize!(master.model)
        @test master.S * JuMP.objective_value(master.model) ≈ ref["objective_value"]

        # a larger cable capacitance moves the moments
        m_c = shhc_moments(hd2; networks = networks_at(hd2, model, (kC_cable = 1.1,)))
        @test m_c.alpha_V != m_nom.alpha_V
    end
end
