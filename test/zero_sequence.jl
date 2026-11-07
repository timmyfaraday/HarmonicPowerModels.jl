################################################################################
# test/zero_sequence.jl                                                        #
# Zero-sequence network of the second-order cone program: the winding rules,  #
# the conversion of a harmonic network, missing data, floating parts, and the #
# stochastic and robust evaluations at a zero-sequence harmonic.              #
################################################################################

@testset "Zero Sequence" begin

    # ── A. winding rules ─────────────────────────────────────────────────────
    @testset "winding pairs" begin
        z, n1, n2 = 0.1 + 0.5im, 0.3 + 0.0im, 0.0 + 0.6im
        @test HPM._zs_combination(:earthed, :earthed, z, n1, n2) == (:series, z + n1 + n2)
        @test HPM._zs_combination(:earthed, :star, z, n1, n2)    == (:series, z + n1)
        @test HPM._zs_combination(:star, :earthed, z, n1, n2)    == (:series, z + n2)
        @test HPM._zs_combination(:earthed, :delta, z, n1, n2)   == (:shunt_fr, z + n1)
        @test HPM._zs_combination(:delta, :earthed, z, n1, n2)   == (:shunt_to, z + n2)
        @test HPM._zs_combination(:star, :delta, z, n1, n2)      == (:shunt_fr, z)
        @test HPM._zs_combination(:delta, :star, z, n1, n2)      == (:shunt_to, z)
        for (a, b) in ((:open, :earthed), (:earthed, :open), (:delta, :delta), (:open, :open),
                       (:delta, :open), (:star, :open))
            @test first(HPM._zs_combination(a, b, z, n1, n2)) === :open
        end
        @test HPM._zs_side('D', 1, false) === :delta
        @test HPM._zs_side('Y', 1, false) === :earthed
        @test HPM._zs_side('Y', 0, false) === :open
        @test HPM._zs_side('D', 0, true)  === :star
    end

    # A small harmonic network, already scaled: an earthed source at bus 1, a line 1–2,
    # a transformer branch 2–3 and a line 3–4.
    function network(vg; gnd1 = 1, gnd2 = 1, re2 = 0.0, star3 = false)
        bus(i) = Dict{String,Any}("bus_i" => i, "star_point" => (i == 3 && star3))
        line(i, f, t) = Dict{String,Any}("index" => i, "f_bus" => f, "t_bus" => t,
            "br_r" => 0.01, "br_x" => 0.1, "g_fr" => 0.0, "b_fr" => 0.02, "g_to" => 0.0,
            "b_to" => 0.02, "br_r0" => 0.03, "br_x0" => 0.3, "b0_fr" => 0.01, "b0_to" => 0.01,
            "tap" => 1.0, "shift" => 0.0)
        xf = Dict{String,Any}("index" => 2, "f_bus" => 2, "t_bus" => 3, "br_r" => 0.0,
            "br_x" => 0.2, "g_fr" => 0.0, "b_fr" => 0.0, "g_to" => 0.0, "b_to" => 0.0,
            "tap" => 1.0, "shift" => 0.0, "vg" => vg, "gnd1" => gnd1, "gnd2" => gnd2,
            "re1" => 0.0, "xe1" => 0.0, "re2" => re2, "xe2" => 0.0)
        return Dict{String,Any}("harmonic_sequence" => :zero,
            "bus" => Dict{String,Any}("$i" => bus(i) for i in 1:4),
            "branch" => Dict{String,Any}("1" => line(1, 1, 2), "2" => xf, "3" => line(3, 3, 4)),
            "gen" => Dict{String,Any}("1" => Dict{String,Any}("gen_bus" => 1, "gsc" => 1.0,
                "bsc" => -10.0, "gsc0" => 2.0, "bsc0" => -20.0)),
            "shunt" => Dict{String,Any}("1" => Dict{String,Any}("shunt_bus" => 4, "gs" => 0.0,
                "bs" => 0.5)),
            "load" => Dict{String,Any}())
    end
    idx = Dict(i => i for i in 1:4)
    yl  = 1 / (0.03 + 0.3im)                      # zero-sequence line series admittance

    # ── B. conversion of a harmonic network ─────────────────────────────────
    @testset "network view" begin
        ntw = network("Yy0"; re2 = 0.1)
        pos = copy(ntw); pos["harmonic_sequence"] = :positive
        @test HPM._zero_sequence_view(pos) === pos             # other sequences untouched
        v = HPM._zero_sequence_view(ntw)
        @test HPM._zero_sequence_view(v) === v                 # idempotent
        @test v["branch"]["1"]["br_x"] == 0.3 && v["branch"]["1"]["b_fr"] == 0.01
        @test v["gen"]["1"]["gsc"] == 2.0 && v["gen"]["1"]["bsc"] == -20.0
        @test isempty(v["shunt"])                              # an unearthed star bank
        @test ntw["branch"]["1"]["br_x"] == 0.1                # the input is not modified

        # YNyn0 with a neutral resistance on side 2: series z + 3·Rn
        Y = HPM._harmonic_admittance(ntw, idx)
        yt = 1 / (0.3 + 0.2im)
        @test Y[2, 3] ≈ -yt && Y[3, 2] ≈ -yt
        @test Y[1, 1] ≈ yl + 0.01im + (2.0 - 20.0im)
        @test Y[3, 3] ≈ yt + yl + 0.01im
        @test Y[4, 4] ≈ yl + 0.01im                            # no capacitor bank

        # YNd11: a shunt at the earthed side, the delta side is blocked
        Y = HPM._harmonic_admittance(network("Yd11"), idx)
        @test Y[2, 3] == 0 && Y[3, 2] == 0
        @test Y[2, 2] ≈ yl + 0.01im + 1 / 0.2im
        @test Y[3, 3] ≈ yl + 0.01im                            # bus 3 keeps its line

        # an unearthed star blocks the zero sequence altogether
        Y = HPM._harmonic_admittance(network("Yy0"; gnd1 = 0), idx)
        @test Y[2, 3] == 0 && Y[2, 2] ≈ yl + 0.01im && Y[3, 3] ≈ yl + 0.01im

        # a delta winding to a star point earths the star point (three-winding T model)
        Y = HPM._harmonic_admittance(network("Yd11"; gnd1 = 0, star3 = false), idx)
        @test Y[2, 3] == 0
        ntw = network("Dy11"; gnd2 = 0, star3 = true)
        Y = HPM._harmonic_admittance(ntw, idx)
        @test Y[2, 3] == 0 && Y[3, 3] ≈ yl + 0.01im + 1 / 0.2im && Y[2, 2] ≈ yl + 0.01im

        # an earthed bank is kept
        ntw = network("Yy0"); ntw["shunt"]["1"]["earthed"] = true
        @test HPM._harmonic_admittance(ntw, idx)[4, 4] ≈ yl + 0.01im + 0.5im
    end

    # ── C. missing zero-sequence data are an error, not a fallback ──────────
    @testset "missing data" begin
        ntw = network("Yy0"); delete!(ntw["branch"]["1"], "br_x0")
        @test_throws ErrorException HPM._zero_sequence_view(ntw)
        ntw["branch"]["1"]["zero_sequence_open"] = true
        @test HPM._zero_sequence_view(ntw)["branch"]["1"]["br_x"] == 0.0
        ntw = network("Yy0"); delete!(ntw["gen"]["1"], "gsc0")
        @test_throws ErrorException HPM._zero_sequence_view(ntw)
        ntw["gen"]["1"]["zero_sequence_open"] = true
        @test HPM._zero_sequence_view(ntw)["gen"]["1"]["gsc"] == 0.0
        ntw = network("Yy0"); delete!(ntw["branch"]["2"], "gnd2")
        @test_throws ErrorException HPM._zero_sequence_view(ntw)
        ntw = network("Yd11"); delete!(ntw["branch"]["2"], "gnd2")   # delta side needs none
        @test HPM._zero_sequence_view(ntw) isa Dict
        # a network model's own zero-sequence data are kept
        ntw = network("Yy0"); delete!(ntw["branch"]["1"], "br_x0")
        ntw["branch"]["1"]["sequence_data"] = :zero
        @test HPM._zero_sequence_view(ntw)["branch"]["1"]["br_x"] == 0.1
    end

    # ── D. parts of the network with no path to earth ───────────────────────
    @testset "floating parts" begin
        Y = ComplexF64[2 -1 0 0; -1 2 0 0; 0 0 1 -1; 0 0 -1 1]    # buses 3–4 float
        Y2, fl = HPM._earth_floating(Y)
        @test fl == [[3, 4]]
        @test Y2[1:2, 1:2] == Y[1:2, 1:2] && Y2[3, 3] == 2 && Y2[3, 4] == -1
        Z = inv(Y2)
        @test Z[1:2, 3:4] == zeros(2, 2)                          # no coupling to the rest
        @test isempty(last(HPM._earth_floating(Y[1:2, 1:2])))
    end

    # ── E. the stochastic and robust evaluations at a zero-sequence harmonic ──
    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    for br in values(data["branch"])        # zero-sequence data: three times the series impedance
        br["br_r0"] = 3 * br["br_r"]; br["br_x0"] = 3 * br["br_x"]
        br["b0_fr"] = 0.6 * br["b_fr"]; br["b0_to"] = 0.6 * br["b_to"]
    end
    pp = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli",
        "harmonics" => Dict(3 => (0.0, 30d2r), 5 => (4.75d2r, 40.05d2r)))

    @testset "zero-sequence data equal to positive give the positive network" begin
        same = deepcopy(data)
        for br in values(same["branch"])
            br["br_r0"] = br["br_r"]; br["br_x0"] = br["br_x"]
            br["b0_fr"] = br["b_fr"]; br["b0_to"] = br["b_to"]
        end
        hd  = HPM.replicate(same; H = [1, 3])
        nw  = hd["nw"]["3"]
        @test nw["harmonic_sequence"] === :zero
        pos = copy(nw); pos["harmonic_sequence"] = :positive
        bi  = Dict(b["bus_i"] => i for (i, b) in enumerate(sort(collect(values(nw["bus"])), by = b -> b["bus_i"])))
        @test HPM._harmonic_admittance(nw, bi) == HPM._harmonic_admittance(pos, bi)
    end

    @testset "sHHC_SOC with a zero-sequence harmonic" begin
        hd  = HPM.replicate(deepcopy(data); H = [1, 3, 5])
        sol = HPM.solve_shhc(hd, sHHC_SOC, solver_soc; pce_params = deepcopy(pp),
                             hpf_optimizer = solver_nlp)
        @test sol["termination_status"] == OPTIMAL
        @test sol["solution"][3] > 0.0 && sol["solution"][5] > 0.0

        # no zero-sequence path to earth at all: the units cannot inject the 3rd harmonic
        open = deepcopy(data)
        for br in values(open["branch"]); br["b0_fr"] = 0.0; br["b0_to"] = 0.0; end
        open["branch"]["1"]["zero_sequence_open"] = true        # the source branch
        delete!(open["branch"]["1"], "br_x0")
        hd2  = HPM.replicate(open; H = [1, 3, 5])
        sol2 = HPM.solve_shhc(hd2, sHHC_SOC, solver_soc; pce_params = deepcopy(pp),
                              hpf_optimizer = solver_nlp)
        @test sol2["termination_status"] == OPTIMAL
        @test abs(sol2["solution"][3]) <= 1e-8
        @test sol2["solution"][5] >= sol["solution"][5] - 1e-9    # the 3rd no longer shares the limits
    end

    @testset "robust evaluation at a zero-sequence harmonic" begin
        hd  = HPM.replicate(deepcopy(data); H = [1, 3, 5])
        HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, solver_nlp)
        HPM.update_hdata_with_pce_data!(hd, deepcopy(pp))
        cables = (2, 3)
        len = 1.0e3
        zy(b) = (p, f, seq) -> begin
            br = data["branch"]["$b"]; k = seq === :zero ? 3.0 : 1.0; ky = seq === :zero ? 0.6 : 1.0
            ((k * br["br_r"] * sqrt(f / 50) + im * k * br["br_x"] * f / 50) / len,
             im * ky * p.kC * (br["b_fr"] + br["b_to"]) * f / 50 / len)
        end
        model = ExactPiModel(Dict(b => (zy = zy(b), len = len, zbase = 1.0) for b in cables))
        box   = ParameterBox(:kC => (1.0, 0.8, 1.2))
        ce    = ConstraintEvaluator(hd, model, box)
        g     = Dict(3 => 1e-4, 5 => 1e-4)
        for ξ in ([0.0], [0.7], [-1.0])
            m    = moments_at(hd, model, latent(box, ξ))
            full = [normalized_violation(c, g, m) for c in shhc_constraints(hd, m)]
            fast = [fast_violation(ce, j, ξ, g) for j in eachindex(ce.cons)]
            @test maximum(abs.(fast .- full)) < 1e-10
        end
        # the adjoint gradient agrees with central differences at the zero-sequence harmonic
        j = findfirst(c -> c.kind === :thd, ce.cons)
        v, gr = fast_violation_gradient(ce, j, [0.3], g)
        δ = 1e-5
        fd = (fast_violation(ce, j, [0.3 + δ], g) - fast_violation(ce, j, [0.3 - δ], g)) / (2δ)
        @test isapprox(gr[1], fd; rtol = 1e-4, atol = 1e-8)
    end
end
