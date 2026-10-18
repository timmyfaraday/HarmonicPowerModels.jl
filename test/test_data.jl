################################################################################
# test/test_data.jl                                                            #
# Data-layer invariants of `replicate` and of the harmonic admittance builders. #
#                                                                               #
# These pin model choices that no objective value protects: the industrial HOPF #
# network is the only test case carrying bus shunts, and its harmonic voltages  #
# sit at ~1e-5 p.u., so a wrong shunt model changes nothing a solution-based    #
# assertion would catch. The invariants asserted here are exact instead.        #
#                                                                               #
#   A. `replicate` scales a bus shunt's susceptance with frequency, by sign.    #
#   B. `_harmonic_admittance` (reduced SOCP network) stamps bus shunts and      #
#      branch conductance.                                                     #
#   C. `build_harmonic_Ybus` (kernel path) stamps them with the same scaling.   #
#   D. dHHC_SOC ≡ sHHC_SOC on a shunt-bearing network — the equivalence that    #
#      the stamps in B exist to preserve.                                       #
#   E. `replicate` derives a load's harmonic current rating from its rated      #
#      fundamental current, and an explicit rating in the data wins.            #
################################################################################

@testset "Data layer — harmonic scaling and shunt modelling" begin

    _shunt_path = joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                           "industrial_network_hopf.m")

    # ── A. bus shunt susceptance scales with frequency ────────────────────────
    @testset "A. replicate scales bus shunt susceptance" begin
        data = HPM.parse_file(_shunt_path)
        @test !isempty(data["shunt"])

        bs1 = Dict(i => s["bs"] for (i, s) in data["shunt"])
        gs1 = Dict(i => s["gs"] for (i, s) in data["shunt"])
        @test any(!iszero, values(bs1))          # the case really has shunts

        H = [1, 3, 5, 13]
        hdata = HPM.replicate(deepcopy(data); H = H)

        for h in H, (i, s) in hdata["nw"]["$h"]["shunt"]
            # capacitive (bs ≥ 0) ⇒ B(h) = bs·h ;  gs is frequency-independent
            @test s["bs"] ≈ bs1[i] * h
            @test s["gs"] ≈ gs1[i]
        end

        # inductive shunts scale the other way: B(h) = bs/h
        ind = HPM.parse_file(_shunt_path)
        for (_, s) in ind["shunt"]
            s["bs"] = -0.10
        end
        hdata_i = HPM.replicate(deepcopy(ind); H = H)
        for h in H, (_, s) in hdata_i["nw"]["$h"]["shunt"]
            @test s["bs"] ≈ -0.10 / h
        end
    end

    # ── B. reduced admittance matrix stamps shunts and branch conductance ─────
    @testset "B. _harmonic_admittance stamps shunt and g_fr/g_to" begin
        data  = HPM.parse_file(_shunt_path)
        hdata = HPM.replicate(deepcopy(data); H = [1, 5])
        ntw   = hdata["nw"]["5"]

        bus_ids = sort([b["bus_i"] for (_, b) in ntw["bus"]])
        bus_idx = Dict(id => k for (k, id) in enumerate(bus_ids))

        Y_ref = HPM._harmonic_admittance(ntw, bus_idx)

        # dropping every shunt must change exactly the shunt buses' diagonals,
        # each by the shunt admittance gs + j·bs
        stripped = deepcopy(ntw)
        empty!(stripped["shunt"])
        Y_nosh = HPM._harmonic_admittance(stripped, bus_idx)

        touched = Set{Int}()
        for (_, s) in ntw["shunt"]
            k = bus_idx[s["shunt_bus"]]
            push!(touched, k)
            @test Y_ref[k, k] - Y_nosh[k, k] ≈ s["gs"] + im * s["bs"]
        end
        @test !isempty(touched)
        for k in 1:length(bus_ids), l in 1:length(bus_ids)
            (k == l && k in touched) && continue
            @test Y_ref[k, l] ≈ Y_nosh[k, l]
        end

        # branch series conductance enters the from/to diagonals (behind the tap
        # on the from side), mirroring the QCQP's PCE current balance
        withg = deepcopy(ntw)
        for (_, br) in withg["branch"]
            br["g_fr"] = 0.05
            br["g_to"] = 0.07
        end
        Y_g = HPM._harmonic_admittance(withg, bus_idx)
        dg  = zeros(ComplexF64, length(bus_ids), length(bus_ids))
        for (_, br) in ntw["branch"]
            tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
            dg[bus_idx[br["f_bus"]], bus_idx[br["f_bus"]]] += 0.05 / tm^2
            dg[bus_idx[br["t_bus"]], bus_idx[br["t_bus"]]] += 0.07
        end
        @test Y_g ≈ Y_ref .+ dg
    end

    # ── C. kernel Ybus applies the same shunt model ───────────────────────────
    @testset "C. build_harmonic_Ybus stamps shunts consistently" begin
        data = HPM.parse_file(_shunt_path)
        for h in (5, 13)
            Y = HPM.build_harmonic_Ybus(data, h)
            stripped = deepcopy(data); empty!(stripped["shunt"])
            Y0 = HPM.build_harmonic_Ybus(stripped, h)
            bus_ids = sort([b["bus_i"] for (_, b) in data["bus"]])
            bus_idx = Dict(id => k for (k, id) in enumerate(bus_ids))
            for (_, s) in data["shunt"]
                k = bus_idx[s["shunt_bus"]]
                # same scaling rule as replicate: capacitive ⇒ ·h
                @test Y[k, k] - Y0[k, k] ≈ s["gs"] + im * s["bs"] * h
            end
        end
    end

    # ── D. dHHC ≡ sHHC_SOC with a shunt present ───────────────────────────────
    @testset "D. dHHC_SOC ≡ sHHC_SOC on a shunt-bearing network" begin
        base = HPM.parse_file(joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                                       "two_bus_example_hpf.m"))
        for (_, b) in base["bus"]
            b["ref_angle"] = 0.0
        end
        base["principle"] = "absolute equality"

        # point mass at E[φ] ⇒ zero variance ⇒ the stochastic model must
        # reproduce the deterministic one exactly, which isolates the network
        # model as the only possible source of disagreement
        α̂, β̂ = 4.0414, 4.7762
        Eω   = α̂ / (α̂ + β̂)
        d2r  = π / 180
        θ    = Dict(5 => 4.75d2r + (40.05d2r - 4.75d2r) * Eω,
                    7 => -22.77d2r + (27.51d2r + 22.77d2r) * Eω)
        pce_params = Dict{String,Any}(
            "alpha_hat" => α̂, "beta_hat" => β̂, "delta" => 2, "epsilon" => 0.10,
            "harmonics" => Dict(h => (v, v) for (h, v) in θ))

        # A capacitor bank at bus 2, added BEFORE replicate so the susceptance
        # goes through the frequency scaling of test A: this exercises scaling
        # and stamp together, end to end. The bank is deliberately large — a
        # small one shifts the answer by only ~0.2 %, too little for the "shunt
        # actually matters" guard below to tell a correct stamp from a missing one.
        base["shunt"] = Dict("1" => Dict("shunt_bus" => 2, "gs" => 0.0,
                                         "bs" => 0.30, "status" => 1,
                                         "index" => 1,
                                         "source_id" => Any["bus", 2]))
        build() = HPM.replicate(deepcopy(base); H = [1, 5, 7])

        hd_det = build()
        for (nw, ntw) in hd_det["nw"]
            nh = parse(Int, nw); nh == 1 && continue
            for (_, ld) in ntw["load"]
                ld["ref_angle"] = rad2deg(θ[nh])
            end
        end
        sol_det = HPM.solve_hhc(hd_det, dHHC_SOC, solver_soc, solver_nlp)
        sol_soc = HPM.solve_shhc(build(), sHHC_SOC, solver_soc;
                                 pce_params = pce_params,
                                 hpf_optimizer = solver_nlp)

        @test sol_det["termination_status"] == OPTIMAL
        @test sol_soc["termination_status"] == OPTIMAL
        @test isapprox(sol_soc["objective_value"], sol_det["objective"]; rtol = 1e-6)

        # and the shunt must actually matter, or the test proves nothing
        nosh = deepcopy(base); empty!(nosh["shunt"])
        sol_nosh = HPM.solve_shhc(HPM.replicate(nosh; H = [1, 5, 7]), sHHC_SOC,
                                  solver_soc; pce_params = pce_params,
                                  hpf_optimizer = solver_nlp)
        @test !isapprox(sol_nosh["objective_value"], sol_soc["objective_value"];
                        rtol = 1e-2)
    end

    # ── E. load harmonic current rating ───────────────────────────────────────
    @testset "E. load c_rating derived from rated current" begin
        data  = HPM.parse_file(joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                                        "two_bus_example_hpf.m"))
        pd, qd = data["load"]["1"]["pd"], data["load"]["1"]["qd"]
        vmin   = data["bus"]["2"]["vmin"]

        hdata = HPM.replicate(deepcopy(data); H = [1, 5, 7])
        for nw in ("1", "5", "7")
            # rated fundamental current |S_d| / V_min, NOT the harmonic-scaled power
            @test hdata["nw"][nw]["load"]["1"]["c_rating"] ≈ sqrt(pd^2 + qd^2) / vmin
        end

        # an explicit rating in the input data wins
        override = HPM.parse_file(joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                                           "two_bus_example_hpf.m"))
        override["load"]["1"]["c_rating"] = 0.123
        hd_o = HPM.replicate(deepcopy(override); H = [1, 5, 7])
        for nw in ("1", "5", "7")
            @test hd_o["nw"][nw]["load"]["1"]["c_rating"] == 0.123
        end
    end

end
