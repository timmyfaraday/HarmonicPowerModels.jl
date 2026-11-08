################################################################################
# test/branch_vector_group.jl                                                  #
# Data-layer invariants for transformer vector groups carried by branches.     #
################################################################################

@testset "Data layer — branch vector groups" begin
    case_path = joinpath(HPM.BASE_DIR, "test", "data", "matpower")
    H = collect(1:13)

    # ── A. vector-group phase angles ─────────────────────────────────────────
    @testset "Vector-group angle" begin
        for (vg, angles) in (("Yy0", (0.0, 0.0, 0.0)),
                             ("Yd11", (-330.0, 330.0, 0.0)),
                             ("Dy1", (-30.0, 30.0, 0.0)))
            for (sequence, angle) in zip((:positive, :negative, :zero), angles)
                @test HPM._vector_group_angle(vg, sequence) == angle
            end
        end
    end

    # ── B. branches without vg retain their shifts exactly ───────────────────
    @testset "No vg, no phase change" begin
        data = HPM.parse_file(joinpath(case_path, "30bus_network_hhc.m"))
        @test all(!haskey(br, "vg") for br in values(data["branch"]))
        hd = HPM.replicate(deepcopy(data); H = H)
        for h in H, (id, br) in hd["nw"]["$h"]["branch"]
            @test isequal(br["shift"], data["branch"][id]["shift"])
        end
    end

    # ── C. an xfmr and its unit-tap branch have the same phase stamp ──────────
    @testset "xfmr ≡ branch with vg" begin
        data = HPM.parse_file(joinpath(case_path, "industrial_network_hhc.m"))
        id = first(sort([id for (id, x) in data["xfmr"] if x["vg"] in ("Yd11", "Dy11")]))
        branch_data = deepcopy(data)
        xfmr = pop!(branch_data["xfmr"], id)
        bid = maximum(parse.(Int, collect(keys(branch_data["branch"])))) + 1
        branch_data["branch"]["$bid"] = Dict{String,Any}(
            "index" => bid, "f_bus" => xfmr["f_bus"], "t_bus" => xfmr["t_bus"],
            "br_r" => xfmr["r1"] + xfmr["r2"], "br_x" => xfmr["xsc"],
            "g_fr" => 0.0, "g_to" => 0.0, "b_fr" => 0.0, "b_to" => 0.0,
            "tap" => 1.0, "shift" => 0.0, "vg" => xfmr["vg"],
            "gnd1" => xfmr["gnd1"], "gnd2" => xfmr["gnd2"], "re1" => xfmr["re1"],
            "xe1" => xfmr["xe1"], "re2" => xfmr["re2"], "xe2" => xfmr["xe2"],
            "rate_a" => xfmr["rateA"] / data["baseMVA"], "br_status" => 1)
        bus_ids, bus_idx = HPM._bus_id_map(data)
        hd_xfmr = HPM.replicate(deepcopy(data); H = H)
        hd_branch = HPM.replicate(deepcopy(branch_data); H = H)
        for h in H
            br = hd_branch["nw"]["$h"]["branch"]["$bid"]
            @test br["tap"] == 1.0
            @test br["shift"] == deg2rad(HPM._vector_group_angle(xfmr["vg"], HPM.harmonic_sequence(h)))
            HPM.harmonic_sequence(h) === :zero && continue
            Y_xfmr = HPM._harmonic_admittance(hd_xfmr["nw"]["$h"], bus_idx)
            Y_branch = HPM._harmonic_admittance(hd_branch["nw"]["$h"], bus_idx)
            @test isapprox(Y_xfmr, Y_branch; atol = 1e-12, rtol = 0)
        end

        # The kernel builder uses lossless xfmr records; match that in the branch.
        for x in values(data["xfmr"])
            x["r1"] = x["r2"] = 0.0
        end
        for x in values(branch_data["xfmr"])
            x["r1"] = x["r2"] = 0.0
        end
        branch_data["branch"]["$bid"]["br_r"] = 0.0
        for h in H
            Y_xfmr = HPM.build_harmonic_Ybus(data, h)
            Y_branch = HPM.build_harmonic_Ybus(branch_data, h)
            @test isapprox(Y_xfmr, Y_branch; atol = 1e-12, rtol = 0)
            if HPM.harmonic_sequence(h) === :zero
                # a delta winding blocks the zero sequence between the two sides, but does
                # not disconnect the delta-side bus from the rest of its network
                f, t = bus_idx[xfmr["f_bus"]], bus_idx[xfmr["t_bus"]]
                @test iszero(Y_branch[f, t]) && iszero(Y_branch[t, f])
                delta_bus = uppercase(xfmr["vg"][1]) == 'D' ? f : t
                @test !all(iszero, Y_branch[delta_bus, :])
            end
        end
    end

    # ── D. radial phase shifts cancel when the terminal buses are eliminated ──
    @testset "Radial transformer Kron reduction" begin
        data = HPM.parse_file(joinpath(case_path, "30bus_network_hhc.m"))
        phased = deepcopy(data)
        for (id, endpoints) in (("13", (9, 11)), ("16", (12, 13)))
            br = phased["branch"][id]
            @test (br["f_bus"], br["t_bus"]) == endpoints
            br["vg"] = "Yd11"
        end
        hd = HPM.replicate(deepcopy(data); H = [1, 5, 7])
        hd_phased = HPM.replicate(phased; H = [1, 5, 7])
        bus_ids, bus_idx = HPM._bus_id_map(data)
        elim = [bus_idx[11], bus_idx[13]]
        keep = [bus_idx[id] for id in bus_ids if id ∉ (11, 13)]
        kron(Y) = Y[keep, keep] - Y[keep, elim] * (Y[elim, elim] \ Y[elim, keep])
        for h in (5, 7)
            Y = HPM._harmonic_admittance(hd["nw"]["$h"], bus_idx)
            Y_phased = HPM._harmonic_admittance(hd_phased["nw"]["$h"], bus_idx)
            @test !isapprox(Y, Y_phased; atol = 1e-10, rtol = 0)
            @test isapprox(kron(Y), kron(Y_phased); atol = 1e-10, rtol = 0)
            for id in ("13", "16")
                @test hd_phased["nw"]["$h"]["branch"][id]["tap"] == data["branch"][id]["tap"]
            end
        end
    end

    # ── E. rescaling the fundamental replaces its recorded vector-group phase ─
    @testset "Rescaling replaces the phase" begin
        data = HPM.parse_file(joinpath(case_path, "30bus_network_hhc.m"))
        data["branch"]["13"]["vg"] = "Yd11"
        shift = data["branch"]["13"]["shift"]
        hd = HPM.replicate(deepcopy(data); H = [1, 5, 7])
        bus_ids, bus_idx = HPM._bus_id_map(data)
        for h in (1, 5, 7)
            ntw = HPM.harmonic_network(hd, h; sequence = HPM.harmonic_sequence(h))
            ref = hd["nw"]["$h"]
            @test isapprox(ntw["branch"]["13"]["shift"], ref["branch"]["13"]["shift"];
                           atol = 1e-12, rtol = 0)
            Y = HPM._harmonic_admittance(ntw, bus_idx)
            Y_ref = HPM._harmonic_admittance(ref, bus_idx)
            @test maximum(abs, Y - Y_ref) <= 1e-12
        end
        ntw = HPM.harmonic_network(hd, 5; sequence = :positive)
        @test isapprox(ntw["branch"]["13"]["shift"], shift + deg2rad(-330);
                       atol = 1e-12, rtol = 0)
    end
end
