################################################################################
# test/shhc_independent.jl                                                     #
# Integration tests for the independent-germ stochastic HHC model on a         #
# multi-unit (N=4) network: one Beta germ per harmonic unit, harmonic-unit     #
# designation, and per-(unit,harmonic) angle ranges.                           #
#                                                                               #
# Covers the data layer (update_hdata_with_pce_data!) and both solvers         #
# (sHHC_SOC, sHHC_QCQP) — the two-bus integration test only has N=1, where the  #
# multivariate basis collapses to the univariate one and these paths are not   #
# exercised.                                                                    #
################################################################################

@testset "SHHC — Independent germs (N=4)" begin

    path  = joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                     "industrial_network_hhc.m")
    _data = HPM.parse_file(path)
    for (_, bus) in _data["bus"]
        bus["ref_angle"] = 0.0
    end
    _data["principle"] = "absolute equality"

    d2r = π / 180
    # value type Any so a harmonic entry can be a bare tuple OR a per-unit Dict
    base_harm = Dict{Int,Any}(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r))
    base_params() = Dict{String,Any}(
        "alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.10,
        "harmonics" => deepcopy(base_harm),
    )
    H_TEST = [1, 5, 7]
    fresh() = HPM.replicate(deepcopy(_data); H=H_TEST)

    load_ids = sort(parse.(Int, collect(keys(_data["load"]))))
    @test length(load_ids) == 4          # N=4 harmonic units available

    # ── 1. Data layer — default (all loads harmonic) ─────────────────────────
    @testset "data layer: default multivariate basis" begin
        hd = fresh(); HPM.update_hdata_with_pce_data!(hd, base_params())
        pce = hd["pce"]; pb = pce["basis"]
        @test length(pce["harmonic_units"]) == 4
        @test pb.dim == binomial(4 + 2, 2)            # P = 15
        @test pce["K"] == 0:(pb.dim - 1)
        @test length(pce["K_phys"]) == 1 + 4 * 2       # 1 + Nδ = 9 single-germ
        @test isapprox(pce["psi_norms"][1], 1.0; atol=1e-12)

        # each unit's coefficients are supported only on its own germ dimension
        ind = pb.ind
        for (uid, gd) in pce["germ_dim"]
            cc = pce["cos_phi"][5][uid]
            nz = [k for k in 1:length(cc) if abs(cc[k]) > 1e-12]
            @test all(all(ind[k, u] == 0 for u in 1:size(ind, 2) if u != gd)
                      for k in nz)
        end
    end

    # ── 2. Data layer — harmonic-unit subset + per-unit ranges ───────────────
    @testset "data layer: subset + per-unit ranges" begin
        p = base_params()
        p["harmonic_units"] = [load_ids[1], load_ids[2]]
        p["harmonics"][5]   = Dict(load_ids[1] => (0.0, 30.0),
                                   load_ids[2] => (10.0, 40.0))
        hd = fresh(); HPM.update_hdata_with_pce_data!(hd, p)
        pce = hd["pce"]
        @test length(pce["harmonic_units"]) == 2
        @test pce["basis"].dim == binomial(2 + 2, 2)   # P = 6
        # normal loads carry zero injection
        @test all(iszero, pce["cos_phi"][5][load_ids[3]])
        @test all(iszero, pce["sin_phi"][5][load_ids[4]])
        # per-unit ranges → distinct coefficients for the two harmonic units
        @test pce["cos_phi"][5][load_ids[1]] != pce["cos_phi"][5][load_ids[2]]
    end

    # ── 2b. Per-unit mean angles centred on each unit's ref_angle ─────────────
    @testset "data layer: angle_centre = ref_angle" begin
        α̂, β̂ = 4.0414, 4.7762
        Eω   = α̂ / (α̂ + β̂)
        mean_of(r) = r[1] + (r[2] - r[1]) * Eω

        # Re-parse the case so the file's OWN reference angles survive: the outer
        # testset zeroes them for the deterministic baselines, but this lever is
        # about exactly that data — 0° at the 155 kV buses, −30° one transformer
        # down, −60° two down, the vector-group cascade of the industrial network.
        dat = HPM.parse_file(path)
        dat["principle"] = "absolute equality"
        @test length(unique(bus["ref_angle"] for (_, bus) in dat["bus"])) == 3
        rep() = HPM.replicate(deepcopy(dat); H = H_TEST)

        # shared mean: every unit sits on the Beta mean of the supplied range
        hd0 = rep(); HPM.update_hdata_with_pce_data!(hd0, base_params())
        m0  = [mean_of(hd0["pce"]["angle_ranges"][5][d]) for d in load_ids]
        @test all(isapprox(m, m0[1]; atol = 1e-12) for m in m0)

        # per-unit mean: each unit sits on its own ref_angle, sequence sign applied
        p = base_params(); p["angle_centre"] = "ref_angle"
        hd = rep(); HPM.update_hdata_with_pce_data!(hd, p)
        for d in load_ids
            bus_d = hd["nw"]["5"]["load"]["$d"]["load_bus"]
            want  = deg2rad(hd["nw"]["5"]["load"]["$d"]["ref_angle"])
            @test mean_of(hd["pce"]["angle_ranges"][5][d]) ≈ want
            # h=5 is negative sequence ⇒ the sign of ref_angle is flipped
            @test want ≈ -deg2rad(dat["bus"]["$bus_d"]["ref_angle"])
        end
        # the width of each support is preserved, only its centre moved
        width(r) = r[2] - r[1]
        w0 = [width(hd0["pce"]["angle_ranges"][5][d]) for d in load_ids]
        w1 = [width(hd["pce"]["angle_ranges"][5][d])  for d in load_ids]
        @test all(isapprox(w1[i], w0[i]; atol = 1e-12) for i in eachindex(w0))
        # and units on different voltage levels now carry different mean angles
        m1 = [mean_of(hd["pce"]["angle_ranges"][5][d]) for d in load_ids]
        @test length(unique(round.(m1; digits = 9))) > 1

        # a degenerate support collapses exactly onto ref_angle, which is what
        # makes the point-mass reduction match the deterministic model per unit
        pp = base_params(); pp["angle_centre"] = "ref_angle"
        pp["harmonics"][5] = (0.3, 0.3)
        hdp = rep(); HPM.update_hdata_with_pce_data!(hdp, pp)
        for d in load_ids
            lo, hi = hdp["pce"]["angle_ranges"][5][d]
            @test lo == hi
            @test lo ≈ deg2rad(hdp["nw"]["5"]["load"]["$d"]["ref_angle"])
        end

        # a load with no ref_angle in the data is centred on zero, not skipped
        bare = HPM.parse_file(path)
        for (_, bus) in bare["bus"]; delete!(bus, "ref_angle"); end
        pb2 = base_params(); pb2["angle_centre"] = "ref_angle"
        hdb = HPM.replicate(deepcopy(bare); H = H_TEST)
        HPM.update_hdata_with_pce_data!(hdb, pb2)
        @test all(abs(mean_of(hdb["pce"]["angle_ranges"][5][d])) < 1e-12 for d in load_ids)

        # an unknown mode is rejected rather than silently ignored
        pbad = base_params(); pbad["angle_centre"] = "nonsense"
        @test_throws ErrorException HPM.update_hdata_with_pce_data!(rep(), pbad)
    end

    # ── 3. sHHC_SOC — independent germs, all units host equally ────────────────
    @testset "sHHC_SOC: independent, equal hosting" begin
        sol = HPM.solve_shhc(fresh(), sHHC_SOC, solver_soc;
                             pce_params = base_params(), hpf_optimizer = solver_nlp)
        # Clarabel sits on its tolerance here (Mosek: OPTIMAL, same I_h to 1e-13)
        @test sol["termination_status"] in (OPTIMAL, ALMOST_OPTIMAL)
        for (_, I_h) in sol["solution"]
            @test I_h > 0.0
        end
    end

    # ── 4. sHHC_QCQP — subset designation zeroes normal loads ────────────────
    @testset "sHHC_QCQP: subset excludes normal loads" begin
        p = base_params()
        p["harmonic_units"] = [load_ids[1], load_ids[2]]
        sol = HPM.solve_shhc(fresh(), sHHC_QCQP, solver_nlp; pce_params = p)
        @test sol["termination_status"] in [LOCALLY_SOLVED, OPTIMAL]

        cmd(h, d) = sol["solution"]["nw"]["$h"]["load"]["$d"]["cmd"]
        for h in (5, 7)
            # harmonic units host equally and positively
            @test cmd(h, load_ids[1]) > 0.0
            @test isapprox(cmd(h, load_ids[1]), cmd(h, load_ids[2]); rtol=1e-5)
            # normal loads host nothing
            @test cmd(h, load_ids[3]) ≈ 0.0
            @test cmd(h, load_ids[4]) ≈ 0.0
        end
    end

    # ── 4b. sHHC_QCQP ≡ sHHC_SOC with an off-nominal branch tap ──────────────
    # The QCQP's KCL must orient each arc by the branch's own f_bus (the tap
    # scaling and the from-side shunt apply at the from-end only), and the SOCP's
    # branch current is the series current y_s·(V_f/t − V_t). Without taps both
    # errors vanish, so this case needs a tap ≠ 1; on branch 4 (buses 6–7) the
    # current limit binds, so both are exercised. Totals are compared: the THD
    # constraint couples the harmonics and leaves their split nearly degenerate.
    @testset "sHHC_QCQP ≡ sHHC_SOC: off-nominal tap" begin
        tapped() = begin
            d = deepcopy(_data)
            d["branch"]["4"]["tap"] = 0.9
            HPM.replicate(d; H=H_TEST)
        end
        soc  = HPM.solve_shhc(tapped(), sHHC_SOC, solver_soc;
                              pce_params = base_params(), hpf_optimizer = solver_nlp)
        qcqp = HPM.solve_shhc(tapped(), sHHC_QCQP, solver_nlp; pce_params = base_params())
        @test qcqp["termination_status"] in [LOCALLY_SOLVED, OPTIMAL]
        @test isapprox(sum(qcqp["solution"]["nw"]["$h"]["load"]["$(load_ids[1])"]["cmd"] for h in (5, 7)),
                       sum(soc["solution"][h] for h in (5, 7)); rtol=1e-3)
    end

end
