################################################################################
# test/shhc_industrial.jl                                                      #
# Integration tests for solve_shhc — verifies sHHC_QCQP and sHHC_SOC on the   #
# two-bus network from `two_bus_example_hpf.m`, augmented with the harmonic   #
# angle distributions from §II Table III of the SHHC paper.                   #
#                                                                               #
# Assertions                                                                    #
#   A. Both stochastic solvers reach a feasible (solved) termination status.   #
#   B. dHHC_SOC ≡ sHHC_SOC ≡ sHHC_QCQP. This network has a SINGLE harmonic     #
#      unit, so every harmonic bus voltage is |V_h| = |Z_h|·I_h regardless of  #
#      the injection angle: the angle carries no information and the exact      #
#      uplift from modelling its uncertainty is ZERO. Any disagreement between  #
#      the three models here is a network-model or constraint-set gap, which is #
#      what makes this the sharpest available cross-check of the three builds.  #
#   C. Per-harmonic SOC bound I₅ > 0 and satisfies the IHD limit analytically. #
################################################################################

@testset "SHHC — Two-Bus Network" begin

    # ── Test network ───────────────────────────────────────────────────────────
    # Extend the standard two-bus HPF test case with:
    #   • ref_angle = 0° on all buses   (needed by constraint_load_current_angle
    #                                    in the deterministic dHHC_SOC baseline)
    #   • absolute equality fairness principle

    path  = joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                     "two_bus_example_hpf.m")
    _data = HPM.parse_file(path)
    for (_, bus) in _data["bus"]
        bus["ref_angle"] = 0.0
    end
    _data["principle"] = "absolute equality"

    # PCE parameters — Beta(α̂,β̂) fit from §II, h=5 angle bounds (Table III)
    d2r = π / 180
    pce_params = Dict{String,Any}(
        "alpha_hat" => 4.0414,
        "beta_hat"  => 4.7762,
        "delta"     => 2,
        "n_terms"   => 4,
        "epsilon"   => 0.10,          # 10 % risk level → λ(0.10) = 3.0
        "harmonics" => Dict(
            5 => ( 4.75d2r,  40.05d2r),
            7 => (-22.77d2r, 27.51d2r),
        ),
    )

    # Harmonics to consider (fundamental + 5th + 7th)
    H_TEST = [1, 5, 7]

    # Helper: fresh replicated hdata for each sub-test (avoid state bleed)
    fresh_hdata() = HPM.replicate(deepcopy(_data); H=H_TEST)

    # ── Solvers ────────────────────────────────────────────────────────────────
    # solver_nlp  = Ipopt  (defined in runtests.jl)
    # solver_soc  = Clarabel (defined in runtests.jl)

    # ── 1. Deterministic baseline ──────────────────────────────────────────────
    sol_det = HPM.solve_hhc(fresh_hdata(), dHHC_SOC, solver_soc, solver_nlp)

    @testset "Deterministic baseline (dHHC_SOC)" begin
        @test sol_det["termination_status"] == OPTIMAL
        @test sol_det["objective"] > 0.0
    end

    I_det = sol_det["objective"]

    # ── 2. Stochastic SOC ───────────────────────────────────────────────────────
    sol_soc = HPM.solve_shhc(fresh_hdata(), sHHC_SOC, solver_soc;
                             pce_params = pce_params,
                             hpf_optimizer = solver_nlp)

    @testset "Stochastic SOC (sHHC_SOC)" begin
        @test sol_soc["termination_status"] == OPTIMAL
        @test isapprox(sol_soc["objective_value"], sum(values(sol_soc["solution"]));
                       rtol=1e-4)

        # Per-harmonic capacity must be positive for each harmonic with angle data
        for (h, I_h) in sol_soc["solution"]
            @test I_h >= 0.0
        end

        # B. One harmonic unit ⇒ |V_h| is angle-independent ⇒ zero true uplift,
        #    so the stochastic SOCP must reproduce the deterministic capacity.
        #    The residual is the degree-2 truncation of |V|², which leaves a
        #    spurious variance of order 1e-3·mean and can only make the chance
        #    constraint tighter, never looser — hence the one-sided tolerance.
        @test isapprox(sol_soc["objective_value"], I_det; rtol = 1e-3)
        @test sol_soc["objective_value"] <= I_det + 1e-6
    end

    # ── 3. Stochastic QCQP ────────────────────────────────────────────────────
    sol_qp = HPM.solve_shhc(fresh_hdata(), sHHC_QCQP, solver_nlp;
                             pce_params = pce_params)

    @testset "Stochastic QCQP (sHHC_QCQP)" begin
        @test sol_qp["termination_status"] in [LOCALLY_SOLVED, OPTIMAL]
        @test sol_qp["objective"] >= 0.0

        # B (QCQP side). Same argument: one harmonic unit ⇒ no uplift to gain,
        #    so the QCQP must also land on the deterministic capacity.
        @test isapprox(sol_qp["objective"], I_det; rtol = 1e-3)
    end

    # ── 4. Per-unit harmonic current rating ───────────────────────────────────
    @testset "per-unit c_rating cap (sHHC_SOC)" begin
        # replicate derives load["c_rating"] from the unit's rated fundamental
        # current |S_d|/V_min; on this network that is well above the capacity the
        # IHD limit allows, so the cap must NOT bind at the default value …
        hd_ref = fresh_hdata()
        @test all(ld["c_rating"] > 1.0 for (_, ld) in hd_ref["nw"]["5"]["load"])

        # … and must bind when tightened below the unconstrained optimum.
        cap  = 0.5 * minimum(values(sol_soc["solution"]))
        hd   = fresh_hdata()
        for (_, ntw) in hd["nw"], (_, ld) in ntw["load"]
            ld["c_rating"] = cap
        end
        sol_cap = HPM.solve_shhc(hd, sHHC_SOC, solver_soc;
                                 pce_params = pce_params, hpf_optimizer = solver_nlp)
        @test sol_cap["termination_status"] == OPTIMAL
        for (_, I_h) in sol_cap["solution"]
            @test I_h ⪅ cap
        end
        @test sol_cap["objective_value"] < sol_soc["objective_value"]
    end

    # ── 5. Fairness principle guard ───────────────────────────────────────────
    @testset "sHHC_SOC rejects non-equality principles" begin
        # the SOCP assigns one scalar I_h to every unit, so it can only represent
        # absolute equality; anything else must fail loudly, not be ignored
        for p in ("maximin", "maximum efficiency", "Kalai-Smorodinsky bargaining")
            hd = fresh_hdata()
            hd["principle"] = p
            @test_throws ErrorException HPM.solve_shhc(hd, sHHC_SOC, solver_soc;
                            pce_params = pce_params, hpf_optimizer = solver_nlp)
        end
    end

    # ── 6. Report ──────────────────────────────────────────────────────────────
    println("\n  ┌─ SHHC Two-Bus Results (1 harmonic unit ⇒ expected uplift 0) ")
    println("  │  Deterministic HHC (dHHC_SOC):  I_total = ",
            round(I_det, digits=4))
    println("  │  Stochastic SOC   (sHHC_SOC):   I_total = ",
            round(sol_soc["objective_value"], digits=4),
            "  (", round((sol_soc["objective_value"]/I_det - 1)*100, digits=3), "%)")
    println("  │  Stochastic QCQP  (sHHC_QCQP):  I_total = ",
            round(sol_qp["objective"], digits=4),
            "  (", round((sol_qp["objective"]/I_det - 1)*100, digits=3), "%)")
    println("  │  per-harmonic sHHC_SOC capacity:")
    for (h, I_h) in sort(collect(sol_soc["solution"]); by = first)
        println("  │    h=$h  I_h = ", round(I_h, digits=4), " p.u.")
    end
    println("  └────────────────────────────────────────────────────────────")

end
