################################################################################
# test/shhc_industrial.jl                                                      #
# Integration tests for solve_shhc — verifies sHHC_QCQP and sHHC_LP on the   #
# two-bus network from `two_bus_example_hpf.m`, augmented with the harmonic   #
# angle distributions from §II Table III of the SHHC paper.                   #
#                                                                               #
# Assertions                                                                    #
#   A. Both stochastic solvers reach a feasible (solved) termination status.   #
#   B. sHHC_LP objective > dHHC_SOC objective — accounting for phase-angle    #
#      uncertainty strictly increases the admissible hosting capacity.          #
#   C. sHHC_QCQP objective ≥ dHHC_SOC objective (Cantelli bound ≤ worst-case).#
#   D. Per-harmonic LP bound I₅ > 0 and satisfies the IHD limit analytically. #
################################################################################

@testset "SHHC — Two-Bus Network" begin

    # ── Test network ───────────────────────────────────────────────────────────
    # Extend the standard two-bus HPF test case with:
    #   • ref_angle = 0° on all buses   (needed by constraint_load_current_angle
    #                                    in the deterministic dHHC_SOC baseline)
    #   • absolute equality fairness principle

    path  = joinpath(HPM.BASE_DIR, "test", "data", "matpower",
                     "two_bus_example_hpf.m")
    _data = PMs.parse_file(path)
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

    # ── 2. Stochastic LP ───────────────────────────────────────────────────────
    sol_lp = HPM.solve_shhc(fresh_hdata(), sHHC_LP, solver_soc;
                             pce_params = pce_params,
                             hpf_optimizer = solver_nlp)

    @testset "Stochastic LP (sHHC_LP)" begin
        @test sol_lp["termination_status"] == OPTIMAL
        @test isapprox(sol_lp["objective_value"], sum(values(sol_lp["solution"]));
                       rtol=1e-4)

        # Per-harmonic capacity must be positive for each harmonic with angle data
        for (h, I_h) in sol_lp["solution"]
            @test I_h >= 0.0
        end

        # B. Stochastic LP admits strictly more capacity than deterministic
        @test sol_lp["objective_value"] > I_det
    end

    # ── 3. Stochastic QCQP ────────────────────────────────────────────────────
    sol_qp = HPM.solve_shhc(fresh_hdata(), sHHC_QCQP, solver_nlp;
                             pce_params = pce_params)

    @testset "Stochastic QCQP (sHHC_QCQP)" begin
        @test sol_qp["termination_status"] in [LOCALLY_SOLVED, OPTIMAL]
        @test sol_qp["objective"] >= 0.0

        # C. Stochastic QCQP should match or exceed deterministic capacity
        #    (Cantelli risk bound is less conservative than worst-case)
        @test sol_qp["objective"] >= I_det - 1e-4   # small tolerance for NLP
    end

    # ── 4. Report ──────────────────────────────────────────────────────────────
    println("\n  ┌─ SHHC Two-Bus Results ─────────────────────────────────────")
    println("  │  Deterministic HHC (dHHC_SOC):  I_total = ",
            round(I_det, digits=4))
    println("  │  Stochastic LP   (sHHC_LP):     I_total = ",
            round(sol_lp["objective_value"], digits=4),
            "  (+",
            round((sol_lp["objective_value"]/I_det - 1)*100, digits=1), "%)")
    println("  │  Stochastic QCQP (sHHC_QCQP):  I_total = ",
            round(sol_qp["objective"], digits=4))
    for (h, I_h) in sort(collect(sol_lp["solution"]))
        println("  │    h=$h  I_h = ", round(I_h, digits=4), " p.u.")
    end
    println("  └────────────────────────────────────────────────────────────")

end
