################################################################################
# test/rhhc_minimal.jl                                                         #
# Minimal L–C example (examples/rhhc/data/lc_minimal.m): one inductive branch  #
# from an ideal source to a single unit, its charging as the capacitance; the  #
# resonance sits at h = 6 and the box (L, C ±20 %) moves it over [5, 7.5].    #
# With one unit the angle is irrelevant, so the robust capacity has two        #
# independent references: the analytic IHD bound √limit / max|Z| per harmonic #
# and the sampled program on a dense grid of the box.                          #
################################################################################

@testset "RHHC — Minimal L–C Example" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "lc_minimal.m"))
    H    = [1, 5, 7, 11, 13]
    pp   = Dict{String,Any}("alpha_hat" => 4.0414, "beta_hat" => 4.7762, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "cantelli",
        "harmonics" => Dict(5 => (4.75d2r, 40.05d2r), 7 => (-22.77d2r, 27.51d2r),
                            11 => (79.69d2r, 164.01d2r), 13 => (49.38d2r, 147.67d2r)))
    hd    = HPM.replicate(deepcopy(data); H)
    model = ScaledBranchModel(hd, Dict(1 => :line))
    box   = ParameterBox(:kL_line => (1.0, 0.8, 1.2), :kC_line => (1.0, 0.8, 1.2))
    res   = solve_rhhc(hd, model, box, solver_soc; pce_params = deepcopy(pp),
                       hpf_optimizer = solver_nlp, nbackground = 64, posterior = 256)

    @test isapprox(only(modal_scan(hd, model, nominal_point(box))).h, 6.0; atol = 0.01)
    @test res["signed_off"]
    @test res["gap"] < 1e-3
    @test res["posterior"].max_violation ≤ 1e-6

    # reference 1: per harmonic, the IHD limit alone allows √limit / max_box |Z_11(h)|
    br = hd["nw"]["1"]["branch"]["1"]
    r, x, bt = br["br_r"], br["br_x"], br["b_to"]
    k  = range(0.8, 1.2; length = 201)
    m0 = HPM.shhc_moments(hd)
    for h in H[2:end]
        zmax = maximum(abs(1 / (1 / (r * sqrt(h) + im * x * h * kL) + im * bt * h * kC))
                       for kL in k, kC in k)
        c = only(c for c in shhc_constraints(hd, m0) if c.kind === :ihd && c.key == (h, 1))
        Iref = sqrt(c.limit) / zmax
        @test res["solution"][h] ≤ Iref * (1 + 1e-4)          # never above the IHD bound
        @test res["solution"][h] ≥ Iref * 0.98                # and close to it
    end

    # reference 2: the sampled program on a dense grid (a relaxation of the robust one)
    master = shhc_master(hd, m0, solver_soc)
    for a in range(-1, 1; length = 21), b in range(-1, 1; length = 21)
        add_shhc_constraints!(master, hd, moments_at(hd, model, latent(box, [a, b])))
    end
    JuMP.optimize!(master.model)
    pgrid = master.S * JuMP.objective_value(master.model)
    @test res["lower_bound"] ≤ pgrid * (1 + 1e-6)
    @test (pgrid - res["lower_bound"]) / pgrid < 1e-3
end
