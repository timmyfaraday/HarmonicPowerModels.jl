################################################################################
# test/rhhc_exact_moments.jl                                                   #
# Exact moments (`"moments" => "exact"`) in sHHC_SOC and the robust program:  #
# the network moments against Monte Carlo, the fast evaluation against the    #
# full one, and the cutting-set method end to end.                            #
################################################################################

using Random, Statistics, Distributions

@testset "RHHC — Exact Moments" begin

    d2r  = π / 180
    data = HPM.parse_file(joinpath(HPM.BASE_DIR, "examples", "rhhc", "data", "hello_world.m"))
    data["principle"] = "absolute equality"
    # wide angle supports, beyond 360° at the 13th harmonic, where δ = 2 loses variance
    pp = Dict{String,Any}("alpha_hat" => 4.1612, "beta_hat" => 3.5541, "delta" => 2,
        "n_terms" => 4, "epsilon" => 0.05, "bound" => "gaussian", "moments" => "exact",
        "harmonics" => Dict(5 => (-143.07d2r, -9.30d2r), 7 => (-263.84d2r, -52.49d2r),
                            13 => (-180.49d2r, 479.68d2r)))

    function prepared(params)
        hd = HPM.replicate(deepcopy(data); H = [1, 5, 7, 13])
        HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, solver_nlp)
        HPM.update_hdata_with_pce_data!(hd, deepcopy(params))
        return hd
    end

    @testset "network moments agree with Monte Carlo" begin
        hd = prepared(pp)
        m  = HPM.shhc_moments(hd)
        @test m.K == collect(0:3) && all(m.psi_norms .== 1.0)
        ids  = sort([b["bus_i"] for b in values(hd["nw"]["1"]["bus"])])
        bi   = Dict(id => i for (i, id) in enumerate(ids))
        keep = [bi[id] for id in ids if !get(hd["nw"]["1"]["bus"]["$id"], "is_clean", false)]
        pos  = Dict(k => j for (j, k) in enumerate(keep))
        units = sort(hd["pce"]["harmonic_units"])
        n = bi[3]                                   # the end of the cable feeder
        Random.seed!(11); N = 200_000
        ω = rand(Beta(4.1612, 3.5541), N, length(units))
        Q = zeros(N, 3)
        for (i, h) in enumerate((5, 7, 13))
            Z  = inv(HPM._harmonic_admittance(hd["nw"]["$h"], bi)[keep, keep])
            Y = zeros(ComplexF64, N)
            for (j, d) in enumerate(units)
                ld = hd["nw"]["$h"]["load"]["$d"]
                a, b = hd["pce"]["angle_ranges"][h][d]
                Y .+= Z[pos[n], pos[bi[ld["load_bus"]]]] .* cis.(a .+ (b - a) .* ω[:, j])
            end
            Q[:, i] = abs2.(Y)
        end
        C = cov(Q)
        for (i, h) in enumerate((5, 7, 13))
            c = m.alpha_V[h][n, :]
            @test isapprox(c[1], mean(Q[:, i]); rtol = 1e-2)
            @test isapprox(sum(abs2, c[2:end]), C[i, i]; rtol = 3e-2)
        end
        # the total distortion uses the covariance across the harmonics
        L = reduce(vcat, [m.alpha_V[h][n, 2:end]' for h in (5, 7, 13)])
        @test isapprox(L * L', C; atol = 3e-2 * maximum(abs.(C)))
    end

    @testset "point masses: exact and expansion moments coincide" begin
        pm = deepcopy(pp); pm["harmonics"] = Dict(h => (r[1], r[1]) for (h, r) in pp["harmonics"])
        hx = prepared(pm); pm2 = deepcopy(pm); pm2["moments"] = "pce"; hp = prepared(pm2)
        mx, mp = HPM.shhc_moments(hx), HPM.shhc_moments(hp)
        for h in (5, 7, 13)
            @test isapprox(mx.alpha_V[h][:, 1], mp.alpha_V[h][:, 1]; rtol = 1e-10, atol = 1e-14)
        end
    end

    @testset "the calibrated factor refuses exact moments" begin
        @test_throws ArgumentError HPM.solve_shhc_calibrated(HPM.replicate(deepcopy(data);
            H = [1, 5, 7, 13]), sHHC_SOC, solver_soc; pce_params = deepcopy(pp),
            hpf_optimizer = solver_nlp)
    end

    @testset "robust evaluation with exact moments" begin
        hd     = prepared(pp)
        cables = (2, 3); len = 1.0e3
        zy(b) = (p, f, seq) -> begin
            br = data["branch"]["$b"]
            ((br["br_r"] * sqrt(f / 50) + im * br["br_x"] * f / 50) / len,
             im * p.kC * (br["b_fr"] + br["b_to"]) * f / 50 / len)
        end
        model = ExactPiModel(Dict(b => (zy = zy(b), len = len, zbase = 1.0) for b in cables))
        box   = ParameterBox(:kC => (1.0, 0.8, 1.2))
        ce    = ConstraintEvaluator(hd, model, box)
        @test ce.exact
        g = Dict(5 => 1e-4, 7 => 1e-4, 13 => 5e-5)
        for ξ in ([0.0], [0.6], [-1.0])
            m    = moments_at(hd, model, latent(box, ξ))
            full = [normalized_violation(c, g, m) for c in shhc_constraints(hd, m)]
            fast = [fast_violation(ce, j, ξ, g) for j in eachindex(ce.cons)]
            @test maximum(abs.(fast .- full)) < 1e-9
        end
        j = findfirst(c -> c.kind === :thd, ce.cons)
        v, gr = fast_violation_gradient(ce, j, [0.2], g)
        @test v == fast_violation(ce, j, [0.2], g)
        fd = (fast_violation(ce, j, [0.2 + 1e-4], g) - fast_violation(ce, j, [0.2 - 1e-4], g)) / 2e-4
        @test isapprox(gr[1], fd; rtol = 1e-3, atol = 1e-8)

        res = solve_rhhc(HPM.replicate(deepcopy(data); H = [1, 5, 7, 13]), model, box, solver_soc;
                         pce_params = deepcopy(pp), hpf_optimizer = solver_nlp,
                         nbackground = 16, posterior = 32)
        @test res["signed_off"]
        @test res["lower_bound"] <= res["upper_bound"] + 1e-12 && res["gap"] < 1e-3
        @test res["posterior"].max_violation <= 1e-6
    end
end
