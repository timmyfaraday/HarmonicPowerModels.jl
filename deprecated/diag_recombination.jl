################################################################################
# Diagnostic: why does angle uncertainty REDUCE hosting capacity (industrial)?  #
#                                                                               #
# Investigates three hypotheses:                                                #
#   (1) does dHHC@E[φ] already recombine?  → compare det vs stochastic MEAN     #
#   (2) transformer rotation reference     → report per-unit injection angles   #
#   (3) Cantelli too conservative?         → decompose E[S]+λσ, re-solve Gaussian#
# Plus: timing (det slower than unc) → JIT first-call test.                     #
################################################################################

using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf

include(joinpath(@__DIR__, "pce_params.jl"))
const HPM = HarmonicPowerModels
const PMs = PowerModels

ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver" => "ma27", "hsllib" => HSL_jll.libhsl_path,
                "print_level" => 0, "tol" => 1e-10, "constr_viol_tol" => 1e-10,
                "acceptable_tol" => 1e-10, "acceptable_constr_viol_tol" => 1e-10)
mosek_soc = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

data = PMs.parse_file(joinpath(@__DIR__, "data", "industrial_network_shhc.m"))
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]

const α̂ = PCE_PARAMS["alpha_hat"]; const β̂ = PCE_PARAMS["beta_hat"]
const Eω = α̂ / (α̂ + β̂)
Eφ_rad(h) = (hb = PCE_PARAMS["harmonics"][h]; hb.v_min + (hb.v_max - hb.v_min) * Eω)

function shhc_params(f; bound="cantelli")
    p = deepcopy(PCE_PARAMS)
    p["harmonics"] = Dict{Int,Any}(h => f(h) for h in harm)
    p["bound"] = bound
    return p
end
det_params(;b="cantelli") = shhc_params(h -> (Eφ_rad(h), Eφ_rad(h)); bound=b)
unc_params(;b="cantelli") = shhc_params(h -> (PCE_PARAMS["harmonics"][h].v_min,
                                              PCE_PARAMS["harmonics"][h].v_max); bound=b)

get_I(sol, h) = (s = sol["solution"]; haskey(s, "nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])

# ── solve helpers that RETURN the hdata so we can read pce internals ───────────
function run_soc(p)
    hd = HPM.replicate(deepcopy(data), H = H)
    t = @elapsed sol = HPM.solve_shhc(hd, sHHC_SOC, mosek_soc; pce_params=p,
                                      hpf_optimizer=ipopt_ma27)
    return sol, hd, t
end
function run_dhhc()
    hd = HPM.replicate(deepcopy(data), H = H)
    for h in harm, (_, load) in hd["nw"]["$h"]["load"]
        load["ref_angle"] = rad2deg(Eφ_rad(h))
    end
    sol = HPM.solve_hhc(hd, dHHC_SOC, mosek_soc, ipopt_ma27)
    return sol, hd
end

println("Solving …")
sol_dhhc, _      = run_dhhc()
sol_det, hd_det, t_det = run_soc(det_params())     # 1st solve_soc → cold (JIT)
sol_unc, hd_unc, t_unc = run_soc(unc_params())     # 2nd          → warm
sol_unc2, _,     t_unc2 = run_soc(unc_params())    # 3rd          → warm again
sol_g,   hd_g,   _      = run_soc(unc_params(b="gaussian"))

# bus index map (mirror solver)
fund = hd_unc["nw"]["1"]
bus_ids_sorted = sort([b["bus_i"] for (_,b) in fund["bus"]])
bus_idx = Dict(id => i for (i,id) in enumerate(bus_ids_sorted))
ref_bus_ids = Set(bus_idx[b["bus_i"]] for (_,b) in fund["bus"] if b["bus_type"]==3)

pce   = hd_unc["pce"]
psi   = pce["psi_norms"]; K = pce["K"]
λ_c   = HPM.cantelli_lambda(PCE_PARAMS["epsilon"])
λ_g   = HPM.gaussian_lambda(PCE_PARAMS["epsilon"])

# evaluate aggregate S = Σ_h α_V[h][n,:]·g_h as a PCE; return (mean, std)
function agg_stats(alpha_V, n, g)
    Svec = zeros(length(K))
    for h in harm
        haskey(alpha_V, h) || continue
        Svec .+= alpha_V[h][n, :] .* g[h]
    end
    mean = Svec[1]
    var  = sum(psi[k+1]*Svec[k+1]^2 for k in K if k != 0; init=0.0)
    return mean, sqrt(max(var,0.0))
end

# use the dHHC hosting currents for BOTH (apples-to-apples constraint comparison)
g_dhhc = Dict(h => get_I(sol_dhhc,h)^2 for h in harm)

aV_det = hd_det["pce"]["alpha_V"]
aV_unc = hd_unc["pce"]["alpha_V"]

# ── report ────────────────────────────────────────────────────────────────────
outfile = joinpath(@__DIR__, "output", "diag_recombination.txt")
mkpath(dirname(outfile))
open(outfile, "w") do io
    sep(c='─') = println(io, repeat(c, 92))
    println(io, repeat('═',92))
    println(io, "DIAGNOSTIC — angle uncertainty vs hosting capacity (industrial)")
    println(io, repeat('═',92))
    @printf(io, "ε=%.2f   Cantelli λ=%.4f   Gaussian λ=%.4f   (ratio %.2f×)\n\n",
            PCE_PARAMS["epsilon"], λ_c, λ_g, λ_c/λ_g)

    # (A) phasor shrinkage / mean resultant length per harmonic
    println(io, "(1)+(2) MEAN-LEVEL RECOMBINATION — phasor shrinkage per harmonic")
    println(io, "    r = |E[e^{jφ}]| = sqrt(E[cosφ]²+E[sinφ]²);  r<1 ⇒ mean phasor shrinks.")
    @printf(io, "    %-4s %-12s %-12s %-12s %-12s %-12s\n",
            "h","E[cosφ]","E[sinφ]","r","1-r [%]","σ_φ[deg]")
    for h in harm
        cc = pce["cos_phi"][h][1]; sc = pce["sin_phi"][h][1]   # unit 1, k=0 term = mean
        Ec, Es = cc[1], sc[1]
        r = sqrt(Ec^2 + Es^2)
        # angle std from Beta over the range
        hb = PCE_PARAMS["harmonics"][h]; rng = hb.v_max - hb.v_min
        sω = sqrt(α̂*β̂/((α̂+β̂)^2*(α̂+β̂+1)))
        @printf(io, "    %-4d %-12.5f %-12.5f %-12.5f %-12.4f %-12.3f\n",
                h, Ec, Es, r, (1-r)*100, rad2deg(rng*sω))
    end

    # (B) binding-constraint decomposition at the dHHC hosting currents
    println(io)
    println(io, "(3) CHANCE-CONSTRAINT DECOMPOSITION at fixed I = I_dHHC (THD/RMS-V per bus)")
    println(io, "    S = Σ_h |V_{n,h}|²  ;  det = S(E[φ]) ;  E[S],σ(S) under uncertainty.")
    println(io, "    Recombination ⇔ E[S] < det.  Constraint binds on E[S]+λσ ≤ limit.")
    @printf(io, "    %-5s %-7s %-12s %-12s %-12s %-12s %-12s %-12s %-10s\n",
            "bus","lim","det S","E[S]","E[S]/det%","σ(S)","+λc·σ","+λg·σ","limit")
    for (bi,fb) in sort(collect(fund["bus"]); by=x->x[2]["bus_i"])
        n = bus_idx[fb["bus_i"]]
        n in ref_bus_ids && continue
        v = fb["vm"]
        for (tag, lim2) in (("THD", haskey(fb,"thdmax") ? (fb["thdmax"]*v)^2 : NaN),
                            ("RMS", haskey(fb,"vmaxrms") ? fb["vmaxrms"]^2 - v^2 : NaN))
            isnan(lim2) && continue
            md, _   = agg_stats(aV_det, n, g_dhhc)
            mu, sg  = agg_stats(aV_unc, n, g_dhhc)
            @printf(io, "    %-5s %-7s %-12.4e %-12.4e %-+12.3f %-12.4e %-12.4e %-12.4e %-10.4e\n",
                    fb["bus_i"], tag, md, mu, (mu/md-1)*100, sg,
                    mu+λ_c*sg, mu+λ_g*sg, lim2)
        end
    end

    # (C) objective comparison: dHHC vs SOC(cantelli) vs SOC(gaussian)
    println(io)
    println(io, "(3) OBJECTIVE — total hosting capacity Σ_h I_h [pu]")
    td = sum(get_I(sol_dhhc,h) for h in harm)
    tc = sum(get_I(sol_unc,h)  for h in harm)
    tg = sum(get_I(sol_g,h)    for h in harm)
    @printf(io, "    dHHC@E[φ]          = %.6e\n", td)
    @printf(io, "    sHHC_SOC Cantelli  = %.6e   (%+.3f %% vs dHHC)\n", tc, (tc/td-1)*100)
    @printf(io, "    sHHC_SOC Gaussian  = %.6e   (%+.3f %% vs dHHC)\n", tg, (tg/td-1)*100)
    @printf(io, "    %-4s %-16s %-16s %-16s\n","h","I_dHHC","I_Cantelli","I_Gaussian")
    for h in harm
        @printf(io, "    %-4d %-16.6e %-16.6e %-16.6e\n",
                h, get_I(sol_dhhc,h), get_I(sol_unc,h), get_I(sol_g,h))
    end

    # (D) timing — JIT first-call hypothesis
    println(io)
    println(io, "(timing) sHHC_SOC wall-time, identical problem size (P=$(length(K))):")
    @printf(io, "    1st solve (det, cold/JIT) = %.3f s\n", t_det)
    @printf(io, "    2nd solve (unc, warm)     = %.3f s\n", t_unc)
    @printf(io, "    3rd solve (unc, warm)     = %.3f s\n", t_unc2)
    println(io, "    ⇒ det-vs-unc gap is Julia first-call compilation, not problem difficulty.")
    println(io, repeat('═',92))
end
println("\n✓ wrote ", outfile)
run(`cat $outfile`)
