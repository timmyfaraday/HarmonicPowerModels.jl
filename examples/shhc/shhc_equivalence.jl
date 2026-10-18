################################################################################
# HarmonicPowerModels.jl — shared SHHC equivalence test-case machinery           #
#                                                                               #
# Provides the two equivalence demonstrations used by shhc_industrial.jl and     #
# shhc_ieee30.jl:                                                                #
#   deterministic_equivalence(cfg) — at E[φ], dHHC ≡ sHHC_SOC ≡ sHHC_QCQP in     #
#       objective, injected current phasors and bus voltage phasors (full data    #
#       + phasor CSV for plotting).                                              #
#   stochastic_equivalence(cfg)    — sHHC_SOC ≡ sHHC_QCQP for the empirical,      #
#       Gaussian and Cantelli chance bounds, cross-checked against a Monte-Carlo  #
#       dHHC objective distribution (CSV + overlay refs for plotting).           #
#                                                                               #
# `cfg` is a NamedTuple built by each network script (see those files).          #
################################################################################

using HarmonicPowerModels, PowerModels
using JuMP, Printf, Distributions, Random, Statistics

const HPM = HarmonicPowerModels
const PMs = PowerModels

# ── small helpers ───────────────────────────────────────────────────────────────
_polar(z) = (abs(z), rad2deg(angle(z)))
get_I(sol, h) = (s = sol["solution"]; haskey(s, "nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])
_obj(sol) = get(sol, "objective_value", get(sol, "objective", NaN))
status(sol) = string(get(sol, "termination_status", "N/A"))

# E[ω] and per-harmonic expected angle [rad]
Eω(cfg) = cfg.α̂ / (cfg.α̂ + cfg.β̂)
function Eφ_rad(cfg, h)
    lo, hi = cfg.base_ranges[h]
    return lo + (hi - lo) * Eω(cfg)
end

# pce_params with explicit ranges + bound + (optional) empirical λ*
function _pce(cfg, ranges; bound = "cantelli", emp = nothing)
    p = deepcopy(cfg.pce_params)
    p["harmonics"]      = ranges
    p["harmonic_units"] = cfg.units
    p["bound"]          = bound
    emp === nothing || (p["empirical_lambdas"] = emp)
    return p
end
det_ranges(cfg) = Dict{Int,Any}(h => (Eφ_rad(cfg, h), Eφ_rad(cfg, h)) for h in cfg.harm)
phys_ranges(cfg) = Dict{Int,Any}(h => cfg.base_ranges[h] for h in cfg.harm)

# ── solve wrappers (fresh replicate each call) ──────────────────────────────────
function _rep(cfg)
    hd = HPM.replicate(deepcopy(cfg.data), H = cfg.H)
    # optional angle reference: "global" (network frame) or "local" (relative to
    # each unit's fundamental bus voltage, see update_hdata_with_angle_offsets!)
    hasproperty(cfg, :angle_reference) && (hd["angle_reference"] = cfg.angle_reference)
    # optional post-replicate hook (e.g. to_single_ground!) applied before unit restriction
    hasproperty(cfg, :postrep) && cfg.postrep !== nothing && cfg.postrep(hd)
    if cfg.restrict_units                       # IEEE-30: keep only the first-N units
        uset = Set(cfg.units)
        for h in cfg.harm, k in collect(keys(hd["nw"]["$h"]["load"]))
            parse(Int, k) in uset || delete!(hd["nw"]["$h"]["load"], k)
        end
    end
    return hd
end

function solve_dhhc(cfg, ranges)               # deterministic point-mass angles
    hd = _rep(cfg)
    for h in cfg.harm, (k, l) in hd["nw"]["$h"]["load"]
        d = parse(Int, k)
        haskey(ranges, h) && haskey(ranges[h], d) ?
            (l["ref_angle"] = rad2deg(ranges[h][d][1])) :
            (l["ref_angle"] = rad2deg(Eφ_rad(cfg, h)))
    end
    sol = HPM.solve_hhc(hd, dHHC_SOC, cfg.mosek, cfg.ipopt)
    return sol, hd
end
solve_dhhc_pm(cfg) = solve_dhhc(cfg, Dict(h => Dict(d => (Eφ_rad(cfg,h),) for d in cfg.units) for h in cfg.harm))

function solve_soc(cfg, p)
    hd = _rep(cfg)
    sol = HPM.solve_shhc(hd, sHHC_SOC, cfg.mosek; pce_params = p, hpf_optimizer = cfg.ipopt)
    return sol, hd
end
function solve_qcqp(cfg, p)
    hd = _rep(cfg)
    sol = HPM.solve_shhc(hd, sHHC_QCQP, cfg.ipopt; pce_params = p)
    return sol, hd
end

# ── phasor extraction (native where available) ──────────────────────────────────
# returns (cur::Dict{(unit,h)=>Complex}, vol::Dict{(bus,h)=>Complex})
function phasors_dhhc(cfg, sol)
    cur = Dict{Tuple{Int,Int},ComplexF64}(); vol = Dict{Tuple{Int,Int},ComplexF64}()
    for h in cfg.harm
        nw = sol["solution"]["nw"]["$h"]
        for d in cfg.units
            ld = nw["load"]["$d"]; cur[(d,h)] = complex(ld["crd"], ld["cid"])
        end
        for (bi, b) in nw["bus"]; vol[(parse(Int,bi),h)] = complex(b["vr"], b["vi"]); end
    end
    return cur, vol
end
function phasors_qcqp(cfg, sol)
    cur = Dict{Tuple{Int,Int},ComplexF64}(); vol = Dict{Tuple{Int,Int},ComplexF64}()
    g0(d) = d["0"]
    for h in cfg.harm
        nw = sol["solution"]["nw"]["$h"]
        for d in cfg.units
            ld = nw["load"]["$d"]; cur[(d,h)] = complex(g0(ld["crd_pce"]), g0(ld["cid_pce"]))
        end
        for (bi, b) in nw["bus"]; vol[(parse(Int,bi),h)] = complex(g0(b["vr_pce"]), g0(b["vi_pce"])); end
    end
    return cur, vol
end
# SOC: reconstruct from I_h and the cached voltage-response PCE (k=0 = deterministic)
function phasors_soc(cfg, sol, hd)
    pce = hd["pce"]; ar = pce["alpha_re"]; ai = pce["alpha_im"]
    cph = pce["cos_phi"]; sph = pce["sin_phi"]
    bus_ids = sort([b["bus_i"] for (_,b) in hd["nw"]["1"]["bus"]])
    bidx = Dict(id => i for (i,id) in enumerate(bus_ids))
    cur = Dict{Tuple{Int,Int},ComplexF64}(); vol = Dict{Tuple{Int,Int},ComplexF64}()
    for h in cfg.harm
        Ih = get_I(sol, h)
        for d in cfg.units
            cur[(d,h)] = Ih * complex(cph[h][d][1], sph[h][d][1])
        end
        haskey(ar, h) || continue
        for id in bus_ids
            vol[(id,h)] = Ih * complex(ar[h][bidx[id],1], ai[h][bidx[id],1])
        end
    end
    return cur, vol
end

# ══════════════════════════════════════════════════════════════════════════════
# 1. Deterministic equivalence
# ══════════════════════════════════════════════════════════════════════════════
function deterministic_equivalence(cfg)
    println("[$(cfg.name)] deterministic equivalence at E[φ]…")
    dr = det_ranges(cfg)
    sol_d, hd_d = solve_dhhc_pm(cfg)
    sol_s, hd_s = solve_soc(cfg, _pce(cfg, dr; bound="empirical"))   # bound irrelevant (σ=0)
    sol_q, _    = solve_qcqp(cfg, _pce(cfg, dr; bound="empirical"))

    cur_d, vol_d = phasors_dhhc(cfg, sol_d)
    cur_q, vol_q = phasors_qcqp(cfg, sol_q)
    cur_s, vol_s = phasors_soc(cfg, sol_s, hd_s)

    bus_ids = sort([b["bus_i"] for (_,b) in hd_d["nw"]["1"]["bus"]])
    od, os, oq = _obj(sol_d), _obj(sol_s), _obj(sol_q)

    # ── phasor CSV (for plotting) ──────────────────────────────────────────────
    csv = joinpath(cfg.outdir, "$(cfg.name)_det_phasors.csv")
    open(csv, "w") do io
        println(io, "method,quantity,id,h,re,im,mag,ang_deg")
        for (meth, cur, vol) in (("dHHC",cur_d,vol_d), ("SOC",cur_s,vol_s), ("QCQP",cur_q,vol_q))
            for h in cfg.harm, d in cfg.units
                z = cur[(d,h)]; m,a = _polar(z)
                println(io, "$meth,I,$d,$h,$(real(z)),$(imag(z)),$m,$a")
            end
            for h in cfg.harm, n in bus_ids
                haskey(vol,(n,h)) || continue
                z = vol[(n,h)]; m,a = _polar(z)
                println(io, "$meth,V,$n,$h,$(real(z)),$(imag(z)),$m,$a")
            end
        end
    end

    # ── text report (full data) ────────────────────────────────────────────────
    txt = joinpath(cfg.outdir, "$(cfg.name)_det_equivalence.txt")
    open(txt, "w") do io
        sep(c='═')=println(io,repeat(c,92))
        sep(); println(io,"DETERMINISTIC EQUIVALENCE @E[φ] — $(uppercase(cfg.name))"); sep()
        println(io,"Units=$(length(cfg.units))  Harmonics=$(cfg.harm)")
        println(io,"Statuses: dHHC=$(status(sol_d))  SOC=$(status(sol_s))  QCQP=$(status(sol_q))")
        @printf(io,"\nObjective Σ_units Σ_h I_h [pu]:\n")
        @printf(io,"  dHHC=%.10e\n  SOC =%.10e  (Δ=%+.3e%%)\n  QCQP=%.10e  (Δ=%+.3e%%)\n",
                od, os, (os/od-1)*100, oq, (oq/od-1)*100)

        println(io,"\nInjected current phasors  I_{unit,h}  [mag∠deg]:")
        @printf(io,"  %-4s %-4s | %-22s | %-22s | %-22s\n","unit","h","dHHC","SOC","QCQP")
        for d in cfg.units, h in cfg.harm
            md,ad=_polar(cur_d[(d,h)]); ms,as=_polar(cur_s[(d,h)]); mq,aq=_polar(cur_q[(d,h)])
            @printf(io,"  %-4d %-4d | %10.4e∠%7.2f | %10.4e∠%7.2f | %10.4e∠%7.2f\n",
                    d,h,md,ad,ms,as,mq,aq)
        end
        println(io,"\nBus voltage phasors  V_{bus,h}  [mag∠deg]  (non-zero buses):")
        @printf(io,"  %-4s %-4s | %-22s | %-22s | %-22s\n","bus","h","dHHC","SOC","QCQP")
        maxdev=0.0
        for n in bus_ids, h in cfg.harm
            haskey(vol_d,(n,h)) && haskey(vol_q,(n,h)) || continue
            (abs(vol_d[(n,h)])<1e-12 && abs(vol_q[(n,h)])<1e-12) && continue
            md,ad=_polar(vol_d[(n,h)]); ms,as=_polar(get(vol_s,(n,h),0.0im)); mq,aq=_polar(vol_q[(n,h)])
            @printf(io,"  %-4d %-4d | %10.4e∠%7.2f | %10.4e∠%7.2f | %10.4e∠%7.2f\n",
                    n,h,md,ad,ms,as,mq,aq)
            maxdev=max(maxdev, abs(vol_q[(n,h)]-vol_d[(n,h)])/max(abs(vol_d[(n,h)]),1e-12))
        end
        @printf(io,"\nMax |V| relative deviation QCQP vs dHHC: %.3e %%\n", maxdev*100)
        println(io,"Phasor CSV: $(basename(csv))")
        sep()
    end
    println("  ✓ $(basename(txt)) + $(basename(csv))  | obj dHHC=$(round(od,sigdigits=7)) SOC Δ=$(round((os/od-1)*100,sigdigits=3))% QCQP Δ=$(round((oq/od-1)*100,sigdigits=3))%")
    return (; od, os, oq, txt, csv)
end

# ══════════════════════════════════════════════════════════════════════════════
# 2. Stochastic equivalence (SOC ≡ QCQP for empirical/gaussian/cantelli) + MC
# ══════════════════════════════════════════════════════════════════════════════
function _solve_bound(cfg, bound)
    pr = phys_ranges(cfg)
    if bound == "empirical"
        ss = HPM.solve_shhc_calibrated(_rep(cfg), sHHC_SOC, cfg.mosek;
                 pce_params=_pce(cfg,pr), hpf_optimizer=cfg.ipopt, nsamp=cfg.nsamp_lambda)
        sq = cfg.qcqp_ok ?
             HPM.solve_shhc_calibrated(_rep(cfg), sHHC_QCQP, cfg.ipopt;
                 pce_params=_pce(cfg,pr), nsamp=cfg.nsamp_lambda) : nothing
    else
        ss, _ = solve_soc(cfg, _pce(cfg, pr; bound=bound))
        sq = cfg.qcqp_ok ? solve_qcqp(cfg, _pce(cfg, pr; bound=bound))[1] : nothing
    end
    return ss, sq
end

function stochastic_equivalence(cfg; bounds = ["empirical","gaussian","cantelli"])
    println("[$(cfg.name)] stochastic equivalence (SOC≡QCQP per bound) + MC…")
    # solve all bounds
    res = Dict{String,Any}()
    for b in bounds
        ss, sq = _solve_bound(cfg, b)
        os = sum(get_I(ss,h) for h in cfg.harm)
        oq = sq === nothing ? NaN : sum(get_I(sq,h) for h in cfg.harm)
        res[b] = (soc=os, qcqp=oq, st_s=status(ss), st_q=sq===nothing ? "skipped" : status(sq))
        println("    $b: SOC=$(round(os,sigdigits=6))  QCQP=$(sq===nothing ? "—" : round(oq,sigdigits=6))")
    end

    # deterministic dHHC at E[φ] (≡ sHHC point-mass) — the no-uncertainty reference
    H_Eφ = sum(get_I(solve_dhhc_pm(cfg)[1], h) for h in cfg.harm)

    # ── Monte-Carlo dHHC objective distribution ─────────────────────────────────
    println("    MC ($(cfg.nmc) samples)…")
    Random.seed!(1); Bd = Distributions.Beta(cfg.α̂, cfg.β̂)
    objmc = Vector{Float64}(undef, cfg.nmc)
    for m in 1:cfg.nmc
        ω = Dict(d => rand(Bd) for d in cfg.units)
        rng = Dict(h => Dict(d => (cfg.base_ranges[h][1] + (cfg.base_ranges[h][2]-cfg.base_ranges[h][1])*ω[d],) for d in cfg.units) for h in cfg.harm)
        sol, _ = solve_dhhc(cfg, rng)
        objmc[m] = sum(get_I(sol,h) for h in cfg.harm)
    end
    q(p)=quantile(objmc,p)

    # ── outputs ─────────────────────────────────────────────────────────────────
    mccsv = joinpath(cfg.outdir, "$(cfg.name)_stoch_mc.csv")
    open(mccsv,"w") do io; println(io,"sample,objective"); for m in 1:cfg.nmc; println(io,m,",",objmc[m]); end; end
    refcsv = joinpath(cfg.outdir, "$(cfg.name)_stoch_refs.csv")
    open(refcsv,"w") do io
        println(io,"name,value")
        for (n,v) in (("q05",q(.05)),("median",q(.5)),("q95",q(.95)),("mean",mean(objmc)),
                      ("dHHC_Ephi",H_Eφ))
            println(io,n,",",v)
        end
        for b in bounds
            println(io,"SOC_$b,",res[b].soc)
            isnan(res[b].qcqp) || println(io,"QCQP_$b,",res[b].qcqp)
        end
    end

    txt = joinpath(cfg.outdir, "$(cfg.name)_stoch_equivalence.txt")
    open(txt,"w") do io
        sep(c='═')=println(io,repeat(c,84))
        sep(); println(io,"STOCHASTIC EQUIVALENCE — $(uppercase(cfg.name))"); sep()
        println(io,"Units=$(length(cfg.units))  ε=$(cfg.ε)  MC samples=$(cfg.nmc)")
        @printf(io,"\nsHHC_SOC vs sHHC_QCQP per chance bound  [Σ_h I_h, pu]:\n")
        @printf(io,"  %-10s %-16s %-16s %-12s %-s\n","bound","SOC","QCQP","SOC-QCQP%","status")
        for b in bounds
            r=res[b]; d = isnan(r.qcqp) ? NaN : (r.soc/r.qcqp-1)*100
            @printf(io,"  %-10s %-16.8e %-16s %-12s %s/%s\n", b, r.soc,
                    isnan(r.qcqp) ? "—" : @sprintf("%.8e",r.qcqp),
                    isnan(d) ? "—" : @sprintf("%+.4f",d), r.st_s, r.st_q)
        end
        @printf(io,"\nMonte-Carlo dHHC objective distribution [pu]:\n")
        @printf(io,"  min=%.6e q05=%.6e median=%.6e q95=%.6e max=%.6e mean=%.6e\n",
                minimum(objmc),q(.05),q(.5),q(.95),maximum(objmc),mean(objmc))
        vp(h)=100*count(<(h),objmc)/cfg.nmc
        @printf(io,"\ndeterministic dHHC@E[φ] = %.6e  (MC pct %.1f%%, viol %.1f%%)\n", H_Eφ, vp(H_Eφ), vp(H_Eφ))
        @printf(io,"\nsHHC objective vs MC (q05 = 95%%-reliable target; viol = P[H<obj]):\n")
        for b in bounds
            @printf(io,"  SOC  %-10s = %.6e  (vs q05 %+.2f%%, viol %.1f%%)\n", b, res[b].soc,(res[b].soc/q(.05)-1)*100, vp(res[b].soc))
            isnan(res[b].qcqp) || @printf(io,"  QCQP %-10s = %.6e  (vs q05 %+.2f%%, viol %.1f%%)\n", b, res[b].qcqp,(res[b].qcqp/q(.05)-1)*100, vp(res[b].qcqp))
        end
        println(io,"\nMC CSV: $(basename(mccsv))   refs: $(basename(refcsv))")
        sep()
    end
    println("  ✓ $(basename(txt)) + $(basename(mccsv))")
    return (; res, objmc, txt, mccsv, refcsv)
end
