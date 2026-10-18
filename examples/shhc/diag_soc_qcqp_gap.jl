################################################################################
# Diagnostic: SOC vs QCQP gap on IEEE-30 (uncertain).                            #
# Take each solver's optimal hosting I_h and evaluate it against the SOC-EXACT   #
# chance constraints (analytic reduced-Z |V|²/|I|² PCE). If QCQP's I_h violates  #
# them → QCQP optimistic (under-enforcing). If it satisfies them → SOC over-     #
# conservative. Also dumps per-constraint chance/limit ratios.                   #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Distributions, Random
const HPM = HarmonicPowerModels
const PMs = PowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver" => "ma27", "hsllib" => HSL_jll.libhsl_path,
                "print_level" => 0, "tol" => 1e-10, "constr_viol_tol" => 1e-10,
                "acceptable_tol" => 1e-10, "acceptable_constr_viol_tol" => 1e-10)
mosek_soc = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

data = PMs.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc.m"))
H = [1,5,7,11,13]; harm = [5,7,11,13]
α̂, β̂ = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]
base_ranges = Dict(h => PCE_PARAMS["harmonics"][h] for h in harm)
N = 4
all_units = sort(parse.(Int, collect(keys(data["load"]))))
get_I(sol,h) = (s=sol["solution"]; haskey(s,"nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])

# worst-case angles (align at first non-ref unit bus, for a concrete common case)
hd0 = HPM.replicate(deepcopy(data), H=H)
Z, bus_pos = HPM.reduced_harmonic_impedance(hd0)
ubuses = [hd0["nw"]["1"]["load"]["$d"]["load_bus"] for d in all_units[1:N]]
nstar = first(b for b in ubuses if haskey(bus_pos,b))
ang = HPM.constructive_worst_case_angles(Z, bus_pos, hd0, nstar)
unc = HPM.recentre_angle_ranges(base_ranges, ang, α̂, β̂)
pce_params = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
pce_params["harmonics"] = unc; pce_params["harmonic_units"] = all_units[1:N]; pce_params["bound"]="cantelli"

# ── solve SOC (keep hd for alpha_V/alpha_I cache) and QCQP ──────────────────────
hd_soc = HPM.replicate(deepcopy(data), H=H)
sol_soc = HPM.solve_shhc(hd_soc, sHHC_SOC, mosek_soc; pce_params=deepcopy(pce_params), hpf_optimizer=ipopt_ma27)
hd_q = HPM.replicate(deepcopy(data), H=H)
sol_q = HPM.solve_shhc(hd_q, sHHC_QCQP, ipopt_ma27; pce_params=deepcopy(pce_params))

pce = hd_soc["pce"]
alpha_V = pce["alpha_V"]; alpha_I = pce["alpha_I"]
psi = pce["psi_norms"]; lam = pce["lambda"]; K = pce["K"]
fund = hd_soc["nw"]["1"]
bus_ids = sort([b["bus_i"] for (_,b) in fund["bus"]])
bus_idx = Dict(id=>i for (i,id) in enumerate(bus_ids))
ref_ids = Set(bus_idx[b["bus_i"]] for (_,b) in fund["bus"] if b["bus_type"]==3)
Kh = [k for k in K if k!=0]

# chance value E[S]+λσ for an aggregate with per-harmonic coef vectors
chance(coef, g) = begin
    ES = sum(coef[h][1]*g[h] for h in keys(coef))
    var = sum(psi[k+1]*(sum(coef[h][k+1]*g[h] for h in keys(coef)))^2 for k in Kh)
    ES + lam*sqrt(max(var,0.0))
end

# evaluate ALL SOC-exact constraints at hosting vector g=(h=>I_h^2); return list of (name, ratio)
function ratios(g)
    out = Tuple{String,Float64}[]
    # IHD per (bus,h)
    for h in harm, (bi,bus) in hd_soc["nw"]["$h"]["bus"]
        haskey(bus,"ihdmax") || continue
        n = bus_idx[bus["bus_i"]]; haskey(alpha_V,h) || continue
        A = alpha_V[h]
        Cscal = A[n,1] + lam*sqrt(sum(psi[k+1]*A[n,idx]^2 for (idx,k) in enumerate(K) if k!=0; init=0.0))
        v = fund["bus"][bi]["vm"]; lim = (bus["ihdmax"]*v)^2
        lim>0 && push!(out, ("IHD b$(bus["bus_i"]) h$h", Cscal*g[h]/lim))
    end
    # THD + RMS-V per bus
    for (bi,fb) in fund["bus"]
        n = bus_idx[fb["bus_i"]]; n in ref_ids && continue
        coef = Dict(h=>alpha_V[h][n,:] for h in harm if haskey(alpha_V,h))
        v = fb["vm"]
        haskey(fb,"thdmax")  && (lim=(fb["thdmax"]*v)^2; lim>0 && push!(out,("THD b$(fb["bus_i"])", chance(coef,g)/lim)))
        haskey(fb,"vmaxrms") && (lim=fb["vmaxrms"]^2-v^2; lim>0 && push!(out,("RMSV b$(fb["bus_i"])", chance(coef,g)/lim)))
    end
    # RMS-I per branch
    for (b_str,fbr) in fund["branch"]
        haskey(fbr,"c_rating") || continue
        b=parse(Int,b_str)
        coef = Dict(h=>alpha_I[h][b] for h in harm if haskey(alpha_I,h) && haskey(alpha_I[h],b))
        isempty(coef) && continue
        lim = fbr["c_rating"]^2 - get(fbr,"cm_fr",0.0)^2
        lim>0 && push!(out,("RMSI br$b", chance(coef,g)/lim))
    end
    out
end

g_soc = Dict(h=>get_I(sol_soc,h)^2 for h in harm)
g_q   = Dict(h=>get_I(sol_q,h)^2 for h in harm)

println("="^70)
@printf("IEEE-30 N=%d  worst-case align @bus %d\n", N, nstar)
@printf("Totals: SOC=%.6e  QCQP=%.6e  (QCQP/SOC %+.3f%%)\n",
        sum(get_I(sol_soc,h) for h in harm), sum(get_I(sol_q,h) for h in harm),
        (sum(get_I(sol_q,h) for h in harm)/sum(get_I(sol_soc,h) for h in harm)-1)*100)
println("Per-harmonic I_h:")
for h in harm; @printf("  h%-3d SOC=%.6e  QCQP=%.6e\n", h, get_I(sol_soc,h), get_I(sol_q,h)); end

for (tag, g) in (("SOC solution", g_soc), ("QCQP solution", g_q))
    rs = sort(ratios(g); by=x->-x[2])
    @printf("\n%s — top chance/limit ratios vs SOC-EXACT constraints:\n", tag)
    @printf("  max ratio = %.4f  (>1.0 ⇒ VIOLATES the exact constraint)\n", rs[1][2])
    for (nm,r) in rs[1:min(6,length(rs))]; @printf("    %-16s %.4f\n", nm, r); end
end
println("="^70)

# ── Monte-Carlo ground truth (exact trig + reduced Z) at the binding bus ─────────
# Single germ per unit (comonotone harmonics, matching ρ≥0.97); exact cos/sin.
# Empirical P[|V_{n*,h}|² > limit] should be ≤ ε=0.05 for a correct model.
Random.seed!(1)
nsamp = 400_000
ε = PCE_PARAMS["epsilon"]
units = all_units[1:N]
ubus  = Dict(d => hd0["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)
np3   = bus_pos[nstar]
Beta  = Distributions.Beta(α̂, β̂)

# precompute complex Z[n*,d] per harmonic and the recentred (lo,hi) per (h,unit)
function mc_violation(g)   # g: h => I_h^2 ; returns Dict h => P[violation]
    out = Dict{Int,Float64}()
    for h in harm
        haskey(alpha_V,h) || continue
        bus3_str = "$nstar"
        haskey(hd_soc["nw"]["$h"]["bus"], bus3_str) || continue
        bus = hd_soc["nw"]["$h"]["bus"][bus3_str]
        haskey(bus,"ihdmax") || (out[h]=NaN; continue)
        v1  = fund["bus"][bus3_str]["vm"]
        lim = (bus["ihdmax"]*v1)^2
        Ih  = sqrt(g[h])
        zs  = Dict(d => (haskey(bus_pos,ubus[d]) ? Z[h][np3, bus_pos[ubus[d]]] : 0.0+0im) for d in units)
        los = Dict(d => unc[h][d][1] for d in units)
        wid = Dict(d => unc[h][d][2]-unc[h][d][1] for d in units)
        viol = 0
        for _ in 1:nsamp
            Vr = 0.0; Vi = 0.0
            for d in units
                ω = rand(Beta); φ = los[d] + wid[d]*ω
                zr, zi = real(zs[d]), imag(zs[d])
                c, s = cos(φ), sin(φ)
                Vr += Ih*(zr*c - zi*s); Vi += Ih*(zi*c + zr*s)
            end
            (Vr*Vr + Vi*Vi) > lim && (viol += 1)
        end
        out[h] = viol/nsamp
    end
    out
end

println("\nMonte-Carlo (n=$nsamp) — empirical P[IHD violated] at bus $nstar (target ε=$ε):")
@printf("  %-4s  %-14s  %-14s\n", "h", "P(viol)|SOC I_h", "P(viol)|QCQP I_h")
vs = mc_violation(g_soc); vq = mc_violation(g_q)
for h in harm
    @printf("  %-4d  %-14.4f  %-14.4f\n", h, get(vs,h,NaN), get(vq,h,NaN))
end
println("  (SOC ≈ ε ⇒ SOC correct/tight; QCQP ≫ ε ⇒ QCQP optimistic / under-enforces.)")
println("="^70)
